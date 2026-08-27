#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "geo_wal.h"

#include "geo_index_io.h"
#include "geo_index_persistence.h"

#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#define GEO_WAL_VERSION 1U
#define GEO_WAL_FRAME_VERSION 1U
#define GEO_WAL_FRAME_MAGIC UINT32_C(0x31465747)
#define GEO_WAL_MIN_SEGMENT_SIZE (64U * 1024U)
#define GEO_WAL_MAX_PAYLOAD_SIZE (64U * 1024U * 1024U)
#define GEO_WAL_FILE_PREFIX "wal-"
#define GEO_WAL_FILE_SUFFIX ".gbw"

typedef struct {
    char magic[8];
    uint32_t version;
    uint32_t endian_marker;
    uint32_t header_size;
    uint32_t reserved;
    uint64_t segment_id;
    uint64_t first_sequence;
    uint64_t checksum;
    uint8_t padding[16];
} GeoWalFileHeader;

typedef struct {
    uint32_t magic;
    uint16_t version;
    uint16_t header_size;
    uint32_t payload_size;
    uint32_t entry_count;
    uint64_t first_sequence;
    uint64_t payload_checksum;
    uint64_t header_checksum;
    uint64_t reserved;
} GeoWalFrameHeader;

_Static_assert(sizeof(GeoWalFileHeader) == 64U, "WAL file header layout is persisted");
_Static_assert(offsetof(GeoWalFileHeader, segment_id) == 24U, "WAL segment ID offset is persisted");
_Static_assert(offsetof(GeoWalFileHeader, checksum) == 40U, "WAL file checksum offset is persisted");
_Static_assert(sizeof(GeoWalFrameHeader) == 48U, "WAL frame header layout is persisted");
_Static_assert(offsetof(GeoWalFrameHeader, first_sequence) == 16U, "WAL sequence offset is persisted");
_Static_assert(offsetof(GeoWalFrameHeader, header_checksum) == 32U, "WAL frame checksum offset is persisted");

struct GeoWal {
    char *directory;
    char *active_path;
    int descriptor;
    off_t write_offset;
    size_t segment_size;
    uint64_t first_segment_id;
    uint64_t segment_id;
    uint64_t next_sequence;
};

static const char GEO_WAL_MAGIC[8] = { 'G', 'B', 'W', 'A', 'L', '1', '\0', '\0' };

static bool wal_write_all(int descriptor, const void *data, size_t size, off_t offset)
{
    const unsigned char *bytes = data;

    while (size) {
        ssize_t written = pwrite(descriptor, bytes, size, offset);

        if (written < 0 && errno == EINTR) {
            continue;
        }

        if (written <= 0) {
            return false;
        }

        bytes += (size_t) written;
        size -= (size_t) written;
        offset += written;
    }

    return true;
}

static ssize_t wal_read_some(int descriptor, void *data, size_t size, off_t offset)
{
    for (;;) {
        ssize_t received = pread(descriptor, data, size, offset);

        if (received < 0 && errno == EINTR) {
            continue;
        }

        return received;
    }
}

static char *wal_segment_path(const char *directory, uint64_t segment_id)
{
    const char *separator = directory[strlen(directory) - 1U] == '/' ? "" : "/";
    int length = snprintf(NULL,
                          0,
                          "%s%s" GEO_WAL_FILE_PREFIX "%020" PRIu64 GEO_WAL_FILE_SUFFIX,
                          directory,
                          separator,
                          segment_id);

    if (length < 0) {
        return NULL;
    }

    char *path = malloc((size_t) length + 1U);

    if (!path) {
        return NULL;
    }

    snprintf(path,
             (size_t) length + 1U,
             "%s%s" GEO_WAL_FILE_PREFIX "%020" PRIu64 GEO_WAL_FILE_SUFFIX,
             directory,
             separator,
             segment_id);

    return path;
}

static uint64_t wal_file_header_checksum(const GeoWalFileHeader *header)
{
    GeoWalFileHeader copy = *header;

    copy.checksum = 0;

    return geo_persisted_checksum_update(geo_persisted_checksum_initial(), &copy, sizeof(copy));
}

static uint64_t wal_frame_header_checksum(const GeoWalFrameHeader *header)
{
    GeoWalFrameHeader copy = *header;

    copy.header_checksum = 0;

    return geo_persisted_checksum_update(geo_persisted_checksum_initial(), &copy, sizeof(copy));
}

static bool wal_create_segment(GeoWal *wal, uint64_t segment_id, uint64_t first_sequence)
{
    char *path = wal_segment_path(wal->directory, segment_id);

    if (!path) {
        return false;
    }

    int descriptor = open(path, O_RDWR | O_CREAT | O_EXCL, S_IRUSR | S_IWUSR);
    GeoWalFileHeader header = {
        .version = GEO_WAL_VERSION,
        .endian_marker = GEO_FILE_ENDIAN_MARKER,
        .header_size = sizeof(GeoWalFileHeader),
        .segment_id = segment_id,
        .first_sequence = first_sequence,
    };

    memcpy(header.magic, GEO_WAL_MAGIC, sizeof(header.magic));
    header.checksum = wal_file_header_checksum(&header);

    bool succeeded = descriptor >= 0 &&
                     wal_write_all(descriptor, &header, sizeof(header), 0) &&
                     fsync(descriptor) == 0 &&
                     geo_io_sync_parent_directory(path);

    if (!succeeded) {
        if (descriptor >= 0) {
            close(descriptor);
            unlink(path);
        }

        free(path);

        return false;
    }

    if (wal->descriptor >= 0) {
        close(wal->descriptor);
    }

    free(wal->active_path);
    wal->active_path = path;
    wal->descriptor = descriptor;
    wal->write_offset = sizeof(header);
    wal->segment_id = segment_id;

    return true;
}

static bool wal_parse_segment_name(const char *name, uint64_t *segment_id)
{
    size_t prefix_size = sizeof(GEO_WAL_FILE_PREFIX) - 1U;
    size_t suffix_size = sizeof(GEO_WAL_FILE_SUFFIX) - 1U;
    size_t length = strlen(name);

    if (length != prefix_size + 20U + suffix_size ||
        memcmp(name, GEO_WAL_FILE_PREFIX, prefix_size) != 0 ||
        memcmp(name + length - suffix_size, GEO_WAL_FILE_SUFFIX, suffix_size) != 0) {
        return false;
    }

    uint64_t value = 0;

    for (size_t position = prefix_size; position < prefix_size + 20U; ++position) {
        unsigned digit = (unsigned) (name[position] - '0');

        if (digit > 9U || value > (UINT64_MAX - digit) / 10U) {
            return false;
        }

        value = value * 10U + digit;
    }

    if (!value) {
        return false;
    }

    *segment_id = value;

    return true;
}

static bool wal_find_segments(const char *directory, uint64_t *first_segment, uint64_t *last_segment)
{
    DIR *stream = opendir(directory);

    if (!stream) {
        return false;
    }

    uint64_t first = UINT64_MAX;
    uint64_t last = 0;
    struct dirent *entry;

    errno = 0;

    while ((entry = readdir(stream)) != NULL) {
        uint64_t segment_id;

        if (!wal_parse_segment_name(entry->d_name, &segment_id)) {
            continue;
        }

        if (segment_id < first) {
            first = segment_id;
        }

        if (segment_id > last) {
            last = segment_id;
        }
    }

    bool succeeded = errno == 0 && closedir(stream) == 0;

    if (!succeeded) {
        return false;
    }

    *first_segment = last ? first : 0;
    *last_segment = last;

    return true;
}

static bool wal_replay_segment(GeoWal *wal,
                               uint64_t segment_id,
                               bool last_segment,
                               GeoWalReplayFunction replay,
                               void *replay_context)
{
    char *path = wal_segment_path(wal->directory, segment_id);
    int descriptor = path ? open(path, last_segment ? O_RDWR : O_RDONLY) : -1;
    GeoWalFileHeader file_header;
    ssize_t header_bytes = descriptor >= 0 ? wal_read_some(descriptor, &file_header, sizeof(file_header), 0) : -1;
    bool succeeded = header_bytes == (ssize_t) sizeof(file_header) &&
                     memcmp(file_header.magic, GEO_WAL_MAGIC, sizeof(file_header.magic)) == 0 &&
                     file_header.version == GEO_WAL_VERSION &&
                     file_header.endian_marker == GEO_FILE_ENDIAN_MARKER &&
                     file_header.header_size == sizeof(file_header) &&
                     file_header.reserved == 0 &&
                     file_header.segment_id == segment_id &&
                     (!wal->next_sequence || file_header.first_sequence == wal->next_sequence) &&
                     file_header.checksum == wal_file_header_checksum(&file_header);
    off_t offset = sizeof(file_header);

    if (succeeded && !wal->next_sequence) {
        wal->next_sequence = file_header.first_sequence;
    }

    while (succeeded) {
        GeoWalFrameHeader frame;
        ssize_t received = wal_read_some(descriptor, &frame, sizeof(frame), offset);

        if (received == 0) {
            break;
        }

        if (received != (ssize_t) sizeof(frame)) {
            succeeded = last_segment && ftruncate(descriptor, offset) == 0;
            break;
        }

        bool frame_valid = frame.magic == GEO_WAL_FRAME_MAGIC &&
                           frame.version == GEO_WAL_FRAME_VERSION &&
                           frame.header_size == sizeof(frame) &&
                           frame.payload_size <= GEO_WAL_MAX_PAYLOAD_SIZE &&
                           frame.entry_count > 0 &&
                           frame.first_sequence == wal->next_sequence &&
                           frame.entry_count <= UINT64_MAX - frame.first_sequence &&
                           frame.reserved == 0 &&
                           frame.header_checksum == wal_frame_header_checksum(&frame);

        if (!frame_valid) {
            succeeded = false;
            break;
        }

        void *payload = frame.payload_size ? malloc(frame.payload_size) : NULL;
        ssize_t payload_bytes = frame.payload_size
                                    ? wal_read_some(descriptor, payload, frame.payload_size, offset + (off_t) sizeof(frame))
                                    : 0;

        if ((frame.payload_size && !payload) || payload_bytes != (ssize_t) frame.payload_size) {
            free(payload);
            succeeded = last_segment && ftruncate(descriptor, offset) == 0;
            break;
        }

        uint64_t payload_checksum = geo_persisted_checksum_update(geo_persisted_checksum_initial(),
                                                                   payload,
                                                                   frame.payload_size);

        if (payload_checksum != frame.payload_checksum ||
            (replay && !replay(replay_context,
                               frame.first_sequence,
                               frame.entry_count,
                               payload,
                               frame.payload_size))) {
            free(payload);
            succeeded = false;
            break;
        }

        free(payload);
        wal->next_sequence += frame.entry_count;
        offset += (off_t) sizeof(frame) + frame.payload_size;
    }

    if (succeeded && last_segment) {
        wal->descriptor = descriptor;
        wal->active_path = path;
        wal->write_offset = offset;
        wal->segment_id = segment_id;
        descriptor = -1;
        path = NULL;
    }

    if (descriptor >= 0) {
        close(descriptor);
    }

    free(path);

    return succeeded;
}

GeoWal *geo_wal_open(const char *directory,
                     size_t segment_size,
                     GeoWalReplayFunction replay,
                     void *replay_context,
                     uint64_t *next_sequence)
{
    if (!directory || !directory[0] || segment_size < GEO_WAL_MIN_SEGMENT_SIZE || !next_sequence) {
        return NULL;
    }

    struct stat status;

    if (stat(directory, &status) != 0 || !S_ISDIR(status.st_mode)) {
        return NULL;
    }

    GeoWal *wal = calloc(1, sizeof(*wal));

    if (!wal) {
        return NULL;
    }

    wal->directory = strdup(directory);
    wal->descriptor = -1;
    wal->segment_size = segment_size;
    wal->next_sequence = 0;
    uint64_t first_segment = 0;
    uint64_t last_segment = 0;
    bool succeeded = wal->directory && wal_find_segments(directory, &first_segment, &last_segment);

    if (succeeded && last_segment) {
        wal->first_segment_id = first_segment;

        for (uint64_t segment = first_segment; succeeded && segment <= last_segment; ++segment) {
            succeeded = wal_replay_segment(wal,
                                           segment,
                                           segment == last_segment,
                                           replay,
                                           replay_context);

            if (segment == UINT64_MAX) {
                break;
            }
        }
    } else if (succeeded) {
        wal->first_segment_id = 1U;
        wal->next_sequence = 1U;
        succeeded = wal_create_segment(wal, 1U, wal->next_sequence);
    }

    if (!succeeded) {
        geo_wal_close(wal);

        return NULL;
    }

    *next_sequence = wal->next_sequence;

    return wal;
}

void geo_wal_close(GeoWal *wal)
{
    if (!wal) {
        return;
    }

    if (wal->descriptor >= 0) {
        close(wal->descriptor);
    }

    free(wal->active_path);
    free(wal->directory);
    free(wal);
}

bool geo_wal_append(GeoWal *wal,
                    uint64_t first_sequence,
                    uint32_t entry_count,
                    const void *payload,
                    uint32_t payload_size)
{
    if (!wal || wal->descriptor < 0 || !entry_count || (!payload && payload_size) ||
        payload_size > GEO_WAL_MAX_PAYLOAD_SIZE || first_sequence != wal->next_sequence ||
        entry_count > UINT64_MAX - first_sequence) {
        return false;
    }

    size_t frame_size = sizeof(GeoWalFrameHeader) + (size_t) payload_size;

    if (wal->write_offset > (off_t) sizeof(GeoWalFileHeader) &&
        frame_size > wal->segment_size - (size_t) wal->write_offset) {
        if (wal->segment_id == UINT64_MAX || fdatasync(wal->descriptor) != 0 ||
            !wal_create_segment(wal, wal->segment_id + 1U, first_sequence)) {
            return false;
        }
    }

    GeoWalFrameHeader frame = {
        .magic = GEO_WAL_FRAME_MAGIC,
        .version = GEO_WAL_FRAME_VERSION,
        .header_size = sizeof(GeoWalFrameHeader),
        .payload_size = payload_size,
        .entry_count = entry_count,
        .first_sequence = first_sequence,
        .payload_checksum = geo_persisted_checksum_update(geo_persisted_checksum_initial(), payload, payload_size),
    };

    frame.header_checksum = wal_frame_header_checksum(&frame);
    off_t payload_offset = wal->write_offset + (off_t) sizeof(frame);
    bool succeeded = wal_write_all(wal->descriptor, &frame, sizeof(frame), wal->write_offset) &&
                     (!payload_size || wal_write_all(wal->descriptor, payload, payload_size, payload_offset));

    if (succeeded) {
        wal->write_offset = payload_offset + payload_size;
        wal->next_sequence += entry_count;
    }

    return succeeded;
}

bool geo_wal_sync(GeoWal *wal)
{
    return wal && wal->descriptor >= 0 && fdatasync(wal->descriptor) == 0;
}

bool geo_wal_checkpoint(GeoWal *wal, uint64_t next_sequence)
{
    if (!wal || wal->descriptor < 0 || next_sequence != wal->next_sequence || wal->segment_id == UINT64_MAX ||
        fdatasync(wal->descriptor) != 0) {
        return false;
    }

    uint64_t previous_first_segment = wal->first_segment_id;
    uint64_t previous_last_segment = wal->segment_id;

    if (!wal_create_segment(wal, previous_last_segment + 1U, next_sequence)) {
        return false;
    }

    wal->first_segment_id = wal->segment_id;
    bool succeeded = true;

    for (uint64_t segment = previous_first_segment; segment <= previous_last_segment; ++segment) {
        char *path = wal_segment_path(wal->directory, segment);

        if (!path || unlink(path) != 0) {
            succeeded = false;
        }

        free(path);

        if (!succeeded) {
            break;
        }

        if (segment == UINT64_MAX) {
            break;
        }
    }

    return geo_io_sync_parent_directory(wal->active_path) && succeeded;
}
