#ifndef GEO_INDEX_PERSISTENCE_H
#define GEO_INDEX_PERSISTENCE_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#define GEO_FILE_LEGACY_VERSION 1U
#define GEO_FILE_PREFIX_VERSION 2U
#define GEO_FILE_DENSITY_VERSION 3U
#define GEO_FILE_VERSION 5U
#define GEO_FILE_ENDIAN_MARKER UINT32_C(0x01020304)
#define GEO_PERSISTED_PREFIX_MIN_RECORDS UINT64_C(4096)
#define GEO_PERSISTED_PREFIX_MIN_BITS 8U
#define GEO_PERSISTED_PREFIX_MAX_BITS 16U
#define GEO_PERSISTED_PREFIX_TARGET_RECORDS UINT64_C(64)
#define GEO_PERSISTED_DENSITY_VERSION 1U
#define GEO_PERSISTED_DENSITY_MAX_LEVELS 6U

typedef struct {
    uint64_t range_min;
    uint64_t begin;
    uint64_t end;
} GeoPersistedDensityCell;

typedef struct {
    uint32_t version;
    uint32_t cell_size;
    uint32_t level_count;
    uint32_t reserved;
    uint64_t cell_count;
    uint64_t level_first[GEO_PERSISTED_DENSITY_MAX_LEVELS + 1U];
    uint8_t level_prefix_bits[GEO_PERSISTED_DENSITY_MAX_LEVELS];
    uint8_t padding[8];
} GeoPersistedDensityHeader;

_Static_assert(sizeof(GeoPersistedDensityCell) == 24, "GeoPersistedDensityCell layout is persisted");
_Static_assert(offsetof(GeoPersistedDensityCell, begin) == 8, "GeoPersistedDensityCell.begin offset is persisted");
_Static_assert(sizeof(GeoPersistedDensityHeader) == 96, "GeoPersistedDensityHeader layout is persisted");
_Static_assert(offsetof(GeoPersistedDensityHeader, level_first) == 24, "GeoPersistedDensityHeader.level_first offset is persisted");

typedef struct {
    char magic[8];
    uint32_t version;
    uint32_t record_size;
    uint32_t endian_marker;
    uint32_t prefix_bits;
    uint64_t count;
    uint64_t records_checksum;
    uint64_t prefix_checksum;
    uint64_t density_checksum;
    uint64_t density_bytes;
} GeoFileHeader;

_Static_assert(sizeof(GeoFileHeader) == 64, "GeoFileHeader layout is part of the mmap format");
_Static_assert(offsetof(GeoFileHeader, count) == 24, "GeoFileHeader.count offset is part of the mmap format");
_Static_assert(offsetof(GeoFileHeader, records_checksum) == 32,
               "GeoFileHeader.records_checksum offset is part of the mmap format");

static const char GEO_FILE_MAGIC[8] = { 'G', 'E', 'O', 'B', 'L', 'T', '1', '\0' };

static inline uint64_t geo_persisted_checksum_initial(void)
{
    return UINT64_C(1469598103934665603);
}

/*
 * Persisted sections and intermediate write chunks are multiples of eight bytes. The polynomial form remains identical across
 * streaming block boundaries. Four-word folding shortens the dependency chain without requiring architecture-specific CRC support.
 */
static inline uint64_t geo_persisted_checksum_update(uint64_t checksum, const void *data, size_t size)
{
    const uint8_t *bytes = data;
    const uint64_t prime = UINT64_C(0x9e3779b185ebca87);
    const uint64_t prime_squared = prime * prime;
    const uint64_t prime_cubed = prime_squared * prime;
    const uint64_t prime_fourth = prime_squared * prime_squared;

    while (size >= 4U * sizeof(uint64_t)) {
        uint64_t word0;
        uint64_t word1;
        uint64_t word2;
        uint64_t word3;

        memcpy(&word0, bytes, sizeof(word0));
        memcpy(&word1, bytes + sizeof(uint64_t), sizeof(word1));
        memcpy(&word2, bytes + 2U * sizeof(uint64_t), sizeof(word2));
        memcpy(&word3, bytes + 3U * sizeof(uint64_t), sizeof(word3));

        checksum = checksum * prime_fourth +
                   word0 * prime_cubed +
                   word1 * prime_squared +
                   word2 * prime +
                   word3;
        bytes += 4U * sizeof(uint64_t);
        size -= 4U * sizeof(uint64_t);
    }

    while (size >= sizeof(uint64_t)) {
        uint64_t word;

        memcpy(&word, bytes, sizeof(word));
        checksum = checksum * prime + word;
        bytes += sizeof(word);
        size -= sizeof(word);
    }

    for (size_t i = 0; i < size; ++i) {
        checksum ^= bytes[i];
        checksum *= UINT64_C(1099511628211);
    }

    return checksum;
}

/*
 * Combines an existing checksum with a suffix checksum evaluated from zero. Persisted record sections are always composed of complete
 * uint64_t words, so independent Morton partitions can checksum locally and combine without rereading the final output.
 */
static inline uint64_t geo_persisted_checksum_combine_aligned(uint64_t checksum,
                                                              uint64_t zero_seed_suffix_checksum,
                                                              size_t suffix_bytes)
{
    const uint64_t prime = UINT64_C(0x9e3779b185ebca87);
    uint64_t factor = prime;
    uint64_t power = UINT64_C(1);
    size_t suffix_words = suffix_bytes / sizeof(uint64_t);

    while (suffix_words) {
        if (suffix_words & 1U) {
            power *= factor;
        }

        factor *= factor;
        suffix_words >>= 1U;
    }

    return checksum * power + zero_seed_suffix_checksum;
}

static inline uint8_t geo_persisted_prefix_bits(uint64_t count)
{
    if (count < GEO_PERSISTED_PREFIX_MIN_RECORDS) {
        return 0;
    }

    uint8_t bits = GEO_PERSISTED_PREFIX_MIN_BITS;

    while (bits < GEO_PERSISTED_PREFIX_MAX_BITS &&
           (UINT64_C(1) << bits) * GEO_PERSISTED_PREFIX_TARGET_RECORDS < count) {
        bits++;
    }

    return bits;
}

#endif
