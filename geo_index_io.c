#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "geo_index_io.h"
#include "geo_index_persistence.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <unistd.h>
#define GEO_ATOMIC_IO_SUPPORTED 1
#else
#define GEO_ATOMIC_IO_SUPPORTED 0
#endif

static char *geo_io_copy_parent_path(const char *path)
{
    size_t length = strlen(path);
    char *directory = malloc(length + 1U);

    if (!directory) {
        return NULL;
    }

    memcpy(directory, path, length + 1U);

    char *separator = strrchr(directory, '/');

    if (!separator) {
        directory[0] = '.';
        directory[1] = '\0';
    } else if (separator == directory) {
        separator[1] = '\0';
    } else {
        *separator = '\0';
    }

    return directory;
}

FILE *geo_io_create_atomic_file(const char *target_path, char **temporary_path)
{
#if GEO_ATOMIC_IO_SUPPORTED
    static const char suffix[] = ".tmp.XXXXXX";

    if (!target_path || !target_path[0] || !temporary_path) {
        return NULL;
    }

    size_t target_length = strlen(target_path);

    if (target_length > SIZE_MAX - sizeof(suffix)) {
        return NULL;
    }

    char *path = malloc(target_length + sizeof(suffix));

    if (!path) {
        return NULL;
    }

    memcpy(path, target_path, target_length);
    memcpy(path + target_length, suffix, sizeof(suffix));

    int descriptor = mkstemp(path);

    if (descriptor < 0) {
        free(path);

        return NULL;
    }

    FILE *file = fdopen(descriptor, "w+b");

    if (!file) {
        close(descriptor);
        unlink(path);
        free(path);

        return NULL;
    }

    *temporary_path = path;

    return file;
#else
    (void) target_path;
    (void) temporary_path;

    return NULL;
#endif
}

bool geo_io_sync_parent_directory(const char *path)
{
#if GEO_ATOMIC_IO_SUPPORTED
    if (!path || !path[0]) {
        return false;
    }

    char *directory = geo_io_copy_parent_path(path);

    if (!directory) {
        return false;
    }

    int descriptor = open(directory, O_RDONLY);
    bool succeeded = descriptor >= 0 && fsync(descriptor) == 0;

    if (descriptor >= 0 && close(descriptor) != 0) {
        succeeded = false;
    }

    free(directory);

    return succeeded;
#else
    (void) path;

    return false;
#endif
}

bool geo_io_publish_atomic_file(FILE *file, const char *temporary_path, const char *target_path)
{
#if GEO_ATOMIC_IO_SUPPORTED
    if (!file || !temporary_path || !target_path) {
        return false;
    }

    bool succeeded = fflush(file) == 0 && fsync(fileno(file)) == 0;

    if (fclose(file) != 0) {
        succeeded = false;
    }

    if (succeeded) {
        succeeded = rename(temporary_path, target_path) == 0;
    }

    if (succeeded) {
        succeeded = geo_io_sync_parent_directory(target_path);
    }

    if (!succeeded) {
        unlink(temporary_path);
    }

    return succeeded;
#else
    (void) file;
    (void) temporary_path;
    (void) target_path;

    return false;
#endif
}

void geo_io_discard_atomic_file(FILE *file, const char *temporary_path)
{
    if (file) {
        (void) fclose(file);
    }

#if GEO_ATOMIC_IO_SUPPORTED
    if (temporary_path) {
        (void) unlink(temporary_path);
    }
#else
    (void) temporary_path;
#endif
}

bool geo_io_write_u64_offsets(FILE *file, const size_t *offsets, size_t count, uint64_t *checksum)
{
    if (!file || (!offsets && count) || !checksum) {
        return false;
    }

    uint64_t encoded[256];
    uint64_t value_checksum = geo_persisted_checksum_initial();

    for (size_t position = 0; position < count;) {
        size_t chunk_count = count - position;

        if (chunk_count > sizeof(encoded) / sizeof(encoded[0])) {
            chunk_count = sizeof(encoded) / sizeof(encoded[0]);
        }

        for (size_t i = 0; i < chunk_count; ++i) {
            encoded[i] = (uint64_t) offsets[position + i];
        }

        size_t bytes = chunk_count * sizeof(*encoded);

        if (fwrite(encoded, 1, bytes, file) != bytes) {
            return false;
        }

        value_checksum = geo_persisted_checksum_update(value_checksum, encoded, bytes);
        position += chunk_count;
    }

    *checksum = value_checksum;

    return true;
}
