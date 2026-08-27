#ifndef GEO_INDEX_IO_H
#define GEO_INDEX_IO_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

// Creates a temporary file beside target_path. The caller owns both the FILE and returned path.
FILE *geo_io_create_atomic_file(const char *target_path, char **temporary_path);

// Flushes and closes file, atomically replaces target_path, and synchronizes its parent directory.
// The caller still owns temporary_path and must free it after this call.
bool geo_io_publish_atomic_file(FILE *file, const char *temporary_path, const char *target_path);

// Closes an unpublished temporary file and removes its directory entry.
void geo_io_discard_atomic_file(FILE *file, const char *temporary_path);

// Makes a prior create, rename, or unlink durable when supported by the host filesystem.
bool geo_io_sync_parent_directory(const char *path);

// Persists native offsets as fixed-width uint64_t values and returns the checksum of the encoded bytes.
bool geo_io_write_u64_offsets(FILE *file, const size_t *offsets, size_t count, uint64_t *checksum);

#endif
