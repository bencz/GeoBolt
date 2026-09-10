#ifndef GEO_ROCKS_BRIDGE_H
#define GEO_ROCKS_BRIDGE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GeoRocksDatabase GeoRocksDatabase;
typedef struct GeoRocksBatch GeoRocksBatch;
typedef struct GeoRocksSnapshot GeoRocksSnapshot;
typedef struct GeoRocksIterator GeoRocksIterator;
typedef struct GeoRocksMultiGet GeoRocksMultiGet;

typedef enum {
    GEO_ROCKS_CF_CATALOG = 0,
    GEO_ROCKS_CF_OBJECTS,
    GEO_ROCKS_CF_SPATIAL_DELTA,
    GEO_ROCKS_CF_IDEMPOTENCY,
    GEO_ROCKS_CF_SECONDARY_INDEX,
    GEO_ROCKS_CF_COUNT,
} GeoRocksColumnFamily;

typedef enum {
    GEO_ROCKS_OK = 0,
    GEO_ROCKS_NOT_FOUND,
    GEO_ROCKS_INVALID_ARGUMENT,
    GEO_ROCKS_IO_ERROR,
    GEO_ROCKS_CORRUPTION,
    GEO_ROCKS_BUSY,
    GEO_ROCKS_OUT_OF_MEMORY,
    GEO_ROCKS_INTERNAL_ERROR,
} GeoRocksStatusCode;

typedef struct {
    GeoRocksStatusCode code;
    char message[256];
} GeoRocksStatus;

typedef struct {
    size_t block_cache_bytes;
    size_t write_buffer_bytes;
    int background_jobs;
    bool create_if_missing;
} GeoRocksConfig;

typedef struct {
    void *data;
    size_t size;
} GeoRocksBuffer;

GeoRocksConfig geo_rocks_default_config(void);
GeoRocksDatabase *geo_rocks_open(const char *directory, const GeoRocksConfig *config, GeoRocksStatus *status);
void geo_rocks_close(GeoRocksDatabase *database);

GeoRocksBatch *geo_rocks_batch_create(GeoRocksDatabase *database, size_t reserved_bytes, GeoRocksStatus *status);
void geo_rocks_batch_destroy(GeoRocksBatch *batch);
bool geo_rocks_batch_put(GeoRocksBatch *batch,
                         GeoRocksColumnFamily column_family,
                         const void *key,
                         size_t key_size,
                         const void *value,
                         size_t value_size,
                         GeoRocksStatus *status);
bool geo_rocks_batch_delete(GeoRocksBatch *batch,
                            GeoRocksColumnFamily column_family,
                            const void *key,
                            size_t key_size,
                            GeoRocksStatus *status);
bool geo_rocks_batch_delete_range(GeoRocksBatch *batch,
                                  GeoRocksColumnFamily column_family,
                                  const void *begin_key,
                                  size_t begin_key_size,
                                  const void *end_key,
                                  size_t end_key_size,
                                  GeoRocksStatus *status);
bool geo_rocks_write(GeoRocksDatabase *database, GeoRocksBatch *batch, bool synchronize, GeoRocksStatus *status);

bool geo_rocks_get(GeoRocksDatabase *database,
                   GeoRocksColumnFamily column_family,
                   const GeoRocksSnapshot *snapshot,
                   const void *key,
                   size_t key_size,
                   GeoRocksBuffer *value,
                   GeoRocksStatus *status);
void geo_rocks_buffer_release(GeoRocksBuffer *buffer);

GeoRocksMultiGet *geo_rocks_multi_get_create(size_t capacity, GeoRocksStatus *status);
void geo_rocks_multi_get_destroy(GeoRocksMultiGet *multi_get);
void geo_rocks_multi_get_release(GeoRocksMultiGet *multi_get);
bool geo_rocks_multi_get_u64_be(GeoRocksDatabase *database,
                                GeoRocksColumnFamily column_family,
                                const GeoRocksSnapshot *snapshot,
                                const uint64_t *keys,
                                size_t key_count,
                                GeoRocksMultiGet *multi_get,
                                GeoRocksStatus *status);
/*
 * Result values are borrowed pins valid only until the next execute call, explicit release, or workspace destruction.
 * A MultiGet workspace is reusable but not safe for concurrent calls.
 */
GeoRocksStatusCode geo_rocks_multi_get_result(const GeoRocksMultiGet *multi_get,
                                              size_t index,
                                              const void **value,
                                              size_t *value_size);

GeoRocksSnapshot *geo_rocks_snapshot_create(GeoRocksDatabase *database, GeoRocksStatus *status);
void geo_rocks_snapshot_destroy(GeoRocksSnapshot *snapshot);
uint64_t geo_rocks_snapshot_sequence(const GeoRocksSnapshot *snapshot);

GeoRocksIterator *geo_rocks_iterator_create(GeoRocksDatabase *database,
                                            GeoRocksColumnFamily column_family,
                                            const GeoRocksSnapshot *snapshot,
                                            GeoRocksStatus *status);
void geo_rocks_iterator_destroy(GeoRocksIterator *iterator);
void geo_rocks_iterator_seek_first(GeoRocksIterator *iterator);
void geo_rocks_iterator_seek(GeoRocksIterator *iterator, const void *key, size_t key_size);
void geo_rocks_iterator_next(GeoRocksIterator *iterator);
bool geo_rocks_iterator_valid(const GeoRocksIterator *iterator);
const void *geo_rocks_iterator_key(const GeoRocksIterator *iterator, size_t *size);
const void *geo_rocks_iterator_value(const GeoRocksIterator *iterator, size_t *size);
bool geo_rocks_iterator_status(const GeoRocksIterator *iterator, GeoRocksStatus *status);

bool geo_rocks_flush(GeoRocksDatabase *database, bool wait, GeoRocksStatus *status);
bool geo_rocks_get_property_u64(GeoRocksDatabase *database,
                                GeoRocksColumnFamily column_family,
                                const char *property,
                                uint64_t *value,
                                GeoRocksStatus *status);

#ifdef __cplusplus
}
#endif

#endif
