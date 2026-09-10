#ifndef GEO_DB_FORMAT_H
#define GEO_DB_FORMAT_H

#include "geobolt/geobolt.h"
#include "geobolt/geodoc.h"
#include "geo_rocks_bridge.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define GEO_DB_FORMAT_VERSION 4U
#define GEO_DB_OBJECT_HEADER_SIZE 32U
#define GEO_DB_SPATIAL_DELTA_SIZE 32U
#define GEO_DB_IDEMPOTENCY_VALUE_SIZE 48U

typedef struct {
    uint64_t fingerprint_low;
    uint64_t fingerprint_high;
    uint64_t expires_at_seconds;
    uint64_t first_sequence;
    uint32_t operation_count;
} GeoDbIdempotencyRecord;

typedef struct {
    uint64_t next_sequence;
    uint64_t committed_operations;
    uint64_t spatial_applied_sequence;
    uint64_t next_server_id;
    uint64_t next_secondary_index_id;
} GeoDbCatalog;

typedef struct {
    uint64_t object_id;
    uint64_t sequence;
    uint64_t morton_code;
    uint32_t operation;
} GeoDbSpatialDelta;

void geo_db_store_u64_be(unsigned char output[8], uint64_t value);
uint64_t geo_db_load_u64_be(const unsigned char input[8]);

bool geo_db_catalog_load(GeoRocksDatabase *rocks,
                         bool create_if_missing,
                         GeoDbCatalog *catalog,
                         GeoDatabaseStatus *status);
bool geo_db_catalog_put(GeoRocksBatch *batch, const GeoDbCatalog *catalog, GeoDatabaseStatus *status);

// Internal serialization boundary: the caller owns an immutable, already validated GeoDoc view.
bool geo_db_object_put(GeoRocksBatch *batch,
                       uint64_t object_id,
                       uint64_t sequence,
                       uint64_t morton_code,
                       GeoDocView document,
                       GeoDatabaseStatus *status);
bool geo_db_object_delete(GeoRocksBatch *batch, uint64_t object_id, GeoDatabaseStatus *status);
bool geo_db_object_decode(const void *value,
                          size_t value_size,
                          uint64_t *sequence,
                          uint64_t *morton_code,
                          const void **document,
                          size_t *document_size);
bool geo_db_object_decode_view(const void *value,
                              size_t value_size,
                              uint64_t *sequence,
                              uint64_t *morton_code,
                              GeoDocView *document);

bool geo_db_spatial_delta_put(GeoRocksBatch *batch,
                              const GeoDbSpatialDelta *delta,
                              GeoDatabaseStatus *status);
bool geo_db_spatial_delta_decode(const void *key,
                                 size_t key_size,
                                 const void *value,
                                 size_t value_size,
                                 GeoDbSpatialDelta *delta);

bool geo_db_idempotency_get(GeoRocksDatabase *rocks,
                            const void *key,
                            size_t key_size,
                            GeoDbIdempotencyRecord *record,
                            bool *found,
                            GeoDatabaseStatus *status);
bool geo_db_idempotency_put(GeoRocksBatch *batch,
                            const void *key,
                            size_t key_size,
                            const GeoDbIdempotencyRecord *record,
                            GeoDatabaseStatus *status);

GeoDatabaseStatus geo_db_status_from_rocks(const GeoRocksStatus *status);

#endif
