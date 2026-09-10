#ifndef GEO_SECONDARY_INDEX_H
#define GEO_SECONDARY_INDEX_H

#include "geobolt/geobolt.h"
#include "geo_db_format.h"
#include "geo_rocks_bridge.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct GeoSecondaryCatalog GeoSecondaryCatalog;

typedef struct {
    int64_t *entry_deltas;
    int64_t *histogram_deltas;
    size_t index_count;
    size_t histogram_bin_count;
} GeoSecondaryWriteState;

typedef struct {
    const GeoDatabaseIndexPredicate *predicate;
    GeoDatabaseIndexType type;
    char json_pointer[GEO_DATABASE_MAX_INDEX_POINTER_SIZE + 1U];
} GeoSecondaryPreparedPredicate;

GeoSecondaryCatalog *geo_secondary_catalog_open(GeoRocksDatabase *rocks,
                                                uint64_t next_secondary_index_id,
                                                GeoDatabaseStatus *status);
void geo_secondary_catalog_destroy(GeoSecondaryCatalog *catalog);

bool geo_secondary_catalog_has_indexes(const GeoSecondaryCatalog *catalog);
GeoDatabaseStatus geo_secondary_catalog_recover(GeoRocksDatabase *rocks, GeoSecondaryCatalog *catalog);
GeoDatabaseStatus geo_secondary_catalog_create(GeoRocksDatabase *rocks,
                                               GeoSecondaryCatalog *catalog,
                                               GeoDbCatalog *database_catalog,
                                               const char *name,
                                               const char *json_pointer,
                                               GeoDatabaseIndexType type);
GeoDatabaseStatus geo_secondary_catalog_drop(GeoRocksDatabase *rocks,
                                             GeoSecondaryCatalog *catalog,
                                             const char *name);
GeoDatabaseStatus geo_secondary_catalog_list(const GeoSecondaryCatalog *catalog,
                                             GeoDatabaseIndexInfo *indexes,
                                             size_t capacity,
                                             size_t *count);

GeoDatabaseStatus geo_secondary_write_begin(const GeoSecondaryCatalog *catalog, GeoSecondaryWriteState *state);
// Views must already be validated and remain immutable for the duration of the write.
GeoDatabaseStatus geo_secondary_write_object(const GeoSecondaryCatalog *catalog,
                                             GeoRocksBatch *batch,
                                             GeoSecondaryWriteState *state,
                                             uint64_t object_id,
                                             GeoDocView old_document,
                                             GeoDocView new_document);
GeoDatabaseStatus geo_secondary_write_statistics(const GeoSecondaryCatalog *catalog,
                                                 GeoRocksBatch *batch,
                                                 const GeoSecondaryWriteState *state);
bool geo_secondary_write_commit(GeoSecondaryCatalog *catalog, const GeoSecondaryWriteState *state);
void geo_secondary_write_destroy(GeoSecondaryWriteState *state);

GeoDatabaseStatus geo_secondary_query(const GeoSecondaryCatalog *catalog,
                                      GeoRocksDatabase *rocks,
                                      const GeoDatabaseIndexPredicate *predicate,
                                      GeoIdResult *result,
                                      GeoDatabaseIndexQueryStats *stats);
GeoDatabaseStatus geo_secondary_query_filtered(const GeoSecondaryCatalog *catalog,
                                               GeoRocksDatabase *rocks,
                                               const GeoDatabaseIndexPredicate *predicate,
                                               GeoIdResult *result,
                                               GeoDatabaseIndexQueryStats *stats,
                                               bool (*filter)(uint64_t object_id, void *context),
                                               void *filter_context);
GeoDatabaseStatus geo_secondary_estimate(const GeoSecondaryCatalog *catalog,
                                         const GeoDatabaseIndexPredicate *predicate,
                                         uint64_t *estimated_matches);
GeoDatabaseStatus geo_secondary_prepare_predicates(const GeoSecondaryCatalog *catalog,
                                                   const GeoDatabaseIndexPredicate *predicates,
                                                   size_t predicate_count,
                                                   GeoSecondaryPreparedPredicate *prepared);
GeoDatabaseStatus geo_secondary_document_matches(const GeoSecondaryPreparedPredicate *predicates,
                                                 size_t predicate_count,
                                                 GeoDocView document,
                                                 bool *matches);

#endif
