#include "geo_query_planner.h"

#include "geo_database_internal.h"
#include "geo_index_private.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define GEO_QUERY_MINIMUM_HASH_CAPACITY 16U
#define GEO_QUERY_CANONICAL_READ_BATCH 256U
#define GEO_QUERY_SECONDARY_SCAN_COST 1U
#define GEO_QUERY_SECONDARY_LOOKUP_COST 24U
#define GEO_QUERY_SPATIAL_CANDIDATE_COST 12U

struct GeoDatabaseQueryWorkspace {
    GeoIdResult *predicate_results;
    uint64_t *membership;
    uint32_t *membership_epochs;
    GeoRocksMultiGet *canonical_reads;
    GeoSecondaryPreparedPredicate *prepared_predicates;
    GeoRadiusQueryPlan spatial_plan;
    GeoObjectCacheView canonical_views[GEO_QUERY_CANONICAL_READ_BATCH];
    uint64_t canonical_ids[GEO_QUERY_CANONICAL_READ_BATCH];
    uint64_t canonical_miss_ids[GEO_QUERY_CANONICAL_READ_BATCH];
    size_t canonical_positions[GEO_QUERY_CANONICAL_READ_BATCH];
    uint64_t predicate_estimates[GEO_DATABASE_MAX_QUERY_PREDICATES];
    size_t predicate_order[GEO_DATABASE_MAX_QUERY_PREDICATES];
    size_t predicate_capacity;
    size_t membership_capacity;
    uint32_t membership_epoch;
};

static uint64_t query_saturating_add(uint64_t first, uint64_t second)
{
    return second > UINT64_MAX - first ? UINT64_MAX : first + second;
}

static uint64_t query_saturating_multiply(uint64_t value, uint64_t multiplier)
{
    return value && multiplier > UINT64_MAX / value ? UINT64_MAX : value * multiplier;
}

static uint64_t query_hash_id(uint64_t value)
{
    value ^= value >> 30U;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27U;
    value *= UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

static bool query_membership_reserve(GeoDatabaseQueryWorkspace *workspace, size_t candidate_count)
{
    if (candidate_count > SIZE_MAX / 2U) {
        return false;
    }

    size_t required = candidate_count * 2U;
    size_t capacity = GEO_QUERY_MINIMUM_HASH_CAPACITY;

    while (capacity < required) {
        if (capacity > SIZE_MAX / 2U) {
            return false;
        }

        capacity *= 2U;
    }

    if (capacity > workspace->membership_capacity) {
        if (capacity > SIZE_MAX / sizeof(*workspace->membership) ||
            capacity > SIZE_MAX / sizeof(*workspace->membership_epochs)) {
            return false;
        }

        uint64_t *membership = malloc(capacity * sizeof(*membership));
        uint32_t *epochs = calloc(capacity, sizeof(*epochs));

        if (!membership || !epochs) {
            free(epochs);
            free(membership);
            return false;
        }

        free(workspace->membership_epochs);
        free(workspace->membership);
        workspace->membership = membership;
        workspace->membership_epochs = epochs;
        workspace->membership_capacity = capacity;
        workspace->membership_epoch = 0U;
    }

    workspace->membership_epoch++;

    if (!workspace->membership_epoch) {
        memset(workspace->membership_epochs,
               0,
               workspace->membership_capacity * sizeof(*workspace->membership_epochs));
        workspace->membership_epoch = 1U;
    }

    return true;
}

static void query_membership_insert(GeoDatabaseQueryWorkspace *workspace, uint64_t object_id)
{
    size_t mask = workspace->membership_capacity - 1U;
    size_t slot = (size_t) query_hash_id(object_id) & mask;

    while (workspace->membership_epochs[slot] == workspace->membership_epoch &&
           workspace->membership[slot] != object_id) {
        slot = (slot + 1U) & mask;
    }

    workspace->membership[slot] = object_id;
    workspace->membership_epochs[slot] = workspace->membership_epoch;
}

static bool query_membership_contains(const GeoDatabaseQueryWorkspace *workspace, uint64_t object_id)
{
    size_t mask = workspace->membership_capacity - 1U;
    size_t slot = (size_t) query_hash_id(object_id) & mask;

    while (workspace->membership_epochs[slot] == workspace->membership_epoch) {
        if (workspace->membership[slot] == object_id) {
            return true;
        }

        slot = (slot + 1U) & mask;
    }

    return false;
}

static bool query_membership_build(GeoDatabaseQueryWorkspace *workspace, const GeoIdResult *candidates)
{
    if (!query_membership_reserve(workspace, candidates->count)) {
        return false;
    }

    for (size_t index = 0U; index < candidates->count; ++index) {
        if (!candidates->ids[index]) {
            return false;
        }

        query_membership_insert(workspace, candidates->ids[index]);
    }

    return true;
}

static bool query_membership_filter_id(uint64_t object_id, void *context)
{
    return query_membership_contains(context, object_id);
}

static GeoDatabaseStatus query_order_predicates(const GeoDatabase *database,
                                                const GeoDatabaseRadiusQuery *query,
                                                GeoDatabaseQueryWorkspace *workspace)
{
    for (size_t index = 0U; index < query->predicate_count; ++index) {
        uint64_t estimate;
        GeoDatabaseStatus status = geo_secondary_estimate(database->secondary_indexes,
                                                          query->predicates + index,
                                                          &estimate);

        if (status != GEO_DATABASE_OK) {
            return status;
        }

        size_t insertion = index;

        while (insertion && estimate < workspace->predicate_estimates[insertion - 1U]) {
            workspace->predicate_estimates[insertion] = workspace->predicate_estimates[insertion - 1U];
            workspace->predicate_order[insertion] = workspace->predicate_order[insertion - 1U];
            insertion--;
        }

        workspace->predicate_estimates[insertion] = estimate;
        workspace->predicate_order[insertion] = index;
    }

    return GEO_DATABASE_OK;
}

static GeoDatabaseStatus query_secondary_intersection(const GeoDatabase *database,
                                                      const GeoDatabaseRadiusQuery *query,
                                                      GeoDatabaseQueryWorkspace *workspace,
                                                      GeoIdResult **candidates,
                                                      uint64_t *entries_scanned)
{
    *entries_scanned = 0U;
    GeoDatabaseStatus status = GEO_DATABASE_OK;

    GeoIdResult *driver = NULL;

    for (size_t rank = 0U; rank < query->predicate_count; ++rank) {
        GeoDatabaseIndexQueryStats index_stats;
        GeoIdResult *partial = workspace->predicate_results + rank;
        size_t predicate_index = workspace->predicate_order[rank];

        if (!driver) {
            status = geo_secondary_query(database->secondary_indexes,
                                         database->rocks,
                                         query->predicates + predicate_index,
                                         partial,
                                         &index_stats);
        } else {
            status = geo_secondary_query_filtered(database->secondary_indexes,
                                                  database->rocks,
                                                  query->predicates + predicate_index,
                                                  partial,
                                                  &index_stats,
                                                  query_membership_filter_id,
                                                  workspace);
        }

        if (status != GEO_DATABASE_OK) {
            return status;
        }

        *entries_scanned = query_saturating_add(*entries_scanned, index_stats.scanned_entries);
        driver = partial;

        if (!driver->count) {
            *candidates = driver;
            return GEO_DATABASE_OK;
        }

        if (rank + 1U < query->predicate_count && !query_membership_build(workspace, driver)) {
            return GEO_DATABASE_OUT_OF_MEMORY;
        }
    }

    *candidates = driver;
    return GEO_DATABASE_OK;
}

static bool query_longitude_matches(double longitude, double minimum, double maximum)
{
    return minimum <= maximum ? longitude >= minimum && longitude <= maximum
                              : longitude >= minimum || longitude <= maximum;
}

static bool query_point_matches(const GeoDatabaseRadiusQuery *query,
                                double minimum_latitude,
                                double maximum_latitude,
                                double minimum_longitude,
                                double maximum_longitude,
                                uint64_t morton_code)
{
    GeoPoint point = geo_decode(morton_code);

    if (point.lat < minimum_latitude || point.lat > maximum_latitude ||
        !query_longitude_matches(point.lng, minimum_longitude, maximum_longitude)) {
        return false;
    }

    return geo_haversine_km(query->latitude, query->longitude, point.lat, point.lng) <= query->radius_km;
}

static GeoDatabaseStatus query_load_canonical_batch(GeoDatabase *database,
                                                    const uint64_t *object_ids,
                                                    const GeoRecord *expected_records,
                                                    size_t object_count,
                                                    GeoDatabaseQueryWorkspace *workspace,
                                                    GeoDatabaseQueryStats *stats)
{
    size_t miss_count = 0U;

    for (size_t index = 0U; index < object_count; ++index) {
        GeoObjectCacheView cached;

        if (geo_object_cache_find(database->object_cache, object_ids[index], &cached)) {
            if (expected_records && cached.morton_code != expected_records[index].z) {
                return GEO_DATABASE_CORRUPTION;
            }

            workspace->canonical_views[index] = cached;
            continue;
        }

        workspace->canonical_miss_ids[miss_count] = object_ids[index];
        workspace->canonical_positions[miss_count] = index;
        miss_count++;
    }

    if (stats) {
        stats->object_lookups = query_saturating_add(stats->object_lookups, miss_count);
    }
    if (!miss_count) {
        return GEO_DATABASE_OK;
    }

    GeoRocksStatus rocks_status;

    if (!geo_rocks_multi_get_u64_be(database->rocks,
                                    GEO_ROCKS_CF_OBJECTS,
                                    NULL,
                                    workspace->canonical_miss_ids,
                                    miss_count,
                                    workspace->canonical_reads,
                                    &rocks_status)) {
        return geo_db_status_from_rocks(&rocks_status);
    }

    for (size_t miss = 0U; miss < miss_count; ++miss) {
        const void *encoded;
        size_t encoded_size;
        GeoRocksStatusCode result_status = geo_rocks_multi_get_result(workspace->canonical_reads,
                                                                      miss,
                                                                      &encoded,
                                                                      &encoded_size);

        if (result_status == GEO_ROCKS_NOT_FOUND) {
            return GEO_DATABASE_CORRUPTION;
        }
        if (result_status != GEO_ROCKS_OK) {
            GeoRocksStatus result_error = { .code = result_status };
            return geo_db_status_from_rocks(&result_error);
        }

        size_t position = workspace->canonical_positions[miss];
        GeoObjectCacheView loaded;

        if (!geo_db_object_decode_view(encoded,
                                       encoded_size,
                                       &loaded.sequence,
                                       &loaded.morton_code,
                                       &loaded.document) ||
            (expected_records && loaded.morton_code != expected_records[position].z)) {
            return GEO_DATABASE_CORRUPTION;
        }

        workspace->canonical_views[position] = loaded;

        if (database->object_cache) {
            GeoObjectCacheEntry *entry = geo_object_cache_entry_create(object_ids[position],
                                                                        loaded.sequence,
                                                                        loaded.morton_code,
                                                                        loaded.document);
            geo_object_cache_fill(database->object_cache, entry);
        }
    }

    return GEO_DATABASE_OK;
}

static GeoDatabaseStatus query_from_secondary(GeoDatabase *database,
                                              const GeoDatabaseRadiusQuery *query,
                                              const GeoIdResult *candidates,
                                              GeoDatabaseQueryWorkspace *workspace,
                                              GeoSearchResult *result,
                                              GeoDatabaseQueryStats *stats)
{
    double minimum_latitude;
    double maximum_latitude;
    double minimum_longitude;
    double maximum_longitude;

    geo_bounding_box(query->latitude,
                     query->longitude,
                     query->radius_km,
                     &minimum_latitude,
                     &maximum_latitude,
                     &minimum_longitude,
                     &maximum_longitude);

    for (size_t first = 0U; first < candidates->count;) {
        size_t remaining = candidates->count - first;
        size_t batch_size = remaining < GEO_QUERY_CANONICAL_READ_BATCH ? remaining : GEO_QUERY_CANONICAL_READ_BATCH;
        GeoDatabaseStatus status = query_load_canonical_batch(database,
                                                              candidates->ids + first,
                                                              NULL,
                                                              batch_size,
                                                              workspace,
                                                              stats);

        if (status != GEO_DATABASE_OK) {
            return status;
        }
        if (stats) {
            stats->spatial.records_scanned = query_saturating_add(stats->spatial.records_scanned, batch_size);
        }

        for (size_t index = 0U; index < batch_size; ++index) {
            uint64_t morton_code = workspace->canonical_views[index].morton_code;

            if (!query_point_matches(query,
                                     minimum_latitude,
                                     maximum_latitude,
                                     minimum_longitude,
                                     maximum_longitude,
                                     morton_code)) {
                continue;
            }

            GeoRecord record = { .id = candidates->ids[first + index], .z = morton_code };

            if (!geo_result_add(result, &record)) {
                return GEO_DATABASE_OUT_OF_MEMORY;
            }

            if (stats) {
                stats->spatial.records_matched++;
            }
        }

        first += batch_size;
    }

    return GEO_DATABASE_OK;
}

static GeoDatabaseStatus query_from_spatial_documents(GeoDatabase *database,
                                                      const GeoDatabaseRadiusQuery *query,
                                                      GeoDatabaseQueryWorkspace *workspace,
                                                      GeoSearchResult *result,
                                                      GeoDatabaseQueryStats *stats)
{
    GeoDatabaseStatus status = geo_secondary_prepare_predicates(database->secondary_indexes,
                                                                query->predicates,
                                                                query->predicate_count,
                                                                workspace->prepared_predicates);

    if (status != GEO_DATABASE_OK) {
        return status;
    }

    GeoSearchStats *spatial_stats = stats ? &stats->spatial : NULL;
    bool searched = geo_segment_set_search_radius_plan_reuse(database->segments,
                                                             &workspace->spatial_plan,
                                                             result,
                                                             spatial_stats) &&
                    geo_spatial_memtable_search_radius(database->spatial_memtable,
                                                       query->latitude,
                                                       query->longitude,
                                                       query->radius_km,
                                                       result,
                                                       NULL,
                                                       NULL,
                                                       NULL,
                                                       spatial_stats);

    if (!searched) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    size_t matched_count = 0U;

    for (size_t first = 0U; first < result->count;) {
        size_t remaining = result->count - first;
        size_t batch_size = remaining < GEO_QUERY_CANONICAL_READ_BATCH ? remaining : GEO_QUERY_CANONICAL_READ_BATCH;

        for (size_t index = 0U; index < batch_size; ++index) {
            workspace->canonical_ids[index] = result->results[first + index].id;
        }

        status = query_load_canonical_batch(database,
                                            workspace->canonical_ids,
                                            result->results + first,
                                            batch_size,
                                            workspace,
                                            stats);

        if (status != GEO_DATABASE_OK) {
            return status;
        }

        for (size_t index = 0U; index < batch_size; ++index) {
            bool matches;

            status = geo_secondary_document_matches(workspace->prepared_predicates,
                                                    query->predicate_count,
                                                    workspace->canonical_views[index].document,
                                                    &matches);

            if (status != GEO_DATABASE_OK) {
                return status;
            }
            if (matches) {
                result->results[matched_count++] = result->results[first + index];
            }
        }

        first += batch_size;
    }

    result->count = matched_count;

    if (stats) {
        stats->metadata_candidates = matched_count;
        stats->spatial.records_matched = matched_count;
    }

    return GEO_DATABASE_OK;
}

GeoDatabaseQueryWorkspace *geo_database_query_workspace_create(size_t predicate_capacity)
{
    if (!predicate_capacity || predicate_capacity > GEO_DATABASE_MAX_QUERY_PREDICATES) {
        return NULL;
    }

    GeoDatabaseQueryWorkspace *workspace = calloc(1U, sizeof(*workspace));

    if (!workspace) {
        return NULL;
    }

    workspace->predicate_results = calloc(predicate_capacity, sizeof(*workspace->predicate_results));
    workspace->prepared_predicates = calloc(predicate_capacity, sizeof(*workspace->prepared_predicates));
    GeoRocksStatus rocks_status;

    workspace->canonical_reads = geo_rocks_multi_get_create(GEO_QUERY_CANONICAL_READ_BATCH, &rocks_status);

    if (!workspace->predicate_results || !workspace->prepared_predicates || !workspace->canonical_reads) {
        geo_rocks_multi_get_destroy(workspace->canonical_reads);
        free(workspace->prepared_predicates);
        free(workspace->predicate_results);
        free(workspace);
        return NULL;
    }

    workspace->predicate_capacity = predicate_capacity;
    return workspace;
}

void geo_database_query_workspace_destroy(GeoDatabaseQueryWorkspace *workspace)
{
    if (!workspace) {
        return;
    }

    for (size_t index = 0U; index < workspace->predicate_capacity; ++index) {
        free(workspace->predicate_results[index].ids);
    }

    geo_rocks_multi_get_destroy(workspace->canonical_reads);
    free(workspace->membership);
    free(workspace->membership_epochs);
    free(workspace->prepared_predicates);
    free(workspace->predicate_results);
    free(workspace);
}

GeoDatabaseStatus geo_database_query_radius_reuse(const GeoDatabase *database,
                                                  const GeoDatabaseRadiusQuery *query,
                                                  GeoDatabaseQueryWorkspace *workspace,
                                                  GeoSearchResult *result,
                                                  GeoDatabaseQueryStats *stats)
{
    if (!database || !query || !workspace || !result || !query->predicates || !query->predicate_count ||
        query->predicate_count > workspace->predicate_capacity ||
        query->predicate_count > GEO_DATABASE_MAX_QUERY_PREDICATES ||
        !geo_is_valid_point(query->latitude, query->longitude) || !isfinite(query->radius_km) || query->radius_km < 0.0) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoDatabase *mutable_database = (GeoDatabase *) database;

    if (!geo_visibility_gate_read_lock(&mutable_database->visibility_gate)) {
        return GEO_DATABASE_FAILED_STATE;
    }

    double start = stats ? geo_get_time_ms() : 0.0;
    GeoDatabaseStatus status = atomic_load_explicit(&database->failed, memory_order_acquire)
                                   ? GEO_DATABASE_FAILED_STATE
                                   : GEO_DATABASE_OK;
    GeoIdResult *candidates = NULL;
    uint64_t secondary_entries_scanned = 0U;
    uint64_t estimated_spatial_candidates = 0U;

    geo_result_clear(result);

    if (stats) {
        *stats = (GeoDatabaseQueryStats) { .predicate_count = (uint32_t) query->predicate_count };
    }

    if (status == GEO_DATABASE_OK) {
        status = query_order_predicates(database, query, workspace);
    }

    if (status == GEO_DATABASE_OK &&
        !geo_segment_set_prepare_radius_query(database->segments,
                                              query->latitude,
                                              query->longitude,
                                              query->radius_km,
                                              sizeof(GeoRecord),
                                              &workspace->spatial_plan)) {
        status = GEO_DATABASE_INVALID_ARGUMENT;
    }

    if (status == GEO_DATABASE_OK) {
        estimated_spatial_candidates = workspace->spatial_plan.candidate_records;
    }

    if (status == GEO_DATABASE_OK) {
        estimated_spatial_candidates = query_saturating_add(
            estimated_spatial_candidates,
            geo_spatial_memtable_active_operations(database->spatial_memtable));
        estimated_spatial_candidates = query_saturating_add(
            estimated_spatial_candidates,
            geo_spatial_memtable_frozen_operations(database->spatial_memtable));

        uint64_t estimated_secondary_scans = 0U;

        for (size_t index = 0U; index < query->predicate_count; ++index) {
            estimated_secondary_scans = query_saturating_add(estimated_secondary_scans,
                                                             workspace->predicate_estimates[index]);
        }

        uint64_t estimated_metadata_candidates = workspace->predicate_estimates[0];
        uint64_t secondary_cost = query_saturating_add(
            query_saturating_multiply(estimated_secondary_scans, GEO_QUERY_SECONDARY_SCAN_COST),
            query_saturating_multiply(estimated_metadata_candidates, GEO_QUERY_SECONDARY_LOOKUP_COST));
        uint64_t spatial_cost = query_saturating_multiply(estimated_spatial_candidates,
                                                          GEO_QUERY_SPATIAL_CANDIDATE_COST);
        bool secondary_driven = secondary_cost <= spatial_cost;

        if (secondary_driven) {
            status = query_secondary_intersection(database,
                                                  query,
                                                  workspace,
                                                  &candidates,
                                                  &secondary_entries_scanned);
        }

        if (stats) {
            stats->secondary_entries_scanned = secondary_entries_scanned;
            stats->metadata_candidates = secondary_driven && candidates ? candidates->count : 0U;
            stats->estimated_spatial_candidates = estimated_spatial_candidates;
            stats->plan = secondary_driven ? GEO_DATABASE_QUERY_PLAN_SECONDARY : GEO_DATABASE_QUERY_PLAN_SPATIAL;
        }

        if (status == GEO_DATABASE_OK && secondary_driven && candidates && candidates->count) {
            status = query_from_secondary(mutable_database, query, candidates, workspace, result, stats);
        }
        if (status == GEO_DATABASE_OK && !secondary_driven) {
            status = query_from_spatial_documents(mutable_database, query, workspace, result, stats);
        }
    }

    if (stats) {
        stats->spatial.search_time_ms = geo_get_time_ms() - start;
    }

    /* PinnableSlice values must never outlive the RocksDB instance that supplied them. */
    geo_rocks_multi_get_release(workspace->canonical_reads);

    if (!geo_visibility_gate_read_unlock(&mutable_database->visibility_gate)) {
        status = GEO_DATABASE_FAILED_STATE;
    }

    if (status != GEO_DATABASE_OK) {
        geo_result_clear(result);
    }

    return status;
}
