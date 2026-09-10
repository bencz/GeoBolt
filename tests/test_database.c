#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "geobolt/geobolt.h"
#include "geobolt/geodoc.h"
#include "geo_db_format.h"
#include "geo_rocks_bridge.h"
#include "test_support.h"

#include <errno.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define TEST_ASSERT(condition, message)                                                                                   \
    do {                                                                                                                  \
        if (!(condition)) {                                                                                               \
            fprintf(stderr, "FAIL: %s (%s:%d)\n", message, __FILE__, __LINE__);                                         \
            succeeded = false;                                                                                            \
            goto cleanup;                                                                                                 \
        }                                                                                                                 \
    } while (0)

static GeoDatabaseConfig test_database_config(void)
{
    GeoDatabaseConfig config = geo_database_default_config();

    config.block_cache_bytes = 8U * 1024U * 1024U;
    config.object_cache_bytes = 4U * 1024U * 1024U;
    config.write_buffer_bytes = 4U * 1024U * 1024U;
    config.max_active_segments = 4U;
    config.background_jobs = 2;
    return config;
}

static bool result_contains(const GeoSearchResult *result, uint64_t object_id, uint64_t morton_code)
{
    for (size_t index = 0; index < result->count; ++index) {
        if (result->results[index].id == object_id && result->results[index].z == morton_code) {
            return true;
        }
    }

    return false;
}

static bool id_result_contains(const GeoIdResult *result, uint64_t object_id)
{
    for (size_t index = 0U; index < result->count; ++index) {
        if (result->ids[index] == object_id) {
            return true;
        }
    }

    return false;
}

static bool build_secondary_document(bool available,
                                     int64_t score,
                                     uint64_t rank,
                                     double price,
                                     int64_t updated_at,
                                     const char *label,
                                     const unsigned char token[2],
                                     GeoDocBuffer *document)
{
    GeoDocBuilder *builder = geo_doc_builder_create();
    bool succeeded = builder &&
                     geo_doc_builder_add_bool(builder, 0U, "available", available, NULL) == GEO_DOC_OK &&
                     geo_doc_builder_add_int64(builder, 0U, "score", score, NULL) == GEO_DOC_OK &&
                     geo_doc_builder_add_uint64(builder, 0U, "rank", rank, NULL) == GEO_DOC_OK &&
                     geo_doc_builder_add_double(builder, 0U, "price", price, NULL) == GEO_DOC_OK &&
                     geo_doc_builder_add_int64(builder, 0U, "updated_at", updated_at, NULL) == GEO_DOC_OK &&
                     geo_doc_builder_add_string(builder, 0U, "label", label, strlen(label), NULL) == GEO_DOC_OK &&
                     geo_doc_builder_add_bytes(builder, 0U, "token", token, 2U, NULL) == GEO_DOC_OK &&
                     geo_doc_builder_finish(builder, document) == GEO_DOC_OK;

    geo_doc_builder_destroy(builder);
    return succeeded;
}

static bool remove_derived_manifest(const char *directory)
{
    static const char *const names[] = {
        "spatial-manifest.gbm",
        "spatial-manifest.gbm.mutations",
        "spatial-manifest.gbm.mutations.checkpoint",
    };
    char path[512];

    for (size_t index = 0; index < sizeof(names) / sizeof(names[0]); ++index) {
        int length = snprintf(path, sizeof(path), "%s/%s", directory, names[index]);

        if (length <= 0 || (size_t) length >= sizeof(path) || (unlink(path) != 0 && errno != ENOENT)) {
            return false;
        }
    }

    return true;
}

static bool rewrite_spatial_applied_sequence(const char *directory, uint64_t spatial_applied_sequence)
{
    char rocks_path[512];
    int length = snprintf(rocks_path, sizeof(rocks_path), "%s/rocksdb", directory);

    if (length <= 0 || (size_t) length >= sizeof(rocks_path)) {
        return false;
    }

    GeoRocksConfig config = geo_rocks_default_config();
    GeoRocksStatus rocks_status;
    GeoRocksDatabase *rocks = geo_rocks_open(rocks_path, &config, &rocks_status);

    if (!rocks) {
        return false;
    }

    GeoDbCatalog catalog;
    GeoDatabaseStatus database_status = GEO_DATABASE_OK;
    bool succeeded = geo_db_catalog_load(rocks, false, &catalog, &database_status) &&
                     spatial_applied_sequence <= catalog.next_sequence;
    GeoRocksBatch *batch = NULL;

    if (succeeded) {
        catalog.spatial_applied_sequence = spatial_applied_sequence;
        batch = geo_rocks_batch_create(rocks, 128U, &rocks_status);
        succeeded = batch && geo_db_catalog_put(batch, &catalog, &database_status) &&
                    geo_rocks_write(rocks, batch, true, &rocks_status);
    }

    geo_rocks_batch_destroy(batch);
    geo_rocks_close(rocks);
    return succeeded;
}

static bool count_spatial_deltas(const char *directory, size_t *count)
{
    char rocks_path[512];
    int length = snprintf(rocks_path, sizeof(rocks_path), "%s/rocksdb", directory);

    if (!count || length <= 0 || (size_t) length >= sizeof(rocks_path)) {
        return false;
    }

    GeoRocksConfig config = geo_rocks_default_config();
    config.create_if_missing = false;
    GeoRocksStatus status;
    GeoRocksDatabase *rocks = geo_rocks_open(rocks_path, &config, &status);

    if (!rocks) {
        return false;
    }

    GeoRocksIterator *iterator = geo_rocks_iterator_create(rocks, GEO_ROCKS_CF_SPATIAL_DELTA, NULL, &status);
    size_t delta_count = 0U;

    if (iterator) {
        geo_rocks_iterator_seek_first(iterator);

        while (geo_rocks_iterator_valid(iterator)) {
            delta_count++;
            geo_rocks_iterator_next(iterator);
        }
    }

    bool succeeded = iterator && geo_rocks_iterator_status(iterator, &status);
    geo_rocks_iterator_destroy(iterator);
    geo_rocks_close(rocks);

    if (succeeded) {
        *count = delta_count;
    }

    return succeeded;
}

static bool test_database_spatial_watermark_crash_recovery(void)
{
    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-database-watermark-XXXXXX";
    char *directory = mkdtemp(directory_template);
    char manifest_path[512];
    GeoDatabase *database = NULL;
    GeoSegmentSet *segments = NULL;
    GeoSearchResult *result = NULL;

    TEST_ASSERT(directory != NULL, "watermark recovery directory must be created");
    int length = snprintf(manifest_path, sizeof(manifest_path), "%s/spatial-manifest.gbm", directory);
    TEST_ASSERT(length > 0 && (size_t) length < sizeof(manifest_path), "derived manifest path must fit");

    GeoDatabaseConfig config = test_database_config();
    GeoDatabaseStatus status;
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "watermark recovery database must open");

    GeoDatabaseMutation mutations[] = {
        { .object_id = 11U, .morton_code = geo_encode(-23.55, -46.63), .operation = GEO_DATABASE_UPSERT },
        { .object_id = 22U, .morton_code = geo_encode(-22.90, -43.17), .operation = GEO_DATABASE_UPSERT },
    };
    TEST_ASSERT(geo_database_write(database, mutations, 2U) == GEO_DATABASE_OK,
                "watermark recovery seed batch must commit");

    mutations[0].morton_code = geo_encode(51.50, -0.12);
    mutations[1] = (GeoDatabaseMutation) { .object_id = 22U, .operation = GEO_DATABASE_DELETE };
    TEST_ASSERT(geo_database_write(database, mutations, 2U) == GEO_DATABASE_OK,
                "watermark recovery mixed update/delete must commit");
    geo_database_close(database);
    database = NULL;

    TEST_ASSERT(rewrite_spatial_applied_sequence(directory, 1U),
                "test must reproduce the crash window before the catalog watermark update");

    config.create_if_missing = false;
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK,
                "manifest-ahead recovery must reconcile the catalog without replaying published deltas");

    GeoDatabaseStats stats;
    TEST_ASSERT(geo_database_get_stats(database, &stats) && stats.next_sequence == 5U && stats.recovered_operations == 0U,
                "manifest-ahead recovery must retain the exact canonical sequence without duplicate replay");

    result = geo_result_create(4U);
    TEST_ASSERT(result != NULL && geo_database_search_radius_reuse(database, 0.0, 0.0, 25000.0, result, NULL),
                "manifest-ahead recovery query must succeed");
    TEST_ASSERT(result->count == 1U && result_contains(result, 11U, mutations[0].morton_code),
                "manifest-ahead recovery must preserve the final update/delete visibility");
    geo_database_close(database);
    database = NULL;

    TEST_ASSERT(remove_derived_manifest(directory), "derived manifest must be replaceable with an older valid generation");
    segments = geo_segment_set_create(manifest_path);
    TEST_ASSERT(segments != NULL, "older valid derived manifest must be created");
    geo_segment_set_destroy(segments);
    segments = NULL;

    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK,
                "catalog-ahead recovery must rebuild the derived index from canonical objects");
    TEST_ASSERT(geo_database_search_radius_reuse(database, 0.0, 0.0, 25000.0, result, NULL) &&
                    result->count == 1U && result_contains(result, 11U, mutations[0].morton_code),
                "catalog-ahead rebuild must reproduce exact canonical visibility");
    geo_database_close(database);
    database = NULL;

    segments = geo_segment_set_open(manifest_path);
    TEST_ASSERT(segments != NULL && geo_segment_set_durable_watermark(segments) == 5U,
                "rebuild must publish the canonical next sequence into the manifest before updating the catalog");

cleanup:
    geo_result_destroy(result);
    geo_segment_set_destroy(segments);
    geo_database_close(database);

    if (directory && !geobolt_test_remove_tree(directory)) {
        fprintf(stderr, "FAIL: watermark recovery database tree could not be removed\n");
        succeeded = false;
    }

    return succeeded;
}

static bool test_database_objects_metadata_reopen_and_rebuild(void)
{
    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-database-test-XXXXXX";
    char *directory = mkdtemp(directory_template);
    GeoDatabase *database = NULL;
    GeoSearchResult *result = NULL;
    GeoDocBuilder *builder = NULL;
    GeoDocBuffer document = { 0 };
    GeoDatabaseObject object = { 0 };

    TEST_ASSERT(directory != NULL, "temporary database directory must be created");

    builder = geo_doc_builder_create();
    TEST_ASSERT(builder != NULL, "metadata builder must be created");
    TEST_ASSERT(geo_doc_builder_add_string(builder, 0U, "kind", "vehicle", 7U, NULL) == GEO_DOC_OK,
                "generic object kind must be encoded");
    TEST_ASSERT(geo_doc_builder_add_string(builder, 0U, "driver_name", "Ana", 3U, NULL) == GEO_DOC_OK,
                "metadata string must be encoded");
    TEST_ASSERT(geo_doc_builder_add_uint64(builder, 0U, "driver_id", 998877U, NULL) == GEO_DOC_OK,
                "metadata integer must be encoded");
    TEST_ASSERT(geo_doc_builder_finish(builder, &document) == GEO_DOC_OK, "metadata document must serialize");

    GeoDatabaseConfig config = test_database_config();
    GeoDatabaseStatus status;
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "new RocksDB-backed database must open");

    GeoDatabaseObjectMutation initial[] = {
        {
            .object_id = 101U,
            .morton_code = geo_encode(-23.5505, -46.6333),
            .document = document.data,
            .document_size = document.size,
            .operation = GEO_DATABASE_UPSERT,
        },
        { .object_id = 202U, .morton_code = geo_encode(-22.9068, -43.1729), .operation = GEO_DATABASE_UPSERT },
        { .object_id = 303U, .morton_code = geo_encode(40.7128, -74.0060), .operation = GEO_DATABASE_UPSERT },
    };

    TEST_ASSERT(geo_database_write_objects(database, initial, 3U) == GEO_DATABASE_OK,
                "objects, metadata, spatial deltas, and catalog must commit atomically");
    TEST_ASSERT(geo_database_get(database, 101U, &object) == GEO_DATABASE_OK, "committed object must be readable by ID");
    TEST_ASSERT(object.object_id == 101U && object.sequence == 1U && object.morton_code == initial[0].morton_code,
                "object identity, sequence, and Morton coordinate must round-trip");
    TEST_ASSERT(object.document_size == document.size && memcmp(object.document, document.data, document.size) == 0,
                "GeoDoc metadata must round-trip byte-for-byte");
    geo_database_object_release(&object);

    GeoDatabaseMutation update_and_delete[] = {
        { .object_id = 101U, .morton_code = geo_encode(51.5074, -0.1278), .operation = GEO_DATABASE_UPSERT },
        { .object_id = 202U, .operation = GEO_DATABASE_DELETE },
    };

    TEST_ASSERT(geo_database_write(database, update_and_delete, 2U) == GEO_DATABASE_OK,
                "mixed update and delete must commit");
    TEST_ASSERT(geo_database_get(database, 202U, &object) == GEO_DATABASE_NOT_FOUND,
                "deleted object must not remain in canonical storage");

    result = geo_result_create(8U);
    TEST_ASSERT(result != NULL, "reusable result must allocate");
    TEST_ASSERT(geo_database_search_radius_reuse(database, 0.0, 0.0, 25000.0, result, NULL),
                "global spatial query must succeed");
    TEST_ASSERT(result->count == 2U && result_contains(result, 101U, update_and_delete[0].morton_code),
                "spatial visibility must implement last-write-wins and deletion");

    GeoDatabaseMutation duplicate_ids[] = {
        { .object_id = 404U, .morton_code = geo_encode(1.0, 1.0), .operation = GEO_DATABASE_UPSERT },
        { .object_id = 404U, .operation = GEO_DATABASE_DELETE },
    };
    TEST_ASSERT(geo_database_write(database, duplicate_ids, 2U) == GEO_DATABASE_INVALID_ARGUMENT,
                "duplicate IDs inside one atomic batch must be rejected");

    uint64_t generated_id_a = 0U;
    uint64_t generated_id_b = 0U;
    TEST_ASSERT(geo_database_insert_generated_document(database,
                                                       35.6762,
                                                       139.6503,
                                                       document.data,
                                                       document.size,
                                                       &generated_id_a) == GEO_DATABASE_OK,
                "database-generated object ID must commit with metadata");
    TEST_ASSERT(geo_database_insert_generated_document(database,
                                                       -33.8688,
                                                       151.2093,
                                                       NULL,
                                                       0U,
                                                       &generated_id_b) == GEO_DATABASE_OK,
                "successive database-generated object ID must commit");
    TEST_ASSERT(generated_id_a >= (UINT64_C(1) << 63U) && generated_id_b == generated_id_a + 1U,
                "database-generated IDs must use the durable reserved monotonic range");

    static const char idempotency_key[] = "ingest-partition-7-offset-991";
    GeoDatabaseObjectMutation idempotent_mutation = {
        .object_id = 505U,
        .morton_code = geo_encode(37.7749, -122.4194),
        .operation = GEO_DATABASE_UPSERT,
    };
    bool replayed;
    TEST_ASSERT(geo_database_write_objects_idempotent(database,
                                                      &idempotent_mutation,
                                                      1U,
                                                      idempotency_key,
                                                      sizeof(idempotency_key) - 1U,
                                                      3600U,
                                                      &replayed) == GEO_DATABASE_OK &&
                    !replayed,
                "first idempotent write must commit normally");
    TEST_ASSERT(geo_database_write_objects_idempotent(database,
                                                      &idempotent_mutation,
                                                      1U,
                                                      idempotency_key,
                                                      sizeof(idempotency_key) - 1U,
                                                      3600U,
                                                      &replayed) == GEO_DATABASE_OK &&
                    replayed,
                "same idempotency key and payload must replay without another mutation");

    GeoDatabaseObjectMutation conflicting_retry = idempotent_mutation;
    conflicting_retry.morton_code = geo_encode(37.0, -122.0);
    TEST_ASSERT(geo_database_write_objects_idempotent(database,
                                                      &conflicting_retry,
                                                      1U,
                                                      idempotency_key,
                                                      sizeof(idempotency_key) - 1U,
                                                      3600U,
                                                      &replayed) == GEO_DATABASE_IDEMPOTENCY_CONFLICT,
                "reusing an idempotency key for different content must be rejected");

    GeoDatabaseStats stats;
    TEST_ASSERT(geo_database_get_stats(database, &stats), "database statistics must be readable");
    TEST_ASSERT(stats.committed_operations == 8U && stats.next_sequence == 9U &&
                    stats.active_spatial_operations == 8U && stats.frozen_spatial_operations == 0U,
                "durable catalog and active spatial memtable must track exact committed progress");
    TEST_ASSERT(geo_database_checkpoint(database) == GEO_DATABASE_OK, "checkpoint must flush canonical and derived state");
    TEST_ASSERT(geo_database_get_stats(database, &stats) &&
                    stats.active_spatial_operations == 0U && stats.frozen_spatial_operations == 0U &&
                    stats.active_segments >= 1U,
                "checkpoint must drain both spatial memtable generations into durable Morton segments");

    geo_database_close(database);
    database = NULL;
    size_t spatial_delta_count = SIZE_MAX;
    TEST_ASSERT(count_spatial_deltas(directory, &spatial_delta_count) && spatial_delta_count == 0U,
                "checkpoint must atomically prune the applied spatial-delta prefix");
    TEST_ASSERT(remove_derived_manifest(directory), "derived spatial manifest must be removable for rebuild validation");

    config.create_if_missing = false;
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK,
                "canonical objects must rebuild a missing derived spatial index automatically");
    TEST_ASSERT(geo_database_search_radius_reuse(database, 0.0, 0.0, 25000.0, result, NULL),
                "query after automatic rebuild must succeed");
    TEST_ASSERT(result->count == 5U && result_contains(result, 101U, update_and_delete[0].morton_code),
                "automatic rebuild must preserve exact last-write-wins state");
    TEST_ASSERT(geo_database_get(database, 101U, &object) == GEO_DATABASE_OK && object.document_size == 0U,
                "metadata-less update must replace the complete prior object document");
    geo_database_object_release(&object);

    TEST_ASSERT(geo_database_write_objects_idempotent(database,
                                                      &idempotent_mutation,
                                                      1U,
                                                      idempotency_key,
                                                      sizeof(idempotency_key) - 1U,
                                                      3600U,
                                                      &replayed) == GEO_DATABASE_OK &&
                    replayed,
                "idempotency decisions must survive close and reopen");

    uint64_t generated_id_after_reopen = 0U;
    TEST_ASSERT(geo_database_insert_generated_document(database,
                                                       0.0,
                                                       0.0,
                                                       NULL,
                                                       0U,
                                                       &generated_id_after_reopen) == GEO_DATABASE_OK &&
                    generated_id_after_reopen == generated_id_b + 1U,
                "database-generated ID allocator must resume monotonically after reopen");

cleanup:
    geo_database_object_release(&object);
    geo_doc_buffer_release(&document);
    geo_doc_builder_destroy(builder);
    geo_result_destroy(result);
    geo_database_close(database);

    if (directory && !geobolt_test_remove_tree(directory)) {
        fprintf(stderr, "FAIL: temporary database tree could not be removed\n");
        succeeded = false;
    }

    return succeeded;
}

static bool test_typed_secondary_indexes_transactional_reopen(void)
{
    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-secondary-index-XXXXXX";
    char *directory = mkdtemp(directory_template);
    GeoDatabase *database = NULL;
    GeoIdResult *result = NULL;
    GeoSearchResult *spatial_result = NULL;
    GeoDatabaseQueryWorkspace *query_workspace = NULL;
    GeoDocBuffer documents[4] = { 0 };
    GeoDocBuffer invalid_document = { 0 };
    const unsigned char tokens[4][2] = {
        { 0U, 1U },
        { 0U, 2U },
        { 1U, 0U },
        { 2U, 0U },
    };

    TEST_ASSERT(directory != NULL, "secondary-index database directory must be created");
    TEST_ASSERT(build_secondary_document(true, -10, 5U, 1.5, 1000, "alpha", tokens[0], documents),
                "first typed metadata document must serialize");
    TEST_ASSERT(build_secondary_document(false, 0, 10U, 2.5, 2000, "beta", tokens[1], documents + 1U),
                "second typed metadata document must serialize");
    TEST_ASSERT(build_secondary_document(true, 20, 15U, 3.5, 3000, "gamma", tokens[2], documents + 2U),
                "third typed metadata document must serialize");

    GeoDatabaseConfig config = test_database_config();
    GeoDatabaseStatus status;
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "secondary-index database must open");

    GeoDatabaseObjectMutation mutations[3];

    for (size_t index = 0U; index < 3U; ++index) {
        mutations[index] = (GeoDatabaseObjectMutation) {
            .object_id = index + 1U,
            .morton_code = geo_encode(-23.5 + (double) index, -46.6 + (double) index),
            .document = documents[index].data,
            .document_size = documents[index].size,
            .operation = GEO_DATABASE_UPSERT,
        };
    }

    TEST_ASSERT(geo_database_write_objects(database, mutations, 3U) == GEO_DATABASE_OK,
                "canonical objects must exist before online index construction");
    TEST_ASSERT(geo_database_create_index(database, "available_idx", "/available", GEO_DATABASE_INDEX_BOOL) == GEO_DATABASE_OK,
                "boolean secondary index must build from canonical objects");
    TEST_ASSERT(geo_database_create_index(database, "score_idx", "/score", GEO_DATABASE_INDEX_INT64) == GEO_DATABASE_OK,
                "signed secondary index must build from canonical objects");
    TEST_ASSERT(geo_database_create_index(database, "rank_idx", "/rank", GEO_DATABASE_INDEX_UINT64) == GEO_DATABASE_OK,
                "unsigned secondary index must build from canonical objects");
    TEST_ASSERT(geo_database_create_index(database, "price_idx", "/price", GEO_DATABASE_INDEX_DOUBLE) == GEO_DATABASE_OK,
                "floating-point secondary index must build from canonical objects");
    TEST_ASSERT(geo_database_create_index(database, "updated_idx", "/updated_at", GEO_DATABASE_INDEX_DATETIME) == GEO_DATABASE_OK,
                "datetime secondary index must build from canonical objects");
    TEST_ASSERT(geo_database_create_index(database, "label_idx", "/label", GEO_DATABASE_INDEX_STRING) == GEO_DATABASE_OK,
                "string secondary index must build from canonical objects");
    TEST_ASSERT(geo_database_create_index(database, "token_idx", "/token", GEO_DATABASE_INDEX_BYTES) == GEO_DATABASE_OK,
                "binary secondary index must build from canonical objects");
    TEST_ASSERT(geo_database_checkpoint(database) == GEO_DATABASE_OK,
                "typed-index fixture must publish an immutable spatial generation for density-aware planning");

    GeoDatabaseIndexInfo index_info[7];
    size_t index_count = 0U;
    TEST_ASSERT(geo_database_list_indexes(database, index_info, 7U, &index_count) == GEO_DATABASE_OK && index_count == 7U,
                "typed index catalog must list every ready definition");

    for (size_t index = 0U; index < index_count; ++index) {
        TEST_ASSERT(index_info[index].entry_count == 3U, "initial typed-index cardinality must be exact");
    }

    result = geo_id_result_create(4U);
    TEST_ASSERT(result != NULL, "secondary query result must allocate");
    spatial_result = geo_result_create(4U);
    query_workspace = geo_database_query_workspace_create(8U);
    TEST_ASSERT(spatial_result != NULL && query_workspace != NULL,
                "reusable conjunctive-query result and workspace must allocate");

    GeoDatabaseIndexPredicate predicate = {
        .index_name = "available_idx",
        .operation = GEO_DATABASE_INDEX_EQUAL,
        .lower = { .type = GEO_DATABASE_INDEX_BOOL, .as.boolean = true },
    };
    GeoDatabaseIndexQueryStats query_stats;
    TEST_ASSERT(geo_database_query_index_reuse(database, &predicate, result, &query_stats) == GEO_DATABASE_OK &&
                    result->count == 2U && id_result_contains(result, 1U) && id_result_contains(result, 3U) &&
                    query_stats.scanned_entries == 2U,
                "boolean equality must seek only the matching key range");

    GeoDatabaseIndexPredicate conjunctive_predicates[2] = {
        predicate,
        {
            .index_name = "score_idx",
            .operation = GEO_DATABASE_INDEX_LESS_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_INT64, .as.signed_integer = 0 },
        },
    };
    GeoDatabaseRadiusQuery radius_query = {
        .predicates = conjunctive_predicates,
        .predicate_count = 2U,
        .latitude = -23.5,
        .longitude = -46.6,
        .radius_km = 25000.0,
    };
    GeoDatabaseQueryStats planner_stats;
    TEST_ASSERT(geo_database_query_radius_reuse(database,
                                                &radius_query,
                                                query_workspace,
                                                spatial_result,
                                                &planner_stats) == GEO_DATABASE_OK &&
                    spatial_result->count == 1U && spatial_result->results[0].id == 1U &&
                    planner_stats.plan == GEO_DATABASE_QUERY_PLAN_SPATIAL && planner_stats.metadata_candidates == 1U &&
                    planner_stats.secondary_entries_scanned == 0U && planner_stats.object_lookups == 0U,
                "small spatial input must evaluate the complete metadata conjunction from the canonical object cache");

    geo_database_query_workspace_destroy(query_workspace);
    query_workspace = geo_database_query_workspace_create(8U);
    TEST_ASSERT(query_workspace != NULL,
                "a fresh workspace must support a spatial-driven query without prior membership state");

    radius_query.predicate_count = 1U;
    radius_query.radius_km = 10.0;
    TEST_ASSERT(geo_database_query_radius_reuse(database,
                                                &radius_query,
                                                query_workspace,
                                                spatial_result,
                                                &planner_stats) == GEO_DATABASE_OK &&
                    spatial_result->count == 1U && spatial_result->results[0].id == 1U &&
                    planner_stats.plan == GEO_DATABASE_QUERY_PLAN_SPATIAL && planner_stats.secondary_entries_scanned == 0U &&
                    planner_stats.object_lookups == 0U,
                "tight radius must evaluate broad metadata directly on cached canonical spatial candidates");

    GeoDatabaseIndexPredicate all_type_predicates[8] = {
        conjunctive_predicates[0],
        conjunctive_predicates[1],
        {
            .index_name = "rank_idx",
            .operation = GEO_DATABASE_INDEX_GREATER_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_UINT64, .as.unsigned_integer = 5U },
        },
        {
            .index_name = "price_idx",
            .operation = GEO_DATABASE_INDEX_GREATER_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_DOUBLE, .as.floating_point = 1.0 },
        },
        {
            .index_name = "price_idx",
            .operation = GEO_DATABASE_INDEX_LESS_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_DOUBLE, .as.floating_point = 2.0 },
        },
        {
            .index_name = "updated_idx",
            .operation = GEO_DATABASE_INDEX_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_DATETIME, .as.datetime = 1000 },
        },
        {
            .index_name = "label_idx",
            .operation = GEO_DATABASE_INDEX_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_STRING, .as.bytes = { .data = "alpha", .size = 5U } },
        },
        {
            .index_name = "token_idx",
            .operation = GEO_DATABASE_INDEX_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_BYTES, .as.bytes = { .data = tokens[0], .size = 2U } },
        },
    };

    radius_query.predicates = all_type_predicates;
    radius_query.predicate_count = 8U;
    TEST_ASSERT(geo_database_query_radius_reuse(database,
                                                &radius_query,
                                                query_workspace,
                                                spatial_result,
                                                &planner_stats) == GEO_DATABASE_OK &&
                    spatial_result->count == 1U && spatial_result->results[0].id == 1U &&
                    planner_stats.plan == GEO_DATABASE_QUERY_PLAN_SPATIAL && planner_stats.secondary_entries_scanned == 0U,
                "spatial-driven canonical filtering must preserve every typed comparison operation");

    all_type_predicates[7].lower.as.bytes.data = tokens[1];
    TEST_ASSERT(geo_database_query_radius_reuse(database,
                                                &radius_query,
                                                query_workspace,
                                                spatial_result,
                                                &planner_stats) == GEO_DATABASE_OK &&
                    spatial_result->count == 0U && planner_stats.plan == GEO_DATABASE_QUERY_PLAN_SPATIAL &&
                    planner_stats.secondary_entries_scanned == 0U,
                "a nonmatching final binary predicate must reject the spatial candidate");

    predicate = (GeoDatabaseIndexPredicate) {
        .index_name = "score_idx",
        .operation = GEO_DATABASE_INDEX_BETWEEN,
        .lower = { .type = GEO_DATABASE_INDEX_INT64, .as.signed_integer = -5 },
        .upper = { .type = GEO_DATABASE_INDEX_INT64, .as.signed_integer = 25 },
    };
    TEST_ASSERT(geo_database_query_index_reuse(database, &predicate, result, NULL) == GEO_DATABASE_OK && result->count == 2U &&
                    id_result_contains(result, 2U) && id_result_contains(result, 3U),
                "signed range must preserve negative-to-positive ordering");

    predicate = (GeoDatabaseIndexPredicate) {
        .index_name = "price_idx",
        .operation = GEO_DATABASE_INDEX_LESS_EQUAL,
        .lower = { .type = GEO_DATABASE_INDEX_DOUBLE, .as.floating_point = 2.5 },
    };
    TEST_ASSERT(geo_database_query_index_reuse(database, &predicate, result, NULL) == GEO_DATABASE_OK && result->count == 2U &&
                    id_result_contains(result, 1U) && id_result_contains(result, 2U),
                "floating-point range must preserve IEEE numeric ordering");

    predicate = (GeoDatabaseIndexPredicate) {
        .index_name = "label_idx",
        .operation = GEO_DATABASE_INDEX_GREATER_EQUAL,
        .lower = { .type = GEO_DATABASE_INDEX_STRING, .as.bytes = { .data = "beta", .size = 4U } },
    };
    TEST_ASSERT(geo_database_query_index_reuse(database, &predicate, result, NULL) == GEO_DATABASE_OK && result->count == 2U &&
                    id_result_contains(result, 2U) && id_result_contains(result, 3U),
                "string range must use preserved lexicographic ordering");

    predicate = (GeoDatabaseIndexPredicate) {
        .index_name = "token_idx",
        .operation = GEO_DATABASE_INDEX_EQUAL,
        .lower = { .type = GEO_DATABASE_INDEX_BYTES, .as.bytes = { .data = tokens[1], .size = 2U } },
    };
    TEST_ASSERT(geo_database_query_index_reuse(database, &predicate, result, NULL) == GEO_DATABASE_OK &&
                    result->count == 1U && result->ids[0] == 2U,
                "binary equality must handle embedded zero bytes without prefix ambiguity");

    TEST_ASSERT(build_secondary_document(true, 100, 30U, 4.5, 4000, "delta", tokens[3], documents + 3U),
                "updated typed metadata document must serialize");
    mutations[0] = (GeoDatabaseObjectMutation) {
        .object_id = 2U,
        .morton_code = geo_encode(-21.5, -44.6),
        .document = documents[3].data,
        .document_size = documents[3].size,
        .operation = GEO_DATABASE_UPDATE,
    };
    TEST_ASSERT(geo_database_write_objects(database, mutations, 1U) == GEO_DATABASE_OK,
                "object update must move every changed secondary key atomically");

    GeoDatabaseIndexPredicate updated_predicates[2] = {
        {
            .index_name = "available_idx",
            .operation = GEO_DATABASE_INDEX_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_BOOL, .as.boolean = true },
        },
        {
            .index_name = "label_idx",
            .operation = GEO_DATABASE_INDEX_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_STRING, .as.bytes = { .data = "delta", .size = 5U } },
        },
    };
    radius_query = (GeoDatabaseRadiusQuery) {
        .predicates = updated_predicates,
        .predicate_count = 2U,
        .latitude = -21.5,
        .longitude = -44.6,
        .radius_km = 1.0,
    };
    TEST_ASSERT(geo_database_query_radius_reuse(database,
                                                &radius_query,
                                                query_workspace,
                                                spatial_result,
                                                &planner_stats) == GEO_DATABASE_OK &&
                    spatial_result->count == 1U && spatial_result->results[0].id == 2U &&
                    planner_stats.object_lookups == 0U,
                "object update must publish new coordinates and metadata into the canonical cache atomically");

    predicate = (GeoDatabaseIndexPredicate) {
        .index_name = "rank_idx",
        .operation = GEO_DATABASE_INDEX_EQUAL,
        .lower = { .type = GEO_DATABASE_INDEX_UINT64, .as.unsigned_integer = 10U },
    };
    TEST_ASSERT(geo_database_query_index_reuse(database, &predicate, result, NULL) == GEO_DATABASE_OK && result->count == 0U,
                "updated object must disappear from its old unsigned key");

    TEST_ASSERT(geo_database_remove(database, 1U) == GEO_DATABASE_OK,
                "object deletion must remove every secondary key in the canonical transaction");
    TEST_ASSERT(geo_database_list_indexes(database, index_info, 7U, &index_count) == GEO_DATABASE_OK,
                "typed index catalog must remain readable after update and delete");

    for (size_t index = 0U; index < index_count; ++index) {
        TEST_ASSERT(index_info[index].entry_count == 2U, "typed-index cardinality must track deletion exactly");
    }

    GeoDocBuilder *invalid_builder = geo_doc_builder_create();
    TEST_ASSERT(invalid_builder != NULL &&
                    geo_doc_builder_add_string(invalid_builder, 0U, "score", "wrong", 5U, NULL) == GEO_DOC_OK &&
                    geo_doc_builder_finish(invalid_builder, &invalid_document) == GEO_DOC_OK,
                "wrong-type document must serialize before schema validation");
    geo_doc_builder_destroy(invalid_builder);
    invalid_builder = NULL;
    mutations[0].object_id = 3U;
    mutations[0].document = invalid_document.data;
    mutations[0].document_size = invalid_document.size;
    TEST_ASSERT(geo_database_write_objects(database, mutations, 1U) == GEO_DATABASE_INVALID_ARGUMENT,
                "present metadata with an indexed type mismatch must reject the entire write");

    geo_database_close(database);
    database = NULL;
    config.create_if_missing = false;
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "database with typed secondary indexes must reopen");

    predicate = (GeoDatabaseIndexPredicate) {
        .index_name = "updated_idx",
        .operation = GEO_DATABASE_INDEX_GREATER,
        .lower = { .type = GEO_DATABASE_INDEX_DATETIME, .as.datetime = 2500 },
    };
    TEST_ASSERT(geo_database_query_index_reuse(database, &predicate, result, NULL) == GEO_DATABASE_OK && result->count == 2U &&
                    id_result_contains(result, 2U) && id_result_contains(result, 3U),
                "datetime ordering and transactional maintenance must survive reopen");

    radius_query = (GeoDatabaseRadiusQuery) {
        .predicates = updated_predicates,
        .predicate_count = 2U,
        .latitude = -21.5,
        .longitude = -44.6,
        .radius_km = 1.0,
    };
    TEST_ASSERT(geo_database_query_radius_reuse(database,
                                                &radius_query,
                                                query_workspace,
                                                spatial_result,
                                                &planner_stats) == GEO_DATABASE_OK,
                "the first canonical query after reopen must succeed");
    TEST_ASSERT(spatial_result->count == 1U && spatial_result->results[0].id == 2U,
                "the first canonical query after reopen must return the updated durable object");
    TEST_ASSERT(planner_stats.object_lookups == 2U,
                "the first canonical query after reopen must load and validate both durable spatial candidates");

    TEST_ASSERT(geo_database_query_radius_reuse(database,
                                                &radius_query,
                                                query_workspace,
                                                spatial_result,
                                                &planner_stats) == GEO_DATABASE_OK,
                "a repeated canonical query after reopen must succeed");
    TEST_ASSERT(spatial_result->count == 1U && spatial_result->results[0].id == 2U,
                "a repeated canonical query after reopen must preserve the exact durable result");
    TEST_ASSERT(planner_stats.object_lookups == 0U,
                "a repeated query after reopen must reuse the validated canonical cache entry");

    GeoDatabaseIndexPredicate reordered_predicates[2] = {
        {
            .index_name = "available_idx",
            .operation = GEO_DATABASE_INDEX_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_BOOL, .as.boolean = true },
        },
        {
            .index_name = "label_idx",
            .operation = GEO_DATABASE_INDEX_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_STRING, .as.bytes = { .data = "aardvark", .size = 8U } },
        },
    };
    radius_query = (GeoDatabaseRadiusQuery) {
        .predicates = reordered_predicates,
        .predicate_count = 2U,
        .latitude = 0.0,
        .longitude = 0.0,
        .radius_km = 25000.0,
    };
    TEST_ASSERT(geo_database_query_radius_reuse(database,
                                                &radius_query,
                                                query_workspace,
                                                spatial_result,
                                                &planner_stats) == GEO_DATABASE_OK &&
                    spatial_result->count == 0U && planner_stats.metadata_candidates == 0U &&
                    planner_stats.secondary_entries_scanned == 0U,
                "persisted selectivity statistics must run an empty selective predicate before a broad predicate");

    TEST_ASSERT(geo_database_drop_index(database, "score_idx") == GEO_DATABASE_OK,
                "dropping a secondary index must remove its catalog definition and key range atomically");
    predicate.index_name = "score_idx";
    predicate.lower.type = GEO_DATABASE_INDEX_INT64;
    predicate.lower.as.signed_integer = 0;
    TEST_ASSERT(geo_database_query_index_reuse(database, &predicate, result, NULL) == GEO_DATABASE_NOT_FOUND,
                "a dropped secondary index must not remain queryable");

    TEST_ASSERT(geo_database_remove(database, 2U) == GEO_DATABASE_OK &&
                    geo_database_remove(database, 3U) == GEO_DATABASE_OK,
                "deleting the final indexed objects must commit zeroed histogram counters");

    geo_database_close(database);
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK,
                "an index retaining static bins with zero current entries must reopen");
    TEST_ASSERT(geo_database_list_indexes(database, index_info, 7U, &index_count) == GEO_DATABASE_OK && index_count == 6U,
                "all remaining typed-index definitions must survive zero-cardinality reopen");

    for (size_t index = 0U; index < index_count; ++index) {
        TEST_ASSERT(index_info[index].entry_count == 0U,
                    "zero-cardinality reopen must preserve exact per-index counts");
    }

    predicate = (GeoDatabaseIndexPredicate) {
        .index_name = "label_idx",
        .operation = GEO_DATABASE_INDEX_EQUAL,
        .lower = { .type = GEO_DATABASE_INDEX_STRING, .as.bytes = { .data = "gamma", .size = 5U } },
    };
    TEST_ASSERT(geo_database_query_index_reuse(database, &predicate, result, NULL) == GEO_DATABASE_OK && result->count == 0U,
                "zeroed persisted histogram bins must never retain stale query results");

cleanup:
    geo_database_query_workspace_destroy(query_workspace);
    geo_result_destroy(spatial_result);
    geo_id_result_destroy(result);
    geo_database_close(database);
    geo_doc_buffer_release(&invalid_document);

    for (size_t index = 0U; index < 4U; ++index) {
        geo_doc_buffer_release(documents + index);
    }

    if (directory && !geobolt_test_remove_tree(directory)) {
        fprintf(stderr, "FAIL: secondary-index database tree could not be removed\n");
        succeeded = false;
    }

    return succeeded;
}

static bool test_query_planner_chunked_canonical_reads(void)
{
    enum {
        MATCHING_OBJECTS = 257,
        TOTAL_OBJECTS = 1025,
    };

    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-query-multiget-XXXXXX";
    char *directory = mkdtemp(directory_template);
    GeoDatabase *database = NULL;
    GeoDatabaseQueryWorkspace *workspace = NULL;
    GeoSearchResult *result = NULL;
    GeoDatabaseObjectMutation *mutations = NULL;
    GeoDocBuilder *builder = NULL;
    GeoDocBuffer matching_document = { 0 };
    GeoDocBuffer filler_document = { 0 };

    TEST_ASSERT(directory != NULL, "chunked-query database directory must be created");

    builder = geo_doc_builder_create();
    TEST_ASSERT(builder != NULL &&
                    geo_doc_builder_add_string(builder, 0U, "group", "matching", 8U, NULL) == GEO_DOC_OK &&
                    geo_doc_builder_finish(builder, &matching_document) == GEO_DOC_OK,
                "matching planner document must serialize");
    geo_doc_builder_destroy(builder);
    builder = geo_doc_builder_create();
    TEST_ASSERT(builder != NULL &&
                    geo_doc_builder_add_string(builder, 0U, "group", "filler", 6U, NULL) == GEO_DOC_OK &&
                    geo_doc_builder_finish(builder, &filler_document) == GEO_DOC_OK,
                "filler planner document must serialize");
    geo_doc_builder_destroy(builder);
    builder = NULL;

    GeoDatabaseConfig config = test_database_config();
    config.object_cache_bytes = 0U;
    config.group_commit_max_operations = TOTAL_OBJECTS;
    config.spatial_memtable_max_operations = TOTAL_OBJECTS * 2U;
    GeoDatabaseStatus status;

    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "chunked-query database must open");

    mutations = calloc(TOTAL_OBJECTS, sizeof(*mutations));
    TEST_ASSERT(mutations != NULL, "chunked-query mutations must allocate");

    for (size_t index = 0U; index < TOTAL_OBJECTS; ++index) {
        bool matching = index < MATCHING_OBJECTS;

        mutations[index] = (GeoDatabaseObjectMutation) {
            .object_id = index + 1U,
            .morton_code = geo_encode(-80.0 + (double) index * 0.1, -170.0 + (double) index * 0.2),
            .document = matching ? matching_document.data : filler_document.data,
            .document_size = matching ? matching_document.size : filler_document.size,
            .operation = GEO_DATABASE_UPSERT,
        };
    }

    TEST_ASSERT(geo_database_write_objects(database, mutations, TOTAL_OBJECTS) == GEO_DATABASE_OK,
                "planner dataset must commit atomically");
    TEST_ASSERT(geo_database_create_index(database, "group_idx", "/group", GEO_DATABASE_INDEX_STRING) == GEO_DATABASE_OK,
                "planner string index must build over the complete dataset");

    workspace = geo_database_query_workspace_create(1U);
    result = geo_result_create(MATCHING_OBJECTS);
    TEST_ASSERT(workspace != NULL && result != NULL, "chunked-query reusable storage must allocate");

    GeoDatabaseIndexPredicate predicate = {
        .index_name = "group_idx",
        .operation = GEO_DATABASE_INDEX_EQUAL,
        .lower = { .type = GEO_DATABASE_INDEX_STRING, .as.bytes = { .data = "matching", .size = 8U } },
    };
    GeoDatabaseRadiusQuery query = {
        .predicates = &predicate,
        .predicate_count = 1U,
        .latitude = 0.0,
        .longitude = 0.0,
        .radius_km = 25000.0,
    };
    GeoDatabaseQueryStats query_stats;

    TEST_ASSERT(geo_database_query_radius_reuse(database, &query, workspace, result, &query_stats) == GEO_DATABASE_OK,
                "secondary-driven query must execute across a full MultiGet chunk and tail");
    TEST_ASSERT(query_stats.plan == GEO_DATABASE_QUERY_PLAN_SECONDARY &&
                    query_stats.metadata_candidates == MATCHING_OBJECTS &&
                    query_stats.object_lookups == MATCHING_OBJECTS && result->count == MATCHING_OBJECTS,
                "chunked canonical reads must preserve exact planner cardinality and lookup accounting");

    for (size_t index = 0U; index < result->count; ++index) {
        TEST_ASSERT(result->results[index].id > 0U && result->results[index].id <= MATCHING_OBJECTS,
                    "chunked canonical reads must not leak filler objects across the metadata conjunction");
    }

cleanup:
    geo_doc_builder_destroy(builder);
    geo_doc_buffer_release(&filler_document);
    geo_doc_buffer_release(&matching_document);
    free(mutations);
    geo_result_destroy(result);
    geo_database_query_workspace_destroy(workspace);
    geo_database_close(database);

    if (directory && !geobolt_test_remove_tree(directory)) {
        fprintf(stderr, "FAIL: chunked-query database tree could not be removed\n");
        succeeded = false;
    }

    return succeeded;
}

static bool test_secondary_histogram_corruption_rejected(void)
{
    static const unsigned char histogram_metadata_prefix[] = {
        'i', 'n', 'd', 'e', 'x', '-', 'h', 'i', 's', 't', '-', 'm', 'e', 't', 'a', '/',
    };

    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-histogram-corruption-XXXXXX";
    char rocks_path[512];
    char *directory = mkdtemp(directory_template);
    GeoDatabase *database = NULL;
    GeoDocBuilder *builder = NULL;
    GeoDocBuffer document = { 0 };
    GeoRocksDatabase *rocks = NULL;
    GeoRocksBatch *batch = NULL;
    GeoRocksBuffer metadata = { 0 };

    TEST_ASSERT(directory != NULL, "histogram-corruption database directory must be created");
    int path_length = snprintf(rocks_path, sizeof(rocks_path), "%s/rocksdb", directory);
    TEST_ASSERT(path_length > 0 && (size_t) path_length < sizeof(rocks_path),
                "histogram-corruption RocksDB path must fit");

    builder = geo_doc_builder_create();
    TEST_ASSERT(builder != NULL &&
                    geo_doc_builder_add_string(builder, 0U, "label", "alpha", 5U, NULL) == GEO_DOC_OK &&
                    geo_doc_builder_finish(builder, &document) == GEO_DOC_OK,
                "histogram-corruption indexed document must serialize");

    GeoDatabaseConfig database_config = test_database_config();
    GeoDatabaseStatus database_status;
    database = geo_database_open(directory, &database_config, &database_status);
    TEST_ASSERT(database != NULL && database_status == GEO_DATABASE_OK,
                "histogram-corruption database must open");
    TEST_ASSERT(geo_database_upsert_document(database,
                                             1U,
                                             0.0,
                                             0.0,
                                             document.data,
                                             document.size) == GEO_DATABASE_OK &&
                    geo_database_create_index(database,
                                              "label_idx",
                                              "/label",
                                              GEO_DATABASE_INDEX_STRING) == GEO_DATABASE_OK,
                "histogram-corruption fixture must persist a nonempty histogram");
    geo_database_close(database);
    database = NULL;

    GeoRocksConfig rocks_config = geo_rocks_default_config();
    rocks_config.create_if_missing = false;
    GeoRocksStatus rocks_status;

    rocks = geo_rocks_open(rocks_path, &rocks_config, &rocks_status);
    TEST_ASSERT(rocks != NULL, "histogram-corruption canonical RocksDB must open directly");

    unsigned char metadata_key[sizeof(histogram_metadata_prefix) + 8U];

    memcpy(metadata_key, histogram_metadata_prefix, sizeof(histogram_metadata_prefix));
    geo_db_store_u64_be(metadata_key + sizeof(histogram_metadata_prefix), 1U);
    TEST_ASSERT(geo_rocks_get(rocks,
                              GEO_ROCKS_CF_CATALOG,
                              NULL,
                              metadata_key,
                              sizeof(metadata_key),
                              &metadata,
                              &rocks_status) &&
                    metadata.size > 40U,
                "histogram-corruption metadata record must be readable");

    ((unsigned char *) metadata.data)[metadata.size - 1U] ^= UINT8_C(0x01);
    batch = geo_rocks_batch_create(rocks, metadata.size + sizeof(metadata_key), &rocks_status);
    TEST_ASSERT(batch &&
                    geo_rocks_batch_put(batch,
                                        GEO_ROCKS_CF_CATALOG,
                                        metadata_key,
                                        sizeof(metadata_key),
                                        metadata.data,
                                        metadata.size,
                                        &rocks_status) &&
                    geo_rocks_write(rocks, batch, true, &rocks_status),
                "histogram-corruption fixture must durably alter the checksummed metadata");

    geo_rocks_batch_destroy(batch);
    batch = NULL;
    geo_rocks_buffer_release(&metadata);
    geo_rocks_close(rocks);
    rocks = NULL;

    database_config.create_if_missing = false;
    database = geo_database_open(directory, &database_config, &database_status);
    TEST_ASSERT(database == NULL && database_status == GEO_DATABASE_CORRUPTION,
                "database open must reject corrupt histogram boundaries before catalog publication");

cleanup:
    geo_rocks_batch_destroy(batch);
    geo_rocks_buffer_release(&metadata);
    geo_rocks_close(rocks);
    geo_database_close(database);
    geo_doc_buffer_release(&document);
    geo_doc_builder_destroy(builder);

    if (directory && !geobolt_test_remove_tree(directory)) {
        fprintf(stderr, "FAIL: histogram-corruption database tree could not be removed\n");
        succeeded = false;
    }

    return succeeded;
}

typedef struct {
    const GeoDatabase *database;
    atomic_bool *start;
    atomic_bool *stop;
    atomic_bool *failed;
    size_t expected_count;
    size_t queries;
} DatabaseReaderContext;

static void *database_reader_main(void *argument)
{
    DatabaseReaderContext *context = argument;
    GeoDatabaseQueryWorkspace *workspace = geo_database_query_workspace_create(1U);
    GeoSearchResult *result = geo_result_create(context->expected_count);
    GeoDatabaseIndexPredicate predicate = {
        .index_name = "active_idx",
        .operation = GEO_DATABASE_INDEX_EQUAL,
        .lower = { .type = GEO_DATABASE_INDEX_BOOL, .as.boolean = true },
    };
    GeoDatabaseRadiusQuery query = {
        .predicates = &predicate,
        .predicate_count = 1U,
        .latitude = 0.0,
        .longitude = 0.0,
        .radius_km = 25000.0,
    };

    if (!workspace || !result) {
        atomic_store_explicit(context->failed, true, memory_order_release);
        geo_result_destroy(result);
        geo_database_query_workspace_destroy(workspace);
        return NULL;
    }

    while (!atomic_load_explicit(context->start, memory_order_acquire)) {
    }

    while (!atomic_load_explicit(context->stop, memory_order_acquire)) {
        if (geo_database_query_radius_reuse(context->database,
                                            &query,
                                            workspace,
                                            result,
                                            NULL) != GEO_DATABASE_OK ||
            result->count != context->expected_count) {
            atomic_store_explicit(context->failed, true, memory_order_release);
            break;
        }

        context->queries++;
    }

    geo_result_destroy(result);
    geo_database_query_workspace_destroy(workspace);
    return NULL;
}

static bool test_concurrent_insert_delete_batches_and_queries(void)
{
    enum {
        OBJECT_COUNT = 256,
        READER_COUNT = 3,
        WRITE_ROUNDS = 12,
    };

    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-database-concurrent-XXXXXX";
    char *directory = mkdtemp(directory_template);
    GeoDatabase *database = NULL;
    GeoDatabaseObjectMutation *mutations = NULL;
    GeoDocBuilder *builder = NULL;
    GeoDocBuffer document = { 0 };
    pthread_t readers[READER_COUNT];
    DatabaseReaderContext contexts[READER_COUNT];
    size_t readers_started = 0U;
    atomic_bool start = ATOMIC_VAR_INIT(false);
    atomic_bool stop = ATOMIC_VAR_INIT(false);
    atomic_bool failed = ATOMIC_VAR_INIT(false);

    TEST_ASSERT(directory != NULL, "concurrent test directory must be created");

    builder = geo_doc_builder_create();
    TEST_ASSERT(builder &&
                    geo_doc_builder_add_bool(builder, 0U, "active", true, NULL) == GEO_DOC_OK &&
                    geo_doc_builder_finish(builder, &document) == GEO_DOC_OK,
                "concurrent canonical-cache metadata must serialize");

    GeoDatabaseConfig config = test_database_config();
    config.group_commit_max_operations = OBJECT_COUNT * 2U;
    config.spatial_memtable_max_operations = OBJECT_COUNT * 2U;
    GeoDatabaseStatus status;
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "concurrent database must open");

    mutations = calloc(OBJECT_COUNT * 2U, sizeof(*mutations));
    TEST_ASSERT(mutations != NULL, "concurrent mutation batch must allocate");

    for (size_t object_index = 0; object_index < OBJECT_COUNT; ++object_index) {
        mutations[object_index] = (GeoDatabaseObjectMutation) {
            .object_id = object_index + 1U,
            .morton_code = geo_encode(-40.0 + (double) object_index * 0.25, -120.0 + (double) object_index * 0.5),
            .document = document.data,
            .document_size = document.size,
            .operation = GEO_DATABASE_UPSERT,
        };
    }

    TEST_ASSERT(geo_database_write_objects(database, mutations, OBJECT_COUNT) == GEO_DATABASE_OK,
                "initial concurrent batch must commit");
    TEST_ASSERT(geo_database_create_index(database, "active_idx", "/active", GEO_DATABASE_INDEX_BOOL) == GEO_DATABASE_OK,
                "concurrent cache queries require a production boolean secondary index");

    geo_database_close(database);
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK,
                "concurrent database must reopen with an empty process-local cache");

    for (size_t reader = 0; reader < READER_COUNT; ++reader) {
        contexts[reader] = (DatabaseReaderContext) {
            .database = database,
            .start = &start,
            .stop = &stop,
            .failed = &failed,
            .expected_count = OBJECT_COUNT,
        };
        TEST_ASSERT(pthread_create(readers + reader, NULL, database_reader_main, contexts + reader) == 0,
                    "database query reader must start");
        readers_started++;
    }

    atomic_store_explicit(&start, true, memory_order_release);

    for (size_t round = 0; round < WRITE_ROUNDS; ++round) {
        uint64_t old_base = round % 2U ? 1001U : 1U;
        uint64_t new_base = round % 2U ? 1U : 1001U;

        for (size_t object_index = 0; object_index < OBJECT_COUNT; ++object_index) {
            mutations[object_index] = (GeoDatabaseObjectMutation) {
                .object_id = old_base + object_index,
                .operation = GEO_DATABASE_DELETE,
            };
            mutations[OBJECT_COUNT + object_index] = (GeoDatabaseObjectMutation) {
                .object_id = new_base + object_index,
                .morton_code = geo_encode(-60.0 + (double) object_index * 0.25,
                                          -170.0 + (double) ((object_index * 7U + round * 29U) % 680U) * 0.5),
                .document = document.data,
                .document_size = document.size,
                .operation = GEO_DATABASE_UPSERT,
            };
        }

        TEST_ASSERT(geo_database_write_objects(database, mutations, OBJECT_COUNT * 2U) == GEO_DATABASE_OK,
                    "simultaneous randomized insert/delete replacement batch must commit");
    }

    atomic_store_explicit(&stop, true, memory_order_release);

    for (size_t reader = 0; reader < readers_started; ++reader) {
        TEST_ASSERT(pthread_join(readers[reader], NULL) == 0, "database query reader must join");
    }
    readers_started = 0U;

    TEST_ASSERT(!atomic_load_explicit(&failed, memory_order_acquire),
                "readers must never observe a partially published insert/delete batch");

    size_t total_queries = 0U;

    for (size_t reader = 0; reader < READER_COUNT; ++reader) {
        total_queries += contexts[reader].queries;
    }
    TEST_ASSERT(total_queries > 0U, "concurrent readers must execute queries during writes");

cleanup:
    atomic_store_explicit(&stop, true, memory_order_release);
    atomic_store_explicit(&start, true, memory_order_release);

    for (size_t reader = 0; reader < readers_started; ++reader) {
        (void) pthread_join(readers[reader], NULL);
    }

    geo_doc_buffer_release(&document);
    geo_doc_builder_destroy(builder);
    free(mutations);
    geo_database_close(database);

    if (directory && !geobolt_test_remove_tree(directory)) {
        fprintf(stderr, "FAIL: concurrent database tree could not be removed\n");
        succeeded = false;
    }

    return succeeded;
}

typedef struct {
    GeoDatabase *database;
    atomic_bool *start;
    uint64_t object_id;
    GeoDatabaseStatus status;
} GroupCommitWriterContext;

static void *group_commit_writer_main(void *argument)
{
    GroupCommitWriterContext *context = argument;

    while (!atomic_load_explicit(context->start, memory_order_acquire)) {
    }

    context->status = geo_database_upsert(context->database,
                                          context->object_id,
                                          -20.0 + (double) context->object_id * 0.01,
                                          -40.0 + (double) context->object_id * 0.01);
    return NULL;
}

static bool test_group_commit_coalesces_concurrent_writers(void)
{
    enum {
        WRITER_COUNT = 24,
    };

    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-database-group-XXXXXX";
    char *directory = mkdtemp(directory_template);
    GeoDatabase *database = NULL;
    pthread_t writers[WRITER_COUNT];
    GroupCommitWriterContext contexts[WRITER_COUNT];
    size_t writers_started = 0U;
    atomic_bool start = ATOMIC_VAR_INIT(false);

    TEST_ASSERT(directory != NULL, "group-commit test directory must be created");

    GeoDatabaseConfig config = test_database_config();
    config.group_commit_delay_us = 20000U;
    config.group_commit_max_operations = 16U;
    config.spatial_memtable_max_operations = 16U;
    config.max_active_segments = 64U;
    GeoDatabaseStatus status;
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "group-commit database must open");

    for (size_t writer = 0; writer < WRITER_COUNT; ++writer) {
        contexts[writer] = (GroupCommitWriterContext) {
            .database = database,
            .start = &start,
            .object_id = writer + 1U,
            .status = GEO_DATABASE_FAILED_STATE,
        };
        TEST_ASSERT(pthread_create(writers + writer, NULL, group_commit_writer_main, contexts + writer) == 0,
                    "group-commit writer must start");
        writers_started++;
    }

    atomic_store_explicit(&start, true, memory_order_release);

    for (size_t writer = 0; writer < writers_started; ++writer) {
        TEST_ASSERT(pthread_join(writers[writer], NULL) == 0, "group-commit writer must join");
        TEST_ASSERT(contexts[writer].status == GEO_DATABASE_OK, "every coalesced write must report its durable outcome");
    }
    writers_started = 0U;

    GeoDatabaseStats stats;
    size_t count = 0U;
    TEST_ASSERT(geo_database_get_stats(database, &stats), "group-commit statistics must be readable");
    TEST_ASSERT(stats.committed_operations == WRITER_COUNT && stats.active_segments <= 4U,
                "concurrent unitary writes must coalesce into far fewer spatial publications");
    TEST_ASSERT(geo_database_search_radius_count(database, 0.0, 0.0, 25000.0, &count, NULL) && count == WRITER_COUNT,
                "coalesced commit must publish every object exactly once");

cleanup:
    atomic_store_explicit(&start, true, memory_order_release);

    for (size_t writer = 0; writer < writers_started; ++writer) {
        (void) pthread_join(writers[writer], NULL);
    }

    geo_database_close(database);

    if (directory && !geobolt_test_remove_tree(directory)) {
        fprintf(stderr, "FAIL: group-commit database tree could not be removed\n");
        succeeded = false;
    }

    return succeeded;
}

int main(void)
{
    printf("========================================\n");
    printf("GEOBOLT DATABASE TEST SUITE\n");
    printf("========================================\n\n");

    if (!test_database_objects_metadata_reopen_and_rebuild()) {
        return EXIT_FAILURE;
    }
    printf("[PASS] RocksDB objects, GeoDoc metadata, reopen, delete, and derived-index rebuild\n");

    if (!test_database_spatial_watermark_crash_recovery()) {
        return EXIT_FAILURE;
    }
    printf("[PASS] spatial watermark crash windows and canonical rebuild recovery\n");

    if (!test_typed_secondary_indexes_transactional_reopen()) {
        return EXIT_FAILURE;
    }
    printf("[PASS] typed indexes and cost-planned geo conjunctions build, query, reopen, and drop transactionally\n");

    if (!test_query_planner_chunked_canonical_reads()) {
        return EXIT_FAILURE;
    }
    printf("[PASS] secondary-driven planner reuses pinned MultiGet chunks and preserves exact tails\n");

    if (!test_secondary_histogram_corruption_rejected()) {
        return EXIT_FAILURE;
    }
    printf("[PASS] checksummed histogram metadata rejects corruption before catalog publication\n");

    if (!test_concurrent_insert_delete_batches_and_queries()) {
        return EXIT_FAILURE;
    }
    printf("[PASS] concurrent atomic insert/delete publication and snapshot queries\n");

    if (!test_group_commit_coalesces_concurrent_writers()) {
        return EXIT_FAILURE;
    }
    printf("[PASS] persistent group-commit coordinator coalesces concurrent unitary writers\n");
    return EXIT_SUCCESS;
}
