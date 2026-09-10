#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "benchmark_common.h"
#include "geobolt/geobolt.h"
#include "geobolt/geodoc.h"
#include "test_support.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

static double elapsed_seconds(const struct timespec *begin, const struct timespec *end)
{
    return (double) (end->tv_sec - begin->tv_sec) + (double) (end->tv_nsec - begin->tv_nsec) / 1000000000.0;
}

int main(int argc, char **argv)
{
    size_t object_count = 100000U;
    size_t query_count = 50U;
    size_t object_cache_megabytes = 128U;

    if ((argc > 1 && !geo_benchmark_parse_size(argv[1], &object_count)) ||
        (argc > 2 && !geo_benchmark_parse_size(argv[2], &query_count)) ||
        (argc > 3 && !geo_benchmark_parse_size(argv[3], &object_cache_megabytes)) ||
        argc > 4 || !object_count || !query_count || object_count > SIZE_MAX / sizeof(GeoDatabaseObjectMutation) ||
        object_cache_megabytes > SIZE_MAX / (1024U * 1024U)) {
        fprintf(stderr, "usage: %s [objects] [queries] [object-cache-mb]\n", argv[0]);
        return EXIT_FAILURE;
    }

    char directory_template[] = "/tmp/geobolt-planner-benchmark-XXXXXX";
    char *directory = mkdtemp(directory_template);
    const char *failure_stage = "temporary directory creation";
    bool succeeded = directory != NULL;
    GeoDatabase *database = NULL;
    GeoDatabaseQueryWorkspace *workspace = NULL;
    GeoSearchResult *result = NULL;
    GeoDatabaseObjectMutation *mutations = NULL;
    GeoDocBuilder *builder = NULL;
    GeoDocBuffer document = { 0 };
    uint64_t result_checksum = 0U;
    uint64_t scanned_entries = 0U;
    uint64_t object_lookups = 0U;
    uint64_t ranges_checked = 0U;
    uint64_t estimated_spatial_candidates = 0U;
    size_t result_count = 0U;

    if (succeeded) {
        builder = geo_doc_builder_create();
        failure_stage = "metadata serialization";
        succeeded = builder &&
                    geo_doc_builder_add_bool(builder, 0U, "available", true, NULL) == GEO_DOC_OK &&
                    geo_doc_builder_finish(builder, &document) == GEO_DOC_OK;
    }

    GeoDatabaseConfig config = geo_database_default_config();

    config.object_cache_bytes = object_cache_megabytes * 1024U * 1024U;
    config.group_commit_max_operations = object_count;
    config.spatial_memtable_max_operations = object_count;

    if (succeeded) {
        GeoDatabaseStatus status;

        failure_stage = "database open";
        database = geo_database_open(directory, &config, &status);
        succeeded = database && status == GEO_DATABASE_OK;
    }

    if (succeeded) {
        failure_stage = "mutation allocation";
        mutations = malloc(object_count * sizeof(*mutations));
        succeeded = mutations != NULL;
    }

    size_t grid_width = 1U;

    while (grid_width <= object_count / grid_width) {
        grid_width++;
    }

    for (size_t object = 0U; succeeded && object < object_count; ++object) {
        size_t row = object / grid_width;
        size_t column = object % grid_width;
        double latitude = -24.0 + 2.0 * (double) row / (double) grid_width;
        double longitude = -47.0 + 2.0 * (double) column / (double) grid_width;

        mutations[object] = (GeoDatabaseObjectMutation) {
            .object_id = object + 1U,
            .morton_code = geo_encode(latitude, longitude),
            .document = document.data,
            .document_size = document.size,
            .operation = GEO_DATABASE_UPSERT,
        };
    }

    if (succeeded) {
        failure_stage = "bulk load and index build";
        succeeded = geo_database_write_objects(database, mutations, object_count) == GEO_DATABASE_OK &&
                    geo_database_create_index(database,
                                              "available_idx",
                                              "/available",
                                              GEO_DATABASE_INDEX_BOOL) == GEO_DATABASE_OK &&
                    geo_database_checkpoint(database) == GEO_DATABASE_OK;
    }

    if (succeeded) {
        failure_stage = "query workspace allocation";
        workspace = geo_database_query_workspace_create(1U);
        result = geo_result_create(64U);
        succeeded = workspace && result;
    }

    GeoDatabaseIndexPredicate predicate = {
        .index_name = "available_idx",
        .operation = GEO_DATABASE_INDEX_EQUAL,
        .lower = { .type = GEO_DATABASE_INDEX_BOOL, .as.boolean = true },
    };
    GeoDatabaseRadiusQuery query = {
        .predicates = &predicate,
        .predicate_count = 1U,
        .latitude = -23.5,
        .longitude = -46.5,
        .radius_km = 1.0,
    };
    GeoDatabaseQueryStats stats;

    if (succeeded) {
        failure_stage = "warmup query";
        succeeded = geo_database_query_radius_reuse(database, &query, workspace, result, &stats) == GEO_DATABASE_OK &&
                    stats.plan == GEO_DATABASE_QUERY_PLAN_SPATIAL;
    }

    struct timespec begin = { 0 };
    struct timespec end = { 0 };

    if (succeeded) {
        clock_gettime(CLOCK_MONOTONIC, &begin);
    }

    for (size_t iteration = 0U; succeeded && iteration < query_count; ++iteration) {
        query.latitude = -23.9 + 1.8 * (double) ((iteration * 37U) % 997U) / 997.0;
        query.longitude = -46.9 + 1.8 * (double) ((iteration * 71U) % 991U) / 991.0;
        failure_stage = "measured query";
        succeeded = geo_database_query_radius_reuse(database, &query, workspace, result, &stats) == GEO_DATABASE_OK &&
                    stats.plan == GEO_DATABASE_QUERY_PLAN_SPATIAL;

        if (succeeded) {
            scanned_entries += stats.secondary_entries_scanned;
            object_lookups += stats.object_lookups;
            ranges_checked += stats.spatial.ranges_checked;
            estimated_spatial_candidates += stats.estimated_spatial_candidates;
            result_count += result->count;

            for (size_t index = 0U; index < result->count; ++index) {
                result_checksum ^= result->results[index].id * UINT64_C(0x9e3779b97f4a7c15);
            }
        }
    }

    if (succeeded) {
        clock_gettime(CLOCK_MONOTONIC, &end);
    }

    double seconds = succeeded ? elapsed_seconds(&begin, &end) : 0.0;

    free(mutations);
    geo_result_destroy(result);
    geo_database_query_workspace_destroy(workspace);
    geo_database_close(database);
    geo_doc_buffer_release(&document);
    geo_doc_builder_destroy(builder);

    if (directory && !geobolt_test_remove_tree(directory)) {
        failure_stage = "temporary directory cleanup";
        succeeded = false;
    }

    if (!succeeded) {
        fprintf(stderr, "query-planner benchmark failed during %s\n", failure_stage);
        return EXIT_FAILURE;
    }

    geo_benchmark_print_environment();
    printf("objects=%zu queries=%zu object_cache_mb=%zu radius_km=1 metadata_selectivity=1.0\n",
           object_count,
           query_count,
           object_cache_megabytes);
    printf("seconds=%.6f queries_per_second=%.2f mean_ms=%.6f\n",
           seconds,
           (double) query_count / seconds,
           seconds * 1000.0 / (double) query_count);
    printf("secondary_entries_scanned=%llu object_lookups=%llu result_count=%zu checksum=%llu\n",
           (unsigned long long) scanned_entries,
           (unsigned long long) object_lookups,
           result_count,
           (unsigned long long) result_checksum);
    printf("mean_ranges_checked=%.2f mean_estimated_spatial_candidates=%.2f\n",
           (double) ranges_checked / (double) query_count,
           (double) estimated_spatial_candidates / (double) query_count);
    return EXIT_SUCCESS;
}
