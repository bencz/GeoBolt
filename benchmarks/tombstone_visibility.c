#include "benchmark_common.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

typedef struct {
    double elapsed_ms;
    uint64_t result_checksum;
    uint64_t records_scanned;
} QueryMeasurement;

static bool measure_global_counts(GeoSegmentSet *set, size_t iterations, QueryMeasurement *measurement)
{
    size_t warmup_count;
    GeoSearchStats warmup_stats;

    if (!geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &warmup_count, &warmup_stats)) {
        return false;
    }

    uint64_t result_checksum = 0;
    uint64_t records_scanned = 0;
    double start = geo_get_time_ms();

    for (size_t iteration = 0; iteration < iterations; ++iteration) {
        size_t count = 0;
        GeoSearchStats stats;

        if (!geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &count, &stats)) {
            return false;
        }

        result_checksum += count;
        records_scanned += stats.records_scanned;
    }

    measurement->elapsed_ms = geo_get_time_ms() - start;
    measurement->result_checksum = result_checksum;
    measurement->records_scanned = records_scanned;

    return true;
}

int main(int argc, char **argv)
{
    size_t record_count = 1000000U;
    size_t iterations = 20U;
    size_t removal_count = 1U;

    if ((argc > 1 && !geo_benchmark_parse_size(argv[1], &record_count)) ||
        (argc > 2 && !geo_benchmark_parse_size(argv[2], &iterations)) ||
        (argc > 3 && !geo_benchmark_parse_size(argv[3], &removal_count)) ||
        argc > 4 ||
        !record_count ||
        !iterations ||
        !removal_count ||
        removal_count > record_count ||
        removal_count > SIZE_MAX / sizeof(uint64_t)) {
        fprintf(stderr, "usage: %s [records] [iterations] [removals]\n", argv[0]);

        return 2;
    }

    char index_path[160];
    char compacted_path[160];
    char clean_rewrite_path[160];
    char manifest_path[160];
    char mutation_path[192];

    snprintf(index_path, sizeof(index_path), "/tmp/geobolt-tombstone-index-%ld.bin", (long) getpid());
    snprintf(compacted_path, sizeof(compacted_path), "/tmp/geobolt-tombstone-compacted-%ld.bin", (long) getpid());
    snprintf(clean_rewrite_path, sizeof(clean_rewrite_path), "/tmp/geobolt-tombstone-clean-%ld.bin", (long) getpid());
    snprintf(manifest_path, sizeof(manifest_path), "/tmp/geobolt-tombstone-manifest-%ld.bin", (long) getpid());
    snprintf(mutation_path, sizeof(mutation_path), "%s.mutations", manifest_path);
    remove(index_path);
    remove(compacted_path);
    remove(clean_rewrite_path);
    remove(manifest_path);
    remove(mutation_path);

    GeoIndex *index = geo_index_create(record_count);
    uint64_t *removed_ids = removal_count ? malloc(removal_count * sizeof(*removed_ids)) : NULL;
    bool succeeded = index && removed_ids;

    for (size_t i = 0; succeeded && i < removal_count; ++i) {
        removed_ids[i] = UINT64_C(1000000) + i * record_count / removal_count;
    }

    for (size_t i = 0; succeeded && i < record_count; ++i) {
        double latitude = (double) (i % 180000U) * 0.001 - 90.0;
        double longitude = (double) ((i * 17U) % 360000U) * 0.001 - 180.0;

        succeeded = geo_index_add(index, UINT64_C(1000000) + i, latitude, longitude);
    }

    if (succeeded) {
        succeeded = geo_index_build(index) && geo_index_save(index, index_path);
    }

    geo_index_destroy(index);

    GeoSegmentSet *set = succeeded ? geo_segment_set_create(manifest_path) : NULL;

    if (succeeded) {
        succeeded = set && geo_segment_set_add_file(set, index_path);
    }

    QueryMeasurement baseline;
    QueryMeasurement tombstone;
    QueryMeasurement compacted;
    GeoSegmentCompactionStats compaction_stats;
    GeoSegmentCompactionStats clean_rewrite_stats;
    double compaction_wall_ms = 0.0;
    double clean_rewrite_wall_ms = 0.0;

    if (succeeded) {
        succeeded = measure_global_counts(set, iterations, &baseline) &&
                    geo_segment_set_remove_ids(set, removed_ids, removal_count) &&
                    measure_global_counts(set, iterations, &tombstone);
    }

    if (succeeded) {
        double compaction_start = geo_get_time_ms();

        succeeded = geo_segment_set_compact(set, compacted_path, &compaction_stats);
        compaction_wall_ms = geo_get_time_ms() - compaction_start;
    }

    if (succeeded) {
        succeeded = measure_global_counts(set, iterations, &compacted);
    }

    if (succeeded) {
        double clean_rewrite_start = geo_get_time_ms();

        succeeded = geo_segment_set_compact(set, clean_rewrite_path, &clean_rewrite_stats);
        clean_rewrite_wall_ms = geo_get_time_ms() - clean_rewrite_start;
    }

    if (succeeded) {
        geo_benchmark_print_environment();
        printf("records=%zu iterations=%zu removals=%zu\n", record_count, iterations, removal_count);
        printf("baseline_ms=%.3f checksum=%" PRIu64 " scanned=%" PRIu64 "\n",
               baseline.elapsed_ms,
               baseline.result_checksum,
               baseline.records_scanned);
        printf("tombstone_ms=%.3f checksum=%" PRIu64 " scanned=%" PRIu64 "\n",
               tombstone.elapsed_ms,
               tombstone.result_checksum,
               tombstone.records_scanned);
        printf("compacted_ms=%.3f checksum=%" PRIu64 " scanned=%" PRIu64 "\n",
               compacted.elapsed_ms,
               compacted.result_checksum,
               compacted.records_scanned);
        printf("compaction_ms=%.3f compaction_wall_ms=%.3f input_segments=%zu records_written=%" PRIu64 "\n",
               compaction_stats.merge_time_ms,
               compaction_wall_ms,
               compaction_stats.input_segments,
               compaction_stats.records_written);
        printf("clean_rewrite_ms=%.3f clean_rewrite_wall_ms=%.3f records_written=%" PRIu64 "\n",
               clean_rewrite_stats.merge_time_ms,
               clean_rewrite_wall_ms,
               clean_rewrite_stats.records_written);
    }

    geo_segment_set_destroy(set);
    free(removed_ids);
    remove(index_path);
    remove(compacted_path);
    remove(clean_rewrite_path);
    remove(manifest_path);
    remove(mutation_path);

    return succeeded ? 0 : 1;
}
