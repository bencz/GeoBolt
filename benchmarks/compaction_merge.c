#include "benchmark_common.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define COMPACTION_BENCHMARK_BATCH_RECORDS 8192U
#define COMPACTION_BENCHMARK_PATH_BYTES 192U

typedef struct {
    double milliseconds;
    size_t workers;
    size_t partitions;
    uint64_t records;
} CompactionMeasurement;

static uint64_t compaction_random_u64(uint64_t *state)
{
    uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));

    value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);

    return value ^ (value >> 31U);
}

static double compaction_random_unit(uint64_t *state)
{
    return (double) (compaction_random_u64(state) >> 11U) * (1.0 / 9007199254740992.0);
}

static bool compaction_write_source(const char *path,
                                    size_t segment_index,
                                    size_t record_count,
                                    size_t build_workers,
                                    bool hotspot)
{
    GeoIndex *index = geo_index_create(record_count);
    GeoRecord *batch = malloc(COMPACTION_BENCHMARK_BATCH_RECORDS * sizeof(*batch));
    uint64_t random_state = UINT64_C(0x6a09e667f3bcc909) ^ segment_index;
    size_t position = 0;
    bool succeeded = index && batch;

    while (succeeded && position < record_count) {
        size_t count = record_count - position;

        if (count > COMPACTION_BENCHMARK_BATCH_RECORDS) {
            count = COMPACTION_BENCHMARK_BATCH_RECORDS;
        }

        for (size_t i = 0; i < count; ++i) {
            bool hot_record = hotspot && compaction_random_unit(&random_state) < 0.8;
            double latitude = hot_record
                                  ? -23.5505 + (compaction_random_unit(&random_state) - 0.5) * 0.04
                                  : compaction_random_unit(&random_state) * 180.0 - 90.0;
            double longitude = hot_record
                                   ? -46.6333 + (compaction_random_unit(&random_state) - 0.5) * 0.04
                                   : compaction_random_unit(&random_state) * 360.0 - 180.0;

            batch[i] = (GeoRecord) {
                .id = UINT64_C(1000000000) + segment_index * record_count + position + i,
                .z = geo_encode(latitude, longitude),
            };
        }

        succeeded = geo_index_add_records(index, batch, count);
        position += succeeded ? count : 0;
    }

    if (succeeded) {
        succeeded = geo_index_build_parallel(index, build_workers) && geo_index_save(index, path);
    }

    free(batch);
    geo_index_destroy(index);

    return succeeded;
}

static bool compaction_measure(char paths[][COMPACTION_BENCHMARK_PATH_BYTES],
                               size_t segment_count,
                               size_t requested_workers,
                               size_t sequence,
                               CompactionMeasurement *measurement)
{
    char manifest_path[COMPACTION_BENCHMARK_PATH_BYTES];
    char mutation_path[COMPACTION_BENCHMARK_PATH_BYTES + 32U];
    char checkpoint_path[COMPACTION_BENCHMARK_PATH_BYTES + 48U];
    char output_path[COMPACTION_BENCHMARK_PATH_BYTES];
    long process_id = (long) getpid();

    snprintf(manifest_path,
             sizeof(manifest_path),
             "/tmp/geobolt-compaction-bench-manifest-%ld-%zu.bin",
             process_id,
             sequence);
    snprintf(mutation_path, sizeof(mutation_path), "%s.mutations", manifest_path);
    snprintf(checkpoint_path, sizeof(checkpoint_path), "%s.mutations.checkpoint", manifest_path);
    snprintf(output_path,
             sizeof(output_path),
             "/tmp/geobolt-compaction-bench-output-%ld-%zu.bin",
             process_id,
             sequence);
    remove(checkpoint_path);
    remove(mutation_path);
    remove(manifest_path);
    remove(output_path);

    GeoSegmentSet *set = geo_segment_set_create(manifest_path);
    bool succeeded = set != NULL;

    for (size_t segment = 0; succeeded && segment < segment_count; ++segment) {
        succeeded = geo_segment_set_add_file(set, paths[segment]);
    }

    GeoSegmentCompactionStats stats = { 0 };

    if (succeeded) {
        succeeded = geo_segment_set_compact_with_workers(set, output_path, requested_workers, &stats);
    }

    size_t visible_count = 0;

    if (succeeded) {
        succeeded = geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &visible_count, NULL) &&
                    visible_count == stats.records_written;
    }

    GeoIndex *persisted = succeeded ? geo_index_open_mmap(output_path) : NULL;

    succeeded = succeeded && persisted;
    geo_index_destroy(persisted);

    if (succeeded) {
        measurement->milliseconds = stats.merge_time_ms;
        measurement->workers = stats.worker_count;
        measurement->partitions = stats.partition_count;
        measurement->records = stats.records_written;
    }

    geo_segment_set_destroy(set);
    remove(checkpoint_path);
    remove(mutation_path);
    remove(manifest_path);
    remove(output_path);

    return succeeded;
}

int main(int argc, char **argv)
{
    size_t segment_count = 8U;
    size_t records_per_segment = 1000000U;
    size_t parallel_workers = 8U;
    size_t rounds = 3U;
    const char *distribution = "uniform";

    if ((argc > 1 && !geo_benchmark_parse_size(argv[1], &segment_count)) ||
        (argc > 2 && !geo_benchmark_parse_size(argv[2], &records_per_segment)) ||
        (argc > 3 && !geo_benchmark_parse_size(argv[3], &parallel_workers)) ||
        (argc > 4 && !geo_benchmark_parse_size(argv[4], &rounds)) ||
        argc > 6 ||
        segment_count < 2U ||
        !records_per_segment ||
        parallel_workers < 2U ||
        parallel_workers > 16U ||
        !rounds ||
        segment_count > SIZE_MAX / sizeof(char[COMPACTION_BENCHMARK_PATH_BYTES]) ||
        records_per_segment > SIZE_MAX / segment_count ||
        records_per_segment > (UINT64_MAX - UINT64_C(1000000000)) / segment_count) {
        fprintf(stderr,
                "usage: %s [segments>=2] [records-per-segment] [workers=2..16] [rounds] [uniform|hotspot]\n",
                argv[0]);

        return 2;
    }

    if (argc > 5) {
        distribution = argv[5];
    }

    bool hotspot = strcmp(distribution, "hotspot") == 0;

    if (!hotspot && strcmp(distribution, "uniform") != 0) {
        fprintf(stderr, "distribution must be 'uniform' or 'hotspot'\n");

        return 2;
    }

    char (*paths)[COMPACTION_BENCHMARK_PATH_BYTES] = calloc(segment_count, sizeof(*paths));
    bool succeeded = paths != NULL;
    long process_id = (long) getpid();

    geo_benchmark_print_environment();
    printf("workload=parallel_morton_compaction segments=%zu records_per_segment=%zu records=%zu workers=%zu rounds=%zu"
           " distribution=%s\n",
           segment_count,
           records_per_segment,
           segment_count * records_per_segment,
           parallel_workers,
           rounds,
           distribution);

    for (size_t segment = 0; succeeded && segment < segment_count; ++segment) {
        snprintf(paths[segment],
                 sizeof(paths[segment]),
                 "/tmp/geobolt-compaction-bench-source-%ld-%zu.bin",
                 process_id,
                 segment);
        remove(paths[segment]);
        succeeded = compaction_write_source(paths[segment], segment, records_per_segment, parallel_workers, hotspot);
    }

    double serial_total_ms = 0.0;
    double parallel_total_ms = 0.0;
    uint64_t expected_records = (uint64_t) segment_count * records_per_segment;

    for (size_t round = 0; succeeded && round < rounds; ++round) {
        CompactionMeasurement first;
        CompactionMeasurement second;
        bool parallel_first = (round & 1U) != 0;
        size_t first_workers = parallel_first ? parallel_workers : 1U;
        size_t second_workers = parallel_first ? 1U : parallel_workers;

        succeeded = compaction_measure(paths, segment_count, first_workers, round * 2U, &first) &&
                    compaction_measure(paths, segment_count, second_workers, round * 2U + 1U, &second) &&
                    first.records == expected_records &&
                    second.records == expected_records;

        if (!succeeded) {
            break;
        }

        CompactionMeasurement serial = parallel_first ? second : first;
        CompactionMeasurement parallel = parallel_first ? first : second;

        serial_total_ms += serial.milliseconds;
        parallel_total_ms += parallel.milliseconds;
        printf("round=%zu serial_ms=%.3f parallel_ms=%.3f workers=%zu partitions=%zu speedup=%.3fx\n",
               round + 1U,
               serial.milliseconds,
               parallel.milliseconds,
               parallel.workers,
               parallel.partitions,
               serial.milliseconds / parallel.milliseconds);
    }

    if (succeeded) {
        double serial_average_ms = serial_total_ms / (double) rounds;
        double parallel_average_ms = parallel_total_ms / (double) rounds;

        printf("average serial_ms=%.3f parallel_ms=%.3f speedup=%.3fx status=PASS\n",
               serial_average_ms,
               parallel_average_ms,
               serial_average_ms / parallel_average_ms);
    }

    for (size_t segment = 0; segment < segment_count; ++segment) {
        remove(paths[segment]);
    }

    free(paths);

    if (!succeeded) {
        fprintf(stderr, "compaction benchmark failed\n");
    }

    return succeeded ? 0 : 1;
}
