#include "benchmark_common.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char **argv)
{
    size_t record_count = 3000000;
    size_t chunk_capacity = 100000;
    size_t sort_threads = 4;
    size_t rounds = 5;
    const char *distribution = "uniform";

    if ((argc > 1 && !geo_benchmark_parse_size(argv[1], &record_count)) ||
        (argc > 2 && !geo_benchmark_parse_size(argv[2], &chunk_capacity)) ||
        (argc > 3 && !geo_benchmark_parse_size(argv[3], &sort_threads)) ||
        (argc > 4 && !geo_benchmark_parse_size(argv[4], &rounds)) ||
        argc > 6 ||
        !record_count ||
        !chunk_capacity ||
        !rounds ||
        record_count > SIZE_MAX / sizeof(GeoRecord)) {
        fprintf(stderr, "usage: %s [records] [chunk-records] [sort-threads] [rounds] [uniform|hotspot]\n", argv[0]);

        return 2;
    }

    if (argc > 5) {
        distribution = argv[5];

        if (strcmp(distribution, "uniform") != 0 && strcmp(distribution, "hotspot") != 0) {
            fprintf(stderr, "distribution must be uniform or hotspot\n");

            return 2;
        }
    }

    GeoRecord *records = malloc(record_count * sizeof(*records));

    if (!records) {
        return 1;
    }

    uint64_t state = UINT64_C(0x243f6a8885a308d3);
    uint64_t input_checksum = 0;
    bool hotspot = strcmp(distribution, "hotspot") == 0;

    for (size_t i = 0; i < record_count; ++i) {
        state = state * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);

        uint64_t morton = state;

        if (hotspot && i % 5U != 0) {
            morton = UINT64_C(0x6a5c000000000000) | (state & UINT64_C(0x0000ffffffffffff));
        }

        records[i] = (GeoRecord) {
            .id = i,
            .z = morton,
        };
        input_checksum ^= records[i].z + records[i].id * UINT64_C(0x9e3779b97f4a7c15);
    }

    geo_benchmark_print_environment();
    printf("records=%zu chunk_records=%zu sort_threads=%zu rounds=%zu distribution=%s input_checksum=%" PRIu64 "\n",
           record_count,
           chunk_capacity,
           sort_threads,
           rounds,
           distribution,
           input_checksum);

    bool succeeded = true;

    for (size_t round = 0; succeeded && round < rounds; ++round) {
        char path[192];

        snprintf(path, sizeof(path), "/tmp/geobolt-stream-benchmark-%ld-%zu.bin", (long) getpid(), round);
        remove(path);

        double start = geo_get_time_ms();
        GeoStreamBuilder *builder = geo_stream_builder_create_parallel(path, "/tmp", chunk_capacity, sort_threads);
        GeoStreamBuildStats stats;

        succeeded = builder &&
                    geo_stream_builder_add_records(builder, records, record_count) &&
                    geo_stream_builder_finish(builder, &stats);

        double total_ms = geo_get_time_ms() - start;
        GeoIndex *mapped = succeeded ? geo_index_open_mmap(path) : NULL;
        size_t mapped_count = 0;

        succeeded = mapped &&
                    geo_search_radius_count(mapped, 0.0, 0.0, 25000.0, &mapped_count, NULL) &&
                    mapped_count == record_count;

        if (succeeded) {
            printf("round=%zu total_ms=%.3f chunk_ms=%.3f intermediate_merge_ms=%.3f final_merge_ms=%.3f "
                   "runs=%zu intermediate_merges=%zu\n",
                   round,
                   total_ms,
                   stats.chunk_build_time_ms,
                   stats.intermediate_merge_time_ms,
                   stats.merge_time_ms,
                   stats.runs_created,
                   stats.intermediate_merges);
        }

        geo_index_destroy(mapped);
        geo_stream_builder_destroy(builder);
        remove(path);
    }

    free(records);

    return succeeded ? 0 : 1;
}
