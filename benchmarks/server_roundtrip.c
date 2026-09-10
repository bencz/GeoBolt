#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "benchmark_common.h"
#include "geobolt/client.h"
#include "geobolt/geodoc.h"
#include "geobolt/server.h"
#include "test_support.h"

#include <pthread.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef enum {
    SERVER_BENCHMARK_PING,
    SERVER_BENCHMARK_RADIUS_COUNT,
    SERVER_BENCHMARK_UPSERT,
} ServerBenchmarkOperation;

typedef struct {
    GeoServer *server;
    GeoServerStatus status;
} ServerThreadContext;

typedef struct {
    uint16_t port;
    size_t operation_count;
    uint64_t first_object_id;
    uint64_t random_seed;
    uint64_t checksum;
    ServerBenchmarkOperation operation;
    atomic_size_t *ready;
    atomic_bool *start;
    bool succeeded;
    const GeoDocBuffer *documents;
} ClientThreadContext;

static const struct {
    const char *name;
    const char *pointer;
    GeoDatabaseIndexType type;
} metadata_indexes[] = {
    { "available_idx", "/available", GEO_DATABASE_INDEX_BOOL },
    { "score_idx", "/score", GEO_DATABASE_INDEX_INT64 },
    { "rank_idx", "/rank", GEO_DATABASE_INDEX_UINT64 },
    { "price_idx", "/price", GEO_DATABASE_INDEX_DOUBLE },
    { "updated_idx", "/updated", GEO_DATABASE_INDEX_DATETIME },
    { "label_idx", "/label", GEO_DATABASE_INDEX_STRING },
    { "token_idx", "/token", GEO_DATABASE_INDEX_BYTES },
};

static bool build_metadata_documents(GeoDocBuffer documents[2])
{
    unsigned char extra[1024];
    memset(extra, 0x5a, sizeof(extra));

    for (size_t version = 0U; version < 2U; ++version) {
        GeoDocBuilder *builder = geo_doc_builder_create();
        unsigned char token[2] = { 0U, (unsigned char) version };
        bool succeeded = builder &&
                         geo_doc_builder_add_bool(builder, 0U, "available", version != 0U, NULL) == GEO_DOC_OK &&
                         geo_doc_builder_add_int64(builder, 0U, "score", (int64_t) version, NULL) == GEO_DOC_OK &&
                         geo_doc_builder_add_uint64(builder, 0U, "rank", version, NULL) == GEO_DOC_OK &&
                         geo_doc_builder_add_double(builder, 0U, "price", (double) version, NULL) == GEO_DOC_OK &&
                         geo_doc_builder_add_int64(builder, 0U, "updated", (int64_t) version, NULL) == GEO_DOC_OK &&
                         geo_doc_builder_add_string(builder, 0U, "label", version ? "new" : "old", 3U, NULL) == GEO_DOC_OK &&
                         geo_doc_builder_add_bytes(builder, 0U, "token", token, sizeof(token), NULL) == GEO_DOC_OK &&
                         geo_doc_builder_add_bytes(builder, 0U, "extra", extra, sizeof(extra), NULL) == GEO_DOC_OK &&
                         geo_doc_builder_finish(builder, documents + version) == GEO_DOC_OK;
        geo_doc_builder_destroy(builder);

        if (!succeeded) {
            return false;
        }
    }

    return true;
}

static bool verify_metadata_indexes(GeoClient *client, uint64_t expected, uint64_t *checksum)
{
    GeoDatabaseIndexValue values[] = {
        { .type = GEO_DATABASE_INDEX_BOOL, .as.boolean = true },
        { .type = GEO_DATABASE_INDEX_INT64, .as.signed_integer = 1 },
        { .type = GEO_DATABASE_INDEX_UINT64, .as.unsigned_integer = 1U },
        { .type = GEO_DATABASE_INDEX_DOUBLE, .as.floating_point = 1.0 },
        { .type = GEO_DATABASE_INDEX_DATETIME, .as.datetime = 1 },
        { .type = GEO_DATABASE_INDEX_STRING, .as.bytes = { .data = "new", .size = 3U } },
        { .type = GEO_DATABASE_INDEX_BYTES, .as.bytes = { .data = "\0\1", .size = 2U } },
    };
    GeoIdResult result = { 0 };
    bool succeeded = true;

    for (size_t index = 0U; succeeded && index < sizeof(values) / sizeof(*values); ++index) {
        GeoDatabaseIndexPredicate predicate = {
            .index_name = metadata_indexes[index].name,
            .operation = GEO_DATABASE_INDEX_EQUAL,
            .lower = values[index],
        };
        succeeded = geo_client_query_index_reuse(client, &predicate, &result, NULL) == GEO_CLIENT_OK && result.count == expected;

        for (size_t record = 0U; succeeded && record < result.count; ++record) {
            *checksum = *checksum * UINT64_C(0x9e3779b97f4a7c15) + result.ids[record];
        }
    }

    free(result.ids);
    return succeeded;
}

static double elapsed_seconds(const struct timespec *begin, const struct timespec *end)
{
    return (double) (end->tv_sec - begin->tv_sec) + (double) (end->tv_nsec - begin->tv_nsec) / 1000000000.0;
}

static uint64_t random_u64(uint64_t *state)
{
    uint64_t value = (*state += UINT64_C(0x9e3779b97f4a7c15));

    value = (value ^ (value >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27U)) * UINT64_C(0x94d049bb133111eb);

    return value ^ (value >> 31U);
}

static double random_unit(uint64_t *state)
{
    return (double) (random_u64(state) >> 11U) * (1.0 / 9007199254740992.0);
}

static void *server_thread_main(void *argument)
{
    ServerThreadContext *context = argument;

    context->status = geo_server_run(context->server);

    return NULL;
}

static void *client_thread_main(void *argument)
{
    ClientThreadContext *context = argument;
    GeoClientConfig config = geo_client_default_config("127.0.0.1", context->port);
    GeoClientStatus status;
    GeoClient *client = geo_client_connect(&config, &status);

    uint64_t warmup_count;

    context->succeeded = client && status == GEO_CLIENT_OK &&
                         geo_client_ping(client, NULL, 0U) == GEO_CLIENT_OK &&
                         geo_client_radius_count(client, 0.0, 0.0, 25.0, &warmup_count) == GEO_CLIENT_OK;

    for (size_t warmup = 0U; context->succeeded && context->operation == SERVER_BENCHMARK_UPSERT && warmup < 32U; ++warmup) {
        const GeoDocBuffer *document = context->documents ? context->documents + warmup % 2U : NULL;

        context->succeeded = geo_client_upsert(client, context->first_object_id, 0.0, 0.0,
                                               document ? document->data : NULL,
                                               document ? document->size : 0U) == GEO_CLIENT_OK;
    }

    atomic_fetch_add_explicit(context->ready, 1U, memory_order_release);

    while (!atomic_load_explicit(context->start, memory_order_acquire)) {
    }

    uint64_t random_state = context->random_seed;

    for (size_t operation = 0; context->succeeded && operation < context->operation_count; ++operation) {
        if (context->operation == SERVER_BENCHMARK_PING) {
            context->succeeded = geo_client_ping(client, NULL, 0) == GEO_CLIENT_OK;
            continue;
        }

        double latitude = random_unit(&random_state) * 170.0 - 85.0;
        double longitude = random_unit(&random_state) * 360.0 - 180.0;

        if (context->operation == SERVER_BENCHMARK_UPSERT) {
            uint64_t object_id = context->first_object_id + (context->documents ? operation / 2U : operation);
            const GeoDocBuffer *document = context->documents ? context->documents + operation % 2U : NULL;

            context->succeeded = geo_client_upsert(client, object_id, latitude, longitude,
                                                   document ? document->data : NULL,
                                                   document ? document->size : 0U) == GEO_CLIENT_OK;
            continue;
        }

        uint64_t count;

        context->succeeded = geo_client_radius_count(client, latitude, longitude, 25.0, &count) == GEO_CLIENT_OK;

        if (context->succeeded) {
            context->checksum = context->checksum * UINT64_C(0x9e3779b97f4a7c15) + count;
        }
    }

    geo_client_close(client);

    return NULL;
}

static bool run_measurement(uint16_t port,
                            size_t client_count,
                            size_t operations_per_client,
                            ServerBenchmarkOperation operation,
                            uint64_t first_object_id,
                            const GeoDocBuffer *documents,
                            double *seconds,
                            uint64_t *checksum)
{
    pthread_t *threads = calloc(client_count, sizeof(*threads));
    ClientThreadContext *contexts = calloc(client_count, sizeof(*contexts));
    atomic_size_t ready = ATOMIC_VAR_INIT(0U);
    atomic_bool start = ATOMIC_VAR_INIT(false);
    size_t threads_started = 0;
    bool succeeded = threads && contexts;

    for (size_t client = 0; succeeded && client < client_count; ++client) {
        contexts[client] = (ClientThreadContext) {
            .port = port,
            .operation_count = operations_per_client,
            .first_object_id = first_object_id + client * operations_per_client,
            .documents = documents,
            .random_seed = UINT64_C(0x243f6a8885a308d3) + client,
            .operation = operation,
            .ready = &ready,
            .start = &start,
        };

        if (pthread_create(threads + client, NULL, client_thread_main, contexts + client) != 0) {
            succeeded = false;
            break;
        }

        threads_started++;
    }

    while (succeeded && atomic_load_explicit(&ready, memory_order_acquire) != client_count) {
    }

    struct timespec begin;
    struct timespec end;

    clock_gettime(CLOCK_MONOTONIC, &begin);
    atomic_store_explicit(&start, true, memory_order_release);

    for (size_t client = 0; client < threads_started; ++client) {
        if (pthread_join(threads[client], NULL) != 0 || !contexts[client].succeeded) {
            succeeded = false;
        }

        *checksum ^= contexts[client].checksum;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    *seconds = elapsed_seconds(&begin, &end);
    free(contexts);
    free(threads);

    return succeeded;
}

int main(int argc, char **argv)
{
    size_t client_count = 4U;
    size_t operations_per_client = 1000U;
    size_t indexed_objects = 10000U;
    bool metadata = argc == 5 && strcmp(argv[4], "metadata") == 0;

    if ((argc > 1 && !geo_benchmark_parse_size(argv[1], &client_count)) ||
        (argc > 2 && !geo_benchmark_parse_size(argv[2], &operations_per_client)) ||
        (argc > 3 && !geo_benchmark_parse_size(argv[3], &indexed_objects)) ||
        argc > 5 || (argc == 5 && !metadata) || (metadata && operations_per_client % 2U) ||
        !client_count || client_count > SIZE_MAX / 4U || !operations_per_client ||
        !indexed_objects || indexed_objects > SIZE_MAX / sizeof(GeoDatabaseMutation) ||
        operations_per_client > (UINT64_MAX - indexed_objects) / client_count) {
        fprintf(stderr, "usage: %s [clients] [operations-per-client] [indexed-objects] [metadata]\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *directory_base = getenv("GEOBOLT_BENCHMARK_DIRECTORY");

    if (!directory_base || !directory_base[0]) {
        directory_base = "/tmp";
    }

    char directory_template[4096];
    int path_length = snprintf(directory_template, sizeof(directory_template), "%s/geobolt-server-benchmark-XXXXXX", directory_base);
    char *directory = path_length > 0 && (size_t) path_length < sizeof(directory_template) ? mkdtemp(directory_template) : NULL;
    bool succeeded = directory != NULL;
    const char *failure_stage = "temporary directory creation";
    GeoServer *server = NULL;
    GeoClient *setup_client = NULL;
    pthread_t server_thread;
    bool server_thread_started = false;
    ServerThreadContext server_context = { 0 };
    GeoDocBuffer documents[2] = { 0 };

    if (succeeded && metadata) {
        failure_stage = "metadata document construction";
        succeeded = build_metadata_documents(documents);
    }

    if (succeeded) {
        failure_stage = "server creation";
        GeoServerConfig server_config = geo_server_default_config(directory);

        server_config.port = 0;
        server_config.worker_threads = client_count < 16U ? client_count : 16U;
        server_config.work_queue_capacity = client_count * 4U;
        GeoServerStatus status;

        server = geo_server_create(&server_config, &status);
        succeeded = server && status == GEO_SERVER_OK;
    }

    if (succeeded) {
        failure_stage = "server thread startup";
        server_context.server = server;
        succeeded = pthread_create(&server_thread, NULL, server_thread_main, &server_context) == 0;
        server_thread_started = succeeded;
    }

    GeoDatabaseMutation *mutations = succeeded ? malloc(indexed_objects * sizeof(*mutations)) : NULL;

    if (succeeded && !mutations) {
        failure_stage = "setup mutation allocation";
        succeeded = false;
    }

    uint64_t random_state = UINT64_C(0x13198a2e03707344);

    for (size_t object = 0; succeeded && object < indexed_objects; ++object) {
        double latitude = random_unit(&random_state) * 180.0 - 90.0;
        double longitude = random_unit(&random_state) * 360.0 - 180.0;

        mutations[object] = (GeoDatabaseMutation) {
            .object_id = object + 1U,
            .morton_code = geo_encode(latitude, longitude),
            .operation = GEO_DATABASE_UPSERT,
        };
    }

    if (succeeded) {
        failure_stage = "initial bulk load";
        GeoClientConfig config = geo_client_default_config("127.0.0.1", geo_server_port(server));
        GeoClientStatus status;

        setup_client = geo_client_connect(&config, &status);
        succeeded = setup_client && status == GEO_CLIENT_OK &&
                    geo_client_write(setup_client, mutations, indexed_objects) == GEO_CLIENT_OK;
    }

    free(mutations);

    for (size_t index = 0U; succeeded && metadata && index < sizeof(metadata_indexes) / sizeof(*metadata_indexes); ++index) {
        failure_stage = "typed index creation";
        succeeded = geo_client_create_index(setup_client, metadata_indexes[index].name,
                                            metadata_indexes[index].pointer, metadata_indexes[index].type) == GEO_CLIENT_OK;
    }

    geo_client_close(setup_client);

    double ping_seconds = 0.0;
    double query_seconds = 0.0;
    double upsert_seconds = 0.0;
    uint64_t checksum = 0U;

    if (succeeded) {
        failure_stage = "ping measurement";
        succeeded = run_measurement(geo_server_port(server),
                                    client_count,
                                    operations_per_client,
                                    SERVER_BENCHMARK_PING,
                                    0U,
                                    NULL,
                                    &ping_seconds,
                                    &checksum);
    }

    if (succeeded) {
        failure_stage = "radius-count measurement";
        succeeded = run_measurement(geo_server_port(server),
                                    client_count,
                                    operations_per_client,
                                    SERVER_BENCHMARK_RADIUS_COUNT,
                                    0U,
                                    NULL,
                                    &query_seconds,
                                    &checksum);
    }

    if (succeeded) {
        failure_stage = "durable-upsert measurement";
        succeeded = run_measurement(geo_server_port(server),
                                    client_count,
                                    operations_per_client,
                                    SERVER_BENCHMARK_UPSERT,
                                    indexed_objects + 1U,
                                    metadata ? documents : NULL,
                                    &upsert_seconds,
                                    &checksum);
    }

    uint64_t visible_objects = 0U;
    uint64_t metadata_checksum = 0U;

    if (succeeded) {
        failure_stage = "final canonical visibility validation";
        GeoClientConfig config = geo_client_default_config("127.0.0.1", geo_server_port(server));
        GeoClientStatus status;
        GeoClient *verifier = geo_client_connect(&config, &status);

        succeeded = verifier && status == GEO_CLIENT_OK &&
                    geo_client_radius_count(verifier, 0.0, 0.0, 25000.0, &visible_objects) == GEO_CLIENT_OK &&
                    visible_objects == indexed_objects + (uint64_t) client_count * (operations_per_client / (metadata ? 2U : 1U));

        if (succeeded && metadata) {
            failure_stage = "typed index update verification";
            succeeded = verify_metadata_indexes(verifier, (uint64_t) client_count * (operations_per_client / 2U), &metadata_checksum);
        }
        geo_client_close(verifier);
    }

    GeoServerStats stats = { 0 };

    if (server) {
        (void) geo_server_get_stats(server, &stats);
        geo_server_request_stop(server);
    }

    if (server_thread_started) {
        bool shutdown_succeeded = pthread_join(server_thread, NULL) == 0 && server_context.status == GEO_SERVER_OK;

        if (succeeded && !shutdown_succeeded) {
            failure_stage = "server shutdown";
        }
        succeeded = shutdown_succeeded && succeeded;
    }

    geo_server_destroy(server);
    geo_doc_buffer_release(documents);
    geo_doc_buffer_release(documents + 1U);

    if (directory && !geobolt_test_remove_tree(directory)) {
        if (succeeded) {
            failure_stage = "temporary directory cleanup";
        }
        succeeded = false;
    }

    if (!succeeded) {
        fprintf(stderr, "server benchmark failed during %s\n", failure_stage);
        return EXIT_FAILURE;
    }

    double total_operations = (double) client_count * (double) operations_per_client;

    geo_benchmark_print_environment();
    printf("data_directory_base=%s\n", directory_base);
    printf("clients=%zu operations_per_client=%zu indexed_objects=%zu\n",
           client_count,
           operations_per_client,
           indexed_objects);
    printf("query_checksum=%llu visible_objects=%llu\n", (unsigned long long) checksum, (unsigned long long) visible_objects);
    printf("metadata=%s metadata_checksum=%llu\n", metadata ? "seven-typed-indexes" : "none",
           (unsigned long long) metadata_checksum);
    printf("ping_ops_per_second=%.0f ping_aggregate_us_per_op=%.3f\n",
           total_operations / ping_seconds,
           ping_seconds * 1000000.0 / total_operations);
    printf("radius_count_ops_per_second=%.0f radius_count_aggregate_us_per_op=%.3f\n",
           total_operations / query_seconds,
           query_seconds * 1000000.0 / total_operations);
    printf("durable_upsert_ops_per_second=%.0f durable_upsert_aggregate_us_per_op=%.3f\n",
           total_operations / upsert_seconds,
           upsert_seconds * 1000000.0 / total_operations);
    printf("completed_requests=%llu backpressure_rejections=%llu peak_payload_bytes=%llu\n",
           (unsigned long long) stats.completed_requests,
           (unsigned long long) stats.backpressure_rejections,
           (unsigned long long) stats.peak_inflight_payload_bytes);

    return EXIT_SUCCESS;
}
