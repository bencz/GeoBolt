#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "benchmark_common.h"
#include "geobolt/client.h"
#include "geobolt/server.h"
#include "test_support.h"

#include <pthread.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

typedef enum {
    SERVER_BENCHMARK_PING,
    SERVER_BENCHMARK_RADIUS_COUNT,
} ServerBenchmarkOperation;

typedef struct {
    GeoServer *server;
    GeoServerStatus status;
} ServerThreadContext;

typedef struct {
    uint16_t port;
    size_t operation_count;
    ServerBenchmarkOperation operation;
    atomic_size_t *ready;
    atomic_bool *start;
    bool succeeded;
} ClientThreadContext;

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

    context->succeeded = client && status == GEO_CLIENT_OK;
    atomic_fetch_add_explicit(context->ready, 1U, memory_order_release);

    while (!atomic_load_explicit(context->start, memory_order_acquire)) {
    }

    uint64_t random_state = UINT64_C(0x243f6a8885a308d3) ^ (uintptr_t) context;

    for (size_t operation = 0; context->succeeded && operation < context->operation_count; ++operation) {
        if (context->operation == SERVER_BENCHMARK_PING) {
            context->succeeded = geo_client_ping(client, NULL, 0) == GEO_CLIENT_OK;
            continue;
        }

        double latitude = random_unit(&random_state) * 170.0 - 85.0;
        double longitude = random_unit(&random_state) * 360.0 - 180.0;
        uint64_t count;

        context->succeeded = geo_client_radius_count(client, latitude, longitude, 25.0, &count) == GEO_CLIENT_OK;
    }

    geo_client_close(client);

    return NULL;
}

static bool run_measurement(uint16_t port,
                            size_t client_count,
                            size_t operations_per_client,
                            ServerBenchmarkOperation operation,
                            double *seconds)
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

    if ((argc > 1 && !geo_benchmark_parse_size(argv[1], &client_count)) ||
        (argc > 2 && !geo_benchmark_parse_size(argv[2], &operations_per_client)) ||
        (argc > 3 && !geo_benchmark_parse_size(argv[3], &indexed_objects)) ||
        argc > 4 || !client_count || client_count > SIZE_MAX / 4U || !operations_per_client ||
        !indexed_objects || indexed_objects > SIZE_MAX / sizeof(GeoDatabaseMutation)) {
        fprintf(stderr, "usage: %s [clients] [operations-per-client] [indexed-objects]\n", argv[0]);
        return EXIT_FAILURE;
    }

    char directory_template[] = "/tmp/geobolt-server-benchmark-XXXXXX";
    char *directory = mkdtemp(directory_template);
    bool succeeded = directory != NULL;
    GeoServer *server = NULL;
    GeoClient *setup_client = NULL;
    pthread_t server_thread;
    bool server_thread_started = false;
    ServerThreadContext server_context = { 0 };

    if (succeeded) {
        GeoServerConfig server_config = geo_server_default_config(directory);

        server_config.port = 0;
        server_config.worker_threads = client_count < 16U ? client_count : 16U;
        server_config.work_queue_capacity = client_count * 4U;
        GeoServerStatus status;

        server = geo_server_create(&server_config, &status);
        succeeded = server && status == GEO_SERVER_OK;
    }

    if (succeeded) {
        server_context.server = server;
        succeeded = pthread_create(&server_thread, NULL, server_thread_main, &server_context) == 0;
        server_thread_started = succeeded;
    }

    GeoDatabaseMutation *mutations = succeeded ? malloc(indexed_objects * sizeof(*mutations)) : NULL;

    if (succeeded && !mutations) {
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
        GeoClientConfig config = geo_client_default_config("127.0.0.1", geo_server_port(server));
        GeoClientStatus status;

        setup_client = geo_client_connect(&config, &status);
        succeeded = setup_client && status == GEO_CLIENT_OK &&
                    geo_client_write(setup_client, mutations, indexed_objects) == GEO_CLIENT_OK;
    }

    free(mutations);
    geo_client_close(setup_client);

    double ping_seconds = 0.0;
    double query_seconds = 0.0;

    if (succeeded) {
        succeeded = run_measurement(geo_server_port(server),
                                    client_count,
                                    operations_per_client,
                                    SERVER_BENCHMARK_PING,
                                    &ping_seconds) &&
                    run_measurement(geo_server_port(server),
                                    client_count,
                                    operations_per_client,
                                    SERVER_BENCHMARK_RADIUS_COUNT,
                                    &query_seconds);
    }

    GeoServerStats stats = { 0 };

    if (server) {
        (void) geo_server_get_stats(server, &stats);
        geo_server_request_stop(server);
    }

    if (server_thread_started) {
        succeeded = pthread_join(server_thread, NULL) == 0 && server_context.status == GEO_SERVER_OK && succeeded;
    }

    geo_server_destroy(server);

    if (directory && !geobolt_test_remove_tree(directory)) {
        succeeded = false;
    }

    if (!succeeded) {
        fprintf(stderr, "server benchmark failed\n");
        return EXIT_FAILURE;
    }

    double total_operations = (double) client_count * (double) operations_per_client;

    geo_benchmark_print_environment();
    printf("clients=%zu operations_per_client=%zu indexed_objects=%zu\n",
           client_count,
           operations_per_client,
           indexed_objects);
    printf("ping_ops_per_second=%.0f ping_mean_us=%.3f\n",
           total_operations / ping_seconds,
           ping_seconds * 1000000.0 / total_operations);
    printf("radius_count_ops_per_second=%.0f radius_count_mean_us=%.3f\n",
           total_operations / query_seconds,
           query_seconds * 1000000.0 / total_operations);
    printf("completed_requests=%llu backpressure_rejections=%llu peak_payload_bytes=%llu\n",
           (unsigned long long) stats.completed_requests,
           (unsigned long long) stats.backpressure_rejections,
           (unsigned long long) stats.peak_inflight_payload_bytes);

    return EXIT_SUCCESS;
}
