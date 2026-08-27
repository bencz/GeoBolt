#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "geobolt/client.h"
#include "geobolt/server.h"
#include "geo_job_queue.h"
#include "geo_protocol.h"
#include "test_support.h"

#include <arpa/inet.h>
#include <pthread.h>
#include <sched.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <unistd.h>

#define SERVER_TEST_ASSERT(condition, message)                                                                             \
    do {                                                                                                                   \
        if (!(condition)) {                                                                                                \
            fprintf(stderr, "FAIL: %s (%s:%d)\n", message, __FILE__, __LINE__);                                          \
            succeeded = false;                                                                                             \
            goto cleanup;                                                                                                  \
        }                                                                                                                  \
    } while (0)

typedef struct {
    GeoServer *server;
    GeoServerStatus status;
} ServerThreadContext;

static void *server_thread_main(void *argument)
{
    ServerThreadContext *context = argument;

    context->status = geo_server_run(context->server);

    return NULL;
}

typedef struct {
    atomic_bool *release;
    atomic_uint *started;
} BlockingJobContext;

static void blocking_job(void *argument)
{
    BlockingJobContext *context = argument;

    atomic_fetch_add_explicit(context->started, 1U, memory_order_release);

    while (!atomic_load_explicit(context->release, memory_order_acquire)) {
    }
}

static bool test_bounded_job_queue_backpressure(void)
{
    bool succeeded = true;
    GeoJobQueue *queue = geo_job_queue_create(1U, 1U);
    atomic_bool release = ATOMIC_VAR_INIT(false);
    atomic_uint started = ATOMIC_VAR_INIT(0U);
    BlockingJobContext context = {
        .release = &release,
        .started = &started,
    };

    SERVER_TEST_ASSERT(queue != NULL, "bounded job queue must be created");
    SERVER_TEST_ASSERT(geo_job_queue_try_submit(queue, blocking_job, &context), "running job must be accepted");

    while (atomic_load_explicit(&started, memory_order_acquire) == 0U) {
    }

    SERVER_TEST_ASSERT(geo_job_queue_try_submit(queue, blocking_job, &context), "one queued job must fit the configured capacity");
    SERVER_TEST_ASSERT(!geo_job_queue_try_submit(queue, blocking_job, &context), "a saturated queue must reject without blocking");

cleanup:
    atomic_store_explicit(&release, true, memory_order_release);

    if (queue) {
        if (!geo_job_queue_wait_idle(queue)) {
            succeeded = false;
        }

        geo_job_queue_destroy(queue);
    }

    return succeeded;
}

typedef struct {
    uint16_t port;
    uint64_t first_id;
    size_t object_count;
    bool succeeded;
} ClientWorkerContext;

static void *client_worker_main(void *argument)
{
    ClientWorkerContext *context = argument;
    GeoClientConfig config = geo_client_default_config("127.0.0.1", context->port);

    config.authentication_token = "correct horse battery staple";
    GeoClientStatus status;
    GeoClient *client = geo_client_connect(&config, &status);
    GeoDatabaseMutation *mutations = calloc(context->object_count, sizeof(*mutations));

    context->succeeded = client && mutations && status == GEO_CLIENT_OK;

    for (size_t object = 0; context->succeeded && object < context->object_count; ++object) {
        mutations[object] = (GeoDatabaseMutation) {
            .object_id = context->first_id + object,
            .morton_code = geo_encode(-45.0 + (double) object * 0.01, -90.0 + (double) object * 0.02),
            .operation = GEO_DATABASE_UPSERT,
        };
    }

    if (context->succeeded) {
        context->succeeded = geo_client_ping(client, "concurrent", 10U) == GEO_CLIENT_OK &&
                             geo_client_write(client, mutations, context->object_count) == GEO_CLIENT_OK;
    }

    free(mutations);
    geo_client_close(client);

    return NULL;
}

typedef struct {
    GeoClient *client;
    atomic_bool *start;
    atomic_bool *failed;
    size_t worker_index;
} SharedClientWorkerContext;

static void *shared_client_worker_main(void *argument)
{
    SharedClientWorkerContext *context = argument;

    while (!atomic_load_explicit(context->start, memory_order_acquire)) {
    }

    for (size_t request = 0; request < 64U; ++request) {
        unsigned char payload[16];

        memset(payload, (int) (context->worker_index + 1U), sizeof(payload));
        geo_protocol_store_u64(payload, request);

        if (geo_client_ping(context->client, payload, sizeof(payload)) != GEO_CLIENT_OK) {
            atomic_store_explicit(context->failed, true, memory_order_release);
            break;
        }
    }

    return NULL;
}

static bool test_shared_client_serialization(GeoClient *client)
{
    enum {
        SHARED_WORKERS = 4,
    };

    pthread_t workers[SHARED_WORKERS];
    SharedClientWorkerContext contexts[SHARED_WORKERS];
    atomic_bool start = ATOMIC_VAR_INIT(false);
    atomic_bool failed = ATOMIC_VAR_INIT(false);
    size_t workers_started = 0;

    for (size_t worker = 0; worker < SHARED_WORKERS; ++worker) {
        contexts[worker] = (SharedClientWorkerContext) {
            .client = client,
            .start = &start,
            .failed = &failed,
            .worker_index = worker,
        };

        if (pthread_create(workers + worker, NULL, shared_client_worker_main, contexts + worker) != 0) {
            atomic_store_explicit(&failed, true, memory_order_release);
            break;
        }

        workers_started++;
    }

    atomic_store_explicit(&start, true, memory_order_release);

    for (size_t worker = 0; worker < workers_started; ++worker) {
        if (pthread_join(workers[worker], NULL) != 0) {
            atomic_store_explicit(&failed, true, memory_order_release);
        }
    }

    return workers_started == SHARED_WORKERS && !atomic_load_explicit(&failed, memory_order_acquire);
}

static bool send_corrupted_header(uint16_t port)
{
    int descriptor = socket(AF_INET, SOCK_STREAM | SOCK_CLOEXEC, 0);
    struct sockaddr_in address = {
        .sin_family = AF_INET,
        .sin_port = htons(port),
        .sin_addr.s_addr = htonl(INADDR_LOOPBACK),
    };

    if (descriptor < 0 || connect(descriptor, (const struct sockaddr *) &address, sizeof(address)) != 0) {
        if (descriptor >= 0) {
            (void) close(descriptor);
        }

        return false;
    }

    GeoProtocolHeader header = {
        .request_id = 99U,
        .payload_checksum = geo_protocol_checksum(NULL, 0),
        .opcode = GEO_PROTOCOL_PING,
    };
    unsigned char encoded[GEO_PROTOCOL_HEADER_SIZE];

    geo_protocol_encode_header(encoded, &header);
    encoded[40] ^= UINT8_C(0x80);
    ssize_t written = send(descriptor, encoded, sizeof(encoded), MSG_NOSIGNAL);
    char response;
    ssize_t received = written == (ssize_t) sizeof(encoded) ? recv(descriptor, &response, sizeof(response), 0) : -1;

    (void) close(descriptor);

    return written == (ssize_t) sizeof(encoded) && received == 0;
}

static int connect_raw_client(uint16_t port)
{
    int descriptor = socket(AF_INET, SOCK_STREAM | SOCK_CLOEXEC, 0);
    struct sockaddr_in address = {
        .sin_family = AF_INET,
        .sin_port = htons(port),
        .sin_addr.s_addr = htonl(INADDR_LOOPBACK),
    };

    if (descriptor < 0 || connect(descriptor, (const struct sockaddr *) &address, sizeof(address)) != 0) {
        if (descriptor >= 0) {
            (void) close(descriptor);
        }

        return -1;
    }

    return descriptor;
}

static bool wait_for_active_connections(GeoServer *server, uint64_t expected)
{
    for (size_t attempt = 0; attempt < 1000000U; ++attempt) {
        GeoServerStats stats;

        if (!geo_server_get_stats(server, &stats)) {
            return false;
        }

        if (stats.active_connections == expected) {
            return true;
        }

        sched_yield();
    }

    return false;
}

static bool test_connection_limit_backpressure(GeoServer *server, uint16_t port, size_t max_connections)
{
    if (max_connections < 2U || max_connections > 64U) {
        return false;
    }

    int retained_clients[64];
    size_t retained_count = 0;
    bool succeeded = true;

    while (retained_count + 1U < max_connections) {
        int descriptor = connect_raw_client(port);

        if (descriptor < 0) {
            succeeded = false;
            break;
        }

        retained_clients[retained_count++] = descriptor;
    }

    if (succeeded) {
        succeeded = wait_for_active_connections(server, max_connections);
    }

    int rejected_client = succeeded ? connect_raw_client(port) : -1;

    if (rejected_client < 0) {
        succeeded = false;
    } else {
        struct timeval timeout = {
            .tv_sec = 1,
        };
        unsigned char ignored;

        succeeded = setsockopt(rejected_client, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(timeout)) == 0 &&
                    recv(rejected_client, &ignored, sizeof(ignored), 0) == 0;
        (void) close(rejected_client);
    }

    for (size_t index = 0; index < retained_count; ++index) {
        (void) close(retained_clients[index]);
    }

    return succeeded && wait_for_active_connections(server, 1U);
}

static bool test_payload_memory_backpressure(GeoServer *server, uint16_t port, uint32_t payload_size)
{
    int reserving_client = connect_raw_client(port);

    if (reserving_client < 0) {
        return false;
    }

    GeoProtocolHeader request = {
        .request_id = 700U,
        .payload_checksum = geo_protocol_checksum(NULL, 0),
        .opcode = GEO_PROTOCOL_PING,
        .payload_size = payload_size,
    };
    unsigned char encoded[GEO_PROTOCOL_HEADER_SIZE];

    geo_protocol_encode_header(encoded, &request);

    if (send(reserving_client, encoded, sizeof(encoded), MSG_NOSIGNAL) != (ssize_t) sizeof(encoded)) {
        (void) close(reserving_client);
        return false;
    }

    GeoServerStats stats;
    size_t attempts = 0;

    do {
        if (!geo_server_get_stats(server, &stats)) {
            (void) close(reserving_client);
            return false;
        }

        if (stats.current_inflight_payload_bytes == payload_size) {
            break;
        }

        sched_yield();
        attempts++;
    } while (attempts < 1000000U);

    if (stats.current_inflight_payload_bytes != payload_size) {
        (void) close(reserving_client);
        return false;
    }

    int rejected_client = connect_raw_client(port);

    request.request_id++;
    geo_protocol_encode_header(encoded, &request);

    bool succeeded = rejected_client >= 0 &&
                     send(rejected_client, encoded, sizeof(encoded), MSG_NOSIGNAL) == (ssize_t) sizeof(encoded) &&
                     recv(rejected_client, encoded, sizeof(encoded), MSG_WAITALL) == (ssize_t) sizeof(encoded);
    GeoProtocolHeader response;

    succeeded = succeeded && geo_protocol_decode_header(encoded, &response) &&
                response.flags == GEO_PROTOCOL_RESPONSE_FLAG &&
                response.request_id == request.request_id &&
                response.status == GEO_PROTOCOL_STATUS_BUSY;

    if (rejected_client >= 0) {
        (void) close(rejected_client);
    }

    (void) close(reserving_client);

    return succeeded;
}

static bool test_server_end_to_end(void)
{
    enum {
        INITIAL_OBJECTS = 128,
        CLIENT_THREADS = 6,
        OBJECTS_PER_CLIENT = 32,
    };

    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-server-test-XXXXXX";
    char *directory = mkdtemp(directory_template);
    GeoServer *server = NULL;
    GeoClient *client = NULL;
    pthread_t server_thread;
    bool server_thread_started = false;
    ServerThreadContext server_context = { 0 };
    pthread_t client_threads[CLIENT_THREADS];
    ClientWorkerContext client_contexts[CLIENT_THREADS];
    size_t clients_started = 0;

    SERVER_TEST_ASSERT(directory != NULL, "server database directory must be created");

    GeoServerConfig server_config = geo_server_default_config(directory);

    server_config.port = 0;
    server_config.authentication_token = "correct horse battery staple";
    server_config.worker_threads = 4U;
    server_config.work_queue_capacity = 32U;
    server_config.max_connections = 8U;
    server_config.max_frame_size = 1024U * 1024U;
    server_config.max_inflight_payload_bytes = server_config.max_frame_size;

    GeoServerStatus server_status;
    server = geo_server_create(&server_config, &server_status);

    SERVER_TEST_ASSERT(server && server_status == GEO_SERVER_OK && geo_server_port(server) != 0,
                       "server must bind an ephemeral loopback port");

    server_context.server = server;
    SERVER_TEST_ASSERT(pthread_create(&server_thread, NULL, server_thread_main, &server_context) == 0,
                       "server reactor thread must start");
    server_thread_started = true;

    GeoClientConfig bad_config = geo_client_default_config("127.0.0.1", geo_server_port(server));

    bad_config.authentication_token = "wrong token";
    GeoClientStatus client_status;
    GeoClient *bad_client = geo_client_connect(&bad_config, &client_status);

    SERVER_TEST_ASSERT(!bad_client && client_status == GEO_CLIENT_UNAUTHORIZED, "invalid credentials must be rejected");

    GeoClientConfig missing_token_config = geo_client_default_config("127.0.0.1", geo_server_port(server));
    GeoClient *missing_token_client = geo_client_connect(&missing_token_config, &client_status);

    SERVER_TEST_ASSERT(!missing_token_client && client_status == GEO_CLIENT_UNAUTHORIZED,
                       "a driver connection without the required token must fail during its authentication handshake");

    GeoClientConfig client_config = geo_client_default_config("127.0.0.1", geo_server_port(server));

    client_config.authentication_token = "correct horse battery staple";
    client = geo_client_connect(&client_config, &client_status);
    SERVER_TEST_ASSERT(client && client_status == GEO_CLIENT_OK, "authenticated client must connect");
    SERVER_TEST_ASSERT(geo_client_ping(client, "GeoBolt", 7U) == GEO_CLIENT_OK, "ping payload must round-trip exactly");
    SERVER_TEST_ASSERT(test_shared_client_serialization(client),
                       "one shared driver connection must serialize concurrent request/response transactions safely");

    GeoDatabaseMutation initial[INITIAL_OBJECTS];

    for (size_t object = 0; object < INITIAL_OBJECTS; ++object) {
        initial[object] = (GeoDatabaseMutation) {
            .object_id = object + 1U,
            .morton_code = geo_encode(-30.0 + (double) object * 0.1, -60.0 + (double) object * 0.2),
            .operation = GEO_DATABASE_UPSERT,
        };
    }

    SERVER_TEST_ASSERT(geo_client_write(client, initial, INITIAL_OBJECTS) == GEO_CLIENT_OK, "network batch must commit");

    uint64_t count;
    SERVER_TEST_ASSERT(geo_client_radius_count(client, 0.0, 0.0, 25000.0, &count) == GEO_CLIENT_OK && count == INITIAL_OBJECTS,
                       "network radius count must observe the committed batch");

    uint16_t port = geo_server_port(server);

    SERVER_TEST_ASSERT(test_connection_limit_backpressure(server, port, server_config.max_connections),
                       "connection admission must reject sockets beyond the configured live-connection limit");
    SERVER_TEST_ASSERT(test_payload_memory_backpressure(server, port, 768U * 1024U),
                       "global payload budget must reject memory pressure before allocating another frame");

    for (size_t worker = 0; worker < CLIENT_THREADS; ++worker) {
        client_contexts[worker] = (ClientWorkerContext) {
            .port = port,
            .first_id = UINT64_C(10000) + worker * OBJECTS_PER_CLIENT,
            .object_count = OBJECTS_PER_CLIENT,
        };
        SERVER_TEST_ASSERT(pthread_create(client_threads + worker, NULL, client_worker_main, client_contexts + worker) == 0,
                           "concurrent client must start");
        clients_started++;
    }

    for (size_t worker = 0; worker < clients_started; ++worker) {
        SERVER_TEST_ASSERT(pthread_join(client_threads[worker], NULL) == 0, "concurrent client must join");
        SERVER_TEST_ASSERT(client_contexts[worker].succeeded, "concurrent client commands must succeed");
    }

    clients_started = 0;
    uint64_t expected_count = INITIAL_OBJECTS + CLIENT_THREADS * OBJECTS_PER_CLIENT;
    SERVER_TEST_ASSERT(geo_client_radius_count(client, 0.0, 0.0, 25000.0, &count) == GEO_CLIENT_OK && count == expected_count,
                       "all independent connections must publish visible objects");

    GeoDatabaseStats database_stats;
    SERVER_TEST_ASSERT(geo_client_get_database_stats(client, &database_stats) == GEO_CLIENT_OK &&
                           database_stats.committed_operations == expected_count,
                       "database stats must cross the wire without native-layout dependence");
    SERVER_TEST_ASSERT(geo_client_checkpoint(client) == GEO_CLIENT_OK, "remote checkpoint must complete");
    SERVER_TEST_ASSERT(send_corrupted_header(port), "corrupted header checksum must close the connection");

    GeoServerStats stats;
    SERVER_TEST_ASSERT(geo_server_get_stats(server, &stats) &&
                           stats.authentication_failures == 2U &&
                           stats.protocol_errors >= 1U &&
                           stats.rejected_connections >= 1U,
                       "server monitoring must account for auth and protocol failures");

    geo_client_close(client);
    client = NULL;
    geo_server_request_stop(server);
    SERVER_TEST_ASSERT(pthread_join(server_thread, NULL) == 0, "server reactor must stop cleanly");
    server_thread_started = false;
    SERVER_TEST_ASSERT(server_context.status == GEO_SERVER_OK, "normal stop must not report a reactor failure");

    geo_server_destroy(server);
    server = NULL;

    server_config.create_database_if_missing = false;
    server = geo_server_create(&server_config, &server_status);
    SERVER_TEST_ASSERT(server && server_status == GEO_SERVER_OK, "server must reopen the existing durable database");

    server_context = (ServerThreadContext) { .server = server };
    SERVER_TEST_ASSERT(pthread_create(&server_thread, NULL, server_thread_main, &server_context) == 0,
                       "reopened server reactor must start");
    server_thread_started = true;

    client_config.port = geo_server_port(server);
    client = geo_client_connect(&client_config, &client_status);
    SERVER_TEST_ASSERT(client && client_status == GEO_CLIENT_OK, "client must reconnect after server restart");
    SERVER_TEST_ASSERT(geo_client_radius_count(client, 0.0, 0.0, 25000.0, &count) == GEO_CLIENT_OK && count == expected_count,
                       "server restart must preserve committed spatial state");

cleanup:
    for (size_t worker = 0; worker < clients_started; ++worker) {
        (void) pthread_join(client_threads[worker], NULL);
    }

    geo_client_close(client);

    if (server_thread_started) {
        geo_server_request_stop(server);
        (void) pthread_join(server_thread, NULL);
    }

    geo_server_destroy(server);

    if (directory && !geobolt_test_remove_tree(directory)) {
        fprintf(stderr, "FAIL: server test directory could not be removed\n");
        succeeded = false;
    }

    return succeeded;
}

int main(void)
{
    printf("========================================\n");
    printf("GEOBOLT SERVER TEST SUITE\n");
    printf("========================================\n\n");

    if (!test_server_end_to_end()) {
        return EXIT_FAILURE;
    }

    printf("[PASS] auth, protocol, shared/multiple clients, admission limits, shutdown, and restart\n");

    if (!test_bounded_job_queue_backpressure()) {
        return EXIT_FAILURE;
    }

    printf("[PASS] bounded queue backpressure\n");

    return EXIT_SUCCESS;
}
