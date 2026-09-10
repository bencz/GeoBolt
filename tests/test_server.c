#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "geobolt/client.h"
#include "geobolt/geodoc.h"
#include "geobolt/server.h"
#include "geo_job_queue.h"
#include "geo_protocol.h"
#include "test_support.h"

#include <arpa/inet.h>
#include <errno.h>
#include <inttypes.h>
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
    bool disconnected;
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

        GeoClientStatus expected = context->disconnected ? GEO_CLIENT_IO_ERROR : GEO_CLIENT_OK;

        if (geo_client_ping(context->client, payload, sizeof(payload)) != expected) {
            atomic_store_explicit(context->failed, true, memory_order_release);
            break;
        }
    }

    return NULL;
}

static bool test_shared_client_serialization(GeoClient *client, bool disconnected)
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
            .disconnected = disconnected,
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

static bool transfer_socket_bytes(int descriptor, void *buffer, size_t size, bool sending)
{
    unsigned char *bytes = buffer;

    while (size) {
        ssize_t transferred = sending ? send(descriptor, bytes, size, MSG_NOSIGNAL) : recv(descriptor, bytes, size, 0);

        if (transferred < 0 && errno == EINTR) {
            continue;
        }
        if (transferred <= 0) {
            return false;
        }

        bytes += (size_t) transferred;
        size -= (size_t) transferred;
    }

    return true;
}

static bool send_raw_frame(int descriptor, uint16_t opcode, uint64_t request_id, void *payload, uint32_t payload_size)
{
    GeoProtocolHeader header = {
        .request_id = request_id,
        .payload_checksum = geo_protocol_checksum(payload, payload_size),
        .opcode = opcode,
        .payload_size = payload_size,
    };
    unsigned char encoded[GEO_PROTOCOL_HEADER_SIZE];

    geo_protocol_encode_header(encoded, &header);

    return transfer_socket_bytes(descriptor, encoded, sizeof(encoded), true) &&
           transfer_socket_bytes(descriptor, payload, payload_size, true);
}

static bool receive_raw_response(int descriptor,
                                  uint16_t opcode,
                                  uint64_t request_id,
                                  GeoProtocolStatus status,
                                  const void *payload,
                                  uint32_t payload_size)
{
    unsigned char encoded[GEO_PROTOCOL_HEADER_SIZE];
    GeoProtocolHeader response;

    if (!transfer_socket_bytes(descriptor, encoded, sizeof(encoded), false) ||
        !geo_protocol_decode_header(encoded, &response) ||
        response.request_id != request_id || response.opcode != opcode ||
        response.flags != GEO_PROTOCOL_RESPONSE_FLAG || response.status != (uint32_t) status ||
        response.payload_size != payload_size || response.payload_checksum != geo_protocol_checksum(payload, payload_size)) {
        return false;
    }

    const unsigned char *expected = payload;
    size_t remaining = payload_size;

    while (remaining) {
        unsigned char received[4096];
        size_t chunk_size = remaining < sizeof(received) ? remaining : sizeof(received);

        if (!transfer_socket_bytes(descriptor, received, chunk_size, false) || memcmp(received, expected, chunk_size) != 0) {
            return false;
        }

        expected += chunk_size;
        remaining -= chunk_size;
    }

    return true;
}

static bool test_incomplete_and_malformed_writes(uint16_t port)
{
    unsigned char frame[GEO_PROTOCOL_HEADER_SIZE + sizeof(uint32_t) + GEO_PROTOCOL_OBJECT_OPERATION_SIZE] = { 0 };
    unsigned char *payload = frame + GEO_PROTOCOL_HEADER_SIZE;

    geo_protocol_store_u32(payload, 1U);
    geo_protocol_store_u64(payload + 4U, 123U);
    geo_protocol_store_u64(payload + 12U, geo_encode(-23.5, -46.5));
    geo_protocol_store_u32(payload + 20U, GEO_DATABASE_INSERT);
    GeoProtocolHeader header = {
        .request_id = 333U,
        .opcode = GEO_PROTOCOL_WRITE_OBJECTS,
        .payload_size = sizeof(frame) - GEO_PROTOCOL_HEADER_SIZE,
        .payload_checksum = geo_protocol_checksum(payload, sizeof(frame) - GEO_PROTOCOL_HEADER_SIZE),
    };
    geo_protocol_encode_header(frame, &header);
    const size_t cutoffs[] = { 1U, GEO_PROTOCOL_HEADER_SIZE - 1U, GEO_PROTOCOL_HEADER_SIZE, sizeof(frame) - 1U, sizeof(frame) };

    for (size_t index = 0U; index < sizeof(cutoffs) / sizeof(cutoffs[0]); ++index) {
        int descriptor = connect_raw_client(port);
        struct timeval timeout = { .tv_sec = 10 };
        bool succeeded = descriptor >= 0 &&
                         setsockopt(descriptor, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(timeout)) == 0 &&
                         transfer_socket_bytes(descriptor, frame, cutoffs[index], true) && shutdown(descriptor, SHUT_WR) == 0;

        if (succeeded && cutoffs[index] == sizeof(frame)) {
            succeeded = receive_raw_response(descriptor, header.opcode, header.request_id, GEO_PROTOCOL_STATUS_OK, NULL, 0U);
        }

        unsigned char unexpected;
        succeeded = succeeded && recv(descriptor, &unexpected, 1U, 0) == 0;

        if (descriptor >= 0 && close(descriptor) != 0) {
            succeeded = false;
        }
        if (!succeeded) {
            return false;
        }
    }

    /* These declarations must be rejected before allocating an operation array. */
    int descriptor = connect_raw_client(port);
    struct timeval timeout = { .tv_sec = 10 };
    bool succeeded = descriptor >= 0 &&
                     setsockopt(descriptor, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(timeout)) == 0;
    const uint32_t counts[] = { 2U, UINT32_MAX };

    for (size_t index = 0U; succeeded && index < sizeof(counts) / sizeof(counts[0]); ++index) {
        geo_protocol_store_u32(payload, counts[index]);
        succeeded = send_raw_frame(descriptor, header.opcode, header.request_id, payload, header.payload_size) &&
                    receive_raw_response(descriptor, header.opcode, header.request_id,
                                         GEO_PROTOCOL_STATUS_INVALID_REQUEST, NULL, 0U);
    }

    if (descriptor >= 0 && close(descriptor) != 0) {
        succeeded = false;
    }

    GeoClientConfig config = geo_client_default_config("127.0.0.1", port);
    GeoClientStatus status;
    GeoClient *client = geo_client_connect(&config, &status);
    GeoDatabaseStats stats;

    succeeded = succeeded && client && status == GEO_CLIENT_OK &&
                geo_client_get_database_stats(client, &stats) == GEO_CLIENT_OK && stats.committed_operations == 1U;
    geo_client_close(client);
    return succeeded;
}

static bool test_large_response_after_half_close(uint16_t port)
{
    const uint32_t payload_size = 8U * 1024U * 1024U;
    unsigned char *payload = malloc(payload_size);
    int descriptor = connect_raw_client(port);
    struct timeval timeout = { .tv_sec = 10 };
    int receive_buffer = 128 * 1024;
    bool succeeded = payload && descriptor >= 0 &&
                     setsockopt(descriptor, SOL_SOCKET, SO_RCVBUF, &receive_buffer, sizeof(receive_buffer)) == 0 &&
                     setsockopt(descriptor, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(timeout)) == 0 &&
                     setsockopt(descriptor, SOL_SOCKET, SO_SNDTIMEO, &timeout, sizeof(timeout)) == 0;

    if (payload) {
        for (size_t index = 0U; index < payload_size; ++index) {
            payload[index] = (unsigned char) (index * 37U + index / 4096U);
        }
    }

    succeeded = succeeded && send_raw_frame(descriptor, GEO_PROTOCOL_PING, 444U, payload, payload_size) &&
                shutdown(descriptor, SHUT_WR) == 0;

    /* Hold the large response unread while another connection completes work. The reactor
     * must return from a blocked send, serve this client, then resume the response on EPOLLOUT. */
    GeoClientConfig config = geo_client_default_config("127.0.0.1", port);
    GeoClientStatus status;
    GeoClient *probe = geo_client_connect(&config, &status);

    succeeded = succeeded && probe && status == GEO_CLIENT_OK && geo_client_ping(probe, NULL, 0U) == GEO_CLIENT_OK;
    geo_client_close(probe);

    succeeded = succeeded && receive_raw_response(descriptor, GEO_PROTOCOL_PING, 444U,
                                                  GEO_PROTOCOL_STATUS_OK, payload, payload_size);

    unsigned char unexpected;
    succeeded = succeeded && recv(descriptor, &unexpected, 1U, 0) == 0;

    if (descriptor >= 0 && close(descriptor) != 0) {
        succeeded = false;
    }

    free(payload);
    return succeeded;
}

typedef struct {
    uint16_t port;
    size_t worker;
    bool succeeded;
} DisconnectWorkerContext;

static void *disconnect_worker_main(void *argument)
{
    DisconnectWorkerContext *context = argument;
    unsigned char payload[4096];
    const uint32_t sizes[] = { 0U, 1U, sizeof(payload) };

    memset(payload, (int) context->worker + 1, sizeof(payload));
    context->succeeded = true;

    for (size_t round = 0U; context->succeeded && round < 128U; ++round) {
        int descriptor = connect_raw_client(context->port);
        struct timeval timeout = { .tv_sec = 10 };

        context->succeeded = descriptor >= 0 &&
                             setsockopt(descriptor, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(timeout)) == 0 &&
                             setsockopt(descriptor, SOL_SOCKET, SO_SNDTIMEO, &timeout, sizeof(timeout)) == 0;

        for (size_t frame = 0U; context->succeeded && frame < 3U; ++frame) {
            context->succeeded = send_raw_frame(descriptor, GEO_PROTOCOL_PING, round * 3U + frame, payload, sizes[frame]);
        }

        if (context->succeeded && round % 4U == 0U) {
            /* RST races queued/running jobs and other connections' completion events. */
            struct linger reset = { .l_onoff = 1, .l_linger = 0 };

            context->succeeded = setsockopt(descriptor, SOL_SOCKET, SO_LINGER, &reset, sizeof(reset)) == 0;
        } else if (context->succeeded) {
            context->succeeded = shutdown(descriptor, SHUT_WR) == 0;

            for (size_t frame = 0U; context->succeeded && frame < 3U; ++frame) {
                context->succeeded = receive_raw_response(descriptor, GEO_PROTOCOL_PING, round * 3U + frame,
                                                          GEO_PROTOCOL_STATUS_OK, payload, sizes[frame]);
            }

            if (context->succeeded) {
                unsigned char unexpected;

                context->succeeded = recv(descriptor, &unexpected, 1U, 0) == 0;
            }
        }

        if (descriptor >= 0 && close(descriptor) != 0) {
            context->succeeded = false;
        }
    }

    return NULL;
}

static bool wait_for_active_connections(GeoServer *server, uint64_t expected);

static bool test_half_close_and_concurrent_disconnects(void)
{
    enum {
        WORKERS = 8,
    };
    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-disconnect-test-XXXXXX";
    char *directory = mkdtemp(directory_template);
    GeoServer *server = NULL;
    pthread_t reactor;
    bool reactor_started = false;
    ServerThreadContext reactor_context = { 0 };
    pthread_t workers[WORKERS];
    DisconnectWorkerContext contexts[WORKERS];
    size_t workers_started = 0U;

    SERVER_TEST_ASSERT(directory != NULL, "disconnect test directory must be created");
    GeoServerConfig config = geo_server_default_config(directory);
    GeoServerStatus status;

    config.port = 0U;
    server = geo_server_create(&config, &status);
    SERVER_TEST_ASSERT(server && status == GEO_SERVER_OK, "disconnect test server must open");
    reactor_context.server = server;
    SERVER_TEST_ASSERT(pthread_create(&reactor, NULL, server_thread_main, &reactor_context) == 0,
                       "disconnect test reactor must start");
    reactor_started = true;

    SERVER_TEST_ASSERT(test_incomplete_and_malformed_writes(geo_server_port(server)),
                       "only a complete validated write may commit, and its response must survive FIN");
    SERVER_TEST_ASSERT(test_large_response_after_half_close(geo_server_port(server)),
                       "a blocked large response must resume after FIN while other clients progress");

    for (size_t worker = 0U; worker < WORKERS; ++worker) {
        contexts[worker] = (DisconnectWorkerContext) { .port = geo_server_port(server), .worker = worker };
        SERVER_TEST_ASSERT(pthread_create(workers + worker, NULL, disconnect_worker_main, contexts + worker) == 0,
                           "disconnect worker must start");
        workers_started++;
    }

cleanup:
    for (size_t worker = 0U; worker < workers_started; ++worker) {
        if (pthread_join(workers[worker], NULL) != 0 || !contexts[worker].succeeded) {
            fprintf(stderr, "FAIL: worker %zu lost ordered responses after half-close\n", worker);
            succeeded = false;
        }
    }

    if (succeeded) {
        GeoServerStats stats;

        succeeded = wait_for_active_connections(server, 0U) && geo_server_get_stats(server, &stats) &&
                    stats.current_inflight_payload_bytes == 0U;

        if (!succeeded) {
            fprintf(stderr, "FAIL: closed connections must release all payloads before server destruction\n");
        }
    }

    if (reactor_started) {
        geo_server_request_stop(server);

        if (pthread_join(reactor, NULL) != 0 || reactor_context.status != GEO_SERVER_OK) {
            succeeded = false;
        }
    }

    geo_server_destroy(server);

    if (directory && !geobolt_test_remove_tree(directory)) {
        succeeded = false;
    }

    return succeeded;
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
    GeoDocBuilder *metadata_builder = NULL;
    GeoDocBuffer metadata = { 0 };
    GeoDatabaseObject remote_object = { 0 };
    GeoIdResult *index_result = NULL;
    GeoSearchResult *radius_result = NULL;
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

    if (!server || server_status != GEO_SERVER_OK || geo_server_port(server) == 0U) {
        fprintf(stderr, "server creation diagnostic: status=%d port=%u\n",
                server_status,
                server ? geo_server_port(server) : 0U);
    }

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
    SERVER_TEST_ASSERT(test_shared_client_serialization(client, false),
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

    metadata_builder = geo_doc_builder_create();
    SERVER_TEST_ASSERT(metadata_builder != NULL, "remote metadata builder must allocate");
    SERVER_TEST_ASSERT(geo_doc_builder_add_string(metadata_builder, 0U, "kind", "vehicle", 7U, NULL) == GEO_DOC_OK &&
                           geo_doc_builder_add_string(metadata_builder, 0U, "driver", "Ana", 3U, NULL) == GEO_DOC_OK &&
                           geo_doc_builder_finish(metadata_builder, &metadata) == GEO_DOC_OK,
                       "remote GeoDoc metadata must serialize");
    SERVER_TEST_ASSERT(geo_client_insert(client, 9000U, -23.5505, -46.6333, metadata.data, metadata.size) == GEO_CLIENT_OK,
                       "remote conditional insert with metadata must commit");
    SERVER_TEST_ASSERT(geo_client_create_index(client, "kind_idx", "/kind", GEO_DATABASE_INDEX_STRING) == GEO_CLIENT_OK,
                       "remote typed-index creation must build through the production server path");
    index_result = geo_id_result_create(4U);
    radius_result = geo_result_create(4U);
    SERVER_TEST_ASSERT(index_result != NULL && radius_result != NULL, "remote typed and spatial results must allocate");
    GeoDatabaseIndexPredicate index_predicate = {
        .index_name = "kind_idx",
        .operation = GEO_DATABASE_INDEX_EQUAL,
        .lower = { .type = GEO_DATABASE_INDEX_STRING, .as.bytes = { .data = "vehicle", .size = 7U } },
    };
    SERVER_TEST_ASSERT(geo_client_query_index_reuse(client, &index_predicate, index_result, NULL) == GEO_CLIENT_OK &&
                           index_result->count == 1U && index_result->ids[0] == 9000U,
                       "remote typed equality query must return IDs through the bounded wire response");
    GeoDatabaseRadiusQuery remote_radius_query = {
        .predicates = &index_predicate,
        .predicate_count = 1U,
        .latitude = -23.5505,
        .longitude = -46.6333,
        .radius_km = 5.0,
    };
    GeoDatabaseQueryStats remote_query_stats = { 0 };
    GeoClientStatus radius_query_status = geo_client_query_radius_reuse(client,
                                                                        &remote_radius_query,
                                                                        radius_result,
                                                                        &remote_query_stats);

    if (radius_query_status != GEO_CLIENT_OK || radius_result->count != 1U || radius_result->results[0].id != 9000U) {
        fprintf(stderr,
                "radius query diagnostic: status=%d count=%zu metadata=%" PRIu64 " plan=%d\n",
                radius_query_status,
                radius_result->count,
                remote_query_stats.metadata_candidates,
                remote_query_stats.plan);
    }

    SERVER_TEST_ASSERT(radius_query_status == GEO_CLIENT_OK &&
                           radius_result->count == 1U && radius_result->results[0].id == 9000U &&
                           remote_query_stats.metadata_candidates == 1U && remote_query_stats.predicate_count == 1U,
                       "remote conjunctive radius query must execute through the reusable production planner");
    SERVER_TEST_ASSERT(geo_client_insert(client, 9000U, -23.5505, -46.6333, metadata.data, metadata.size) ==
                           GEO_CLIENT_ALREADY_EXISTS,
                       "remote insert must reject an existing generic object ID");
    SERVER_TEST_ASSERT(geo_client_get(client, 9000U, &remote_object) == GEO_CLIENT_OK &&
                           remote_object.document_size == metadata.size &&
                           memcmp(remote_object.document, metadata.data, metadata.size) == 0,
                       "remote GET must return exact location, sequence, and GeoDoc bytes");
    geo_database_object_release(&remote_object);
    SERVER_TEST_ASSERT(geo_client_update(client, 9000U, 51.5074, -0.1278, metadata.data, metadata.size) == GEO_CLIENT_OK,
                       "remote conditional update must replace an existing object");
    SERVER_TEST_ASSERT(geo_client_update(client, 9900U, 0.0, 0.0, metadata.data, metadata.size) == GEO_CLIENT_NOT_FOUND,
                       "remote update must reject a missing object");
    SERVER_TEST_ASSERT(geo_client_delete(client, 9000U) == GEO_CLIENT_OK,
                       "remote delete must remove location and metadata atomically");
    SERVER_TEST_ASSERT(geo_client_get(client, 9000U, &remote_object) == GEO_CLIENT_NOT_FOUND,
                       "remote GET must report a deleted object as not found");
    SERVER_TEST_ASSERT(geo_client_query_index_reuse(client, &index_predicate, index_result, NULL) == GEO_CLIENT_OK &&
                           index_result->count == 0U,
                       "remote delete must remove the typed secondary key in the same canonical commit");

    uint64_t generated_object_id = 0U;
    GeoClientStatus generated_status = geo_client_insert_generated(client,
                                                                   48.8566,
                                                                   2.3522,
                                                                   metadata.data,
                                                                   metadata.size,
                                                                   &generated_object_id);

    SERVER_TEST_ASSERT(generated_status == GEO_CLIENT_OK &&
                           generated_object_id >= (UINT64_C(1) << 63U),
                       "remote generated insert must return a durable generic object ID");
    SERVER_TEST_ASSERT(geo_client_get(client, generated_object_id, &remote_object) == GEO_CLIENT_OK &&
                           remote_object.document_size == metadata.size,
                       "remote GET must resolve a database-generated object ID");
    geo_database_object_release(&remote_object);

    GeoDatabaseObjectMutation idempotent_mutation = {
        .object_id = 9100U,
        .morton_code = geo_encode(41.1579, -8.6291),
        .document = metadata.data,
        .document_size = metadata.size,
        .operation = GEO_DATABASE_UPSERT,
    };
    static const char idempotency_key[] = "remote-ingest-42";
    bool replayed;
    SERVER_TEST_ASSERT(geo_client_write_objects_idempotent(client,
                                                           &idempotent_mutation,
                                                           1U,
                                                           idempotency_key,
                                                           sizeof(idempotency_key) - 1U,
                                                           3600U,
                                                           &replayed) == GEO_CLIENT_OK &&
                           !replayed,
                       "first remote idempotent write must commit");
    SERVER_TEST_ASSERT(geo_client_write_objects_idempotent(client,
                                                           &idempotent_mutation,
                                                           1U,
                                                           idempotency_key,
                                                           sizeof(idempotency_key) - 1U,
                                                           3600U,
                                                           &replayed) == GEO_CLIENT_OK &&
                           replayed,
                       "remote retry must return the durable replay decision");
    idempotent_mutation.morton_code = geo_encode(42.0, -8.0);
    SERVER_TEST_ASSERT(geo_client_write_objects_idempotent(client,
                                                           &idempotent_mutation,
                                                           1U,
                                                           idempotency_key,
                                                           sizeof(idempotency_key) - 1U,
                                                           3600U,
                                                           &replayed) == GEO_CLIENT_IDEMPOTENCY_CONFLICT,
                       "remote key reuse with different content must be rejected");

    GeoDatabaseIndexInfo remote_indexes[1];
    size_t remote_index_count = 0U;
    SERVER_TEST_ASSERT(geo_client_list_indexes(client, remote_indexes, 1U, &remote_index_count) == GEO_CLIENT_OK &&
                           remote_index_count == 1U && remote_indexes[0].entry_count == 2U &&
                           strcmp(remote_indexes[0].name, "kind_idx") == 0,
                       "remote index catalog must expose exact post-mutation cardinality");

    uint64_t count;
    SERVER_TEST_ASSERT(geo_client_radius_count(client, 0.0, 0.0, 25000.0, &count) == GEO_CLIENT_OK &&
                           count == INITIAL_OBJECTS + 2U,
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
    uint64_t expected_count = INITIAL_OBJECTS + CLIENT_THREADS * OBJECTS_PER_CLIENT + 2U;
    SERVER_TEST_ASSERT(geo_client_radius_count(client, 0.0, 0.0, 25000.0, &count) == GEO_CLIENT_OK && count == expected_count,
                       "all independent connections must publish visible objects");

    GeoDatabaseStats database_stats;
    SERVER_TEST_ASSERT(geo_client_get_database_stats(client, &database_stats) == GEO_CLIENT_OK &&
                           database_stats.committed_operations == expected_count + 3U,
                       "database stats must cross the wire without native-layout dependence");
    SERVER_TEST_ASSERT(geo_client_checkpoint(client) == GEO_CLIENT_OK, "remote checkpoint must complete");
    SERVER_TEST_ASSERT(send_corrupted_header(port), "corrupted header checksum must close the connection");

    GeoServerStats stats;
    SERVER_TEST_ASSERT(geo_server_get_stats(server, &stats) &&
                           stats.authentication_failures == 2U &&
                           stats.protocol_errors >= 1U &&
                           stats.rejected_connections >= 1U,
                       "server monitoring must account for auth and protocol failures");

    geo_server_request_stop(server);
    SERVER_TEST_ASSERT(pthread_join(server_thread, NULL) == 0, "server reactor must stop cleanly");
    server_thread_started = false;
    SERVER_TEST_ASSERT(server_context.status == GEO_SERVER_OK, "normal stop must not report a reactor failure");

    geo_server_destroy(server);
    server = NULL;

    SERVER_TEST_ASSERT(test_shared_client_serialization(client, true),
                       "shared client must serialize disconnect state and report I/O failure to every caller");
    geo_client_close(client);
    client = NULL;

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
    SERVER_TEST_ASSERT(geo_client_query_index_reuse(client, &index_predicate, index_result, NULL) == GEO_CLIENT_OK &&
                           index_result->count == 2U,
                       "server restart must preserve typed secondary definitions and keys");
    SERVER_TEST_ASSERT(geo_client_drop_index(client, "kind_idx") == GEO_CLIENT_OK,
                       "remote typed-index drop must remove the definition through protocol v4");

cleanup:
    for (size_t worker = 0; worker < clients_started; ++worker) {
        (void) pthread_join(client_threads[worker], NULL);
    }

    geo_client_close(client);
    geo_result_destroy(radius_result);
    geo_id_result_destroy(index_result);
    geo_database_object_release(&remote_object);
    geo_doc_buffer_release(&metadata);
    geo_doc_builder_destroy(metadata_builder);

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

    if (!test_half_close_and_concurrent_disconnects()) {
        return EXIT_FAILURE;
    }

    printf("[PASS] ordered pipelined responses after half-close and concurrent connection resets\n");

    return EXIT_SUCCESS;
}
