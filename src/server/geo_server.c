#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "geo_server_internal.h"

#include <arpa/inet.h>
#include <errno.h>
#include <fcntl.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>
#include <sys/epoll.h>
#include <sys/eventfd.h>
#include <sys/socket.h>
#include <unistd.h>

#define GEO_SERVER_DEFAULT_WORKERS 4U
#define GEO_SERVER_DEFAULT_QUEUE_CAPACITY 1024U
#define GEO_SERVER_DEFAULT_MAX_CONNECTIONS 4096U
#define GEO_SERVER_DEFAULT_MAX_INFLIGHT_PAYLOAD_BYTES (256U * 1024U * 1024U)
#define GEO_SERVER_DEFAULT_BACKLOG 512
#define GEO_SERVER_MAX_EVENTS 256
#define GEO_SERVER_MAX_TOKEN_SIZE 4096U

static char *server_copy_string(const char *value)
{
    size_t size = strlen(value) + 1U;
    char *copy = malloc(size);

    if (copy) {
        memcpy(copy, value, size);
    }

    return copy;
}

static void server_erase_secret(char *secret, size_t size)
{
    volatile unsigned char *bytes = (volatile unsigned char *) secret;

    while (size) {
        *bytes++ = 0;
        size--;
    }
}

static bool server_set_nonblocking(int descriptor)
{
    int flags = fcntl(descriptor, F_GETFL, 0);

    return flags >= 0 && fcntl(descriptor, F_SETFL, flags | O_NONBLOCK) == 0;
}

static bool server_epoll_add(GeoServer *server, int descriptor, uint32_t events, GeoEventSource *source)
{
    struct epoll_event event = {
        .events = events,
        .data.ptr = source,
    };

    return epoll_ctl(server->epoll_descriptor, EPOLL_CTL_ADD, descriptor, &event) == 0;
}

static bool server_epoll_modify(GeoConnection *connection, uint32_t events)
{
    struct epoll_event event = {
        .events = events,
        .data.ptr = &connection->source,
    };

    return epoll_ctl(connection->server->epoll_descriptor, EPOLL_CTL_MOD, connection->descriptor, &event) == 0;
}

static void server_unlink_connection(GeoConnection *connection)
{
    GeoServer *server = connection->server;

    if (connection->previous) {
        connection->previous->next = connection->next;
    } else {
        server->connections = connection->next;
    }

    if (connection->next) {
        connection->next->previous = connection->previous;
    }

    atomic_fetch_sub_explicit(&server->active_connections, 1U, memory_order_relaxed);
}

static void server_release_request_payload(GeoConnection *connection)
{
    if (!connection->request_payload) {
        return;
    }

    atomic_fetch_sub_explicit(&connection->server->inflight_payload_bytes,
                              connection->request_header.payload_size,
                              memory_order_relaxed);
    free(connection->request_payload);
    connection->request_payload = NULL;
}

static void server_free_connection(GeoConnection *connection)
{
    server_unlink_connection(connection);
    free(connection->response);
    server_release_request_payload(connection);
    free(connection);
}

static void server_close_connection(GeoConnection *connection)
{
    if (connection->closing) {
        return;
    }

    connection->closing = true;

    if (connection->descriptor >= 0) {
        (void) epoll_ctl(connection->server->epoll_descriptor, EPOLL_CTL_DEL, connection->descriptor, NULL);
        (void) close(connection->descriptor);
        connection->descriptor = -1;
    }

    if (!connection->job_inflight) {
        server_free_connection(connection);
    }
}

bool geo_server_prepare_response(GeoConnection *connection,
                                 GeoProtocolStatus status,
                                 const void *payload,
                                 uint32_t payload_size)
{
    if (payload_size > connection->server->max_frame_size) {
        return false;
    }

    size_t response_size = GEO_PROTOCOL_HEADER_SIZE + (size_t) payload_size;
    unsigned char *response = malloc(response_size);

    if (!response) {
        return false;
    }

    GeoProtocolHeader header = {
        .request_id = connection->request_header.request_id,
        .payload_checksum = geo_protocol_checksum(payload, payload_size),
        .opcode = connection->request_header.opcode,
        .flags = GEO_PROTOCOL_RESPONSE_FLAG,
        .status = status,
        .payload_size = payload_size,
    };

    geo_protocol_encode_header(response, &header);

    if (payload_size) {
        memcpy(response + GEO_PROTOCOL_HEADER_SIZE, payload, payload_size);
    }

    free(connection->response);
    connection->response = response;
    connection->response_size = response_size;
    connection->response_sent = 0;

    return true;
}

void geo_server_complete_job(GeoConnection *connection)
{
    GeoServer *server = connection->server;

    if (pthread_mutex_lock(&server->completion_lock) != 0) {
        atomic_store_explicit(&server->event_loop_failed, true, memory_order_release);
        geo_server_request_stop(server);
        return;
    }

    connection->completion_next = NULL;

    if (server->completed_tail) {
        server->completed_tail->completion_next = connection;
    } else {
        server->completed_head = connection;
    }

    server->completed_tail = connection;
    pthread_mutex_unlock(&server->completion_lock);

    uint64_t signal_value = 1U;
    ssize_t ignored = write(server->wakeup_descriptor, &signal_value, sizeof(signal_value));

    (void) ignored;
}

static void server_reset_request(GeoConnection *connection)
{
    server_release_request_payload(connection);
    connection->header_received = 0;
    connection->payload_received = 0;
    connection->header_decoded = false;
    memset(&connection->request_header, 0, sizeof(connection->request_header));
}

static bool server_dispatch_request(GeoConnection *connection)
{
    if (!geo_job_queue_try_submit(connection->server->work_queue, geo_server_execute_request, connection)) {
        atomic_fetch_add_explicit(&connection->server->backpressure_rejections, 1U, memory_order_relaxed);

        if (!geo_server_prepare_response(connection, GEO_PROTOCOL_STATUS_BUSY, NULL, 0)) {
            return false;
        }

        server_reset_request(connection);

        return server_epoll_modify(connection, EPOLLOUT | EPOLLRDHUP);
    }

    connection->job_inflight = true;

    return server_epoll_modify(connection, EPOLLRDHUP);
}

static bool server_read_connection(GeoConnection *connection)
{
    for (;;) {
        unsigned char *destination;
        size_t remaining;

        if (!connection->header_decoded) {
            destination = connection->encoded_header + connection->header_received;
            remaining = GEO_PROTOCOL_HEADER_SIZE - connection->header_received;
        } else {
            destination = connection->request_payload + connection->payload_received;
            remaining = connection->request_header.payload_size - connection->payload_received;
        }

        if (!remaining) {
            if (!connection->header_decoded) {
                if (!geo_protocol_decode_header(connection->encoded_header, &connection->request_header) ||
                    connection->request_header.flags != 0 ||
                    connection->request_header.status != 0) {
                    atomic_fetch_add_explicit(&connection->server->protocol_errors, 1U, memory_order_relaxed);
                    return false;
                }

                if (connection->request_header.payload_size > connection->server->max_frame_size) {
                    connection->close_after_response = true;
                    return geo_server_prepare_response(connection, GEO_PROTOCOL_STATUS_FRAME_TOO_LARGE, NULL, 0) &&
                           server_epoll_modify(connection, EPOLLOUT | EPOLLRDHUP);
                }

                connection->header_decoded = true;

                if (connection->request_header.payload_size) {
                    uint64_t payload_size = connection->request_header.payload_size;
                    uint64_t inflight = atomic_load_explicit(&connection->server->inflight_payload_bytes, memory_order_relaxed);

                    if (inflight > connection->server->max_inflight_payload_bytes ||
                        payload_size > connection->server->max_inflight_payload_bytes - inflight) {
                        atomic_fetch_add_explicit(&connection->server->backpressure_rejections, 1U, memory_order_relaxed);
                        connection->close_after_response = true;

                        return geo_server_prepare_response(connection, GEO_PROTOCOL_STATUS_BUSY, NULL, 0) &&
                               server_epoll_modify(connection, EPOLLOUT | EPOLLRDHUP);
                    }

                    connection->request_payload = malloc(connection->request_header.payload_size);

                    if (!connection->request_payload) {
                        return false;
                    }

                    uint64_t updated = atomic_fetch_add_explicit(&connection->server->inflight_payload_bytes,
                                                                  payload_size,
                                                                  memory_order_relaxed) + payload_size;
                    uint64_t peak = atomic_load_explicit(&connection->server->peak_inflight_payload_bytes, memory_order_relaxed);

                    while (updated > peak &&
                           !atomic_compare_exchange_weak_explicit(&connection->server->peak_inflight_payload_bytes,
                                                                  &peak,
                                                                  updated,
                                                                  memory_order_relaxed,
                                                                  memory_order_relaxed)) {
                    }

                    continue;
                }
            }

            if (geo_protocol_checksum(connection->request_payload, connection->request_header.payload_size) !=
                connection->request_header.payload_checksum) {
                atomic_fetch_add_explicit(&connection->server->protocol_errors, 1U, memory_order_relaxed);
                return false;
            }

            return server_dispatch_request(connection);
        }

        ssize_t received = recv(connection->descriptor, destination, remaining, 0);

        if (received > 0) {
            if (connection->header_decoded) {
                connection->payload_received += (size_t) received;
            } else {
                connection->header_received += (size_t) received;
            }

            continue;
        }

        if (received == 0) {
            return false;
        }

        if (errno == EINTR) {
            continue;
        }

        return errno == EAGAIN || errno == EWOULDBLOCK;
    }
}

static bool server_write_connection(GeoConnection *connection)
{
    while (connection->response_sent < connection->response_size) {
        ssize_t written = send(connection->descriptor,
                               connection->response + connection->response_sent,
                               connection->response_size - connection->response_sent,
                               MSG_NOSIGNAL);

        if (written > 0) {
            connection->response_sent += (size_t) written;
            continue;
        }

        if (written < 0 && errno == EINTR) {
            continue;
        }

        return written < 0 && (errno == EAGAIN || errno == EWOULDBLOCK);
    }

    free(connection->response);
    connection->response = NULL;
    connection->response_size = 0;
    connection->response_sent = 0;

    if (connection->close_after_response) {
        return false;
    }

    return server_epoll_modify(connection, EPOLLIN | EPOLLRDHUP);
}

static bool server_process_completions(GeoServer *server)
{
    uint64_t value;

    while (read(server->wakeup_descriptor, &value, sizeof(value)) < 0 && errno == EINTR) {
    }

    if (pthread_mutex_lock(&server->completion_lock) != 0) {
        return false;
    }

    GeoConnection *completed = server->completed_head;

    server->completed_head = NULL;
    server->completed_tail = NULL;
    pthread_mutex_unlock(&server->completion_lock);

    while (completed) {
        GeoConnection *next = completed->completion_next;

        completed->completion_next = NULL;
        completed->job_inflight = false;
        atomic_fetch_add_explicit(&server->completed_requests, 1U, memory_order_relaxed);
        server_reset_request(completed);

        if (completed->closing || !completed->response ||
            !server_epoll_modify(completed, EPOLLOUT | EPOLLRDHUP)) {
            if (completed->closing) {
                server_free_connection(completed);
            } else {
                server_close_connection(completed);
            }
        }

        completed = next;
    }

    return true;
}

static void server_accept_connections(GeoServer *server)
{
    for (;;) {
        int descriptor = accept4(server->listener_descriptor, NULL, NULL, SOCK_NONBLOCK | SOCK_CLOEXEC);

        if (descriptor < 0) {
            if (errno == EINTR) {
                continue;
            }

            return;
        }

        if (atomic_load_explicit(&server->active_connections, memory_order_relaxed) >= server->max_connections) {
            atomic_fetch_add_explicit(&server->rejected_connections, 1U, memory_order_relaxed);
            (void) close(descriptor);
            continue;
        }

        int enabled = 1;

        if (setsockopt(descriptor, IPPROTO_TCP, TCP_NODELAY, &enabled, sizeof(enabled)) != 0 ||
            setsockopt(descriptor, SOL_SOCKET, SO_KEEPALIVE, &enabled, sizeof(enabled)) != 0) {
            atomic_fetch_add_explicit(&server->rejected_connections, 1U, memory_order_relaxed);
            (void) close(descriptor);
            continue;
        }

        GeoConnection *connection = calloc(1, sizeof(*connection));

        if (!connection) {
            atomic_fetch_add_explicit(&server->rejected_connections, 1U, memory_order_relaxed);
            (void) close(descriptor);
            continue;
        }

        connection->source.type = GEO_EVENT_CONNECTION;
        connection->server = server;
        connection->descriptor = descriptor;
        connection->authenticated = server->authentication_token_size == 0;
        connection->next = server->connections;

        if (server->connections) {
            server->connections->previous = connection;
        }

        server->connections = connection;
        atomic_fetch_add_explicit(&server->active_connections, 1U, memory_order_relaxed);

        if (!server_epoll_add(server, descriptor, EPOLLIN | EPOLLRDHUP, &connection->source)) {
            server_free_connection(connection);
            (void) close(descriptor);
            continue;
        }

        atomic_fetch_add_explicit(&server->accepted_connections, 1U, memory_order_relaxed);
    }
}

GeoServerConfig geo_server_default_config(const char *database_directory)
{
    return (GeoServerConfig) {
        .database_directory = database_directory,
        .bind_address = "127.0.0.1",
        .port = 7447U,
        .worker_threads = GEO_SERVER_DEFAULT_WORKERS,
        .work_queue_capacity = GEO_SERVER_DEFAULT_QUEUE_CAPACITY,
        .max_connections = GEO_SERVER_DEFAULT_MAX_CONNECTIONS,
        .max_frame_size = GEO_PROTOCOL_DEFAULT_MAX_FRAME,
        .max_inflight_payload_bytes = GEO_SERVER_DEFAULT_MAX_INFLIGHT_PAYLOAD_BYTES,
        .listen_backlog = GEO_SERVER_DEFAULT_BACKLOG,
        .create_database_if_missing = true,
    };
}

GeoServer *geo_server_create(const GeoServerConfig *config, GeoServerStatus *status)
{
    GeoServerStatus create_status = GEO_SERVER_INVALID_ARGUMENT;

    if (!config || !config->database_directory || !config->database_directory[0] ||
        !config->bind_address || !config->bind_address[0] ||
        !config->worker_threads || !config->work_queue_capacity || !config->max_connections ||
        !config->max_frame_size || config->max_frame_size > GEO_PROTOCOL_ABSOLUTE_MAX_FRAME ||
        config->max_inflight_payload_bytes < config->max_frame_size ||
        config->listen_backlog <= 0 ||
        (config->authentication_token && strlen(config->authentication_token) > GEO_SERVER_MAX_TOKEN_SIZE)) {
        goto complete;
    }

    GeoServer *server = calloc(1, sizeof(*server));

    if (!server) {
        create_status = GEO_SERVER_OUT_OF_MEMORY;
        goto complete;
    }

    server->epoll_descriptor = -1;
    server->listener_descriptor = -1;
    server->wakeup_descriptor = -1;
    server->listener_source.type = GEO_EVENT_LISTENER;
    server->wakeup_source.type = GEO_EVENT_WAKEUP;
    server->max_connections = config->max_connections;
    server->max_frame_size = config->max_frame_size;
    server->max_inflight_payload_bytes = config->max_inflight_payload_bytes;

    if (config->authentication_token && config->authentication_token[0]) {
        server->authentication_token = server_copy_string(config->authentication_token);

        if (!server->authentication_token) {
            create_status = GEO_SERVER_OUT_OF_MEMORY;
            geo_server_destroy(server);
            server = NULL;
            goto complete;
        }

        server->authentication_token_size = strlen(server->authentication_token);
    }

    if (pthread_mutex_init(&server->completion_lock, NULL) != 0) {
        create_status = GEO_SERVER_SYNCHRONIZATION_ERROR;
        geo_server_destroy(server);
        server = NULL;
        goto complete;
    }

    server->completion_lock_initialized = true;
    GeoDatabaseConfig database_config = geo_database_default_config();

    database_config.create_if_missing = config->create_database_if_missing;
    GeoDatabaseStatus database_status;
    server->database = geo_database_open(config->database_directory, &database_config, &database_status);

    if (!server->database) {
        create_status = GEO_SERVER_DATABASE_ERROR;
        geo_server_destroy(server);
        server = NULL;
        goto complete;
    }

    server->work_queue = geo_job_queue_create(config->worker_threads, config->work_queue_capacity);

    if (!server->work_queue) {
        create_status = GEO_SERVER_WORKER_ERROR;
        geo_server_destroy(server);
        server = NULL;
        goto complete;
    }

    server->epoll_descriptor = epoll_create1(EPOLL_CLOEXEC);
    server->wakeup_descriptor = eventfd(0, EFD_NONBLOCK | EFD_CLOEXEC);
    server->listener_descriptor = socket(AF_INET, SOCK_STREAM | SOCK_CLOEXEC, 0);

    if (server->epoll_descriptor < 0 || server->wakeup_descriptor < 0 || server->listener_descriptor < 0) {
        create_status = GEO_SERVER_DESCRIPTOR_ERROR;
        geo_server_destroy(server);
        server = NULL;
        goto complete;
    }

    int reuse_address = 1;
    struct sockaddr_in address = {
        .sin_family = AF_INET,
        .sin_port = htons(config->port),
    };

    if (inet_pton(AF_INET, config->bind_address, &address.sin_addr) != 1 ||
        setsockopt(server->listener_descriptor, SOL_SOCKET, SO_REUSEADDR, &reuse_address, sizeof(reuse_address)) != 0 ||
        !server_set_nonblocking(server->listener_descriptor) ||
        bind(server->listener_descriptor, (const struct sockaddr *) &address, sizeof(address)) != 0 ||
        listen(server->listener_descriptor, config->listen_backlog) != 0) {
        create_status = GEO_SERVER_SOCKET_ERROR;
        geo_server_destroy(server);
        server = NULL;
        goto complete;
    }

    socklen_t address_size = sizeof(address);

    if (getsockname(server->listener_descriptor, (struct sockaddr *) &address, &address_size) != 0 ||
        !server_epoll_add(server, server->listener_descriptor, EPOLLIN, &server->listener_source) ||
        !server_epoll_add(server, server->wakeup_descriptor, EPOLLIN, &server->wakeup_source)) {
        create_status = GEO_SERVER_EVENT_REGISTRATION_ERROR;
        geo_server_destroy(server);
        server = NULL;
        goto complete;
    }

    server->port = ntohs(address.sin_port);
    create_status = GEO_SERVER_OK;

complete:
    if (status) {
        *status = create_status;
    }

    return create_status == GEO_SERVER_OK ? server : NULL;
}

void geo_server_request_stop(GeoServer *server)
{
    if (!server) {
        return;
    }

    atomic_store_explicit(&server->stop_requested, true, memory_order_release);

    if (server->wakeup_descriptor >= 0) {
        uint64_t signal_value = 1U;
        ssize_t ignored = write(server->wakeup_descriptor, &signal_value, sizeof(signal_value));

        (void) ignored;
    }
}

GeoServerStatus geo_server_run(GeoServer *server)
{
    if (!server) {
        return GEO_SERVER_INVALID_ARGUMENT;
    }

    bool expected = false;

    if (!atomic_compare_exchange_strong_explicit(&server->running,
                                                  &expected,
                                                  true,
                                                  memory_order_acq_rel,
                                                  memory_order_acquire)) {
        return GEO_SERVER_ALREADY_RUNNING;
    }

    GeoServerStatus status = GEO_SERVER_OK;
    struct epoll_event events[GEO_SERVER_MAX_EVENTS];

    while (!atomic_load_explicit(&server->stop_requested, memory_order_acquire)) {
        int event_count = epoll_wait(server->epoll_descriptor, events, GEO_SERVER_MAX_EVENTS, -1);

        if (event_count < 0) {
            if (errno == EINTR) {
                continue;
            }

            status = GEO_SERVER_EVENT_LOOP_ERROR;
            break;
        }

        for (int event_index = 0; event_index < event_count; ++event_index) {
            GeoEventSource *source = events[event_index].data.ptr;

            if (source->type == GEO_EVENT_LISTENER) {
                server_accept_connections(server);
                continue;
            }

            if (source->type == GEO_EVENT_WAKEUP) {
                if (!server_process_completions(server)) {
                    status = GEO_SERVER_SYNCHRONIZATION_ERROR;
                    atomic_store_explicit(&server->stop_requested, true, memory_order_release);
                }

                continue;
            }

            GeoConnection *connection = (GeoConnection *) source;
            uint32_t ready = events[event_index].events;
            bool succeeded = true;

            if ((ready & EPOLLIN) && !connection->job_inflight && !connection->response) {
                succeeded = server_read_connection(connection);
            }

            if (succeeded && (ready & EPOLLOUT) && connection->response) {
                succeeded = server_write_connection(connection);
            }

            if (!succeeded || (ready & (EPOLLERR | EPOLLHUP | EPOLLRDHUP))) {
                server_close_connection(connection);
            }
        }
    }

    if (atomic_load_explicit(&server->event_loop_failed, memory_order_acquire)) {
        status = GEO_SERVER_SYNCHRONIZATION_ERROR;
    }

    atomic_store_explicit(&server->running, false, memory_order_release);

    return status;
}

void geo_server_destroy(GeoServer *server)
{
    if (!server) {
        return;
    }

    geo_server_request_stop(server);

    if (server->listener_descriptor >= 0) {
        (void) close(server->listener_descriptor);
    }

    for (GeoConnection *connection = server->connections; connection; connection = connection->next) {
        if (connection->descriptor >= 0) {
            (void) close(connection->descriptor);
            connection->descriptor = -1;
        }

        connection->closing = true;
    }

    geo_job_queue_destroy(server->work_queue);

    while (server->connections) {
        server_free_connection(server->connections);
    }

    if (server->wakeup_descriptor >= 0) {
        (void) close(server->wakeup_descriptor);
    }

    if (server->epoll_descriptor >= 0) {
        (void) close(server->epoll_descriptor);
    }

    geo_database_close(server->database);

    if (server->completion_lock_initialized) {
        (void) pthread_mutex_destroy(&server->completion_lock);
    }

    server_erase_secret(server->authentication_token, server->authentication_token_size);
    free(server->authentication_token);
    free(server);
}

uint16_t geo_server_port(const GeoServer *server)
{
    return server ? server->port : 0U;
}

bool geo_server_get_stats(const GeoServer *server, GeoServerStats *stats)
{
    if (!server || !stats) {
        return false;
    }

    *stats = (GeoServerStats) {
        .accepted_connections = atomic_load_explicit(&server->accepted_connections, memory_order_relaxed),
        .rejected_connections = atomic_load_explicit(&server->rejected_connections, memory_order_relaxed),
        .active_connections = atomic_load_explicit(&server->active_connections, memory_order_relaxed),
        .completed_requests = atomic_load_explicit(&server->completed_requests, memory_order_relaxed),
        .protocol_errors = atomic_load_explicit(&server->protocol_errors, memory_order_relaxed),
        .authentication_failures = atomic_load_explicit(&server->authentication_failures, memory_order_relaxed),
        .backpressure_rejections = atomic_load_explicit(&server->backpressure_rejections, memory_order_relaxed),
        .current_inflight_payload_bytes = atomic_load_explicit(&server->inflight_payload_bytes, memory_order_relaxed),
        .peak_inflight_payload_bytes = atomic_load_explicit(&server->peak_inflight_payload_bytes, memory_order_relaxed),
    };

    return true;
}
