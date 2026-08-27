#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "geobolt/client.h"

#include "geo_protocol.h"

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <netdb.h>
#include <netinet/tcp.h>
#include <poll.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <unistd.h>

#define GEO_CLIENT_DEFAULT_TIMEOUT_MS 5000U

struct GeoClient {
    int descriptor;
    pthread_mutex_t lock;
    uint64_t next_request_id;
    size_t max_frame_size;
    bool lock_initialized;
};

static void client_close_descriptor(GeoClient *client)
{
    if (client->descriptor >= 0) {
        (void) close(client->descriptor);
        client->descriptor = -1;
    }
}

static bool client_connect_with_timeout(int descriptor,
                                        const struct sockaddr *address,
                                        socklen_t address_size,
                                        int timeout_ms)
{
    int flags = fcntl(descriptor, F_GETFL, 0);

    if (flags < 0 || fcntl(descriptor, F_SETFL, flags | O_NONBLOCK) != 0) {
        return false;
    }

    if (connect(descriptor, address, address_size) != 0) {
        if (errno != EINPROGRESS) {
            return false;
        }

        struct pollfd readiness = {
            .fd = descriptor,
            .events = POLLOUT,
        };
        int poll_status;

        do {
            poll_status = poll(&readiness, 1, timeout_ms);
        } while (poll_status < 0 && errno == EINTR);

        int socket_error = 0;
        socklen_t socket_error_size = sizeof(socket_error);

        if (poll_status != 1 || !(readiness.revents & POLLOUT) ||
            getsockopt(descriptor, SOL_SOCKET, SO_ERROR, &socket_error, &socket_error_size) != 0 || socket_error != 0) {
            return false;
        }
    }

    return fcntl(descriptor, F_SETFL, flags) == 0;
}

static bool client_write_frame(int descriptor,
                               const unsigned char header[GEO_PROTOCOL_HEADER_SIZE],
                               const void *payload,
                               size_t payload_size)
{
    struct iovec vectors[2] = {
        { .iov_base = (void *) header, .iov_len = GEO_PROTOCOL_HEADER_SIZE },
        { .iov_base = (void *) payload, .iov_len = payload_size },
    };
    struct msghdr message = {
        .msg_iov = vectors,
        .msg_iovlen = payload_size ? 2U : 1U,
    };

    while (message.msg_iovlen) {
        ssize_t written = sendmsg(descriptor, &message, MSG_NOSIGNAL);

        if (written < 0 && errno == EINTR) {
            continue;
        }

        if (written <= 0) {
            return false;
        }

        size_t consumed = (size_t) written;

        while (message.msg_iovlen && consumed >= message.msg_iov[0].iov_len) {
            consumed -= message.msg_iov[0].iov_len;
            message.msg_iov++;
            message.msg_iovlen--;
        }

        if (message.msg_iovlen && consumed) {
            message.msg_iov[0].iov_base = (unsigned char *) message.msg_iov[0].iov_base + consumed;
            message.msg_iov[0].iov_len -= consumed;
        }
    }

    return true;
}

static bool client_read_all(int descriptor, void *data, size_t size)
{
    unsigned char *bytes = data;

    while (size) {
        ssize_t received = recv(descriptor, bytes, size, 0);

        if (received < 0 && errno == EINTR) {
            continue;
        }

        if (received <= 0) {
            return false;
        }

        bytes += (size_t) received;
        size -= (size_t) received;
    }

    return true;
}

static GeoClientStatus client_map_status(uint32_t status)
{
    switch (status) {
        case GEO_PROTOCOL_STATUS_OK:
            return GEO_CLIENT_OK;
        case GEO_PROTOCOL_STATUS_INVALID_REQUEST:
            return GEO_CLIENT_INVALID_ARGUMENT;
        case GEO_PROTOCOL_STATUS_UNAUTHORIZED:
            return GEO_CLIENT_UNAUTHORIZED;
        case GEO_PROTOCOL_STATUS_BUSY:
            return GEO_CLIENT_BUSY;
        case GEO_PROTOCOL_STATUS_PROTOCOL_ERROR:
        case GEO_PROTOCOL_STATUS_FRAME_TOO_LARGE:
            return GEO_CLIENT_PROTOCOL_ERROR;
        default:
            return GEO_CLIENT_SERVER_ERROR;
    }
}

static GeoClientStatus client_request(GeoClient *client,
                                      uint16_t opcode,
                                      const void *payload,
                                      uint32_t payload_size,
                                      unsigned char **response_payload,
                                      uint32_t *response_size)
{
    if (!client || client->descriptor < 0 || (!payload && payload_size) || payload_size > client->max_frame_size ||
        !response_payload || !response_size || pthread_mutex_lock(&client->lock) != 0) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    *response_payload = NULL;
    *response_size = 0;

    // One mutex intentionally spans the blocking transaction: interleaving frames would corrupt this synchronous connection stream.

    uint64_t request_id = client->next_request_id++;
    GeoProtocolHeader request = {
        .request_id = request_id,
        .payload_checksum = geo_protocol_checksum(payload, payload_size),
        .opcode = opcode,
        .payload_size = payload_size,
    };
    unsigned char encoded_header[GEO_PROTOCOL_HEADER_SIZE];

    geo_protocol_encode_header(encoded_header, &request);
    if (!client_write_frame(client->descriptor, encoded_header, payload, payload_size) ||
        !client_read_all(client->descriptor, encoded_header, sizeof(encoded_header))) {
        client_close_descriptor(client);
        pthread_mutex_unlock(&client->lock);

        return GEO_CLIENT_IO_ERROR;
    }

    GeoProtocolHeader response;

    if (!geo_protocol_decode_header(encoded_header, &response) ||
        response.flags != GEO_PROTOCOL_RESPONSE_FLAG ||
        response.opcode != opcode ||
        response.request_id != request_id ||
        response.payload_size > client->max_frame_size) {
        client_close_descriptor(client);
        pthread_mutex_unlock(&client->lock);

        return GEO_CLIENT_PROTOCOL_ERROR;
    }

    unsigned char *received_payload = response.payload_size ? malloc(response.payload_size) : NULL;

    if (response.payload_size && !received_payload) {
        client_close_descriptor(client);
        pthread_mutex_unlock(&client->lock);

        return GEO_CLIENT_OUT_OF_MEMORY;
    }

    if (response.payload_size && !client_read_all(client->descriptor, received_payload, response.payload_size)) {
        free(received_payload);
        client_close_descriptor(client);
        pthread_mutex_unlock(&client->lock);

        return GEO_CLIENT_IO_ERROR;
    }

    if (geo_protocol_checksum(received_payload, response.payload_size) != response.payload_checksum) {
        free(received_payload);
        client_close_descriptor(client);
        pthread_mutex_unlock(&client->lock);

        return GEO_CLIENT_PROTOCOL_ERROR;
    }

    GeoClientStatus status = client_map_status(response.status);

    *response_payload = received_payload;
    *response_size = response.payload_size;

    pthread_mutex_unlock(&client->lock);

    return status;
}

GeoClientConfig geo_client_default_config(const char *host, uint16_t port)
{
    return (GeoClientConfig) {
        .host = host,
        .port = port,
        .timeout_ms = GEO_CLIENT_DEFAULT_TIMEOUT_MS,
        .max_frame_size = GEO_PROTOCOL_DEFAULT_MAX_FRAME,
    };
}

GeoClient *geo_client_connect(const GeoClientConfig *config, GeoClientStatus *status)
{
    GeoClientStatus connect_status = GEO_CLIENT_INVALID_ARGUMENT;

    size_t authentication_token_size = config && config->authentication_token ? strlen(config->authentication_token) : 0;

    if (!config || !config->host || !config->host[0] || !config->port || !config->timeout_ms ||
        config->timeout_ms > INT_MAX || !config->max_frame_size ||
        config->max_frame_size > GEO_PROTOCOL_ABSOLUTE_MAX_FRAME ||
        authentication_token_size > UINT32_MAX || authentication_token_size > config->max_frame_size) {
        goto complete;
    }

    char service[6];
    int service_length = snprintf(service, sizeof(service), "%u", config->port);
    struct addrinfo hints = {
        .ai_family = AF_UNSPEC,
        .ai_socktype = SOCK_STREAM,
        .ai_protocol = IPPROTO_TCP,
    };
    struct addrinfo *addresses = NULL;

    if (service_length <= 0 || (size_t) service_length >= sizeof(service) ||
        getaddrinfo(config->host, service, &hints, &addresses) != 0) {
        connect_status = GEO_CLIENT_CONNECT_ERROR;
        goto complete;
    }

    int descriptor = -1;

    for (const struct addrinfo *address = addresses; descriptor < 0 && address; address = address->ai_next) {
        int candidate = socket(address->ai_family, address->ai_socktype | SOCK_CLOEXEC, address->ai_protocol);

        if (candidate >= 0 &&
            client_connect_with_timeout(candidate, address->ai_addr, address->ai_addrlen, (int) config->timeout_ms)) {
            descriptor = candidate;
        } else if (candidate >= 0) {
            (void) close(candidate);
        }
    }

    freeaddrinfo(addresses);

    if (descriptor < 0) {
        connect_status = GEO_CLIENT_CONNECT_ERROR;
        goto complete;
    }

    struct timeval timeout = {
        .tv_sec = config->timeout_ms / 1000U,
        .tv_usec = (config->timeout_ms % 1000U) * 1000U,
    };
    int enabled = 1;

    if (setsockopt(descriptor, SOL_SOCKET, SO_RCVTIMEO, &timeout, sizeof(timeout)) != 0 ||
        setsockopt(descriptor, SOL_SOCKET, SO_SNDTIMEO, &timeout, sizeof(timeout)) != 0 ||
        setsockopt(descriptor, IPPROTO_TCP, TCP_NODELAY, &enabled, sizeof(enabled)) != 0) {
        (void) close(descriptor);
        connect_status = GEO_CLIENT_CONNECT_ERROR;
        goto complete;
    }

    GeoClient *client = calloc(1, sizeof(*client));

    if (!client) {
        (void) close(descriptor);
        connect_status = GEO_CLIENT_OUT_OF_MEMORY;
        goto complete;
    }

    client->descriptor = descriptor;
    client->next_request_id = 1U;
    client->max_frame_size = config->max_frame_size;

    if (pthread_mutex_init(&client->lock, NULL) != 0) {
        geo_client_close(client);
        client = NULL;
        connect_status = GEO_CLIENT_CONNECT_ERROR;
        goto complete;
    }

    client->lock_initialized = true;

    unsigned char *response = NULL;
    uint32_t response_size = 0;

    connect_status = client_request(client,
                                    GEO_PROTOCOL_AUTH,
                                    config->authentication_token,
                                    (uint32_t) authentication_token_size,
                                    &response,
                                    &response_size);
    free(response);

    if (connect_status != GEO_CLIENT_OK || response_size != 0) {
        if (connect_status == GEO_CLIENT_OK) {
            connect_status = GEO_CLIENT_PROTOCOL_ERROR;
        }

        geo_client_close(client);
        client = NULL;
    }

complete:
    if (status) {
        *status = connect_status;
    }

    return connect_status == GEO_CLIENT_OK ? client : NULL;
}

void geo_client_close(GeoClient *client)
{
    if (!client) {
        return;
    }

    client_close_descriptor(client);

    if (client->lock_initialized) {
        (void) pthread_mutex_destroy(&client->lock);
    }

    free(client);
}

GeoClientStatus geo_client_ping(GeoClient *client, const void *payload, size_t payload_size)
{
    if (payload_size > UINT32_MAX) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char *response = NULL;
    uint32_t response_size = 0;
    GeoClientStatus status = client_request(client, GEO_PROTOCOL_PING, payload, (uint32_t) payload_size, &response, &response_size);

    if (status == GEO_CLIENT_OK && (response_size != payload_size || (payload_size && memcmp(response, payload, payload_size) != 0))) {
        status = GEO_CLIENT_PROTOCOL_ERROR;
    }

    free(response);

    return status;
}

GeoClientStatus geo_client_write(GeoClient *client, const GeoDatabaseMutation *mutations, size_t mutation_count)
{
    if (!mutations || !mutation_count || mutation_count > (UINT32_MAX - sizeof(uint32_t)) / GEO_PROTOCOL_OPERATION_SIZE) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    size_t payload_size = sizeof(uint32_t) + mutation_count * GEO_PROTOCOL_OPERATION_SIZE;

    if (!client || payload_size > client->max_frame_size) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char *payload = malloc(payload_size);

    if (!payload) {
        return GEO_CLIENT_OUT_OF_MEMORY;
    }

    geo_protocol_store_u32(payload, (uint32_t) mutation_count);
    unsigned char *encoded = payload + sizeof(uint32_t);

    for (size_t operation = 0; operation < mutation_count; ++operation) {
        geo_protocol_store_u64(encoded, mutations[operation].object_id);
        geo_protocol_store_u64(encoded + 8U, mutations[operation].morton_code);
        geo_protocol_store_u32(encoded + 16U, mutations[operation].operation);
        geo_protocol_store_u32(encoded + 20U, mutations[operation].reserved);
        encoded += GEO_PROTOCOL_OPERATION_SIZE;
    }

    unsigned char *response = NULL;
    uint32_t response_size = 0;
    GeoClientStatus status = client_request(client,
                                            GEO_PROTOCOL_WRITE,
                                            payload,
                                            (uint32_t) payload_size,
                                            &response,
                                            &response_size);

    if (status == GEO_CLIENT_OK && response_size != 0) {
        status = GEO_CLIENT_PROTOCOL_ERROR;
    }

    free(response);
    free(payload);

    return status;
}

GeoClientStatus geo_client_radius_count(GeoClient *client,
                                        double latitude,
                                        double longitude,
                                        double radius_km,
                                        uint64_t *count)
{
    if (!count) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char payload[GEO_PROTOCOL_RADIUS_SIZE];

    geo_protocol_store_double(payload, latitude);
    geo_protocol_store_double(payload + 8U, longitude);
    geo_protocol_store_double(payload + 16U, radius_km);

    unsigned char *response = NULL;
    uint32_t response_size = 0;
    GeoClientStatus status = client_request(client, GEO_PROTOCOL_RADIUS_COUNT, payload, sizeof(payload), &response, &response_size);

    if (status == GEO_CLIENT_OK && response_size == sizeof(uint64_t)) {
        *count = geo_protocol_load_u64(response);
    } else if (status == GEO_CLIENT_OK) {
        status = GEO_CLIENT_PROTOCOL_ERROR;
    }

    free(response);

    return status;
}

GeoClientStatus geo_client_get_database_stats(GeoClient *client, GeoDatabaseStats *stats)
{
    if (!stats) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char *response = NULL;
    uint32_t response_size = 0;
    GeoClientStatus status = client_request(client, GEO_PROTOCOL_STATS, NULL, 0, &response, &response_size);

    if (status == GEO_CLIENT_OK && response_size == GEO_PROTOCOL_STATS_SIZE) {
        stats->committed_operations = geo_protocol_load_u64(response);
        stats->recovered_operations = geo_protocol_load_u64(response + 8U);
        stats->checkpoints = geo_protocol_load_u64(response + 16U);
        stats->maintenance_failures = geo_protocol_load_u64(response + 24U);
        stats->next_sequence = geo_protocol_load_u64(response + 32U);
        stats->physical_records = geo_protocol_load_u64(response + 40U);
        uint64_t active_segments = geo_protocol_load_u64(response + 48U);

        if (active_segments > SIZE_MAX) {
            status = GEO_CLIENT_PROTOCOL_ERROR;
        } else {
            stats->active_segments = (size_t) active_segments;
        }
    } else if (status == GEO_CLIENT_OK) {
        status = GEO_CLIENT_PROTOCOL_ERROR;
    }

    free(response);

    return status;
}

GeoClientStatus geo_client_checkpoint(GeoClient *client)
{
    unsigned char *response = NULL;
    uint32_t response_size = 0;
    GeoClientStatus status = client_request(client, GEO_PROTOCOL_CHECKPOINT, NULL, 0, &response, &response_size);

    if (status == GEO_CLIENT_OK && response_size != 0) {
        status = GEO_CLIENT_PROTOCOL_ERROR;
    }

    free(response);

    return status;
}
