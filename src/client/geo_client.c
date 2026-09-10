#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "geobolt/client.h"

#include "geo_protocol.h"

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>
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
        ssize_t received;

#if defined(__clang_analyzer__)
        // The synchronous protocol requires exclusive stream ownership across the bounded receive. Suppress only this call.
        __attribute__((suppress))
#endif
        {
            received = recv(descriptor, bytes, size, 0);
        }

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
        case GEO_PROTOCOL_STATUS_NOT_FOUND:
            return GEO_CLIENT_NOT_FOUND;
        case GEO_PROTOCOL_STATUS_ALREADY_EXISTS:
            return GEO_CLIENT_ALREADY_EXISTS;
        case GEO_PROTOCOL_STATUS_IDEMPOTENCY_CONFLICT:
            return GEO_CLIENT_IDEMPOTENCY_CONFLICT;
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
    if (!client || (!payload && payload_size) || payload_size > client->max_frame_size ||
        !response_payload || !response_size || pthread_mutex_lock(&client->lock) != 0) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    *response_payload = NULL;
    *response_size = 0;

    /* Connection failures change the descriptor under this same mutex. Checking it before
     * locking races another caller closing the stream after an I/O or protocol error. */
    if (client->descriptor < 0) {
        pthread_mutex_unlock(&client->lock);
        return GEO_CLIENT_IO_ERROR;
    }

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

static GeoClientStatus client_encode_object_mutations(GeoClient *client,
                                                      const GeoDatabaseObjectMutation *mutations,
                                                      size_t mutation_count,
                                                      unsigned char **payload,
                                                      size_t *payload_size)
{
    if (!client || !mutations || !mutation_count || mutation_count > UINT32_MAX || !payload || !payload_size) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    size_t encoded_size = sizeof(uint32_t);

    for (size_t index = 0; index < mutation_count; ++index) {
        bool operation_valid = mutations[index].operation >= GEO_DATABASE_UPSERT &&
                               mutations[index].operation <= GEO_DATABASE_UPDATE;

        if (!operation_valid || mutations[index].document_size > UINT32_MAX ||
            (mutations[index].document_size && !mutations[index].document) ||
            encoded_size > SIZE_MAX - GEO_PROTOCOL_OBJECT_OPERATION_SIZE ||
            mutations[index].document_size > SIZE_MAX - encoded_size - GEO_PROTOCOL_OBJECT_OPERATION_SIZE) {
            return GEO_CLIENT_INVALID_ARGUMENT;
        }

        encoded_size += GEO_PROTOCOL_OBJECT_OPERATION_SIZE + mutations[index].document_size;
    }

    if (encoded_size > UINT32_MAX || encoded_size > client->max_frame_size) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char *encoded_payload = malloc(encoded_size);

    if (!encoded_payload) {
        return GEO_CLIENT_OUT_OF_MEMORY;
    }

    geo_protocol_store_u32(encoded_payload, (uint32_t) mutation_count);
    unsigned char *encoded = encoded_payload + sizeof(uint32_t);

    for (size_t index = 0; index < mutation_count; ++index) {
        geo_protocol_store_u64(encoded, mutations[index].object_id);
        geo_protocol_store_u64(encoded + 8U, mutations[index].morton_code);
        geo_protocol_store_u32(encoded + 16U, mutations[index].operation);
        geo_protocol_store_u32(encoded + 20U, mutations[index].reserved);
        geo_protocol_store_u32(encoded + 24U, (uint32_t) mutations[index].document_size);
        geo_protocol_store_u32(encoded + 28U, 0U);

        if (mutations[index].document_size) {
            memcpy(encoded + GEO_PROTOCOL_OBJECT_OPERATION_SIZE,
                   mutations[index].document,
                   mutations[index].document_size);
        }

        encoded += GEO_PROTOCOL_OBJECT_OPERATION_SIZE + mutations[index].document_size;
    }

    *payload = encoded_payload;
    *payload_size = encoded_size;
    return GEO_CLIENT_OK;
}

GeoClientStatus geo_client_write_objects(GeoClient *client,
                                         const GeoDatabaseObjectMutation *mutations,
                                         size_t mutation_count)
{
    unsigned char *payload;
    size_t payload_size;
    GeoClientStatus status = client_encode_object_mutations(client, mutations, mutation_count, &payload, &payload_size);

    if (status != GEO_CLIENT_OK) {
        return status;
    }

    unsigned char *response = NULL;
    uint32_t response_size = 0U;
    status = client_request(client,
                            GEO_PROTOCOL_WRITE_OBJECTS,
                            payload,
                            (uint32_t) payload_size,
                            &response,
                            &response_size);

    if (status == GEO_CLIENT_OK && response_size != 0U) {
        status = GEO_CLIENT_PROTOCOL_ERROR;
    }

    free(response);
    free(payload);
    return status;
}

GeoClientStatus geo_client_write_objects_idempotent(GeoClient *client,
                                                    const GeoDatabaseObjectMutation *mutations,
                                                    size_t mutation_count,
                                                    const void *idempotency_key,
                                                    size_t idempotency_key_size,
                                                    uint32_t retention_seconds,
                                                    bool *replayed)
{
    if (!idempotency_key || !idempotency_key_size || idempotency_key_size > GEO_DATABASE_MAX_IDEMPOTENCY_KEY_SIZE ||
        !retention_seconds || retention_seconds > GEO_DATABASE_MAX_IDEMPOTENCY_RETENTION_SECONDS || !replayed) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    *replayed = false;

    unsigned char *mutation_payload;
    size_t mutation_payload_size;
    GeoClientStatus status = client_encode_object_mutations(client,
                                                            mutations,
                                                            mutation_count,
                                                            &mutation_payload,
                                                            &mutation_payload_size);

    if (status != GEO_CLIENT_OK) {
        return status;
    }

    if (idempotency_key_size > SIZE_MAX - 8U || mutation_payload_size > SIZE_MAX - 8U - idempotency_key_size) {
        free(mutation_payload);
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    size_t payload_size = 8U + idempotency_key_size + mutation_payload_size;

    if (payload_size > UINT32_MAX || payload_size > client->max_frame_size) {
        free(mutation_payload);
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char *payload = malloc(payload_size);

    if (!payload) {
        free(mutation_payload);
        return GEO_CLIENT_OUT_OF_MEMORY;
    }

    geo_protocol_store_u32(payload, (uint32_t) idempotency_key_size);
    geo_protocol_store_u32(payload + 4U, retention_seconds);
    memcpy(payload + 8U, idempotency_key, idempotency_key_size);
    memcpy(payload + 8U + idempotency_key_size, mutation_payload, mutation_payload_size);
    free(mutation_payload);

    unsigned char *response = NULL;
    uint32_t response_size = 0U;
    status = client_request(client,
                            GEO_PROTOCOL_WRITE_OBJECTS_IDEMPOTENT,
                            payload,
                            (uint32_t) payload_size,
                            &response,
                            &response_size);

    if (status == GEO_CLIENT_OK && response_size == sizeof(uint64_t)) {
        *replayed = geo_protocol_load_u64(response) != 0U;
    } else if (status == GEO_CLIENT_OK) {
        status = GEO_CLIENT_PROTOCOL_ERROR;
    }

    free(response);
    free(payload);
    return status;
}

static GeoClientStatus client_write_one_object(GeoClient *client,
                                               uint64_t object_id,
                                               double latitude,
                                               double longitude,
                                               const void *document,
                                               size_t document_size,
                                               GeoDatabaseOperation operation)
{
    if (!isfinite(latitude) || !isfinite(longitude) || latitude < GEO_MIN_LAT || latitude > GEO_MAX_LAT ||
        longitude < GEO_MIN_LNG || longitude > GEO_MAX_LNG) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    GeoDatabaseObjectMutation mutation = {
        .object_id = object_id,
        .morton_code = geo_encode(latitude, longitude),
        .document = document,
        .document_size = document_size,
        .operation = operation,
    };

    return geo_client_write_objects(client, &mutation, 1U);
}

GeoClientStatus geo_client_insert(GeoClient *client,
                                  uint64_t object_id,
                                  double latitude,
                                  double longitude,
                                  const void *document,
                                  size_t document_size)
{
    return client_write_one_object(client,
                                   object_id,
                                   latitude,
                                   longitude,
                                   document,
                                   document_size,
                                   GEO_DATABASE_INSERT);
}

GeoClientStatus geo_client_insert_generated(GeoClient *client,
                                            double latitude,
                                            double longitude,
                                            const void *document,
                                            size_t document_size,
                                            uint64_t *object_id)
{
    if (!client || !object_id || !isfinite(latitude) || !isfinite(longitude) || latitude < GEO_MIN_LAT ||
        latitude > GEO_MAX_LAT || longitude < GEO_MIN_LNG || longitude > GEO_MAX_LNG ||
        (document_size && !document) || document_size > UINT32_MAX ||
        document_size > SIZE_MAX - GEO_PROTOCOL_GENERATED_INSERT_SIZE) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    size_t payload_size = GEO_PROTOCOL_GENERATED_INSERT_SIZE + document_size;

    if (payload_size > client->max_frame_size) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char *payload = malloc(payload_size);

    if (!payload) {
        return GEO_CLIENT_OUT_OF_MEMORY;
    }

    geo_protocol_store_u64(payload, geo_encode(latitude, longitude));
    geo_protocol_store_u32(payload + 8U, (uint32_t) document_size);
    geo_protocol_store_u32(payload + 12U, 0U);

    if (document_size) {
        memcpy(payload + GEO_PROTOCOL_GENERATED_INSERT_SIZE, document, document_size);
    }

    unsigned char *response = NULL;
    uint32_t response_size = 0U;
    GeoClientStatus status = client_request(client,
                                            GEO_PROTOCOL_INSERT_GENERATED,
                                            payload,
                                            (uint32_t) payload_size,
                                            &response,
                                            &response_size);

    if (status == GEO_CLIENT_OK && response_size == sizeof(uint64_t)) {
        *object_id = geo_protocol_load_u64(response);
    } else if (status == GEO_CLIENT_OK) {
        status = GEO_CLIENT_PROTOCOL_ERROR;
    }

    free(response);
    free(payload);
    return status;
}

GeoClientStatus geo_client_upsert(GeoClient *client,
                                  uint64_t object_id,
                                  double latitude,
                                  double longitude,
                                  const void *document,
                                  size_t document_size)
{
    return client_write_one_object(client,
                                   object_id,
                                   latitude,
                                   longitude,
                                   document,
                                   document_size,
                                   GEO_DATABASE_UPSERT);
}

GeoClientStatus geo_client_update(GeoClient *client,
                                  uint64_t object_id,
                                  double latitude,
                                  double longitude,
                                  const void *document,
                                  size_t document_size)
{
    return client_write_one_object(client,
                                   object_id,
                                   latitude,
                                   longitude,
                                   document,
                                   document_size,
                                   GEO_DATABASE_UPDATE);
}

GeoClientStatus geo_client_delete(GeoClient *client, uint64_t object_id)
{
    GeoDatabaseObjectMutation mutation = {
        .object_id = object_id,
        .operation = GEO_DATABASE_DELETE,
    };

    return geo_client_write_objects(client, &mutation, 1U);
}

GeoClientStatus geo_client_get(GeoClient *client, uint64_t object_id, GeoDatabaseObject *object)
{
    if (!client || !object || !object_id) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    *object = (GeoDatabaseObject) { 0 };
    unsigned char request[sizeof(uint64_t)];
    unsigned char *response = NULL;
    uint32_t response_size = 0U;

    geo_protocol_store_u64(request, object_id);
    GeoClientStatus status = client_request(client,
                                            GEO_PROTOCOL_GET_OBJECT,
                                            request,
                                            sizeof(request),
                                            &response,
                                            &response_size);

    if (status != GEO_CLIENT_OK) {
        free(response);
        return status;
    }

    if (response_size < GEO_PROTOCOL_OBJECT_RESPONSE_SIZE) {
        free(response);
        return GEO_CLIENT_PROTOCOL_ERROR;
    }

    uint32_t document_size = geo_protocol_load_u32(response + 24U);

    if (geo_protocol_load_u32(response + 28U) != 0U ||
        document_size != response_size - GEO_PROTOCOL_OBJECT_RESPONSE_SIZE) {
        free(response);
        return GEO_CLIENT_PROTOCOL_ERROR;
    }

    void *document = document_size ? malloc(document_size) : NULL;

    if (document_size && !document) {
        free(response);
        return GEO_CLIENT_OUT_OF_MEMORY;
    }

    if (document_size) {
        memcpy(document, response + GEO_PROTOCOL_OBJECT_RESPONSE_SIZE, document_size);
    }

    *object = (GeoDatabaseObject) {
        .object_id = geo_protocol_load_u64(response),
        .sequence = geo_protocol_load_u64(response + 8U),
        .morton_code = geo_protocol_load_u64(response + 16U),
        .document = document,
        .document_size = document_size,
    };
    free(response);

    if (object->object_id != object_id || !object->sequence) {
        geo_database_object_release(object);
        return GEO_CLIENT_PROTOCOL_ERROR;
    }

    return GEO_CLIENT_OK;
}

GeoClientStatus geo_client_create_index(GeoClient *client,
                                        const char *name,
                                        const char *json_pointer,
                                        GeoDatabaseIndexType type)
{
    if (!client || !name || !name[0] || !json_pointer || type < GEO_DATABASE_INDEX_BOOL || type > GEO_DATABASE_INDEX_BYTES) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    size_t name_size = strlen(name);
    size_t pointer_size = strlen(json_pointer);

    if (name_size > GEO_DATABASE_MAX_INDEX_NAME_SIZE || pointer_size > GEO_DATABASE_MAX_INDEX_POINTER_SIZE ||
        name_size > SIZE_MAX - GEO_PROTOCOL_INDEX_CREATE_HEADER_SIZE ||
        pointer_size > SIZE_MAX - GEO_PROTOCOL_INDEX_CREATE_HEADER_SIZE - name_size) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    size_t payload_size = GEO_PROTOCOL_INDEX_CREATE_HEADER_SIZE + name_size + pointer_size;

    if (payload_size > client->max_frame_size) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char *payload = malloc(payload_size);

    if (!payload) {
        return GEO_CLIENT_OUT_OF_MEMORY;
    }

    geo_protocol_store_u32(payload, (uint32_t) type);
    geo_protocol_store_u32(payload + 4U, (uint32_t) name_size);
    geo_protocol_store_u32(payload + 8U, (uint32_t) pointer_size);
    geo_protocol_store_u32(payload + 12U, 0U);
    memcpy(payload + GEO_PROTOCOL_INDEX_CREATE_HEADER_SIZE, name, name_size);
    memcpy(payload + GEO_PROTOCOL_INDEX_CREATE_HEADER_SIZE + name_size, json_pointer, pointer_size);

    unsigned char *response = NULL;
    uint32_t response_size = 0U;
    GeoClientStatus status = client_request(client,
                                            GEO_PROTOCOL_CREATE_INDEX,
                                            payload,
                                            (uint32_t) payload_size,
                                            &response,
                                            &response_size);

    if (status == GEO_CLIENT_OK && response_size != 0U) {
        status = GEO_CLIENT_PROTOCOL_ERROR;
    }

    free(response);
    free(payload);
    return status;
}

GeoClientStatus geo_client_drop_index(GeoClient *client, const char *name)
{
    if (!client || !name || !name[0]) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    size_t name_size = strlen(name);

    if (name_size > GEO_DATABASE_MAX_INDEX_NAME_SIZE) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char payload[GEO_PROTOCOL_INDEX_DROP_HEADER_SIZE + GEO_DATABASE_MAX_INDEX_NAME_SIZE];

    geo_protocol_store_u32(payload, (uint32_t) name_size);
    geo_protocol_store_u32(payload + 4U, 0U);
    memcpy(payload + GEO_PROTOCOL_INDEX_DROP_HEADER_SIZE, name, name_size);

    unsigned char *response = NULL;
    uint32_t response_size = 0U;
    GeoClientStatus status = client_request(client,
                                            GEO_PROTOCOL_DROP_INDEX,
                                            payload,
                                            (uint32_t) (GEO_PROTOCOL_INDEX_DROP_HEADER_SIZE + name_size),
                                            &response,
                                            &response_size);

    if (status == GEO_CLIENT_OK && response_size != 0U) {
        status = GEO_CLIENT_PROTOCOL_ERROR;
    }

    free(response);
    return status;
}

GeoClientStatus geo_client_list_indexes(GeoClient *client,
                                       GeoDatabaseIndexInfo *indexes,
                                       size_t capacity,
                                       size_t *count)
{
    if (!client || !count || (capacity && !indexes)) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char *response = NULL;
    uint32_t response_size = 0U;
    GeoClientStatus status = client_request(client, GEO_PROTOCOL_LIST_INDEXES, NULL, 0U, &response, &response_size);

    if (status != GEO_CLIENT_OK) {
        free(response);
        return status;
    }
    if (response_size < 8U || geo_protocol_load_u32(response + 4U) != 0U) {
        free(response);
        return GEO_CLIENT_PROTOCOL_ERROR;
    }

    uint32_t encoded_count = geo_protocol_load_u32(response);
    const unsigned char *cursor = response + 8U;
    size_t remaining = response_size - 8U;
    size_t decoded_count = 0U;

    while (decoded_count < encoded_count && remaining >= 32U) {
        uint32_t type = geo_protocol_load_u32(cursor + 16U);
        uint32_t name_size = geo_protocol_load_u32(cursor + 20U);
        uint32_t pointer_size = geo_protocol_load_u32(cursor + 24U);
        size_t entry_size = 32U + (size_t) name_size + pointer_size;

        if (geo_protocol_load_u32(cursor + 28U) != 0U || type < GEO_DATABASE_INDEX_BOOL || type > GEO_DATABASE_INDEX_BYTES ||
            name_size == 0U || name_size > GEO_DATABASE_MAX_INDEX_NAME_SIZE ||
            pointer_size > GEO_DATABASE_MAX_INDEX_POINTER_SIZE || entry_size > remaining) {
            status = GEO_CLIENT_PROTOCOL_ERROR;
            break;
        }

        if (decoded_count < capacity) {
            GeoDatabaseIndexInfo *info = indexes + decoded_count;

            *info = (GeoDatabaseIndexInfo) {
                .index_id = geo_protocol_load_u64(cursor),
                .entry_count = geo_protocol_load_u64(cursor + 8U),
                .type = (GeoDatabaseIndexType) type,
            };
            memcpy(info->name, cursor + 32U, name_size);
            memcpy(info->json_pointer, cursor + 32U + name_size, pointer_size);
        }

        cursor += entry_size;
        remaining -= entry_size;
        decoded_count++;
    }

    if (status == GEO_CLIENT_OK && (decoded_count != encoded_count || remaining != 0U)) {
        status = GEO_CLIENT_PROTOCOL_ERROR;
    }

    *count = encoded_count;
    free(response);

    if (status == GEO_CLIENT_OK && encoded_count > capacity) {
        status = GEO_CLIENT_OUT_OF_MEMORY;
    }
    return status;
}

GeoClientStatus geo_client_query_index_reuse(GeoClient *client,
                                             const GeoDatabaseIndexPredicate *predicate,
                                             GeoIdResult *result,
                                             GeoDatabaseIndexQueryStats *stats)
{
    if (!client || !result) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    size_t payload_size = geo_protocol_index_predicate_size(predicate);

    if (!payload_size || payload_size > client->max_frame_size) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char *payload = malloc(payload_size);

    if (!payload) {
        return GEO_CLIENT_OUT_OF_MEMORY;
    }

    if (!geo_protocol_encode_index_predicate(payload, payload_size, predicate)) {
        free(payload);
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char *response = NULL;
    uint32_t response_size = 0U;
    GeoClientStatus status = client_request(client,
                                            GEO_PROTOCOL_QUERY_INDEX,
                                            payload,
                                            (uint32_t) payload_size,
                                            &response,
                                            &response_size);
    free(payload);

    if (status != GEO_CLIENT_OK) {
        free(response);
        return status;
    }
    if (response_size < GEO_PROTOCOL_INDEX_QUERY_RESPONSE_SIZE || geo_protocol_load_u32(response + 20U) != 0U) {
        free(response);
        return GEO_CLIENT_PROTOCOL_ERROR;
    }

    uint32_t result_count = geo_protocol_load_u32(response + 16U);
    size_t id_bytes = response_size - GEO_PROTOCOL_INDEX_QUERY_RESPONSE_SIZE;

    if (id_bytes % sizeof(uint64_t) != 0U || result_count != id_bytes / sizeof(uint64_t)) {
        free(response);
        return GEO_CLIENT_PROTOCOL_ERROR;
    }

    if (!geo_id_result_reserve(result, result_count)) {
        free(response);
        return GEO_CLIENT_OUT_OF_MEMORY;
    }

    for (size_t index = 0U; index < result_count; ++index) {
        result->ids[index] = geo_protocol_load_u64(response + GEO_PROTOCOL_INDEX_QUERY_RESPONSE_SIZE + index * sizeof(uint64_t));
    }
    result->count = result_count;

    if (stats) {
        stats->scanned_entries = geo_protocol_load_u64(response);
        stats->matched_entries = geo_protocol_load_u64(response + 8U);
    }

    free(response);
    return GEO_CLIENT_OK;
}

GeoClientStatus geo_client_query_radius_reuse(GeoClient *client,
                                              const GeoDatabaseRadiusQuery *query,
                                              GeoSearchResult *result,
                                              GeoDatabaseQueryStats *stats)
{
    if (!client || !query || !result || !query->predicates || !query->predicate_count ||
        query->predicate_count > GEO_DATABASE_MAX_QUERY_PREDICATES ||
        !geo_is_valid_point(query->latitude, query->longitude) || !isfinite(query->radius_km) || query->radius_km < 0.0) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    size_t payload_size = GEO_PROTOCOL_RADIUS_QUERY_HEADER_SIZE;

    for (size_t index = 0U; index < query->predicate_count; ++index) {
        size_t predicate_size = geo_protocol_index_predicate_size(query->predicates + index);

        if (!predicate_size || predicate_size > SIZE_MAX - payload_size) {
            return GEO_CLIENT_INVALID_ARGUMENT;
        }

        payload_size += predicate_size;
    }

    if (payload_size > client->max_frame_size || payload_size > UINT32_MAX) {
        return GEO_CLIENT_INVALID_ARGUMENT;
    }

    unsigned char *payload = malloc(payload_size);

    if (!payload) {
        return GEO_CLIENT_OUT_OF_MEMORY;
    }

    geo_protocol_store_double(payload, query->latitude);
    geo_protocol_store_double(payload + 8U, query->longitude);
    geo_protocol_store_double(payload + 16U, query->radius_km);
    geo_protocol_store_u32(payload + 24U, (uint32_t) query->predicate_count);
    geo_protocol_store_u32(payload + 28U, 0U);

    unsigned char *cursor = payload + GEO_PROTOCOL_RADIUS_QUERY_HEADER_SIZE;

    for (size_t index = 0U; index < query->predicate_count; ++index) {
        size_t predicate_size = geo_protocol_index_predicate_size(query->predicates + index);

        if (!geo_protocol_encode_index_predicate(cursor, predicate_size, query->predicates + index)) {
            free(payload);
            return GEO_CLIENT_INVALID_ARGUMENT;
        }

        cursor += predicate_size;
    }

    unsigned char *response = NULL;
    uint32_t response_size = 0U;
    GeoClientStatus status = client_request(client,
                                            GEO_PROTOCOL_QUERY_RADIUS,
                                            payload,
                                            (uint32_t) payload_size,
                                            &response,
                                            &response_size);
    free(payload);

    if (status != GEO_CLIENT_OK) {
        free(response);
        return status;
    }

    if (response_size < GEO_PROTOCOL_RADIUS_QUERY_RESPONSE_SIZE || geo_protocol_load_u32(response + 76U) != 0U) {
        free(response);
        return GEO_CLIENT_PROTOCOL_ERROR;
    }

    uint32_t plan = geo_protocol_load_u32(response + 64U);
    uint32_t predicate_count = geo_protocol_load_u32(response + 68U);
    uint32_t result_count = geo_protocol_load_u32(response + 72U);
    size_t record_bytes = response_size - GEO_PROTOCOL_RADIUS_QUERY_RESPONSE_SIZE;

    if (plan < GEO_DATABASE_QUERY_PLAN_SPATIAL || plan > GEO_DATABASE_QUERY_PLAN_SECONDARY ||
        predicate_count != query->predicate_count || record_bytes % sizeof(GeoRecord) != 0U ||
        result_count != record_bytes / sizeof(GeoRecord)) {
        free(response);
        return GEO_CLIENT_PROTOCOL_ERROR;
    }

    if (!geo_result_reserve(result, result_count)) {
        free(response);
        return GEO_CLIENT_OUT_OF_MEMORY;
    }

    for (size_t index = 0U; index < result_count; ++index) {
        const unsigned char *encoded = response + GEO_PROTOCOL_RADIUS_QUERY_RESPONSE_SIZE + index * sizeof(GeoRecord);
        uint64_t object_id = geo_protocol_load_u64(encoded);

        if (!object_id) {
            geo_result_clear(result);
            free(response);
            return GEO_CLIENT_PROTOCOL_ERROR;
        }

        result->results[index] = (GeoRecord) {
            .id = object_id,
            .z = geo_protocol_load_u64(encoded + 8U),
        };
    }
    result->count = result_count;

    if (stats) {
        *stats = (GeoDatabaseQueryStats) {
            .spatial = {
                .records_scanned = geo_protocol_load_u64(response + 32U),
                .records_matched = geo_protocol_load_u64(response + 40U),
                .ranges_checked = geo_protocol_load_u64(response + 48U),
                .search_time_ms = geo_protocol_load_double(response + 56U),
            },
            .secondary_entries_scanned = geo_protocol_load_u64(response),
            .metadata_candidates = geo_protocol_load_u64(response + 8U),
            .estimated_spatial_candidates = geo_protocol_load_u64(response + 16U),
            .object_lookups = geo_protocol_load_u64(response + 24U),
            .plan = (GeoDatabaseQueryPlan) plan,
            .predicate_count = predicate_count,
        };
    }

    free(response);
    return GEO_CLIENT_OK;
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
