#include "geo_server_internal.h"

#include <stdlib.h>

static bool server_token_matches(const GeoServer *server, const unsigned char *token, size_t token_size)
{
    size_t maximum_size = token_size > server->authentication_token_size ? token_size : server->authentication_token_size;
    unsigned difference = (unsigned) (token_size ^ server->authentication_token_size);

    for (size_t index = 0; index < maximum_size; ++index) {
        unsigned char received = index < token_size ? token[index] : 0;
        unsigned char expected = index < server->authentication_token_size
                                     ? (unsigned char) server->authentication_token[index]
                                     : 0;

        difference |= received ^ expected;
    }

    return difference == 0;
}

static GeoProtocolStatus server_map_database_status(GeoDatabaseStatus status)
{
    if (status == GEO_DATABASE_OK) {
        return GEO_PROTOCOL_STATUS_OK;
    }

    if (status == GEO_DATABASE_INVALID_ARGUMENT) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    return GEO_PROTOCOL_STATUS_DATABASE_ERROR;
}

static GeoProtocolStatus server_execute_write(GeoConnection *connection)
{
    const unsigned char *payload = connection->request_payload;
    uint32_t payload_size = connection->request_header.payload_size;

    if (payload_size < sizeof(uint32_t)) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    uint32_t operation_count = geo_protocol_load_u32(payload);

#if SIZE_MAX <= UINT32_MAX
    if (operation_count > (SIZE_MAX - sizeof(uint32_t)) / GEO_PROTOCOL_OPERATION_SIZE) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }
#endif

    size_t expected_size = sizeof(uint32_t) + (size_t) operation_count * GEO_PROTOCOL_OPERATION_SIZE;

    if (!operation_count || expected_size != payload_size) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    GeoDatabaseMutation *mutations = malloc((size_t) operation_count * sizeof(*mutations));

    if (!mutations) {
        return GEO_PROTOCOL_STATUS_INTERNAL_ERROR;
    }

    const unsigned char *encoded = payload + sizeof(uint32_t);

    for (uint32_t operation = 0; operation < operation_count; ++operation) {
        mutations[operation] = (GeoDatabaseMutation) {
            .object_id = geo_protocol_load_u64(encoded),
            .morton_code = geo_protocol_load_u64(encoded + 8U),
            .operation = geo_protocol_load_u32(encoded + 16U),
            .reserved = geo_protocol_load_u32(encoded + 20U),
        };
        encoded += GEO_PROTOCOL_OPERATION_SIZE;
    }

    GeoDatabaseStatus database_status = geo_database_write(connection->server->database, mutations, operation_count);

    free(mutations);

    return server_map_database_status(database_status);
}

static GeoProtocolStatus server_execute_radius_count(GeoConnection *connection,
                                                     unsigned char output[sizeof(uint64_t)],
                                                     uint32_t *output_size)
{
    if (connection->request_header.payload_size != GEO_PROTOCOL_RADIUS_SIZE) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    double latitude = geo_protocol_load_double(connection->request_payload);
    double longitude = geo_protocol_load_double(connection->request_payload + 8U);
    double radius_km = geo_protocol_load_double(connection->request_payload + 16U);
    size_t count;

    if (!geo_database_search_radius_count(connection->server->database,
                                          latitude,
                                          longitude,
                                          radius_km,
                                          &count,
                                          NULL)) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    geo_protocol_store_u64(output, (uint64_t) count);
    *output_size = sizeof(uint64_t);

    return GEO_PROTOCOL_STATUS_OK;
}

static GeoProtocolStatus server_execute_stats(GeoConnection *connection,
                                              unsigned char output[GEO_PROTOCOL_STATS_SIZE],
                                              uint32_t *output_size)
{
    if (connection->request_header.payload_size != 0) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    GeoDatabaseStats stats;

    if (!geo_database_get_stats(connection->server->database, &stats)) {
        return GEO_PROTOCOL_STATUS_DATABASE_ERROR;
    }

    geo_protocol_store_u64(output, stats.committed_operations);
    geo_protocol_store_u64(output + 8U, stats.recovered_operations);
    geo_protocol_store_u64(output + 16U, stats.checkpoints);
    geo_protocol_store_u64(output + 24U, stats.maintenance_failures);
    geo_protocol_store_u64(output + 32U, stats.next_sequence);
    geo_protocol_store_u64(output + 40U, stats.physical_records);
    geo_protocol_store_u64(output + 48U, stats.active_segments);
    geo_protocol_store_u64(output + 56U, 0U);
    *output_size = GEO_PROTOCOL_STATS_SIZE;

    return GEO_PROTOCOL_STATUS_OK;
}

void geo_server_execute_request(void *context)
{
    GeoConnection *connection = context;
    GeoServer *server = connection->server;
    unsigned char output[GEO_PROTOCOL_STATS_SIZE];
    uint32_t output_size = 0;
    GeoProtocolStatus status;

    if (!connection->authenticated && connection->request_header.opcode != GEO_PROTOCOL_AUTH) {
        status = GEO_PROTOCOL_STATUS_UNAUTHORIZED;
    } else {
        switch (connection->request_header.opcode) {
            case GEO_PROTOCOL_AUTH:
                if (server_token_matches(server,
                                         connection->request_payload,
                                         connection->request_header.payload_size)) {
                    connection->authenticated = true;
                    status = GEO_PROTOCOL_STATUS_OK;
                } else {
                    atomic_fetch_add_explicit(&server->authentication_failures, 1U, memory_order_relaxed);
                    status = GEO_PROTOCOL_STATUS_UNAUTHORIZED;
                }
                break;

            case GEO_PROTOCOL_PING:
                status = GEO_PROTOCOL_STATUS_OK;
                break;

            case GEO_PROTOCOL_WRITE:
                status = server_execute_write(connection);
                break;

            case GEO_PROTOCOL_RADIUS_COUNT:
                status = server_execute_radius_count(connection, output, &output_size);
                break;

            case GEO_PROTOCOL_STATS:
                status = server_execute_stats(connection, output, &output_size);
                break;

            case GEO_PROTOCOL_CHECKPOINT:
                status = connection->request_header.payload_size == 0
                             ? server_map_database_status(geo_database_checkpoint(server->database))
                             : GEO_PROTOCOL_STATUS_INVALID_REQUEST;
                break;

            default:
                status = GEO_PROTOCOL_STATUS_INVALID_REQUEST;
                break;
        }
    }

    const void *response_payload = output;

    if (connection->request_header.opcode == GEO_PROTOCOL_PING && status == GEO_PROTOCOL_STATUS_OK) {
        response_payload = connection->request_payload;
        output_size = connection->request_header.payload_size;
    }

    if (!geo_server_prepare_response(connection, status, response_payload, output_size)) {
        connection->close_after_response = true;
    }

    geo_server_complete_job(connection);
}
