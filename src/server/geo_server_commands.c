#include "geo_server_internal.h"

#include <stdlib.h>
#include <string.h>

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

    if (status == GEO_DATABASE_NOT_FOUND) {
        return GEO_PROTOCOL_STATUS_NOT_FOUND;
    }

    if (status == GEO_DATABASE_ALREADY_EXISTS) {
        return GEO_PROTOCOL_STATUS_ALREADY_EXISTS;
    }

    if (status == GEO_DATABASE_IDEMPOTENCY_CONFLICT) {
        return GEO_PROTOCOL_STATUS_IDEMPOTENCY_CONFLICT;
    }

    return GEO_PROTOCOL_STATUS_DATABASE_ERROR;
}

static GeoProtocolStatus server_decode_object_mutations(const unsigned char *payload,
                                                        size_t payload_size,
                                                        GeoDatabaseObjectMutation **decoded_mutations,
                                                        uint32_t *decoded_count)
{
    if (!payload || payload_size < sizeof(uint32_t) || !decoded_mutations || !decoded_count) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    uint32_t operation_count = geo_protocol_load_u32(payload);

    if (!operation_count || operation_count > (payload_size - sizeof(uint32_t)) / GEO_PROTOCOL_OBJECT_OPERATION_SIZE) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

#if SIZE_MAX <= UINT32_MAX
    if (operation_count > SIZE_MAX / sizeof(GeoDatabaseObjectMutation)) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }
#endif

    GeoDatabaseObjectMutation *mutations = malloc((size_t) operation_count * sizeof(*mutations));

    if (!mutations) {
        return GEO_PROTOCOL_STATUS_INTERNAL_ERROR;
    }

    const unsigned char *encoded = payload + sizeof(uint32_t);
    size_t remaining = payload_size - sizeof(uint32_t);
    bool valid = true;

    for (uint32_t operation = 0U; valid && operation < operation_count; ++operation) {
        if (remaining < GEO_PROTOCOL_OBJECT_OPERATION_SIZE) {
            valid = false;
            break;
        }

        uint32_t document_size = geo_protocol_load_u32(encoded + 24U);

        if (geo_protocol_load_u32(encoded + 28U) != 0U ||
            document_size > remaining - GEO_PROTOCOL_OBJECT_OPERATION_SIZE) {
            valid = false;
            break;
        }

        mutations[operation] = (GeoDatabaseObjectMutation) {
            .object_id = geo_protocol_load_u64(encoded),
            .morton_code = geo_protocol_load_u64(encoded + 8U),
            .document = document_size ? encoded + GEO_PROTOCOL_OBJECT_OPERATION_SIZE : NULL,
            .document_size = document_size,
            .operation = geo_protocol_load_u32(encoded + 16U),
            .reserved = geo_protocol_load_u32(encoded + 20U),
        };
        size_t encoded_size = GEO_PROTOCOL_OBJECT_OPERATION_SIZE + document_size;
        encoded += encoded_size;
        remaining -= encoded_size;
    }

    if (!valid || remaining != 0U) {
        free(mutations);
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    *decoded_mutations = mutations;
    *decoded_count = operation_count;
    return GEO_PROTOCOL_STATUS_OK;
}

static GeoProtocolStatus server_execute_write_objects(GeoConnection *connection)
{
    GeoDatabaseObjectMutation *mutations;
    uint32_t operation_count;
    GeoProtocolStatus status = server_decode_object_mutations(connection->request_payload,
                                                              connection->request_header.payload_size,
                                                              &mutations,
                                                              &operation_count);

    if (status != GEO_PROTOCOL_STATUS_OK) {
        return status;
    }

    status = server_map_database_status(geo_database_write_objects(connection->server->database,
                                                                   mutations,
                                                                   operation_count));
    free(mutations);
    return status;
}

static GeoProtocolStatus server_execute_write_objects_idempotent(GeoConnection *connection,
                                                                 unsigned char output[sizeof(uint64_t)],
                                                                 uint32_t *output_size)
{
    const unsigned char *payload = connection->request_payload;
    size_t payload_size = connection->request_header.payload_size;

    if (payload_size < 8U) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    uint32_t key_size = geo_protocol_load_u32(payload);
    uint32_t retention_seconds = geo_protocol_load_u32(payload + 4U);

    if (!key_size || key_size > GEO_DATABASE_MAX_IDEMPOTENCY_KEY_SIZE || key_size > payload_size - 8U) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    const unsigned char *key = payload + 8U;
    const unsigned char *mutation_payload = key + key_size;
    size_t mutation_payload_size = payload_size - 8U - key_size;
    GeoDatabaseObjectMutation *mutations;
    uint32_t operation_count;
    GeoProtocolStatus status = server_decode_object_mutations(mutation_payload,
                                                              mutation_payload_size,
                                                              &mutations,
                                                              &operation_count);

    if (status != GEO_PROTOCOL_STATUS_OK) {
        return status;
    }

    bool replayed;
    GeoDatabaseStatus database_status = geo_database_write_objects_idempotent(connection->server->database,
                                                                              mutations,
                                                                              operation_count,
                                                                              key,
                                                                              key_size,
                                                                              retention_seconds,
                                                                              &replayed);
    free(mutations);

    if (database_status != GEO_DATABASE_OK) {
        return server_map_database_status(database_status);
    }

    geo_protocol_store_u64(output, replayed ? 1U : 0U);
    *output_size = sizeof(uint64_t);
    return GEO_PROTOCOL_STATUS_OK;
}

static GeoProtocolStatus server_execute_get_object(GeoConnection *connection,
                                                   unsigned char **output,
                                                   uint32_t *output_size)
{
    if (connection->request_header.payload_size != sizeof(uint64_t)) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    uint64_t object_id = geo_protocol_load_u64(connection->request_payload);
    GeoDatabaseObject object;
    GeoDatabaseStatus database_status = geo_database_get(connection->server->database, object_id, &object);

    if (database_status != GEO_DATABASE_OK) {
        return server_map_database_status(database_status);
    }

    if (object.document_size > UINT32_MAX - GEO_PROTOCOL_OBJECT_RESPONSE_SIZE ||
        object.document_size + GEO_PROTOCOL_OBJECT_RESPONSE_SIZE > connection->server->max_frame_size) {
        geo_database_object_release(&object);
        return GEO_PROTOCOL_STATUS_FRAME_TOO_LARGE;
    }

    size_t response_size = GEO_PROTOCOL_OBJECT_RESPONSE_SIZE + object.document_size;
    unsigned char *response = malloc(response_size);

    if (!response) {
        geo_database_object_release(&object);
        return GEO_PROTOCOL_STATUS_INTERNAL_ERROR;
    }

    geo_protocol_store_u64(response, object.object_id);
    geo_protocol_store_u64(response + 8U, object.sequence);
    geo_protocol_store_u64(response + 16U, object.morton_code);
    geo_protocol_store_u32(response + 24U, (uint32_t) object.document_size);
    geo_protocol_store_u32(response + 28U, 0U);

    if (object.document_size) {
        memcpy(response + GEO_PROTOCOL_OBJECT_RESPONSE_SIZE, object.document, object.document_size);
    }

    geo_database_object_release(&object);
    *output = response;
    *output_size = (uint32_t) response_size;
    return GEO_PROTOCOL_STATUS_OK;
}

static GeoProtocolStatus server_execute_insert_generated(GeoConnection *connection,
                                                         unsigned char output[sizeof(uint64_t)],
                                                         uint32_t *output_size)
{
    size_t payload_size = connection->request_header.payload_size;

    if (payload_size < GEO_PROTOCOL_GENERATED_INSERT_SIZE) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    const unsigned char *payload = connection->request_payload;
    uint32_t document_size = geo_protocol_load_u32(payload + 8U);

    if (geo_protocol_load_u32(payload + 12U) != 0U ||
        document_size != payload_size - GEO_PROTOCOL_GENERATED_INSERT_SIZE) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    uint64_t object_id;
    GeoDatabaseStatus database_status = geo_database_insert_generated_object(connection->server->database,
                                                                             geo_protocol_load_u64(payload),
                                                                             document_size
                                                                                 ? payload + GEO_PROTOCOL_GENERATED_INSERT_SIZE
                                                                                 : NULL,
                                                                             document_size,
                                                                             &object_id);

    if (database_status != GEO_DATABASE_OK) {
        return server_map_database_status(database_status);
    }

    geo_protocol_store_u64(output, object_id);
    *output_size = sizeof(object_id);
    return GEO_PROTOCOL_STATUS_OK;
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

static GeoProtocolStatus server_execute_create_index(GeoConnection *connection)
{
    const unsigned char *payload = connection->request_payload;
    size_t payload_size = connection->request_header.payload_size;

    if (payload_size < GEO_PROTOCOL_INDEX_CREATE_HEADER_SIZE) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    uint32_t type = geo_protocol_load_u32(payload);
    uint32_t name_size = geo_protocol_load_u32(payload + 4U);
    uint32_t pointer_size = geo_protocol_load_u32(payload + 8U);
    size_t variable_size = payload_size - GEO_PROTOCOL_INDEX_CREATE_HEADER_SIZE;

    if (geo_protocol_load_u32(payload + 12U) != 0U || type < GEO_DATABASE_INDEX_BOOL || type > GEO_DATABASE_INDEX_BYTES ||
        name_size == 0U || name_size > GEO_DATABASE_MAX_INDEX_NAME_SIZE ||
        pointer_size > GEO_DATABASE_MAX_INDEX_POINTER_SIZE || name_size > variable_size ||
        pointer_size != variable_size - name_size) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    char name[GEO_DATABASE_MAX_INDEX_NAME_SIZE + 1U];
    char json_pointer[GEO_DATABASE_MAX_INDEX_POINTER_SIZE + 1U];

    memcpy(name, payload + GEO_PROTOCOL_INDEX_CREATE_HEADER_SIZE, name_size);
    name[name_size] = '\0';
    memcpy(json_pointer, payload + GEO_PROTOCOL_INDEX_CREATE_HEADER_SIZE + name_size, pointer_size);
    json_pointer[pointer_size] = '\0';
    return server_map_database_status(geo_database_create_index(connection->server->database,
                                                                name,
                                                                json_pointer,
                                                                (GeoDatabaseIndexType) type));
}

static GeoProtocolStatus server_execute_drop_index(GeoConnection *connection)
{
    const unsigned char *payload = connection->request_payload;
    size_t payload_size = connection->request_header.payload_size;

    if (payload_size < GEO_PROTOCOL_INDEX_DROP_HEADER_SIZE) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    uint32_t name_size = geo_protocol_load_u32(payload);

    if (geo_protocol_load_u32(payload + 4U) != 0U || name_size == 0U || name_size > GEO_DATABASE_MAX_INDEX_NAME_SIZE ||
        GEO_PROTOCOL_INDEX_DROP_HEADER_SIZE + (size_t) name_size != payload_size) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    char name[GEO_DATABASE_MAX_INDEX_NAME_SIZE + 1U];

    memcpy(name, payload + GEO_PROTOCOL_INDEX_DROP_HEADER_SIZE, name_size);
    name[name_size] = '\0';
    return server_map_database_status(geo_database_drop_index(connection->server->database, name));
}

static GeoProtocolStatus server_execute_list_indexes(GeoConnection *connection,
                                                     unsigned char **output,
                                                     uint32_t *output_size)
{
    if (connection->request_header.payload_size != 0U) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    size_t count = 0U;
    GeoDatabaseStatus database_status = geo_database_list_indexes(connection->server->database, NULL, 0U, &count);

    if (database_status != GEO_DATABASE_OK && database_status != GEO_DATABASE_OUT_OF_MEMORY) {
        return server_map_database_status(database_status);
    }
    if (count > UINT32_MAX || count > SIZE_MAX / sizeof(GeoDatabaseIndexInfo)) {
        return GEO_PROTOCOL_STATUS_FRAME_TOO_LARGE;
    }

    GeoDatabaseIndexInfo *indexes = NULL;
    size_t capacity = 0U;

    while (database_status != GEO_DATABASE_OK || capacity < count) {
        if (count > UINT32_MAX || count > SIZE_MAX / sizeof(*indexes)) {
            free(indexes);
            return GEO_PROTOCOL_STATUS_FRAME_TOO_LARGE;
        }

        GeoDatabaseIndexInfo *replacement = count ? realloc(indexes, count * sizeof(*indexes)) : indexes;

        if (count && !replacement) {
            free(indexes);
            return GEO_PROTOCOL_STATUS_INTERNAL_ERROR;
        }

        indexes = replacement;
        capacity = count;
        database_status = geo_database_list_indexes(connection->server->database, indexes, capacity, &count);

        if (database_status != GEO_DATABASE_OK && database_status != GEO_DATABASE_OUT_OF_MEMORY) {
            free(indexes);
            return server_map_database_status(database_status);
        }
    }

    size_t response_size = 8U;

    for (size_t index = 0U; index < count; ++index) {
        size_t name_size = strlen(indexes[index].name);
        size_t pointer_size = strlen(indexes[index].json_pointer);
        size_t entry_size = 32U + name_size + pointer_size;

        if (entry_size > SIZE_MAX - response_size) {
            free(indexes);
            return GEO_PROTOCOL_STATUS_FRAME_TOO_LARGE;
        }
        response_size += entry_size;
    }

    if (response_size > UINT32_MAX || response_size > connection->server->max_frame_size) {
        free(indexes);
        return GEO_PROTOCOL_STATUS_FRAME_TOO_LARGE;
    }

    unsigned char *response = malloc(response_size);

    if (!response) {
        free(indexes);
        return GEO_PROTOCOL_STATUS_INTERNAL_ERROR;
    }

    geo_protocol_store_u32(response, (uint32_t) count);
    geo_protocol_store_u32(response + 4U, 0U);
    unsigned char *cursor = response + 8U;

    for (size_t index = 0U; index < count; ++index) {
        size_t name_size = strlen(indexes[index].name);
        size_t pointer_size = strlen(indexes[index].json_pointer);

        geo_protocol_store_u64(cursor, indexes[index].index_id);
        geo_protocol_store_u64(cursor + 8U, indexes[index].entry_count);
        geo_protocol_store_u32(cursor + 16U, (uint32_t) indexes[index].type);
        geo_protocol_store_u32(cursor + 20U, (uint32_t) name_size);
        geo_protocol_store_u32(cursor + 24U, (uint32_t) pointer_size);
        geo_protocol_store_u32(cursor + 28U, 0U);
        memcpy(cursor + 32U, indexes[index].name, name_size);
        memcpy(cursor + 32U + name_size, indexes[index].json_pointer, pointer_size);
        cursor += 32U + name_size + pointer_size;
    }

    free(indexes);
    *output = response;
    *output_size = (uint32_t) response_size;
    return GEO_PROTOCOL_STATUS_OK;
}

static GeoProtocolStatus server_execute_query_index(GeoConnection *connection,
                                                    unsigned char **output,
                                                    uint32_t *output_size)
{
    const unsigned char *payload = connection->request_payload;
    size_t payload_size = connection->request_header.payload_size;

    char name[GEO_DATABASE_MAX_INDEX_NAME_SIZE + 1U];
    GeoDatabaseIndexPredicate predicate;
    size_t consumed = 0U;

    if (!geo_protocol_decode_index_predicate(payload, payload_size, name, &predicate, &consumed) || consumed != payload_size) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    GeoIdResult *result = geo_id_result_create(64U);

    if (!result) {
        return GEO_PROTOCOL_STATUS_INTERNAL_ERROR;
    }

    GeoDatabaseIndexQueryStats stats;
    GeoDatabaseStatus database_status = geo_database_query_index_reuse(connection->server->database,
                                                                       &predicate,
                                                                       result,
                                                                       &stats);

    if (database_status != GEO_DATABASE_OK) {
        geo_id_result_destroy(result);
        return server_map_database_status(database_status);
    }
    if (connection->server->max_frame_size < GEO_PROTOCOL_INDEX_QUERY_RESPONSE_SIZE || result->count > UINT32_MAX ||
        result->count > (connection->server->max_frame_size - GEO_PROTOCOL_INDEX_QUERY_RESPONSE_SIZE) / sizeof(uint64_t)) {
        geo_id_result_destroy(result);
        return GEO_PROTOCOL_STATUS_FRAME_TOO_LARGE;
    }

    size_t response_size = GEO_PROTOCOL_INDEX_QUERY_RESPONSE_SIZE + result->count * sizeof(uint64_t);
    unsigned char *response = malloc(response_size);

    if (!response) {
        geo_id_result_destroy(result);
        return GEO_PROTOCOL_STATUS_INTERNAL_ERROR;
    }

    geo_protocol_store_u64(response, stats.scanned_entries);
    geo_protocol_store_u64(response + 8U, stats.matched_entries);
    geo_protocol_store_u32(response + 16U, (uint32_t) result->count);
    geo_protocol_store_u32(response + 20U, 0U);

    for (size_t index = 0U; index < result->count; ++index) {
        geo_protocol_store_u64(response + GEO_PROTOCOL_INDEX_QUERY_RESPONSE_SIZE + index * sizeof(uint64_t), result->ids[index]);
    }

    geo_id_result_destroy(result);
    *output = response;
    *output_size = (uint32_t) response_size;
    return GEO_PROTOCOL_STATUS_OK;
}

static GeoProtocolStatus server_execute_query_radius(GeoConnection *connection,
                                                     unsigned char **output,
                                                     uint32_t *output_size)
{
    const unsigned char *payload = connection->request_payload;
    size_t payload_size = connection->request_header.payload_size;

    if (payload_size < GEO_PROTOCOL_RADIUS_QUERY_HEADER_SIZE) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    uint32_t predicate_count = geo_protocol_load_u32(payload + 24U);

    if (!predicate_count || predicate_count > GEO_DATABASE_MAX_QUERY_PREDICATES ||
        geo_protocol_load_u32(payload + 28U) != 0U) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    GeoDatabaseIndexPredicate predicates[GEO_DATABASE_MAX_QUERY_PREDICATES];
    char names[GEO_DATABASE_MAX_QUERY_PREDICATES][GEO_DATABASE_MAX_INDEX_NAME_SIZE + 1U];
    const unsigned char *cursor = payload + GEO_PROTOCOL_RADIUS_QUERY_HEADER_SIZE;
    size_t remaining = payload_size - GEO_PROTOCOL_RADIUS_QUERY_HEADER_SIZE;

    for (size_t index = 0U; index < predicate_count; ++index) {
        size_t consumed = 0U;

        if (!geo_protocol_decode_index_predicate(cursor,
                                                 remaining,
                                                 names[index],
                                                 predicates + index,
                                                 &consumed)) {
            return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
        }

        cursor += consumed;
        remaining -= consumed;
    }

    if (remaining != 0U) {
        return GEO_PROTOCOL_STATUS_INVALID_REQUEST;
    }

    if (!connection->query_workspace) {
        connection->query_workspace = geo_database_query_workspace_create(GEO_DATABASE_MAX_QUERY_PREDICATES);

        if (!connection->query_workspace) {
            return GEO_PROTOCOL_STATUS_INTERNAL_ERROR;
        }
    }

    if (!connection->query_result) {
        connection->query_result = geo_result_create(64U);

        if (!connection->query_result) {
            return GEO_PROTOCOL_STATUS_INTERNAL_ERROR;
        }
    }

    GeoDatabaseRadiusQuery query = {
        .predicates = predicates,
        .predicate_count = predicate_count,
        .latitude = geo_protocol_load_double(payload),
        .longitude = geo_protocol_load_double(payload + 8U),
        .radius_km = geo_protocol_load_double(payload + 16U),
    };
    GeoDatabaseQueryStats stats;
    GeoDatabaseStatus database_status = geo_database_query_radius_reuse(connection->server->database,
                                                                        &query,
                                                                        connection->query_workspace,
                                                                        connection->query_result,
                                                                        &stats);

    if (database_status != GEO_DATABASE_OK) {
        return server_map_database_status(database_status);
    }

    GeoSearchResult *result = connection->query_result;

    if (connection->server->max_frame_size < GEO_PROTOCOL_RADIUS_QUERY_RESPONSE_SIZE || result->count > UINT32_MAX ||
        result->count > (connection->server->max_frame_size - GEO_PROTOCOL_RADIUS_QUERY_RESPONSE_SIZE) / sizeof(GeoRecord)) {
        return GEO_PROTOCOL_STATUS_FRAME_TOO_LARGE;
    }

    size_t response_size = GEO_PROTOCOL_RADIUS_QUERY_RESPONSE_SIZE + result->count * sizeof(GeoRecord);
    unsigned char *response = malloc(response_size);

    if (!response) {
        return GEO_PROTOCOL_STATUS_INTERNAL_ERROR;
    }

    geo_protocol_store_u64(response, stats.secondary_entries_scanned);
    geo_protocol_store_u64(response + 8U, stats.metadata_candidates);
    geo_protocol_store_u64(response + 16U, stats.estimated_spatial_candidates);
    geo_protocol_store_u64(response + 24U, stats.object_lookups);
    geo_protocol_store_u64(response + 32U, stats.spatial.records_scanned);
    geo_protocol_store_u64(response + 40U, stats.spatial.records_matched);
    geo_protocol_store_u64(response + 48U, stats.spatial.ranges_checked);
    geo_protocol_store_double(response + 56U, stats.spatial.search_time_ms);
    geo_protocol_store_u32(response + 64U, (uint32_t) stats.plan);
    geo_protocol_store_u32(response + 68U, stats.predicate_count);
    geo_protocol_store_u32(response + 72U, (uint32_t) result->count);
    geo_protocol_store_u32(response + 76U, 0U);

    for (size_t index = 0U; index < result->count; ++index) {
        unsigned char *encoded = response + GEO_PROTOCOL_RADIUS_QUERY_RESPONSE_SIZE + index * sizeof(GeoRecord);

        geo_protocol_store_u64(encoded, result->results[index].id);
        geo_protocol_store_u64(encoded + 8U, result->results[index].z);
    }

    *output = response;
    *output_size = (uint32_t) response_size;
    return GEO_PROTOCOL_STATUS_OK;
}

void geo_server_execute_request(void *context)
{
    GeoConnection *connection = context;
    GeoServer *server = connection->server;
    unsigned char output[GEO_PROTOCOL_STATS_SIZE];
    unsigned char *dynamic_output = NULL;
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

            case GEO_PROTOCOL_WRITE_OBJECTS:
                status = server_execute_write_objects(connection);
                break;

            case GEO_PROTOCOL_GET_OBJECT:
                status = server_execute_get_object(connection, &dynamic_output, &output_size);
                break;

            case GEO_PROTOCOL_INSERT_GENERATED:
                status = server_execute_insert_generated(connection, output, &output_size);
                break;

            case GEO_PROTOCOL_WRITE_OBJECTS_IDEMPOTENT:
                status = server_execute_write_objects_idempotent(connection, output, &output_size);
                break;

            case GEO_PROTOCOL_CREATE_INDEX:
                status = server_execute_create_index(connection);
                break;

            case GEO_PROTOCOL_DROP_INDEX:
                status = server_execute_drop_index(connection);
                break;

            case GEO_PROTOCOL_LIST_INDEXES:
                status = server_execute_list_indexes(connection, &dynamic_output, &output_size);
                break;

            case GEO_PROTOCOL_QUERY_INDEX:
                status = server_execute_query_index(connection, &dynamic_output, &output_size);
                break;

            case GEO_PROTOCOL_QUERY_RADIUS:
                status = server_execute_query_radius(connection, &dynamic_output, &output_size);
                break;

            default:
                status = GEO_PROTOCOL_STATUS_INVALID_REQUEST;
                break;
        }
    }

    const void *response_payload = dynamic_output ? dynamic_output : output;

    if (connection->request_header.opcode == GEO_PROTOCOL_PING && status == GEO_PROTOCOL_STATUS_OK) {
        response_payload = connection->request_payload;
        output_size = connection->request_header.payload_size;
    }

    if (!geo_server_prepare_response(connection, status, response_payload, output_size)) {
        connection->close_after_response = true;
    }

    free(dynamic_output);

    geo_server_complete_job(connection);
}
