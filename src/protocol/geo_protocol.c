#include "geo_protocol.h"

#include <string.h>

#define GEO_PROTOCOL_HEADER_CHECKSUM_OFFSET 40U

_Static_assert(GEO_PROTOCOL_HEADER_CHECKSUM_OFFSET + sizeof(uint64_t) == GEO_PROTOCOL_HEADER_SIZE,
               "Protocol header checksum must be the final field");

void geo_protocol_store_u32(unsigned char *output, uint32_t value)
{
    output[0] = (unsigned char) (value >> 24U);
    output[1] = (unsigned char) (value >> 16U);
    output[2] = (unsigned char) (value >> 8U);
    output[3] = (unsigned char) value;
}

void geo_protocol_store_u64(unsigned char *output, uint64_t value)
{
    geo_protocol_store_u32(output, (uint32_t) (value >> 32U));
    geo_protocol_store_u32(output + sizeof(uint32_t), (uint32_t) value);
}

uint32_t geo_protocol_load_u32(const unsigned char *input)
{
    return (uint32_t) input[0] << 24U |
           (uint32_t) input[1] << 16U |
           (uint32_t) input[2] << 8U |
           (uint32_t) input[3];
}

uint64_t geo_protocol_load_u64(const unsigned char *input)
{
    return (uint64_t) geo_protocol_load_u32(input) << 32U |
           geo_protocol_load_u32(input + sizeof(uint32_t));
}

void geo_protocol_store_double(unsigned char *output, double value)
{
    uint64_t bits;

    memcpy(&bits, &value, sizeof(bits));
    geo_protocol_store_u64(output, bits);
}

double geo_protocol_load_double(const unsigned char *input)
{
    uint64_t bits = geo_protocol_load_u64(input);
    double value;

    memcpy(&value, &bits, sizeof(value));

    return value;
}

size_t geo_protocol_index_value_size(const GeoDatabaseIndexValue *value)
{
    if (!value) {
        return 0U;
    }

    switch (value->type) {
        case GEO_DATABASE_INDEX_BOOL:
            return 1U;
        case GEO_DATABASE_INDEX_INT64:
        case GEO_DATABASE_INDEX_UINT64:
        case GEO_DATABASE_INDEX_DOUBLE:
        case GEO_DATABASE_INDEX_DATETIME:
            return 8U;
        case GEO_DATABASE_INDEX_STRING:
        case GEO_DATABASE_INDEX_BYTES:
            return value->as.bytes.size <= GEO_DATABASE_MAX_INDEXED_VALUE_SIZE ? value->as.bytes.size : 0U;
    }

    return 0U;
}

bool geo_protocol_encode_index_value(unsigned char *output,
                                     size_t output_size,
                                     const GeoDatabaseIndexValue *value)
{
    size_t expected_size = geo_protocol_index_value_size(value);

    if (!value || !output || output_size != expected_size ||
        ((value->type == GEO_DATABASE_INDEX_STRING || value->type == GEO_DATABASE_INDEX_BYTES) &&
         value->as.bytes.size && !value->as.bytes.data)) {
        return false;
    }

    switch (value->type) {
        case GEO_DATABASE_INDEX_BOOL:
            output[0] = value->as.boolean ? 1U : 0U;
            return true;
        case GEO_DATABASE_INDEX_INT64:
            geo_protocol_store_u64(output, (uint64_t) value->as.signed_integer);
            return true;
        case GEO_DATABASE_INDEX_UINT64:
            geo_protocol_store_u64(output, value->as.unsigned_integer);
            return true;
        case GEO_DATABASE_INDEX_DOUBLE:
            geo_protocol_store_double(output, value->as.floating_point);
            return true;
        case GEO_DATABASE_INDEX_DATETIME:
            geo_protocol_store_u64(output, (uint64_t) value->as.datetime);
            return true;
        case GEO_DATABASE_INDEX_STRING:
        case GEO_DATABASE_INDEX_BYTES:
            if (output_size) {
                memcpy(output, value->as.bytes.data, output_size);
            }
            return true;
    }

    return false;
}

bool geo_protocol_decode_index_value(GeoDatabaseIndexType type,
                                     const void *data,
                                     size_t size,
                                     GeoDatabaseIndexValue *value)
{
    if (!data || !value) {
        return false;
    }

    *value = (GeoDatabaseIndexValue) { .type = type };
    const unsigned char *bytes = data;

    switch (type) {
        case GEO_DATABASE_INDEX_BOOL:
            if (size != 1U || bytes[0] > 1U) {
                return false;
            }
            value->as.boolean = bytes[0] != 0U;
            return true;
        case GEO_DATABASE_INDEX_INT64:
            if (size != 8U) {
                return false;
            }
            value->as.signed_integer = (int64_t) geo_protocol_load_u64(bytes);
            return true;
        case GEO_DATABASE_INDEX_UINT64:
            if (size != 8U) {
                return false;
            }
            value->as.unsigned_integer = geo_protocol_load_u64(bytes);
            return true;
        case GEO_DATABASE_INDEX_DOUBLE:
            if (size != 8U) {
                return false;
            }
            value->as.floating_point = geo_protocol_load_double(bytes);
            return true;
        case GEO_DATABASE_INDEX_DATETIME:
            if (size != 8U) {
                return false;
            }
            value->as.datetime = (int64_t) geo_protocol_load_u64(bytes);
            return true;
        case GEO_DATABASE_INDEX_STRING:
        case GEO_DATABASE_INDEX_BYTES:
            if (size > GEO_DATABASE_MAX_INDEXED_VALUE_SIZE) {
                return false;
            }
            value->as.bytes.data = data;
            value->as.bytes.size = size;
            return true;
    }

    return false;
}

size_t geo_protocol_index_predicate_size(const GeoDatabaseIndexPredicate *predicate)
{
    if (!predicate || !predicate->index_name ||
        predicate->operation < GEO_DATABASE_INDEX_EQUAL || predicate->operation > GEO_DATABASE_INDEX_BETWEEN ||
        predicate->lower.type < GEO_DATABASE_INDEX_BOOL || predicate->lower.type > GEO_DATABASE_INDEX_BYTES) {
        return 0U;
    }

    size_t name_size = strnlen(predicate->index_name, GEO_DATABASE_MAX_INDEX_NAME_SIZE + 1U);
    size_t lower_size = geo_protocol_index_value_size(&predicate->lower);
    size_t upper_size = predicate->operation == GEO_DATABASE_INDEX_BETWEEN
                            ? geo_protocol_index_value_size(&predicate->upper)
                            : 0U;
    bool lower_empty_bytes = (predicate->lower.type == GEO_DATABASE_INDEX_STRING ||
                              predicate->lower.type == GEO_DATABASE_INDEX_BYTES) &&
                             predicate->lower.as.bytes.size == 0U;
    bool upper_empty_bytes = (predicate->upper.type == GEO_DATABASE_INDEX_STRING ||
                              predicate->upper.type == GEO_DATABASE_INDEX_BYTES) &&
                             predicate->upper.as.bytes.size == 0U;

    if (!name_size || name_size > GEO_DATABASE_MAX_INDEX_NAME_SIZE || (!lower_size && !lower_empty_bytes) ||
        (predicate->operation == GEO_DATABASE_INDEX_BETWEEN &&
         (predicate->upper.type != predicate->lower.type || (!upper_size && !upper_empty_bytes)))) {
        return 0U;
    }

    return GEO_PROTOCOL_INDEX_QUERY_HEADER_SIZE + name_size + lower_size + upper_size;
}

bool geo_protocol_encode_index_predicate(unsigned char *output,
                                         size_t output_size,
                                         const GeoDatabaseIndexPredicate *predicate)
{
    size_t expected_size = geo_protocol_index_predicate_size(predicate);

    if (!output || !expected_size || output_size != expected_size) {
        return false;
    }

    size_t name_size = strlen(predicate->index_name);
    size_t lower_size = geo_protocol_index_value_size(&predicate->lower);
    size_t upper_size = predicate->operation == GEO_DATABASE_INDEX_BETWEEN
                            ? geo_protocol_index_value_size(&predicate->upper)
                            : 0U;

    geo_protocol_store_u32(output, (uint32_t) predicate->operation);
    geo_protocol_store_u32(output + 4U, (uint32_t) predicate->lower.type);
    geo_protocol_store_u32(output + 8U, (uint32_t) name_size);
    geo_protocol_store_u32(output + 12U, (uint32_t) lower_size);
    geo_protocol_store_u32(output + 16U, (uint32_t) upper_size);
    geo_protocol_store_u32(output + 20U, 0U);
    memcpy(output + GEO_PROTOCOL_INDEX_QUERY_HEADER_SIZE, predicate->index_name, name_size);

    unsigned char *encoded_values = output + GEO_PROTOCOL_INDEX_QUERY_HEADER_SIZE + name_size;
    return geo_protocol_encode_index_value(encoded_values, lower_size, &predicate->lower) &&
           (predicate->operation != GEO_DATABASE_INDEX_BETWEEN ||
            geo_protocol_encode_index_value(encoded_values + lower_size, upper_size, &predicate->upper));
}

bool geo_protocol_decode_index_predicate(const void *data,
                                         size_t size,
                                         char name[GEO_DATABASE_MAX_INDEX_NAME_SIZE + 1U],
                                         GeoDatabaseIndexPredicate *predicate,
                                         size_t *consumed)
{
    if (!data || size < GEO_PROTOCOL_INDEX_QUERY_HEADER_SIZE || !name || !predicate || !consumed) {
        return false;
    }

    const unsigned char *bytes = data;
    uint32_t operation = geo_protocol_load_u32(bytes);
    uint32_t type = geo_protocol_load_u32(bytes + 4U);
    uint32_t name_size = geo_protocol_load_u32(bytes + 8U);
    uint32_t lower_size = geo_protocol_load_u32(bytes + 12U);
    uint32_t upper_size = geo_protocol_load_u32(bytes + 16U);
    size_t variable_size = size - GEO_PROTOCOL_INDEX_QUERY_HEADER_SIZE;

    if (geo_protocol_load_u32(bytes + 20U) != 0U || operation < GEO_DATABASE_INDEX_EQUAL ||
        operation > GEO_DATABASE_INDEX_BETWEEN || type < GEO_DATABASE_INDEX_BOOL || type > GEO_DATABASE_INDEX_BYTES ||
        name_size == 0U || name_size > GEO_DATABASE_MAX_INDEX_NAME_SIZE || name_size > variable_size ||
        lower_size > variable_size - name_size || upper_size > variable_size - name_size - lower_size ||
        (operation != GEO_DATABASE_INDEX_BETWEEN && upper_size != 0U)) {
        return false;
    }

    const unsigned char *encoded_values = bytes + GEO_PROTOCOL_INDEX_QUERY_HEADER_SIZE + name_size;
    GeoDatabaseIndexPredicate decoded = {
        .index_name = name,
        .operation = (GeoDatabaseIndexOperator) operation,
    };

    memcpy(name, bytes + GEO_PROTOCOL_INDEX_QUERY_HEADER_SIZE, name_size);
    name[name_size] = '\0';

    if (!geo_protocol_decode_index_value((GeoDatabaseIndexType) type, encoded_values, lower_size, &decoded.lower) ||
        (operation == GEO_DATABASE_INDEX_BETWEEN &&
         !geo_protocol_decode_index_value((GeoDatabaseIndexType) type,
                                          encoded_values + lower_size,
                                          upper_size,
                                          &decoded.upper))) {
        return false;
    }

    *predicate = decoded;
    *consumed = GEO_PROTOCOL_INDEX_QUERY_HEADER_SIZE + (size_t) name_size + lower_size + upper_size;
    return true;
}

uint64_t geo_protocol_checksum(const void *data, size_t size)
{
    const unsigned char *bytes = data;
    uint64_t lane0 = UINT64_C(0x9e3779b185ebca87) ^ size;
    uint64_t lane1 = UINT64_C(0xc2b2ae3d27d4eb4f) + size;
    uint64_t lane2 = UINT64_C(0x165667b19e3779f9) ^ (size << 1U);
    uint64_t lane3 = UINT64_C(0x85ebca77c2b2ae63) + (size << 1U);

    while (size >= 32U) {
        lane0 ^= geo_protocol_load_u64(bytes);
        lane1 ^= geo_protocol_load_u64(bytes + 8U);
        lane2 ^= geo_protocol_load_u64(bytes + 16U);
        lane3 ^= geo_protocol_load_u64(bytes + 24U);

        lane0 = (lane0 << 29U | lane0 >> 35U) * UINT64_C(0x9fb21c651e98df25);
        lane1 = (lane1 << 31U | lane1 >> 33U) * UINT64_C(0x9e3779b185ebca87);
        lane2 = (lane2 << 33U | lane2 >> 31U) * UINT64_C(0xc2b2ae3d27d4eb4f);
        lane3 = (lane3 << 35U | lane3 >> 29U) * UINT64_C(0x165667b19e3779f9);

        bytes += 32U;
        size -= 32U;
    }

    uint64_t checksum = lane0 ^
                        (lane1 << 17U | lane1 >> 47U) ^
                        (lane2 << 31U | lane2 >> 33U) ^
                        (lane3 << 47U | lane3 >> 17U);

    while (size >= sizeof(uint64_t)) {
        checksum ^= geo_protocol_load_u64(bytes);
        checksum = (checksum << 27U | checksum >> 37U) * UINT64_C(0x9fb21c651e98df25);
        bytes += sizeof(uint64_t);
        size -= sizeof(uint64_t);
    }

    while (size) {
        checksum ^= *bytes++;
        checksum = (checksum << 11U | checksum >> 53U) * UINT64_C(0x9e3779b185ebca87);
        size--;
    }

    checksum ^= checksum >> 33U;
    checksum *= UINT64_C(0xff51afd7ed558ccd);
    checksum ^= checksum >> 29U;
    checksum *= UINT64_C(0xc4ceb9fe1a85ec53);
    checksum ^= checksum >> 32U;

    return checksum;
}

static uint64_t protocol_header_checksum(const unsigned char encoded[GEO_PROTOCOL_HEADER_SIZE])
{
    unsigned char copy[GEO_PROTOCOL_HEADER_SIZE];

    memcpy(copy, encoded, sizeof(copy));
    memset(copy + GEO_PROTOCOL_HEADER_CHECKSUM_OFFSET, 0, sizeof(uint64_t));

    return geo_protocol_checksum(copy, sizeof(copy));
}

void geo_protocol_encode_header(unsigned char output[GEO_PROTOCOL_HEADER_SIZE], const GeoProtocolHeader *header)
{
    memset(output, 0, GEO_PROTOCOL_HEADER_SIZE);
    geo_protocol_store_u32(output, GEO_PROTOCOL_MAGIC);
    output[4] = (unsigned char) (GEO_PROTOCOL_VERSION >> 8U);
    output[5] = (unsigned char) GEO_PROTOCOL_VERSION;
    output[6] = (unsigned char) (header->opcode >> 8U);
    output[7] = (unsigned char) header->opcode;
    geo_protocol_store_u32(output + 8U, header->flags);
    geo_protocol_store_u32(output + 12U, header->status);
    geo_protocol_store_u64(output + 16U, header->request_id);
    geo_protocol_store_u32(output + 24U, header->payload_size);
    geo_protocol_store_u64(output + 32U, header->payload_checksum);
    geo_protocol_store_u64(output + GEO_PROTOCOL_HEADER_CHECKSUM_OFFSET, protocol_header_checksum(output));
}

bool geo_protocol_decode_header(const unsigned char input[GEO_PROTOCOL_HEADER_SIZE], GeoProtocolHeader *header)
{
    if (!header ||
        geo_protocol_load_u32(input) != GEO_PROTOCOL_MAGIC ||
        ((uint16_t) input[4] << 8U | input[5]) != GEO_PROTOCOL_VERSION ||
        geo_protocol_load_u32(input + 28U) != 0 ||
        geo_protocol_load_u64(input + GEO_PROTOCOL_HEADER_CHECKSUM_OFFSET) != protocol_header_checksum(input)) {
        return false;
    }

    header->opcode = (uint16_t) input[6] << 8U | input[7];
    header->flags = geo_protocol_load_u32(input + 8U);
    header->status = geo_protocol_load_u32(input + 12U);
    header->request_id = geo_protocol_load_u64(input + 16U);
    header->payload_size = geo_protocol_load_u32(input + 24U);
    header->payload_checksum = geo_protocol_load_u64(input + 32U);
    header->reserved = 0;

    return header->opcode >= GEO_PROTOCOL_AUTH && header->opcode <= GEO_PROTOCOL_QUERY_RADIUS;
}
