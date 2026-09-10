#ifndef GEO_PROTOCOL_H
#define GEO_PROTOCOL_H

#include "geobolt/geobolt.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define GEO_PROTOCOL_MAGIC UINT32_C(0x47424c54)
#define GEO_PROTOCOL_VERSION 4U
#define GEO_PROTOCOL_HEADER_SIZE 48U
#define GEO_PROTOCOL_OPERATION_SIZE 24U
#define GEO_PROTOCOL_RADIUS_SIZE 24U
#define GEO_PROTOCOL_STATS_SIZE 64U
#define GEO_PROTOCOL_OBJECT_OPERATION_SIZE 32U
#define GEO_PROTOCOL_OBJECT_RESPONSE_SIZE 32U
#define GEO_PROTOCOL_GENERATED_INSERT_SIZE 16U
#define GEO_PROTOCOL_INDEX_CREATE_HEADER_SIZE 16U
#define GEO_PROTOCOL_INDEX_DROP_HEADER_SIZE 8U
#define GEO_PROTOCOL_INDEX_QUERY_HEADER_SIZE 24U
#define GEO_PROTOCOL_INDEX_QUERY_RESPONSE_SIZE 24U
#define GEO_PROTOCOL_RADIUS_QUERY_HEADER_SIZE 32U
#define GEO_PROTOCOL_RADIUS_QUERY_RESPONSE_SIZE 80U
#define GEO_PROTOCOL_DEFAULT_MAX_FRAME (16U * 1024U * 1024U)
#define GEO_PROTOCOL_ABSOLUTE_MAX_FRAME (64U * 1024U * 1024U)
#define GEO_PROTOCOL_RESPONSE_FLAG UINT32_C(1)

typedef enum {
    GEO_PROTOCOL_AUTH = 1,
    GEO_PROTOCOL_PING = 2,
    GEO_PROTOCOL_WRITE = 3,
    GEO_PROTOCOL_RADIUS_COUNT = 4,
    GEO_PROTOCOL_STATS = 5,
    GEO_PROTOCOL_CHECKPOINT = 6,
    GEO_PROTOCOL_WRITE_OBJECTS = 7,
    GEO_PROTOCOL_GET_OBJECT = 8,
    GEO_PROTOCOL_INSERT_GENERATED = 9,
    GEO_PROTOCOL_WRITE_OBJECTS_IDEMPOTENT = 10,
    GEO_PROTOCOL_CREATE_INDEX = 11,
    GEO_PROTOCOL_DROP_INDEX = 12,
    GEO_PROTOCOL_LIST_INDEXES = 13,
    GEO_PROTOCOL_QUERY_INDEX = 14,
    GEO_PROTOCOL_QUERY_RADIUS = 15,
} GeoProtocolOpcode;

typedef enum {
    GEO_PROTOCOL_STATUS_OK = 0,
    GEO_PROTOCOL_STATUS_INVALID_REQUEST = 1,
    GEO_PROTOCOL_STATUS_UNAUTHORIZED = 2,
    GEO_PROTOCOL_STATUS_BUSY = 3,
    GEO_PROTOCOL_STATUS_INTERNAL_ERROR = 4,
    GEO_PROTOCOL_STATUS_DATABASE_ERROR = 5,
    GEO_PROTOCOL_STATUS_FRAME_TOO_LARGE = 6,
    GEO_PROTOCOL_STATUS_PROTOCOL_ERROR = 7,
    GEO_PROTOCOL_STATUS_NOT_FOUND = 8,
    GEO_PROTOCOL_STATUS_ALREADY_EXISTS = 9,
    GEO_PROTOCOL_STATUS_IDEMPOTENCY_CONFLICT = 10,
} GeoProtocolStatus;

typedef struct {
    uint64_t request_id;
    uint64_t payload_checksum;
    uint32_t flags;
    uint32_t status;
    uint32_t payload_size;
    uint16_t opcode;
    uint16_t reserved;
} GeoProtocolHeader;

_Static_assert(sizeof(GeoProtocolHeader) == 32U, "Decoded protocol headers must remain naturally aligned and compact");

uint64_t geo_protocol_checksum(const void *data, size_t size);
void geo_protocol_encode_header(unsigned char output[GEO_PROTOCOL_HEADER_SIZE], const GeoProtocolHeader *header);
bool geo_protocol_decode_header(const unsigned char input[GEO_PROTOCOL_HEADER_SIZE], GeoProtocolHeader *header);

void geo_protocol_store_u32(unsigned char *output, uint32_t value);
void geo_protocol_store_u64(unsigned char *output, uint64_t value);
uint32_t geo_protocol_load_u32(const unsigned char *input);
uint64_t geo_protocol_load_u64(const unsigned char *input);
void geo_protocol_store_double(unsigned char *output, double value);
double geo_protocol_load_double(const unsigned char *input);

size_t geo_protocol_index_value_size(const GeoDatabaseIndexValue *value);
bool geo_protocol_encode_index_value(unsigned char *output,
                                     size_t output_size,
                                     const GeoDatabaseIndexValue *value);
bool geo_protocol_decode_index_value(GeoDatabaseIndexType type,
                                     const void *data,
                                     size_t size,
                                     GeoDatabaseIndexValue *value);
size_t geo_protocol_index_predicate_size(const GeoDatabaseIndexPredicate *predicate);
bool geo_protocol_encode_index_predicate(unsigned char *output,
                                         size_t output_size,
                                         const GeoDatabaseIndexPredicate *predicate);
bool geo_protocol_decode_index_predicate(const void *data,
                                         size_t size,
                                         char name[GEO_DATABASE_MAX_INDEX_NAME_SIZE + 1U],
                                         GeoDatabaseIndexPredicate *predicate,
                                         size_t *consumed);

#endif
