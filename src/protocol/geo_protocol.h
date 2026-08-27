#ifndef GEO_PROTOCOL_H
#define GEO_PROTOCOL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define GEO_PROTOCOL_MAGIC UINT32_C(0x47424c54)
#define GEO_PROTOCOL_VERSION 1U
#define GEO_PROTOCOL_HEADER_SIZE 48U
#define GEO_PROTOCOL_OPERATION_SIZE 24U
#define GEO_PROTOCOL_RADIUS_SIZE 24U
#define GEO_PROTOCOL_STATS_SIZE 64U
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

#endif
