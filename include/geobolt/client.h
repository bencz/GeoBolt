#ifndef GEOBOLT_CLIENT_H
#define GEOBOLT_CLIENT_H

#include "geobolt/geobolt.h"

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GeoClient GeoClient;

typedef enum {
    GEO_CLIENT_OK = 0,
    GEO_CLIENT_INVALID_ARGUMENT,
    GEO_CLIENT_OUT_OF_MEMORY,
    GEO_CLIENT_CONNECT_ERROR,
    GEO_CLIENT_IO_ERROR,
    GEO_CLIENT_PROTOCOL_ERROR,
    GEO_CLIENT_UNAUTHORIZED,
    GEO_CLIENT_BUSY,
    GEO_CLIENT_SERVER_ERROR,
} GeoClientStatus;

typedef struct {
    const char *host;
    const char *authentication_token;
    size_t max_frame_size;
    uint32_t timeout_ms;
    uint16_t port;
} GeoClientConfig;

GeoClientConfig geo_client_default_config(const char *host, uint16_t port);
GeoClient *geo_client_connect(const GeoClientConfig *config, GeoClientStatus *status);
void geo_client_close(GeoClient *client);

GeoClientStatus geo_client_ping(GeoClient *client, const void *payload, size_t payload_size);
GeoClientStatus geo_client_write(GeoClient *client, const GeoDatabaseMutation *mutations, size_t mutation_count);
GeoClientStatus geo_client_radius_count(GeoClient *client,
                                        double latitude,
                                        double longitude,
                                        double radius_km,
                                        uint64_t *count);
GeoClientStatus geo_client_get_database_stats(GeoClient *client, GeoDatabaseStats *stats);
GeoClientStatus geo_client_checkpoint(GeoClient *client);

#ifdef __cplusplus
}
#endif

#endif
