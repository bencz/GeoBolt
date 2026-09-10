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
    GEO_CLIENT_NOT_FOUND,
    GEO_CLIENT_ALREADY_EXISTS,
    GEO_CLIENT_IDEMPOTENCY_CONFLICT,
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
GeoClientStatus geo_client_write_objects(GeoClient *client,
                                         const GeoDatabaseObjectMutation *mutations,
                                         size_t mutation_count);
GeoClientStatus geo_client_write_objects_idempotent(GeoClient *client,
                                                    const GeoDatabaseObjectMutation *mutations,
                                                    size_t mutation_count,
                                                    const void *idempotency_key,
                                                    size_t idempotency_key_size,
                                                    uint32_t retention_seconds,
                                                    bool *replayed);
GeoClientStatus geo_client_insert(GeoClient *client,
                                  uint64_t object_id,
                                  double latitude,
                                  double longitude,
                                  const void *document,
                                  size_t document_size);
GeoClientStatus geo_client_insert_generated(GeoClient *client,
                                            double latitude,
                                            double longitude,
                                            const void *document,
                                            size_t document_size,
                                            uint64_t *object_id);
GeoClientStatus geo_client_upsert(GeoClient *client,
                                  uint64_t object_id,
                                  double latitude,
                                  double longitude,
                                  const void *document,
                                  size_t document_size);
GeoClientStatus geo_client_update(GeoClient *client,
                                  uint64_t object_id,
                                  double latitude,
                                  double longitude,
                                  const void *document,
                                  size_t document_size);
GeoClientStatus geo_client_delete(GeoClient *client, uint64_t object_id);
GeoClientStatus geo_client_get(GeoClient *client, uint64_t object_id, GeoDatabaseObject *object);
GeoClientStatus geo_client_create_index(GeoClient *client,
                                        const char *name,
                                        const char *json_pointer,
                                        GeoDatabaseIndexType type);
GeoClientStatus geo_client_drop_index(GeoClient *client, const char *name);
GeoClientStatus geo_client_list_indexes(GeoClient *client,
                                       GeoDatabaseIndexInfo *indexes,
                                       size_t capacity,
                                       size_t *count);
GeoClientStatus geo_client_query_index_reuse(GeoClient *client,
                                             const GeoDatabaseIndexPredicate *predicate,
                                             GeoIdResult *result,
                                             GeoDatabaseIndexQueryStats *stats);
GeoClientStatus geo_client_query_radius_reuse(GeoClient *client,
                                              const GeoDatabaseRadiusQuery *query,
                                              GeoSearchResult *result,
                                              GeoDatabaseQueryStats *stats);
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
