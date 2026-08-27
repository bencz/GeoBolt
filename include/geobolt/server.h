#ifndef GEOBOLT_SERVER_H
#define GEOBOLT_SERVER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GeoServer GeoServer;

typedef enum {
    GEO_SERVER_OK = 0,
    GEO_SERVER_INVALID_ARGUMENT,
    GEO_SERVER_OUT_OF_MEMORY,
    GEO_SERVER_DATABASE_ERROR,
    GEO_SERVER_WORKER_ERROR,
    GEO_SERVER_SOCKET_ERROR,
    GEO_SERVER_SYNCHRONIZATION_ERROR,
    GEO_SERVER_DESCRIPTOR_ERROR,
    GEO_SERVER_EVENT_REGISTRATION_ERROR,
    GEO_SERVER_EVENT_LOOP_ERROR,
    GEO_SERVER_ALREADY_RUNNING,
} GeoServerStatus;

typedef struct {
    const char *database_directory;
    const char *bind_address;
    const char *authentication_token;
    size_t worker_threads;
    size_t work_queue_capacity;
    size_t max_connections;
    size_t max_frame_size;
    size_t max_inflight_payload_bytes;
    int listen_backlog;
    uint16_t port;
    bool create_database_if_missing;
} GeoServerConfig;

typedef struct {
    uint64_t accepted_connections;
    uint64_t rejected_connections;
    uint64_t active_connections;
    uint64_t completed_requests;
    uint64_t protocol_errors;
    uint64_t authentication_failures;
    uint64_t backpressure_rejections;
    uint64_t current_inflight_payload_bytes;
    uint64_t peak_inflight_payload_bytes;
} GeoServerStats;

GeoServerConfig geo_server_default_config(const char *database_directory);
GeoServer *geo_server_create(const GeoServerConfig *config, GeoServerStatus *status);

// The caller must stop and join the thread running geo_server_run before destruction.
void geo_server_destroy(GeoServer *server);

// Runs the single reactor in the calling thread. Worker threads execute database commands outside the reactor.
GeoServerStatus geo_server_run(GeoServer *server);
void geo_server_request_stop(GeoServer *server);

uint16_t geo_server_port(const GeoServer *server);
bool geo_server_get_stats(const GeoServer *server, GeoServerStats *stats);

#ifdef __cplusplus
}
#endif

#endif
