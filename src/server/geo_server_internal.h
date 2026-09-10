#ifndef GEO_SERVER_INTERNAL_H
#define GEO_SERVER_INTERNAL_H

#include "geobolt/geobolt.h"
#include "geobolt/server.h"
#include "geo_job_queue.h"
#include "geo_protocol.h"

#include <pthread.h>
#include <stdatomic.h>

typedef enum {
    GEO_EVENT_LISTENER,
    GEO_EVENT_WAKEUP,
    GEO_EVENT_CONNECTION,
} GeoEventType;

typedef struct {
    GeoEventType type;
} GeoEventSource;

typedef struct GeoConnection GeoConnection;

struct GeoConnection {
    GeoEventSource source;
    int descriptor;
    struct GeoServer *server;
    GeoDatabaseQueryWorkspace *query_workspace;
    GeoSearchResult *query_result;
    GeoConnection *previous;
    GeoConnection *next;
    GeoConnection *completion_next;
    unsigned char *request_payload;
    unsigned char *response;
    GeoProtocolHeader request_header;
    size_t header_received;
    size_t payload_received;
    size_t response_size;
    size_t response_sent;
    unsigned char encoded_header[GEO_PROTOCOL_HEADER_SIZE];
    bool header_decoded;
    bool authenticated;
    bool job_inflight;
    bool closing;
    bool close_after_response;
};

_Static_assert(sizeof(GeoConnection) == 192U, "Connection state must remain compact enough for high connection counts");

struct GeoServer {
    GeoEventSource listener_source;
    GeoEventSource wakeup_source;
    GeoDatabase *database;
    GeoJobQueue *work_queue;
    GeoConnection *connections;
    GeoConnection *completed_head;
    GeoConnection *completed_tail;
    GeoConnection *retired_connections;
    pthread_mutex_t completion_lock;
    int epoll_descriptor;
    int listener_descriptor;
    int wakeup_descriptor;
    char *authentication_token;
    size_t authentication_token_size;
    size_t max_connections;
    size_t max_frame_size;
    size_t max_inflight_payload_bytes;
    uint16_t port;
    atomic_bool stop_requested;
    atomic_bool running;
    atomic_bool event_loop_failed;
    atomic_uint_fast64_t accepted_connections;
    atomic_uint_fast64_t rejected_connections;
    atomic_uint_fast64_t active_connections;
    atomic_uint_fast64_t completed_requests;
    atomic_uint_fast64_t protocol_errors;
    atomic_uint_fast64_t authentication_failures;
    atomic_uint_fast64_t backpressure_rejections;
    atomic_uint_fast64_t inflight_payload_bytes;
    atomic_uint_fast64_t peak_inflight_payload_bytes;
    bool completion_lock_initialized;
};

bool geo_server_prepare_response(GeoConnection *connection,
                                 GeoProtocolStatus status,
                                 const void *payload,
                                 uint32_t payload_size);
void geo_server_complete_job(GeoConnection *connection);
void geo_server_execute_request(void *context);

#endif
