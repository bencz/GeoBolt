#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "geobolt/client.h"

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static bool parse_port(const char *text, uint16_t *port)
{
    errno = 0;
    char *end = NULL;
    unsigned long value = strtoul(text, &end, 10);

    if (errno != 0 || !text[0] || !end || *end != '\0' || value == 0 || value > UINT16_MAX) {
        return false;
    }

    *port = (uint16_t) value;

    return true;
}

int main(int argc, char **argv)
{
    uint16_t port;

    if (argc != 3 || !parse_port(argv[2], &port)) {
        fprintf(stderr, "usage: GEOBOLT_TOKEN=secret %s HOST PORT\n", argv[0]);
        return EXIT_FAILURE;
    }

    GeoClientConfig config = geo_client_default_config(argv[1], port);

    config.authentication_token = getenv("GEOBOLT_TOKEN");
    GeoClientStatus status;
    GeoClient *client = geo_client_connect(&config, &status);

    if (!client) {
        fprintf(stderr, "connection failed with client status %d\n", status);
        return EXIT_FAILURE;
    }

    static const char ping_payload[] = "client-sample";

    status = geo_client_ping(client, ping_payload, sizeof(ping_payload) - 1U);

    if (status != GEO_CLIENT_OK) {
        fprintf(stderr, "ping failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    const uint64_t object_id = UINT64_C(123456789);

    status = geo_client_insert(client, object_id, -23.5505, -46.6333, NULL, 0U);

    if (status != GEO_CLIENT_OK) {
        fprintf(stderr, "insert failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    GeoDatabaseObject object;

    status = geo_client_get(client, object_id, &object);

    if (status != GEO_CLIENT_OK) {
        fprintf(stderr, "get failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    GeoPoint location = geo_decode(object.morton_code);
    printf("inserted id=%" PRIu64 " sequence=%" PRIu64 " at %.6f, %.6f\n",
           object.object_id,
           object.sequence,
           location.lat,
           location.lng);
    geo_database_object_release(&object);

    status = geo_client_update(client, object_id, -22.9068, -43.1729, NULL, 0U);

    if (status != GEO_CLIENT_OK) {
        fprintf(stderr, "update failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    status = geo_client_upsert(client, object_id, -22.9070, -43.1731, NULL, 0U);

    if (status != GEO_CLIENT_OK) {
        fprintf(stderr, "upsert failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    GeoDatabaseStats stats;

    status = geo_client_get_database_stats(client, &stats);

    if (status != GEO_CLIENT_OK) {
        fprintf(stderr, "stats failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    printf("connected: committed=%llu physical_records=%llu active_segments=%zu\n",
           (unsigned long long) stats.committed_operations,
           (unsigned long long) stats.physical_records,
           stats.active_segments);

    status = geo_client_delete(client, object_id);

    if (status != GEO_CLIENT_OK) {
        fprintf(stderr, "delete failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    uint64_t generated_object_id;
    status = geo_client_insert_generated(client, 40.7128, -74.0060, NULL, 0U, &generated_object_id);

    if (status != GEO_CLIENT_OK) {
        fprintf(stderr, "generated insert failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    printf("server generated generic object id=%" PRIu64 "\n", generated_object_id);

    status = geo_client_delete(client, generated_object_id);

    if (status != GEO_CLIENT_OK) {
        fprintf(stderr, "generated object delete failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    GeoDatabaseObjectMutation retry_safe_write = {
        .object_id = UINT64_C(987654321),
        .morton_code = geo_encode(51.5074, -0.1278),
        .operation = GEO_DATABASE_UPSERT,
    };
    static const char retry_key[] = "client-usage-import-0001";
    bool replayed;

    status = geo_client_write_objects_idempotent(client,
                                                 &retry_safe_write,
                                                 1U,
                                                 retry_key,
                                                 sizeof(retry_key) - 1U,
                                                 3600U,
                                                 &replayed);

    if (status != GEO_CLIENT_OK || replayed) {
        fprintf(stderr, "first idempotent write failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    status = geo_client_write_objects_idempotent(client,
                                                 &retry_safe_write,
                                                 1U,
                                                 retry_key,
                                                 sizeof(retry_key) - 1U,
                                                 3600U,
                                                 &replayed);

    if (status != GEO_CLIENT_OK || !replayed) {
        fprintf(stderr, "idempotent replay failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    printf("idempotent retry was replayed without another commit\n");
    status = geo_client_delete(client, retry_safe_write.object_id);

    if (status != GEO_CLIENT_OK) {
        fprintf(stderr, "retry-safe object delete failed with client status %d\n", status);
        geo_client_close(client);
        return EXIT_FAILURE;
    }

    geo_client_close(client);

    return EXIT_SUCCESS;
}
