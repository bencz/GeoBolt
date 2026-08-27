#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "geobolt/client.h"

#include <errno.h>
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

    geo_client_close(client);

    return EXIT_SUCCESS;
}
