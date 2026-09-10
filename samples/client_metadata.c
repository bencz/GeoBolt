#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "geobolt/client.h"
#include "geobolt/geodoc.h"

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

    if (errno != 0 || !text[0] || !end || *end != '\0' || value == 0U || value > UINT16_MAX) {
        return false;
    }

    *port = (uint16_t) value;
    return true;
}

static GeoDocStatus build_vehicle_metadata(GeoDocBuilder *builder, GeoDocBuffer *document)
{
    GeoDocStatus status = geo_doc_builder_add_string(builder, 0U, "kind", "vehicle", 7U, NULL);

    if (status == GEO_DOC_OK) {
        status = geo_doc_builder_add_uint64(builder, 0U, "vehicle_id", 4821U, NULL);
    }
    if (status == GEO_DOC_OK) {
        status = geo_doc_builder_add_string(builder, 0U, "plate", "ABC1D23", 7U, NULL);
    }
    if (status == GEO_DOC_OK) {
        status = geo_doc_builder_add_uint64(builder, 0U, "driver_id", 998877U, NULL);
    }
    if (status == GEO_DOC_OK) {
        status = geo_doc_builder_add_string(builder, 0U, "driver_name", "Ana", 3U, NULL);
    }
    if (status == GEO_DOC_OK) {
        status = geo_doc_builder_add_bool(builder, 0U, "online", true, NULL);
    }

    return status == GEO_DOC_OK ? geo_doc_builder_finish(builder, document) : status;
}

int main(int argc, char **argv)
{
    uint16_t port;

    if (argc != 3 || !parse_port(argv[2], &port)) {
        fprintf(stderr, "usage: GEOBOLT_TOKEN=secret %s HOST PORT\n", argv[0]);
        return EXIT_FAILURE;
    }

    GeoDocBuilder *builder = geo_doc_builder_create();
    GeoDocBuffer metadata = { 0 };

    if (!builder || build_vehicle_metadata(builder, &metadata) != GEO_DOC_OK) {
        fputs("failed to build GeoDoc metadata\n", stderr);
        geo_doc_builder_destroy(builder);
        return EXIT_FAILURE;
    }

    GeoClientConfig config = geo_client_default_config(argv[1], port);
    config.authentication_token = getenv("GEOBOLT_TOKEN");
    GeoClientStatus status;
    GeoClient *client = geo_client_connect(&config, &status);

    if (!client) {
        fprintf(stderr, "connection failed with client status %d\n", status);
        geo_doc_buffer_release(&metadata);
        geo_doc_builder_destroy(builder);
        return EXIT_FAILURE;
    }

    const uint64_t object_id = UINT64_C(4821);
    status = geo_client_upsert(client, object_id, -23.5505, -46.6333, metadata.data, metadata.size);

    if (status == GEO_CLIENT_OK) {
        GeoClientStatus index_status = geo_client_create_index(client,
                                                               "online_idx",
                                                               "/online",
                                                               GEO_DATABASE_INDEX_BOOL);

        if (index_status != GEO_CLIENT_OK && index_status != GEO_CLIENT_ALREADY_EXISTS) {
            status = index_status;
        }
    }

    if (status == GEO_CLIENT_OK) {
        GeoDatabaseObject object;
        status = geo_client_get(client, object_id, &object);

        if (status == GEO_CLIENT_OK) {
            GeoDocView view;
            GeoDocValue driver_name;
            const void *name;
            size_t name_size;

            if (geo_doc_open(object.document, object.document_size, &view) == GEO_DOC_OK &&
                geo_doc_find_pointer(view, "/driver_name", &driver_name) == GEO_DOC_OK &&
                geo_doc_value_data(driver_name, &name, &name_size) == GEO_DOC_OK) {
                printf("object=%" PRIu64 " driver=%.*s metadata_bytes=%zu\n",
                       object.object_id,
                       (int) name_size,
                       (const char *) name,
                       object.document_size);
            } else {
                status = GEO_CLIENT_PROTOCOL_ERROR;
            }

            geo_database_object_release(&object);
        }
    }

    if (status == GEO_CLIENT_OK) {
        GeoIdResult *online_objects = geo_id_result_create(16U);
        GeoSearchResult *nearby_online = geo_result_create(16U);
        GeoDatabaseIndexPredicate predicate = {
            .index_name = "online_idx",
            .operation = GEO_DATABASE_INDEX_EQUAL,
            .lower = { .type = GEO_DATABASE_INDEX_BOOL, .as.boolean = true },
        };

        if (!online_objects || !nearby_online) {
            status = GEO_CLIENT_OUT_OF_MEMORY;
        } else {
            status = geo_client_query_index_reuse(client, &predicate, online_objects, NULL);

            if (status == GEO_CLIENT_OK) {
                printf("online objects=%zu\n", online_objects->count);
            }

            GeoDatabaseRadiusQuery query = {
                .predicates = &predicate,
                .predicate_count = 1U,
                .latitude = -23.5505,
                .longitude = -46.6333,
                .radius_km = 5.0,
            };
            GeoDatabaseQueryStats query_stats;

            if (status == GEO_CLIENT_OK) {
                status = geo_client_query_radius_reuse(client, &query, nearby_online, &query_stats);
            }
            if (status == GEO_CLIENT_OK) {
                printf("nearby online objects=%zu plan=%u\n", nearby_online->count, (unsigned) query_stats.plan);
            }
        }

        geo_result_destroy(nearby_online);
        geo_id_result_destroy(online_objects);
    }

    if (status != GEO_CLIENT_OK) {
        fprintf(stderr, "metadata operation failed with client status %d\n", status);
    }

    geo_client_close(client);
    geo_doc_buffer_release(&metadata);
    geo_doc_builder_destroy(builder);
    return status == GEO_CLIENT_OK ? EXIT_SUCCESS : EXIT_FAILURE;
}
