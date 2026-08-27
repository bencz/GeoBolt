#include "geobolt/geo_index.h"

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static bool print_radius_results(const GeoIndex *index,
                                 double search_latitude,
                                 double search_longitude,
                                 double radius_km)
{
    GeoSearchStats stats;
    GeoSearchResult *result = geo_search_radius(index, search_latitude, search_longitude, radius_km, &stats);

    if (!result) {
        return false;
    }

    printf("Searching within %.1f km of (%.7f, %.7f):\n", radius_km, search_latitude, search_longitude);

    for (size_t index_position = 0; index_position < result->count; ++index_position) {
        GeoPoint point = geo_decode(result->results[index_position].z);
        double distance_km = geo_haversine_km(search_latitude, search_longitude, point.lat, point.lng);

        printf("  FOUND id=%" PRIu64 " lat=%.7f lng=%.7f dist=%.3f km\n",
               result->results[index_position].id,
               point.lat,
               point.lng,
               distance_km);
    }

    printf("Stats: scanned=%" PRIu64 " matched=%" PRIu64 " time=%.3f ms\n",
           stats.records_scanned,
           stats.records_matched,
           stats.search_time_ms);

    geo_result_destroy(result);

    return true;
}

static bool print_nearest_results(const GeoIndex *index, double search_latitude, double search_longitude)
{
    GeoSearchResult *result = geo_search_knn(index, search_latitude, search_longitude, 3U, 1000.0, NULL);

    if (!result) {
        return false;
    }

    puts("3 nearest neighbors:");

    for (size_t index_position = 0; index_position < result->count; ++index_position) {
        GeoPoint point = geo_decode(result->results[index_position].z);
        double distance_km = geo_haversine_km(search_latitude, search_longitude, point.lat, point.lng);

        printf("  %zu. id=%" PRIu64 " dist=%.3f km\n",
               index_position + 1U,
               result->results[index_position].id,
               distance_km);
    }

    geo_result_destroy(result);

    return true;
}

static bool print_bounding_box_results(const GeoIndex *index)
{
    GeoSearchResult *result = geo_search_bbox(index, -24.0, -23.0, -47.0, -46.0, NULL);

    if (!result) {
        return false;
    }

    puts("Points in bounding box [-24,-23] x [-47,-46]:");

    for (size_t index_position = 0; index_position < result->count; ++index_position) {
        GeoPoint point = geo_decode(result->results[index_position].z);

        printf("  id=%" PRIu64 " lat=%.7f lng=%.7f\n",
               result->results[index_position].id,
               point.lat,
               point.lng);
    }

    geo_result_destroy(result);

    return true;
}

int main(void)
{
    const double latitude = -23.5614123;
    const double longitude = -46.6558819;
    uint64_t morton_code = geo_encode(latitude, longitude);
    GeoPoint decoded = geo_decode(morton_code);

    puts("=== GeoIndex Demo ===");
    puts("\n--- Encode / Decode ---");
    printf("Encoded : %" PRIu64 "\n", morton_code);
    printf("Original: %.7f %.7f\n", latitude, longitude);
    printf("Decoded : %.7f %.7f\n", decoded.lat, decoded.lng);
    printf("Delta   : %.9f %.9f\n", fabs(latitude - decoded.lat), fabs(longitude - decoded.lng));
    printf("Error   : %.4f meters\n", geo_haversine_m(latitude, longitude, decoded.lat, decoded.lng));

    GeoIndex *index = geo_index_create(100U);

    if (!index) {
        fputs("Failed to create the index\n", stderr);
        return EXIT_FAILURE;
    }

    bool populated = geo_index_add(index, 10U, -23.5505200, -46.6333090) &&
                     geo_index_add(index, 11U, -23.5590000, -46.6400000) &&
                     geo_index_add(index, 12U, -23.5874162, -46.6576336) &&
                     geo_index_add(index, 13U, -23.4542000, -46.5333000) &&
                     geo_index_add(index, 14U, -22.9068000, -43.1729000) &&
                     geo_index_build(index);

    if (!populated) {
        fputs("Failed to populate and build the index\n", stderr);
        geo_index_destroy(index);
        return EXIT_FAILURE;
    }

    const double search_latitude = -23.5505200;
    const double search_longitude = -46.6333090;
    bool succeeded = print_radius_results(index, search_latitude, search_longitude, 5.0) &&
                     print_nearest_results(index, search_latitude, search_longitude) &&
                     print_bounding_box_results(index);

    geo_index_destroy(index);

    if (!succeeded) {
        fputs("A sample query failed\n", stderr);
        return EXIT_FAILURE;
    }

    puts("\n=== Demo Complete ===");

    return EXIT_SUCCESS;
}
