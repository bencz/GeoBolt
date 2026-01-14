#include "geo_index.h"
#include <stdio.h>
#include <inttypes.h>
#include <math.h>

int main(void) {
    printf("=== GeoIndex Demo ===\n\n");

    // =========================================================
    // Test 1: Encode / Decode
    // =========================================================
    printf("--- Test 1: Encode / Decode ---\n");

    double lat = -23.5614123;
    double lng = -46.6558819;

    uint64_t z = geo_encode(lat, lng);
    GeoPoint p = geo_decode(z);

    printf("Encoded : %" PRIu64 "\n", z);
    printf("Original: %.7f %.7f\n", lat, lng);
    printf("Decoded : %.7f %.7f\n", p.lat, p.lng);
    printf("Delta   : %.9f %.9f\n", fabs(lat - p.lat), fabs(lng - p.lng));
    printf("Error   : %.4f meters\n\n", geo_haversine_m(lat, lng, p.lat, p.lng));

    // =========================================================
    // Test 2: Spatial Search
    // =========================================================
    printf("--- Test 2: Spatial Search ---\n");

    GeoIndex *index = geo_index_create(100);
    
    geo_index_add(index, 10, -23.5505200, -46.6333090);  // Sao Paulo center
    geo_index_add(index, 11, -23.5590000, -46.6400000);  // Nearby
    geo_index_add(index, 12, -23.5874162, -46.6576336);  // ~5km away
    geo_index_add(index, 13, -23.4542000, -46.5333000);  // ~15km away
    geo_index_add(index, 14, -22.9068000, -43.1729000);  // Rio de Janeiro

    geo_index_build(index);

    double search_lat = -23.5505200;
    double search_lng = -46.6333090;
    double radius_km = 5.0;

    printf("Searching within %.1f km of (%.7f, %.7f):\n", 
           radius_km, search_lat, search_lng);

    GeoSearchStats stats;
    GeoSearchResult *result = geo_search_radius(index, search_lat, search_lng, radius_km, &stats);

    for (size_t i = 0; i < result->count; i++) {
        GeoPoint pt = geo_decode(result->results[i].z);
        double dist = geo_haversine_km(search_lat, search_lng, pt.lat, pt.lng);
        printf("  FOUND id=%" PRIu64 " lat=%.7f lng=%.7f dist=%.3f km\n",
               result->results[i].id, pt.lat, pt.lng, dist);
    }

    printf("\nStats: scanned=%" PRIu64 " matched=%" PRIu64 " time=%.3f ms\n",
           stats.records_scanned, stats.records_matched, stats.search_time_ms);

    geo_result_destroy(result);

    // =========================================================
    // Test 3: K-Nearest Neighbors
    // =========================================================
    printf("\n--- Test 3: K-Nearest Neighbors ---\n");

    result = geo_search_knn(index, search_lat, search_lng, 3, 1000.0, &stats);

    printf("3 nearest neighbors:\n");
    for (size_t i = 0; i < result->count; i++) {
        GeoPoint pt = geo_decode(result->results[i].z);
        double dist = geo_haversine_km(search_lat, search_lng, pt.lat, pt.lng);
        printf("  %zu. id=%" PRIu64 " dist=%.3f km\n", i + 1, result->results[i].id, dist);
    }

    geo_result_destroy(result);

    // =========================================================
    // Test 4: Bounding Box Search
    // =========================================================
    printf("\n--- Test 4: Bounding Box Search ---\n");

    result = geo_search_bbox(index, -24.0, -23.0, -47.0, -46.0, &stats);

    printf("Points in bounding box [-24,-23] x [-47,-46]:\n");
    for (size_t i = 0; i < result->count; i++) {
        GeoPoint pt = geo_decode(result->results[i].z);
        printf("  id=%" PRIu64 " lat=%.7f lng=%.7f\n",
               result->results[i].id, pt.lat, pt.lng);
    }

    geo_result_destroy(result);
    geo_index_destroy(index);

    printf("\n=== Demo Complete ===\n");
    return 0;
}
