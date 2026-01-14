#include "geo_index.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <inttypes.h>
#include <assert.h>

// =========================================================
// Test Framework
// =========================================================

#define TEST_PASSED "\033[32mPASSED\033[0m"
#define TEST_FAILED "\033[31mFAILED\033[0m"

static int g_tests_run = 0;
static int g_tests_passed = 0;
static int g_tests_failed = 0;

#define ASSERT_TRUE(cond, msg) do { \
    if (!(cond)) { \
        printf("  %s: %s (line %d)\n", TEST_FAILED, msg, __LINE__); \
        return 0; \
    } \
} while(0)

#define ASSERT_FALSE(cond, msg) ASSERT_TRUE(!(cond), msg)
#define ASSERT_EQ(a, b, msg) ASSERT_TRUE((a) == (b), msg)

#define ASSERT_NEAR(a, b, eps, msg) do { \
    double _diff = fabs((double)(a) - (double)(b)); \
    if (_diff > (eps)) { \
        printf("  %s: %s (expected %.10f, got %.10f, diff=%.10e, line %d)\n", \
               TEST_FAILED, msg, (double)(b), (double)(a), _diff, __LINE__); \
        return 0; \
    } \
} while(0)

#define RUN_TEST(test_func) do { \
    g_tests_run++; \
    printf("\n[TEST %d] %s\n", g_tests_run, #test_func); \
    if (test_func()) { \
        printf("  %s\n", TEST_PASSED); \
        g_tests_passed++; \
    } else { \
        g_tests_failed++; \
    } \
} while(0)

// =========================================================
// Test Data - Known Geographic Points
// =========================================================

typedef struct {
    const char *name;
    double lat;
    double lng;
} NamedPoint;

static const NamedPoint KNOWN_CITIES[] = {
    {"Sao Paulo",     -23.5505200, -46.6333090},
    {"Rio de Janeiro",-22.9068467, -43.1728965},
    {"New York",       40.7127753, -74.0059728},
    {"London",         51.5073509,  -0.1277583},
    {"Tokyo",          35.6761919, 139.6503106},
    {"Sydney",        -33.8688197, 151.2092955},
    {"Paris",          48.8566140,   2.3522219},
    {"Moscow",         55.7558260,  37.6172999},
    {"Dubai",          25.2048493,  55.2707828},
    {"Singapore",       1.3521150, 103.8198422},
    {"Cape Town",     -33.9248685,  18.4240553},
    {"Buenos Aires",  -34.6036844, -58.3815591},
    {"Mumbai",         19.0759837,  72.8776559},
    {"Beijing",        39.9041999, 116.4073963},
    {"Los Angeles",    34.0522342,-118.2436849},
};

static const int NUM_KNOWN_CITIES = sizeof(KNOWN_CITIES) / sizeof(KNOWN_CITIES[0]);

typedef struct {
    int city1_idx;
    int city2_idx;
    double distance_km;
} KnownDistance;

static const KnownDistance KNOWN_DISTANCES[] = {
    {0, 1, 357.0},
    {2, 6, 5837.0},
    {3, 7, 2500.0},
    {4, 9, 5312.0},
};

// =========================================================
// Unit Tests - Morton Code
// =========================================================

int test_spread_compact_bits_roundtrip(void) {
    uint32_t test_values[] = {0, 1, 255, 65535, 0xFFFFFFFF, 0x12345678, 0xDEADBEEF};
    
    for (size_t i = 0; i < sizeof(test_values)/sizeof(test_values[0]); i++) {
        uint32_t original = test_values[i];
        uint64_t spread = geo_spread_bits(original);
        uint32_t recovered = geo_compact_bits(spread);
        ASSERT_EQ(original, recovered, "spread/compact roundtrip failed");
    }
    return 1;
}

int test_spread_bits_pattern(void) {
    uint32_t val = 3;
    uint64_t spread = geo_spread_bits(val);
    ASSERT_EQ(spread & 0xF, 0x5ULL, "spread pattern incorrect");
    
    val = 15;
    spread = geo_spread_bits(val);
    ASSERT_EQ(spread & 0xFF, 0x55ULL, "spread pattern incorrect for 0xF");
    return 1;
}

// =========================================================
// Unit Tests - Encode/Decode
// =========================================================

int test_encode_decode_roundtrip(void) {
    for (int i = 0; i < NUM_KNOWN_CITIES; i++) {
        double lat = KNOWN_CITIES[i].lat;
        double lng = KNOWN_CITIES[i].lng;
        
        uint64_t z = geo_encode(lat, lng);
        GeoPoint p = geo_decode(z);
        
        ASSERT_NEAR(p.lat, lat, 0.00005, "latitude decode precision");
        ASSERT_NEAR(p.lng, lng, 0.00005, "longitude decode precision");
    }
    return 1;
}

int test_encode_decode_edge_cases(void) {
    double edge_cases[][2] = {
        {GEO_MIN_LAT, GEO_MIN_LNG},
        {GEO_MAX_LAT, GEO_MAX_LNG},
        {GEO_MIN_LAT, GEO_MAX_LNG},
        {GEO_MAX_LAT, GEO_MIN_LNG},
        {0.0, 0.0},
        {0.0, 180.0},
        {0.0, -180.0},
        {90.0, 0.0},
        {-90.0, 0.0},
    };
    
    for (size_t i = 0; i < sizeof(edge_cases)/sizeof(edge_cases[0]); i++) {
        double lat = edge_cases[i][0];
        double lng = edge_cases[i][1];
        
        uint64_t z = geo_encode(lat, lng);
        GeoPoint p = geo_decode(z);
        
        ASSERT_NEAR(p.lat, lat, 0.0001, "edge case latitude");
        ASSERT_NEAR(p.lng, lng, 0.0001, "edge case longitude");
    }
    return 1;
}

int test_encode_clamping(void) {
    uint64_t z1 = geo_encode(-100.0, -200.0);
    uint64_t z2 = geo_encode(GEO_MIN_LAT, GEO_MIN_LNG);
    
    GeoPoint p1 = geo_decode(z1);
    GeoPoint p2 = geo_decode(z2);
    
    ASSERT_NEAR(p1.lat, p2.lat, 0.0001, "clamping latitude");
    ASSERT_NEAR(p1.lng, p2.lng, 0.0001, "clamping longitude");
    return 1;
}

int test_encode_ordering(void) {
    double base_lat = -23.5505;
    double base_lng = -46.6333;
    
    uint64_t z_base = geo_encode(base_lat, base_lng);
    uint64_t z_near = geo_encode(base_lat + 0.001, base_lng + 0.001);
    uint64_t z_far = geo_encode(base_lat + 10.0, base_lng + 10.0);
    
    uint64_t diff_near = (z_near > z_base) ? (z_near - z_base) : (z_base - z_near);
    uint64_t diff_far = (z_far > z_base) ? (z_far - z_base) : (z_base - z_far);
    
    printf("    Z-diff near: %" PRIu64 ", Z-diff far: %" PRIu64 "\n", diff_near, diff_far);
    return 1;
}

// =========================================================
// Unit Tests - Distance Calculations
// =========================================================

int test_haversine_known_distances(void) {
    for (size_t i = 0; i < sizeof(KNOWN_DISTANCES)/sizeof(KNOWN_DISTANCES[0]); i++) {
        int idx1 = KNOWN_DISTANCES[i].city1_idx;
        int idx2 = KNOWN_DISTANCES[i].city2_idx;
        double expected = KNOWN_DISTANCES[i].distance_km;
        
        double calculated = geo_haversine_km(
            KNOWN_CITIES[idx1].lat, KNOWN_CITIES[idx1].lng,
            KNOWN_CITIES[idx2].lat, KNOWN_CITIES[idx2].lng
        );
        
        double tolerance = expected * 0.05;
        printf("    %s to %s: %.1f km (expected ~%.1f km)\n",
               KNOWN_CITIES[idx1].name, KNOWN_CITIES[idx2].name,
               calculated, expected);
        
        ASSERT_NEAR(calculated, expected, tolerance, "haversine distance");
    }
    return 1;
}

int test_haversine_zero_distance(void) {
    double dist = geo_haversine_km(0.0, 0.0, 0.0, 0.0);
    ASSERT_NEAR(dist, 0.0, 0.0001, "same point distance should be 0");
    
    dist = geo_haversine_km(-23.5505, -46.6333, -23.5505, -46.6333);
    ASSERT_NEAR(dist, 0.0, 0.0001, "same point distance should be 0");
    return 1;
}

int test_haversine_symmetry(void) {
    for (int i = 0; i < NUM_KNOWN_CITIES - 1; i++) {
        double d1 = geo_haversine_km(
            KNOWN_CITIES[i].lat, KNOWN_CITIES[i].lng,
            KNOWN_CITIES[i+1].lat, KNOWN_CITIES[i+1].lng
        );
        double d2 = geo_haversine_km(
            KNOWN_CITIES[i+1].lat, KNOWN_CITIES[i+1].lng,
            KNOWN_CITIES[i].lat, KNOWN_CITIES[i].lng
        );
        ASSERT_NEAR(d1, d2, 0.0001, "haversine should be symmetric");
    }
    return 1;
}

int test_haversine_triangle_inequality(void) {
    for (int i = 0; i < NUM_KNOWN_CITIES - 2; i++) {
        double ab = geo_haversine_km(
            KNOWN_CITIES[i].lat, KNOWN_CITIES[i].lng,
            KNOWN_CITIES[i+1].lat, KNOWN_CITIES[i+1].lng
        );
        double bc = geo_haversine_km(
            KNOWN_CITIES[i+1].lat, KNOWN_CITIES[i+1].lng,
            KNOWN_CITIES[i+2].lat, KNOWN_CITIES[i+2].lng
        );
        double ac = geo_haversine_km(
            KNOWN_CITIES[i].lat, KNOWN_CITIES[i].lng,
            KNOWN_CITIES[i+2].lat, KNOWN_CITIES[i+2].lng
        );
        ASSERT_TRUE(ac <= ab + bc + 0.001, "triangle inequality violated");
    }
    return 1;
}

int test_fast_distance_accuracy(void) {
    double lat1 = -23.5505, lng1 = -46.6333;
    double lat2 = -23.5600, lng2 = -46.6400;
    
    double precise = geo_haversine_km(lat1, lng1, lat2, lng2);
    double fast = geo_fast_distance_km(lat1, lng1, lat2, lng2);
    
    double error = fabs(precise - fast) / precise * 100.0;
    printf("    Precise: %.4f km, Fast: %.4f km, Error: %.2f%%\n", precise, fast, error);
    ASSERT_TRUE(error < 10.0, "fast distance error too high for short distance");
    return 1;
}

// =========================================================
// Unit Tests - Index Operations
// =========================================================

int test_index_create_destroy(void) {
    GeoIndex *index = geo_index_create(100);
    ASSERT_TRUE(index != NULL, "index creation failed");
    ASSERT_EQ(index->count, 0, "new index should be empty");
    ASSERT_EQ(index->capacity, 100, "capacity should match");
    ASSERT_FALSE(index->sorted, "new index should not be sorted");
    geo_index_destroy(index);
    return 1;
}

int test_index_add_single(void) {
    GeoIndex *index = geo_index_create(10);
    ASSERT_TRUE(index != NULL, "index creation failed");
    
    bool added = geo_index_add(index, 1, -23.5505, -46.6333);
    ASSERT_TRUE(added, "add should succeed");
    ASSERT_EQ(index->count, 1, "count should be 1");
    geo_index_destroy(index);
    return 1;
}

int test_index_add_batch(void) {
    GeoIndex *index = geo_index_create(10);
    
    uint64_t ids[5] = {1, 2, 3, 4, 5};
    double lats[5] = {-23.5505, -22.9068, 40.7128, 51.5074, 35.6762};
    double lngs[5] = {-46.6333, -43.1729, -74.0060, -0.1278, 139.6503};
    
    bool added = geo_index_add_batch(index, ids, lats, lngs, 5);
    ASSERT_TRUE(added, "batch add should succeed");
    ASSERT_EQ(index->count, 5, "count should be 5");
    geo_index_destroy(index);
    return 1;
}

int test_index_auto_grow(void) {
    GeoIndex *index = geo_index_create(2);
    
    for (int i = 0; i < 100; i++) {
        bool added = geo_index_add(index, i, (double)(i % 180) - 90.0, (double)(i % 360) - 180.0);
        ASSERT_TRUE(added, "add should succeed with auto-grow");
    }
    
    ASSERT_EQ(index->count, 100, "count should be 100");
    ASSERT_TRUE(index->capacity >= 100, "capacity should have grown");
    geo_index_destroy(index);
    return 1;
}

int test_index_build_sorts(void) {
    GeoIndex *index = geo_index_create(10);
    
    geo_index_add(index, 5, 50.0, 50.0);
    geo_index_add(index, 1, -50.0, -50.0);
    geo_index_add(index, 3, 0.0, 0.0);
    geo_index_add(index, 2, -25.0, -25.0);
    geo_index_add(index, 4, 25.0, 25.0);
    
    ASSERT_FALSE(index->sorted, "should not be sorted before build");
    geo_index_build(index);
    ASSERT_TRUE(index->sorted, "should be sorted after build");
    
    for (size_t i = 1; i < index->count; i++) {
        ASSERT_TRUE(index->records[i-1].z <= index->records[i].z, "records not sorted");
    }
    geo_index_destroy(index);
    return 1;
}

// =========================================================
// Unit Tests - Binary Search
// =========================================================

int test_lower_bound_basic(void) {
    GeoIndex *index = geo_index_create(10);
    
    for (int i = 0; i < 10; i++) {
        geo_index_add(index, i, (double)i * 10.0 - 45.0, (double)i * 20.0 - 90.0);
    }
    geo_index_build(index);
    
    size_t pos = geo_lower_bound(index->records, index->count, index->records[5].z);
    ASSERT_EQ(pos, 5, "lower_bound should find exact match");
    
    pos = geo_lower_bound(index->records, index->count, 0);
    ASSERT_EQ(pos, 0, "lower_bound for 0 should return 0");
    
    pos = geo_lower_bound(index->records, index->count, UINT64_MAX);
    ASSERT_EQ(pos, 10, "lower_bound for MAX should return count");
    
    geo_index_destroy(index);
    return 1;
}

// =========================================================
// Unit Tests - Search Operations
// =========================================================

int test_search_radius_basic(void) {
    GeoIndex *index = geo_index_create(100);
    
    for (int i = 0; i < NUM_KNOWN_CITIES; i++) {
        geo_index_add(index, i, KNOWN_CITIES[i].lat, KNOWN_CITIES[i].lng);
    }
    geo_index_build(index);
    
    GeoSearchStats stats;
    GeoSearchResult *result = geo_search_radius(index, -23.5505, -46.6333, 500.0, &stats);
    
    ASSERT_TRUE(result != NULL, "search should return result");
    ASSERT_TRUE(result->count >= 1, "should find at least Sao Paulo");
    
    printf("    Found %zu results, scanned %" PRIu64 " records in %.3f ms\n",
           result->count, stats.records_scanned, stats.search_time_ms);
    
    geo_result_destroy(result);
    geo_index_destroy(index);
    return 1;
}

int test_search_radius_finds_nearby(void) {
    GeoIndex *index = geo_index_create(100);
    
    double center_lat = -23.5505;
    double center_lng = -46.6333;
    
    geo_index_add(index, 1, center_lat, center_lng);
    geo_index_add(index, 2, center_lat + 0.01, center_lng + 0.01);
    geo_index_add(index, 3, center_lat - 0.01, center_lng - 0.01);
    geo_index_add(index, 4, center_lat + 10.0, center_lng + 10.0);
    
    geo_index_build(index);
    
    GeoSearchResult *result = geo_search_radius(index, center_lat, center_lng, 5.0, NULL);
    
    ASSERT_TRUE(result != NULL, "search should return result");
    ASSERT_EQ(result->count, 3, "should find 3 nearby points");
    
    geo_result_destroy(result);
    geo_index_destroy(index);
    return 1;
}

int test_search_radius_empty_result(void) {
    GeoIndex *index = geo_index_create(100);
    
    geo_index_add(index, 1, 0.0, 0.0);
    geo_index_add(index, 2, 10.0, 10.0);
    geo_index_build(index);
    
    GeoSearchResult *result = geo_search_radius(index, 80.0, 80.0, 1.0, NULL);
    
    ASSERT_TRUE(result != NULL, "search should return result");
    ASSERT_EQ(result->count, 0, "should find no points");
    
    geo_result_destroy(result);
    geo_index_destroy(index);
    return 1;
}

int test_search_bbox_basic(void) {
    GeoIndex *index = geo_index_create(100);
    
    for (int i = 0; i < NUM_KNOWN_CITIES; i++) {
        geo_index_add(index, i, KNOWN_CITIES[i].lat, KNOWN_CITIES[i].lng);
    }
    geo_index_build(index);
    
    GeoSearchResult *result = geo_search_bbox(index, -35.0, -20.0, -60.0, -40.0, NULL);
    
    ASSERT_TRUE(result != NULL, "search should return result");
    printf("    Found %zu results in South America bbox\n", result->count);
    
    geo_result_destroy(result);
    geo_index_destroy(index);
    return 1;
}

int test_search_knn_basic(void) {
    GeoIndex *index = geo_index_create(100);
    
    for (int i = 0; i < NUM_KNOWN_CITIES; i++) {
        geo_index_add(index, i, KNOWN_CITIES[i].lat, KNOWN_CITIES[i].lng);
    }
    geo_index_build(index);
    
    GeoSearchResult *result = geo_search_knn(index, -23.5505, -46.6333, 3, 10000.0, NULL);
    
    ASSERT_TRUE(result != NULL, "search should return result");
    ASSERT_TRUE(result->count <= 3, "should return at most k results");
    
    printf("    Found %zu nearest neighbors\n", result->count);
    
    geo_result_destroy(result);
    geo_index_destroy(index);
    return 1;
}

// =========================================================
// Precision Tests
// =========================================================

int test_precision_at_equator(void) {
    double lat = 0.0;
    double lng = 0.0;
    
    uint64_t z = geo_encode(lat, lng);
    GeoPoint p = geo_decode(z);
    
    double error_m = geo_haversine_m(lat, lng, p.lat, p.lng);
    printf("    Error at equator: %.4f meters\n", error_m);
    ASSERT_TRUE(error_m < 2.0, "precision should be < 2m at equator");
    return 1;
}

int test_precision_at_poles(void) {
    double lat = 89.9999;
    double lng = 0.0;
    
    uint64_t z = geo_encode(lat, lng);
    GeoPoint p = geo_decode(z);
    
    double error_m = geo_haversine_m(lat, lng, p.lat, p.lng);
    printf("    Error near north pole: %.4f meters\n", error_m);
    ASSERT_TRUE(error_m < 5.0, "precision should be < 5m near poles");
    return 1;
}

int test_precision_random_points(void) {
    srand(42);
    double max_error = 0.0;
    double total_error = 0.0;
    int count = 10000;
    
    for (int i = 0; i < count; i++) {
        double lat = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        double lng = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
        
        uint64_t z = geo_encode(lat, lng);
        GeoPoint p = geo_decode(z);
        
        double error_m = geo_haversine_m(lat, lng, p.lat, p.lng);
        total_error += error_m;
        if (error_m > max_error) max_error = error_m;
    }
    
    double avg_error = total_error / count;
    printf("    Random points: avg error = %.4f m, max error = %.4f m\n", avg_error, max_error);
    ASSERT_TRUE(max_error < 5.0, "max error should be < 5m");
    ASSERT_TRUE(avg_error < 2.0, "avg error should be < 2m");
    return 1;
}

// =========================================================
// Performance Tests
// =========================================================

int test_perf_encode_decode(void) {
    int iterations = 1000000;
    
    double start = geo_get_time_ms();
    volatile uint64_t sum = 0;
    
    for (int i = 0; i < iterations; i++) {
        double lat = ((double)(i % 18000) / 100.0) - 90.0;
        double lng = ((double)(i % 36000) / 100.0) - 180.0;
        sum += geo_encode(lat, lng);
    }
    
    double encode_time = geo_get_time_ms() - start;
    
    start = geo_get_time_ms();
    volatile double lat_sum = 0;
    
    for (int i = 0; i < iterations; i++) {
        GeoPoint p = geo_decode((uint64_t)i * 12345);
        lat_sum += p.lat;
    }
    
    double decode_time = geo_get_time_ms() - start;
    
    printf("    Encode: %d ops in %.2f ms (%.2f M ops/sec)\n",
           iterations, encode_time, iterations / encode_time / 1000.0);
    printf("    Decode: %d ops in %.2f ms (%.2f M ops/sec)\n",
           iterations, decode_time, iterations / decode_time / 1000.0);
    
    ASSERT_TRUE(encode_time < 500.0, "encode too slow");
    ASSERT_TRUE(decode_time < 500.0, "decode too slow");
    return 1;
}

int test_perf_index_build(void) {
    int sizes[] = {1000, 10000, 100000, 1000000};
    
    for (size_t s = 0; s < sizeof(sizes)/sizeof(sizes[0]); s++) {
        int n = sizes[s];
        GeoIndex *index = geo_index_create(n);
        
        srand(42);
        for (int i = 0; i < n; i++) {
            double lat = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
            double lng = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
            geo_index_add(index, i, lat, lng);
        }
        
        double start = geo_get_time_ms();
        geo_index_build(index);
        double build_time = geo_get_time_ms() - start;
        
        printf("    Build %d records: %.2f ms\n", n, build_time);
        geo_index_destroy(index);
    }
    return 1;
}

int test_perf_search_radius(void) {
    int n = 1000000;
    GeoIndex *index = geo_index_create(n);
    
    srand(42);
    for (int i = 0; i < n; i++) {
        double lat = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        double lng = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
        geo_index_add(index, i, lat, lng);
    }
    geo_index_build(index);
    
    double radii[] = {1.0, 10.0, 100.0, 1000.0};
    int num_searches = 100;
    
    for (size_t r = 0; r < sizeof(radii)/sizeof(radii[0]); r++) {
        double radius = radii[r];
        double total_time = 0;
        size_t total_results = 0;
        
        srand(123);
        for (int i = 0; i < num_searches; i++) {
            double lat = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
            double lng = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
            
            GeoSearchStats stats;
            GeoSearchResult *result = geo_search_radius(index, lat, lng, radius, &stats);
            
            total_time += stats.search_time_ms;
            total_results += result->count;
            geo_result_destroy(result);
        }
        
        printf("    Radius %.0f km: avg %.3f ms, avg results %.1f\n",
               radius, total_time / num_searches, (double)total_results / num_searches);
    }
    
    geo_index_destroy(index);
    return 1;
}

int test_perf_binary_search(void) {
    int n = 10000000;
    GeoIndex *index = geo_index_create(n);
    
    srand(42);
    for (int i = 0; i < n; i++) {
        double lat = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        double lng = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
        geo_index_add(index, i, lat, lng);
    }
    geo_index_build(index);
    
    int searches = 1000000;
    double start = geo_get_time_ms();
    volatile size_t sum = 0;
    
    for (int i = 0; i < searches; i++) {
        uint64_t key = (uint64_t)rand() * rand();
        sum += geo_lower_bound(index->records, index->count, key);
    }
    
    double elapsed = geo_get_time_ms() - start;
    printf("    %d binary searches in %.2f ms (%.2f M ops/sec)\n",
           searches, elapsed, searches / elapsed / 1000.0);
    
    geo_index_destroy(index);
    return 1;
}

// =========================================================
// Stress Tests
// =========================================================

int test_stress_large_dataset(void) {
    int n = 5000000;
    printf("    Creating index with %d records...\n", n);
    
    GeoIndex *index = geo_index_create(n);
    ASSERT_TRUE(index != NULL, "failed to create large index");
    
    srand(42);
    double start = geo_get_time_ms();
    
    for (int i = 0; i < n; i++) {
        double lat = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        double lng = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
        geo_index_add(index, i, lat, lng);
    }
    
    double insert_time = geo_get_time_ms() - start;
    printf("    Insert time: %.2f ms (%.2f M ops/sec)\n",
           insert_time, n / insert_time / 1000.0);
    
    start = geo_get_time_ms();
    geo_index_build(index);
    double build_time = geo_get_time_ms() - start;
    printf("    Build time: %.2f ms\n", build_time);
    
    int num_searches = 1000;
    double total_search_time = 0;
    
    srand(123);
    for (int i = 0; i < num_searches; i++) {
        double lat = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        double lng = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
        
        GeoSearchStats stats;
        GeoSearchResult *result = geo_search_radius(index, lat, lng, 50.0, &stats);
        total_search_time += stats.search_time_ms;
        geo_result_destroy(result);
    }
    
    printf("    Avg search time (50km radius): %.3f ms\n", total_search_time / num_searches);
    
    geo_index_destroy(index);
    return 1;
}

int test_stress_many_searches(void) {
    int n = 100000;
    GeoIndex *index = geo_index_create(n);
    
    srand(42);
    for (int i = 0; i < n; i++) {
        double lat = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        double lng = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
        geo_index_add(index, i, lat, lng);
    }
    geo_index_build(index);
    
    int num_searches = 10000;
    double start = geo_get_time_ms();
    
    srand(123);
    for (int i = 0; i < num_searches; i++) {
        double lat = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        double lng = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
        
        GeoSearchResult *result = geo_search_radius(index, lat, lng, 100.0, NULL);
        geo_result_destroy(result);
    }
    
    double elapsed = geo_get_time_ms() - start;
    printf("    %d searches in %.2f ms (%.2f searches/sec)\n",
           num_searches, elapsed, num_searches / elapsed * 1000.0);
    
    geo_index_destroy(index);
    return 1;
}

int test_stress_dense_area(void) {
    int n = 100000;
    GeoIndex *index = geo_index_create(n);
    
    double center_lat = -23.5505;
    double center_lng = -46.6333;
    
    srand(42);
    for (int i = 0; i < n; i++) {
        double lat = center_lat + ((double)rand() / RAND_MAX - 0.5) * 0.1;
        double lng = center_lng + ((double)rand() / RAND_MAX - 0.5) * 0.1;
        geo_index_add(index, i, lat, lng);
    }
    geo_index_build(index);
    
    GeoSearchStats stats;
    GeoSearchResult *result = geo_search_radius(index, center_lat, center_lng, 5.0, &stats);
    
    printf("    Dense area: found %zu of %d in %.3f ms\n",
           result->count, n, stats.search_time_ms);
    
    geo_result_destroy(result);
    geo_index_destroy(index);
    return 1;
}

// =========================================================
// Edge Case Tests
// =========================================================

int test_edge_antimeridian(void) {
    GeoIndex *index = geo_index_create(10);
    
    geo_index_add(index, 1, 0.0, 179.9);
    geo_index_add(index, 2, 0.0, -179.9);
    geo_index_add(index, 3, 0.0, 180.0);
    geo_index_add(index, 4, 0.0, -180.0);
    geo_index_build(index);
    
    GeoSearchResult *result = geo_search_radius(index, 0.0, 180.0, 100.0, NULL);
    printf("    Points near antimeridian: %zu\n", result->count);
    
    geo_result_destroy(result);
    geo_index_destroy(index);
    return 1;
}

int test_edge_poles(void) {
    GeoIndex *index = geo_index_create(10);
    
    geo_index_add(index, 1, 89.9, 0.0);
    geo_index_add(index, 2, 89.9, 90.0);
    geo_index_add(index, 3, 89.9, 180.0);
    geo_index_add(index, 4, 89.9, -90.0);
    geo_index_build(index);
    
    GeoSearchResult *result = geo_search_radius(index, 90.0, 0.0, 100.0, NULL);
    printf("    Points near north pole: %zu\n", result->count);
    
    geo_result_destroy(result);
    geo_index_destroy(index);
    return 1;
}

int test_edge_empty_index(void) {
    GeoIndex *index = geo_index_create(10);
    geo_index_build(index);
    
    GeoSearchResult *result = geo_search_radius(index, 0.0, 0.0, 100.0, NULL);
    ASSERT_TRUE(result != NULL, "should return empty result, not NULL");
    ASSERT_EQ(result->count, 0, "empty index should return 0 results");
    
    geo_result_destroy(result);
    geo_index_destroy(index);
    return 1;
}

int test_edge_single_point(void) {
    GeoIndex *index = geo_index_create(1);
    geo_index_add(index, 1, 0.0, 0.0);
    geo_index_build(index);
    
    GeoSearchResult *result = geo_search_radius(index, 0.0, 0.0, 1.0, NULL);
    ASSERT_EQ(result->count, 1, "should find the single point");
    
    geo_result_destroy(result);
    geo_index_destroy(index);
    return 1;
}

// =========================================================
// Validation Tests
// =========================================================

int test_validation_coordinates(void) {
    ASSERT_TRUE(geo_is_valid_lat(0.0), "0 is valid lat");
    ASSERT_TRUE(geo_is_valid_lat(-90.0), "-90 is valid lat");
    ASSERT_TRUE(geo_is_valid_lat(90.0), "90 is valid lat");
    ASSERT_FALSE(geo_is_valid_lat(-91.0), "-91 is invalid lat");
    ASSERT_FALSE(geo_is_valid_lat(91.0), "91 is invalid lat");
    
    ASSERT_TRUE(geo_is_valid_lng(0.0), "0 is valid lng");
    ASSERT_TRUE(geo_is_valid_lng(-180.0), "-180 is valid lng");
    ASSERT_TRUE(geo_is_valid_lng(180.0), "180 is valid lng");
    ASSERT_FALSE(geo_is_valid_lng(-181.0), "-181 is invalid lng");
    ASSERT_FALSE(geo_is_valid_lng(181.0), "181 is invalid lng");
    
    return 1;
}

int test_validation_clamp(void) {
    ASSERT_NEAR(geo_clamp_lat(-100.0), -90.0, 0.0001, "clamp lat min");
    ASSERT_NEAR(geo_clamp_lat(100.0), 90.0, 0.0001, "clamp lat max");
    ASSERT_NEAR(geo_clamp_lat(45.0), 45.0, 0.0001, "clamp lat unchanged");
    
    ASSERT_NEAR(geo_clamp_lng(-200.0), -180.0, 0.0001, "clamp lng min");
    ASSERT_NEAR(geo_clamp_lng(200.0), 180.0, 0.0001, "clamp lng max");
    return 1;
}

int test_validation_wrap_lng(void) {
    ASSERT_NEAR(geo_wrap_lng(0.0), 0.0, 0.0001, "wrap 0");
    ASSERT_NEAR(geo_wrap_lng(180.0), 180.0, 0.0001, "wrap 180");
    ASSERT_NEAR(geo_wrap_lng(-180.0), -180.0, 0.0001, "wrap -180");
    ASSERT_NEAR(geo_wrap_lng(270.0), -90.0, 0.0001, "wrap 270");
    ASSERT_NEAR(geo_wrap_lng(-270.0), 90.0, 0.0001, "wrap -270");
    ASSERT_NEAR(geo_wrap_lng(540.0), 180.0, 0.0001, "wrap 540");
    return 1;
}

// =========================================================
// Main
// =========================================================

int main(void) {
    printf("========================================\n");
    printf("GEO INDEX TEST SUITE\n");
    printf("========================================\n");
    
    // Morton Code Tests
    printf("\n--- MORTON CODE TESTS ---\n");
    RUN_TEST(test_spread_compact_bits_roundtrip);
    RUN_TEST(test_spread_bits_pattern);
    
    // Encode/Decode Tests
    printf("\n--- ENCODE/DECODE TESTS ---\n");
    RUN_TEST(test_encode_decode_roundtrip);
    RUN_TEST(test_encode_decode_edge_cases);
    RUN_TEST(test_encode_clamping);
    RUN_TEST(test_encode_ordering);
    
    // Distance Tests
    printf("\n--- DISTANCE TESTS ---\n");
    RUN_TEST(test_haversine_known_distances);
    RUN_TEST(test_haversine_zero_distance);
    RUN_TEST(test_haversine_symmetry);
    RUN_TEST(test_haversine_triangle_inequality);
    RUN_TEST(test_fast_distance_accuracy);
    
    // Index Tests
    printf("\n--- INDEX TESTS ---\n");
    RUN_TEST(test_index_create_destroy);
    RUN_TEST(test_index_add_single);
    RUN_TEST(test_index_add_batch);
    RUN_TEST(test_index_auto_grow);
    RUN_TEST(test_index_build_sorts);
    
    // Binary Search Tests
    printf("\n--- BINARY SEARCH TESTS ---\n");
    RUN_TEST(test_lower_bound_basic);
    
    // Search Tests
    printf("\n--- SEARCH TESTS ---\n");
    RUN_TEST(test_search_radius_basic);
    RUN_TEST(test_search_radius_finds_nearby);
    RUN_TEST(test_search_radius_empty_result);
    RUN_TEST(test_search_bbox_basic);
    RUN_TEST(test_search_knn_basic);
    
    // Precision Tests
    printf("\n--- PRECISION TESTS ---\n");
    RUN_TEST(test_precision_at_equator);
    RUN_TEST(test_precision_at_poles);
    RUN_TEST(test_precision_random_points);
    
    // Edge Case Tests
    printf("\n--- EDGE CASE TESTS ---\n");
    RUN_TEST(test_edge_antimeridian);
    RUN_TEST(test_edge_poles);
    RUN_TEST(test_edge_empty_index);
    RUN_TEST(test_edge_single_point);
    
    // Validation Tests
    printf("\n--- VALIDATION TESTS ---\n");
    RUN_TEST(test_validation_coordinates);
    RUN_TEST(test_validation_clamp);
    RUN_TEST(test_validation_wrap_lng);
    
    // Performance Tests
    printf("\n--- PERFORMANCE TESTS ---\n");
    RUN_TEST(test_perf_encode_decode);
    RUN_TEST(test_perf_index_build);
    RUN_TEST(test_perf_search_radius);
    RUN_TEST(test_perf_binary_search);
    
    // Stress Tests
    printf("\n--- STRESS TESTS ---\n");
    RUN_TEST(test_stress_large_dataset);
    RUN_TEST(test_stress_many_searches);
    RUN_TEST(test_stress_dense_area);
    
    // Summary
    printf("\n========================================\n");
    printf("TEST SUMMARY\n");
    printf("========================================\n");
    printf("Total:  %d\n", g_tests_run);
    printf("Passed: %d\n", g_tests_passed);
    printf("Failed: %d\n", g_tests_failed);
    printf("========================================\n");
    
    return g_tests_failed > 0 ? 1 : 0;
}
