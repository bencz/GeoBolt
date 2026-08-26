#include "geo_index.h"
#include "geo_index_persistence.h"
#include "geo_index_private.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <inttypes.h>
#include <assert.h>
#include <stdatomic.h>

#if defined(__unix__) || defined(__APPLE__)
#include <pthread.h>
#include <sched.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

// =========================================================
// Benchmark Timing Utilities
// =========================================================

typedef struct {
    double start_time;
    double end_time;
    double elapsed_ms;
    const char *name;
} BenchmarkTimer;

static inline void benchmark_start(BenchmarkTimer *timer, const char *name) {
    timer->name = name;
    timer->start_time = geo_get_time_ms();
}

static inline void benchmark_end(BenchmarkTimer *timer) {
    timer->end_time = geo_get_time_ms();
    timer->elapsed_ms = timer->end_time - timer->start_time;
}

__attribute__((unused))
static inline void benchmark_print(const BenchmarkTimer *timer) {
    printf("    [TIMER] %s: %.3f ms\n", timer->name, timer->elapsed_ms);
}

__attribute__((unused))
static inline void benchmark_print_ops(const BenchmarkTimer *timer, size_t ops) {
    double ops_per_sec = (ops / timer->elapsed_ms) * 1000.0;
    printf("    [TIMER] %s: %.3f ms (%.2f M ops/sec)\n", 
           timer->name, timer->elapsed_ms, ops_per_sec / 1000000.0);
}

// =========================================================
// Test Framework
// =========================================================

#define TEST_PASSED "\033[32mPASSED\033[0m"
#define TEST_FAILED "\033[31mFAILED\033[0m"

static int g_tests_run = 0;
static int g_tests_passed = 0;
static int g_tests_failed = 0;
static const char *g_test_filter = NULL;

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

static void run_test(const char *name, int (*test_function)(void))
{
    if (g_test_filter && !strstr(name, g_test_filter)) {
        return;
    }

    g_tests_run++;
    BenchmarkTimer timer;

    benchmark_start(&timer, name);
    printf("\n[TEST %d] %s\n", g_tests_run, name);

    int result = test_function();

    benchmark_end(&timer);

    if (result) {
        printf("  %s (%.3f ms)\n", TEST_PASSED, timer.elapsed_ms);
        g_tests_passed++;
    } else {
        printf("  %s (%.3f ms)\n", TEST_FAILED, timer.elapsed_ms);
        g_tests_failed++;
    }
}

#define RUN_TEST(test_function) run_test(#test_function, test_function)

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

int test_index_add_batch(void)
{
    GeoIndex *index = geo_index_create(10);

    ASSERT_TRUE(index != NULL, "index creation failed");

    const uint64_t ids[5] = {1, 2, 3, 4, 5};
    const double latitudes[5] = {-23.5505, -22.9068, 40.7128, 51.5074, 35.6762};
    const double longitudes[5] = {-46.6333, -43.1729, -74.0060, -0.1278, 139.6503};

    bool added = geo_index_add_batch(index, ids, latitudes, longitudes, 5);
    ASSERT_TRUE(added, "batch add should succeed");
    ASSERT_EQ(index->count, 5, "count should be 5");

    const uint64_t invalid_ids[9] = {6, 7, 8, 9, 10, 11, 12, 13, 14};
    const double invalid_latitudes[9] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, NAN};
    const double invalid_longitudes[9] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};

    added = geo_index_add_batch(index, invalid_ids, invalid_latitudes, invalid_longitudes, 9);
    ASSERT_FALSE(added, "batch with an invalid tail point should fail");
    ASSERT_EQ(index->count, 5, "a rejected batch must not partially modify the index");

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

int test_index_persistence_roundtrip(void)
{
#if defined(__unix__) || defined(__APPLE__)
    char path[128];

    snprintf(path, sizeof(path), "/tmp/geobolt-index-%ld.bin", (long)getpid());

    GeoIndex *index = geo_index_create(4);
    GeoIndex *mapped_index = NULL;
    GeoSearchResult *result = NULL;
    bool succeeded = index != NULL;

    if (succeeded) {
        succeeded = geo_index_add(index, 10, 0.0, 179.9) &&
                    geo_index_add(index, 11, 0.0, -179.9) &&
                    geo_index_add(index, 12, -23.5505, -46.6333);
    }

    if (succeeded) {
        geo_index_build(index);
        succeeded = geo_index_save(index, path);
    }

    if (succeeded) {
        mapped_index = geo_index_open_mmap(path);
        succeeded = mapped_index != NULL &&
                    geo_index_is_read_only(mapped_index) &&
                    mapped_index->count == index->count;
    }

    if (succeeded) {
        result = geo_search_radius(mapped_index, 0.0, 180.0, 100.0, NULL);
        succeeded = result != NULL && result->count == 2;
    }

    geo_result_destroy(result);
    geo_index_destroy(mapped_index);
    geo_index_destroy(index);
    remove(path);

    return succeeded;
#else
    return 1;
#endif
}

int test_persisted_checksum_is_independent_of_write_blocks(void)
{
    uint64_t words[257];
    uint64_t state = UINT64_C(0x6a09e667f3bcc909);

    for (size_t i = 0; i < sizeof(words) / sizeof(words[0]); ++i) {
        state = state * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
        words[i] = state;
    }

    uint64_t contiguous = geo_persisted_checksum_update(geo_persisted_checksum_initial(), words, sizeof(words));
    uint64_t streamed = geo_persisted_checksum_initial();
    const size_t block_words[] = { 1, 7, 2, 19, 4, 31 };
    size_t offset = 0;
    size_t block_index = 0;

    while (offset < sizeof(words) / sizeof(words[0])) {
        size_t count = block_words[block_index % (sizeof(block_words) / sizeof(block_words[0]))];
        size_t remaining = sizeof(words) / sizeof(words[0]) - offset;

        if (count > remaining) {
            count = remaining;
        }

        streamed = geo_persisted_checksum_update(streamed, words + offset, count * sizeof(*words));
        offset += count;
        block_index++;
    }

    ASSERT_EQ(streamed, contiguous, "streamed checksum must not depend on aligned write block sizes");

    uint64_t combined = geo_persisted_checksum_initial();
    offset = 0;
    block_index = 0;

    while (offset < sizeof(words) / sizeof(words[0])) {
        size_t count = block_words[block_index % (sizeof(block_words) / sizeof(block_words[0]))];
        size_t remaining = sizeof(words) / sizeof(words[0]) - offset;

        if (count > remaining) {
            count = remaining;
        }

        size_t bytes = count * sizeof(*words);
        uint64_t suffix = geo_persisted_checksum_update(0, words + offset, bytes);

        combined = geo_persisted_checksum_combine_aligned(combined, suffix, bytes);
        offset += count;
        block_index++;
    }

    ASSERT_EQ(combined, contiguous, "independent aligned partition checksums must combine to the contiguous checksum");

    return 1;
}

int test_index_persistence_rejects_corruption(void)
{
#if defined(__unix__) || defined(__APPLE__)
    char path[128];

    snprintf(path, sizeof(path), "/tmp/geobolt-corruption-%ld.bin", (long) getpid());
    remove(path);

    GeoIndex *index = geo_index_create(128);
    bool succeeded = index != NULL;

    for (size_t i = 0; succeeded && i < 128; ++i) {
        succeeded = geo_index_add(index,
                                  UINT64_C(9000000) + i,
                                  (double) i * 0.5 - 32.0,
                                  (double) i * 1.25 - 80.0);
    }

    if (succeeded) {
        succeeded = geo_index_build(index) && geo_index_save(index, path);
    }

    FILE *file = succeeded ? fopen(path, "r+b") : NULL;
    unsigned char byte = 0;

    if (succeeded) {
        succeeded = file != NULL &&
                    fseeko(file, (off_t) sizeof(GeoFileHeader) + 3, SEEK_SET) == 0 &&
                    fread(&byte, 1, 1, file) == 1 &&
                    fseeko(file, -1, SEEK_CUR) == 0;
    }

    if (succeeded) {
        byte ^= UINT8_C(0x80);
        succeeded = fwrite(&byte, 1, 1, file) == 1 && fflush(file) == 0;
    }

    if (file && fclose(file) != 0) {
        succeeded = false;
    }

    GeoIndex *corrupted = succeeded ? geo_index_open_mmap(path) : NULL;

    if (succeeded) {
        succeeded = corrupted == NULL && geo_index_save(index, path);
    }

    geo_index_destroy(corrupted);

    struct stat status;

    if (succeeded) {
        succeeded = stat(path, &status) == 0 &&
                    status.st_size > 1 &&
                    truncate(path, status.st_size - 1) == 0;
    }

    GeoIndex *truncated = succeeded ? geo_index_open_mmap(path) : NULL;

    if (succeeded) {
        succeeded = truncated == NULL;
    }

    geo_index_destroy(truncated);
    geo_index_destroy(index);
    remove(path);

    return succeeded;
#else
    return 1;
#endif
}

int test_density_metadata_mmap_roundtrip(void)
{
#if defined(__unix__) || defined(__APPLE__)
    const size_t record_count = 8192;
    const size_t clustered_count = 6144;
    const double center_latitude = -23.5505;
    const double center_longitude = -46.6333;
    const double radius_km = 0.5;
    char path[128];

    snprintf(path, sizeof(path), "/tmp/geobolt-density-%ld.bin", (long) getpid());

    GeoIndex *index = geo_index_create(record_count);
    GeoIndex *mapped_index = NULL;
    bool succeeded = index != NULL;

    for (size_t i = 0; succeeded && i < clustered_count; ++i) {
        double latitude_offset = (double) ((i * 104729U) % 10000U) / 100000.0 - 0.05;
        double longitude_offset = (double) ((i * 130363U) % 10000U) / 100000.0 - 0.05;

        succeeded = geo_index_add(index,
                                  i,
                                  center_latitude + latitude_offset,
                                  center_longitude + longitude_offset);
    }

    for (size_t i = clustered_count; succeeded && i < record_count; ++i) {
        double latitude = (double) ((i * 104729U) % 180000U) / 1000.0 - 90.0;
        double longitude = (double) ((i * 130363U) % 360000U) / 1000.0 - 180.0;

        succeeded = geo_index_add(index, i, latitude, longitude);
    }

    if (succeeded) {
        geo_index_build(index);
        succeeded = index->maximum_density_refinement > 0;
    }

    size_t expected_count = 0;

    for (size_t i = 0; succeeded && i < index->count; ++i) {
        GeoPoint point = geo_decode(index->records[i].z);

        if (geo_haversine_km(center_latitude, center_longitude, point.lat, point.lng) <= radius_km) {
            expected_count++;
        }
    }

    size_t actual_count = 0;

    if (succeeded) {
        succeeded = geo_search_radius_count(index,
                                            center_latitude,
                                            center_longitude,
                                            radius_km,
                                            &actual_count,
                                            NULL) &&
                    actual_count == expected_count &&
                    geo_index_save(index, path);
    }

    if (succeeded) {
        mapped_index = geo_index_open_mmap(path);
        succeeded = mapped_index != NULL &&
                    mapped_index->maximum_density_refinement == index->maximum_density_refinement;
    }

    if (succeeded) {
        succeeded = geo_search_radius_count(mapped_index,
                                            center_latitude,
                                            center_longitude,
                                            radius_km,
                                            &actual_count,
                                            NULL) &&
                    actual_count == expected_count;
    }

    geo_index_destroy(mapped_index);
    geo_index_destroy(index);
    remove(path);

    return succeeded;
#else
    return 1;
#endif
}

static int compare_uint64_values(const void *first_value, const void *second_value)
{
    uint64_t first = *(const uint64_t *) first_value;
    uint64_t second = *(const uint64_t *) second_value;

    return (first > second) - (first < second);
}

int test_batch_executor_ids_match_exact_queries(void)
{
#if defined(__unix__) || defined(__APPLE__)
    const size_t record_count = 50000;
    const size_t query_count = 257;
    GeoIndex *index = geo_index_create(record_count);
    GeoIndex *shards[2] = {
        geo_index_create(record_count / 2),
        geo_index_create(record_count / 2),
    };
    double *latitudes = malloc(query_count * sizeof(*latitudes));
    double *longitudes = malloc(query_count * sizeof(*longitudes));
    double *radii = malloc(query_count * sizeof(*radii));
    GeoBatchIdResult *batch_result = geo_batch_id_result_create(query_count, 1024);
    GeoIdResult *exact_result = geo_id_result_create(128);
    GeoQueryExecutorConfig config = {
        .thread_count = 4,
        .scheduling_chunk = 7,
        .pin_workers = false,
    };
    GeoQueryExecutor *executor = NULL;
    bool succeeded = index && shards[0] && shards[1] && latitudes && longitudes && radii && batch_result && exact_result;

    for (size_t record = 0; succeeded && record < record_count; ++record) {
        double latitude = (double) ((record * 104729U) % 180000U) / 1000.0 - 90.0;
        double longitude = (double) ((record * 130363U) % 360000U) / 1000.0 - 180.0;

        succeeded = geo_index_add(index, record + 1000U, latitude, longitude) &&
                    geo_index_add(shards[record & 1U], record + 1000U, latitude, longitude);
    }

    if (succeeded) {
        geo_index_build(index);
        geo_index_build(shards[0]);
        geo_index_build(shards[1]);

        for (size_t query = 0; query < query_count; ++query) {
            latitudes[query] = (double) ((query * 8191U) % 180000U) / 1000.0 - 90.0;
            longitudes[query] = (double) ((query * 131071U) % 360000U) / 1000.0 - 180.0;
            radii[query] = 1.0 + (double) (query % 17U) * 20.0;
        }

        // Exercise the cost planner's copy strategy as part of the same end-to-end batch.
        latitudes[query_count - 1] = 0.0;
        longitudes[query_count - 1] = 0.0;
        radii[query_count - 1] = M_PI * GEO_EARTH_RADIUS_KM;

        executor = geo_query_executor_create(index, &config);
        succeeded = executor &&
                    geo_query_executor_search_radius_ids(executor,
                                                         latitudes,
                                                         longitudes,
                                                         radii,
                                                         query_count,
                                                         batch_result);
    }

    if (succeeded) {
        succeeded = batch_result->query_count == query_count &&
                    batch_result->offsets[0] == 0 &&
                    batch_result->offsets[query_count] == batch_result->id_count;
    }

    for (size_t query = 0; succeeded && query < query_count; ++query) {
        size_t begin = batch_result->offsets[query];
        size_t end = batch_result->offsets[query + 1];

        succeeded = begin <= end &&
                    geo_search_radius_ids_reuse(index,
                                                latitudes[query],
                                                longitudes[query],
                                                radii[query],
                                                exact_result,
                                                NULL) &&
                    exact_result->count == end - begin &&
                    (!exact_result->count ||
                     memcmp(exact_result->ids,
                            batch_result->ids + begin,
                            exact_result->count * sizeof(*exact_result->ids)) == 0);
    }

    geo_query_executor_destroy(executor);
    executor = NULL;

    if (succeeded) {
        executor = geo_query_executor_create_sharded((const GeoIndex *const *) shards, 2, &config);
        succeeded = executor &&
                    geo_query_executor_search_radius_ids(executor,
                                                         latitudes,
                                                         longitudes,
                                                         radii,
                                                         query_count,
                                                         batch_result);
    }

    for (size_t query = 0; succeeded && query < query_count; ++query) {
        size_t begin = batch_result->offsets[query];
        size_t end = batch_result->offsets[query + 1];
        size_t batch_count = end - begin;

        succeeded = geo_search_radius_ids_reuse(index,
                                                latitudes[query],
                                                longitudes[query],
                                                radii[query],
                                                exact_result,
                                                NULL) &&
                    exact_result->count == batch_count;

        if (succeeded && batch_count > 1) {
            qsort(exact_result->ids, exact_result->count, sizeof(*exact_result->ids), compare_uint64_values);
            qsort(batch_result->ids + begin, batch_count, sizeof(*batch_result->ids), compare_uint64_values);
        }

        succeeded = succeeded &&
                    (!batch_count ||
                     memcmp(exact_result->ids,
                            batch_result->ids + begin,
                            batch_count * sizeof(*exact_result->ids)) == 0);
    }

    geo_query_executor_destroy(executor);
    geo_id_result_destroy(exact_result);
    geo_batch_id_result_destroy(batch_result);
    free(radii);
    free(longitudes);
    free(latitudes);
    geo_index_destroy(shards[1]);
    geo_index_destroy(shards[0]);
    geo_index_destroy(index);

    return succeeded;
#else
    return 1;
#endif
}

int test_stream_builder_multi_run_roundtrip(void)
{
#if defined(__unix__) || defined(__APPLE__)
    const size_t point_count = 4096;
    const size_t clustered_count = 3072;
    const double center_latitude = -23.5505;
    const double center_longitude = -46.6333;
    const double radius_km = 0.5;
    char path[128];

    snprintf(path, sizeof(path), "/tmp/geobolt-stream-%ld.bin", (long) getpid());

    uint64_t *ids = malloc(point_count * sizeof(*ids));
    double *latitudes = malloc(point_count * sizeof(*latitudes));
    double *longitudes = malloc(point_count * sizeof(*longitudes));
    GeoRecord *encoded_records = malloc(point_count * sizeof(*encoded_records));
    const size_t chunk_capacity = 4;
    const size_t expected_runs = (point_count + chunk_capacity - 1) / chunk_capacity;
    GeoStreamBuilder *builder = geo_stream_builder_create(path, "/tmp", chunk_capacity);
    GeoIndex *mapped = NULL;
    GeoSearchResult *result = NULL;
    bool succeeded = ids && latitudes && longitudes && encoded_records && builder;

    for (size_t i = 0; succeeded && i < point_count; ++i) {
        ids[i] = UINT64_C(1000000) + i;

        if (i < clustered_count) {
            latitudes[i] = center_latitude + (double) ((i * 104729U) % 1000U) / 1000000.0 - 0.0005;
            longitudes[i] = center_longitude + (double) ((i * 130363U) % 1000U) / 1000000.0 - 0.0005;
        } else {
            latitudes[i] = (double) ((i * 104729U) % 180000U) / 1000.0 - 90.0;
            longitudes[i] = (double) ((i * 130363U) % 360000U) / 1000.0 - 180.0;
        }

        encoded_records[i] = (GeoRecord) {
            .id = ids[i],
            .z = geo_encode(latitudes[i], longitudes[i]),
        };
    }

    if (succeeded) {
        succeeded = geo_stream_builder_add_batch(builder, ids, latitudes, longitudes, 3333) &&
                    geo_stream_builder_add_records(builder, encoded_records + 3333, point_count - 3333);
    }

    GeoStreamBuildStats build_stats;

    if (succeeded) {
        succeeded = geo_stream_builder_finish(builder, &build_stats) &&
                    build_stats.records_written == point_count &&
                    build_stats.runs_created == expected_runs &&
                    build_stats.intermediate_merges == 33 &&
                    build_stats.peak_open_runs == 64;
    }

    if (succeeded) {
        mapped = geo_index_open_mmap(path);
        succeeded = mapped &&
                    mapped->count == point_count &&
                    mapped->prefix_bits > 0 &&
                    mapped->prefix_offsets != NULL &&
                    mapped->maximum_density_refinement > 0;
    }

    if (succeeded) {
        for (size_t i = 1; i < mapped->count; ++i) {
            if (mapped->records[i - 1].z > mapped->records[i].z) {
                succeeded = false;
                break;
            }
        }
    }

    size_t expected_count = 0;

    for (size_t i = 0; succeeded && i < mapped->count; ++i) {
        GeoPoint point = geo_decode(mapped->records[i].z);

        if (geo_haversine_km(center_latitude, center_longitude, point.lat, point.lng) <= radius_km) {
            expected_count++;
        }
    }

    if (succeeded) {
        result = geo_search_radius(mapped, center_latitude, center_longitude, radius_km, NULL);
        succeeded = result && result->count == expected_count;
    }

    geo_result_destroy(result);
    geo_index_destroy(mapped);
    geo_stream_builder_destroy(builder);
    free(encoded_records);
    free(longitudes);
    free(latitudes);
    free(ids);
    remove(path);

    return succeeded;
#else
    return 1;
#endif
}

int test_stream_builder_single_run_direct_copy_roundtrip(void)
{
#if defined(__unix__) || defined(__APPLE__)
    const size_t record_count = 32768;
    char path[160];
    long process_id = (long) getpid();

    snprintf(path, sizeof(path), "/tmp/geobolt-stream-single-run-%ld.bin", process_id);
    remove(path);

    GeoRecord *records = malloc(record_count * sizeof(*records));
    uint64_t state = UINT64_C(0xbb67ae8584caa73b);
    uint64_t expected_id_sum = 0;
    uint64_t expected_z_xor = 0;
    bool succeeded = records != NULL;

    for (size_t i = 0; succeeded && i < record_count; ++i) {
        state = state * UINT64_C(2862933555777941757) + UINT64_C(3037000493);
        records[i] = (GeoRecord) {
            .id = UINT64_C(610000000) + i,
            .z = state,
        };
        expected_id_sum += records[i].id;
        expected_z_xor ^= records[i].z;
    }

    GeoStreamBuilder *builder = succeeded ? geo_stream_builder_create_parallel(path, "/tmp", record_count, 2) : NULL;
    GeoStreamBuildStats stats;

    succeeded = builder &&
                geo_stream_builder_add_records(builder, records, record_count) &&
                geo_stream_builder_finish(builder, &stats) &&
                stats.records_written == record_count &&
                stats.runs_created == 0U &&
                stats.intermediate_merges == 0U;

    GeoIndex *mapped = succeeded ? geo_index_open_mmap(path) : NULL;
    uint64_t actual_id_sum = 0;
    uint64_t actual_z_xor = 0;

    succeeded = mapped && mapped->count == record_count;

    for (size_t i = 0; succeeded && i < mapped->count; ++i) {
        if (i && mapped->records[i - 1].z > mapped->records[i].z) {
            succeeded = false;
            break;
        }

        actual_id_sum += mapped->records[i].id;
        actual_z_xor ^= mapped->records[i].z;
    }

    size_t visible_count = 0;

    succeeded = succeeded &&
                actual_id_sum == expected_id_sum &&
                actual_z_xor == expected_z_xor &&
                geo_search_radius_count(mapped, 0.0, 0.0, 25000.0, &visible_count, NULL) &&
                visible_count == record_count;

    geo_index_destroy(mapped);
    geo_stream_builder_destroy(builder);
    free(records);
    remove(path);

    return succeeded;
#else
    return 1;
#endif
}

int test_parallel_radix_build_preserves_records(void)
{
    const size_t record_count = 600000;
    GeoRecord *records = malloc(record_count * sizeof(*records));
    GeoIndex *index = geo_index_create(record_count);
    uint64_t state = UINT64_C(0x9e3779b97f4a7c15);
    uint64_t expected_id_sum = 0;
    uint64_t expected_z_xor = 0;
    bool succeeded = records && index;

    for (size_t i = 0; succeeded && i < record_count; ++i) {
        state = state * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);

        records[i] = (GeoRecord) {
            .id = UINT64_C(10000000) + i,
            .z = i % 17 ? state : UINT64_C(0x123456789abcdef0),
        };
        expected_id_sum += records[i].id;
        expected_z_xor ^= records[i].z;
    }

    if (succeeded) {
        succeeded = geo_index_add_records(index, records, record_count) &&
                    geo_index_build_parallel(index, 4) &&
                    index->sorted &&
                    index->prefix_bits > 0;
    }

    uint64_t actual_id_sum = 0;
    uint64_t actual_z_xor = 0;
    uint64_t previous_duplicate_id = 0;
    bool found_duplicate = false;

    for (size_t i = 0; succeeded && i < index->count; ++i) {
        if (i && index->records[i - 1].z > index->records[i].z) {
            succeeded = false;
            break;
        }

        actual_id_sum += index->records[i].id;
        actual_z_xor ^= index->records[i].z;

        if (index->records[i].z == UINT64_C(0x123456789abcdef0)) {
            if (found_duplicate && index->records[i].id < previous_duplicate_id) {
                succeeded = false;
                break;
            }

            previous_duplicate_id = index->records[i].id;
            found_duplicate = true;
        }
    }

    if (succeeded) {
        succeeded = found_duplicate && actual_id_sum == expected_id_sum && actual_z_xor == expected_z_xor;
    }

    geo_index_destroy(index);
    free(records);

    return succeeded;
}

int test_parallel_radix_sorter_reuses_workers_and_scratch(void)
{
    const size_t record_count = 300000;
    const size_t rounds = 3;
    GeoParallelSorter *sorter = geo_parallel_sorter_create(record_count, 4);
    GeoIndex *index = geo_index_create(record_count);
    GeoRecord *records = malloc(record_count * sizeof(*records));
    bool succeeded = sorter && index && records;

    for (size_t round = 0; succeeded && round < rounds; ++round) {
        uint64_t state = UINT64_C(0x6a09e667f3bcc909) ^ round;
        uint64_t expected_id_sum = 0;
        uint64_t expected_z_xor = 0;

        for (size_t i = 0; i < record_count; ++i) {
            state = state * UINT64_C(2862933555777941757) + UINT64_C(3037000493);

            records[i] = (GeoRecord) {
                .id = round * record_count + i,
                .z = i % 29 ? state : UINT64_C(0x4f1bbcdc676f2b5d),
            };
            expected_id_sum += records[i].id;
            expected_z_xor ^= records[i].z;
        }

        succeeded = geo_index_add_records(index, records, record_count) &&
                    geo_index_sort_transient_with_sorter(index, sorter);

        uint64_t actual_id_sum = 0;
        uint64_t actual_z_xor = 0;

        for (size_t i = 0; succeeded && i < index->count; ++i) {
            if (i && index->records[i - 1].z > index->records[i].z) {
                succeeded = false;
                break;
            }

            actual_id_sum += index->records[i].id;
            actual_z_xor ^= index->records[i].z;
        }

        succeeded = succeeded && actual_id_sum == expected_id_sum && actual_z_xor == expected_z_xor;
        geo_index_clear(index);
    }

    free(records);
    geo_index_destroy(index);
    geo_parallel_sorter_destroy(sorter);

    return succeeded;
}

#if defined(__unix__) || defined(__APPLE__)
typedef struct {
    GeoSegmentSet *set;
    size_t expected_count;
    atomic_bool ready;
    atomic_bool stop;
    atomic_bool failed;
} SegmentQueryWorker;

static void *segment_query_during_compaction(void *argument)
{
    SegmentQueryWorker *worker = argument;
    size_t count;

    if (!geo_segment_set_search_radius_count(worker->set, -23.5505, -46.6333, 5.0, &count, NULL) ||
        count != worker->expected_count) {
        atomic_store_explicit(&worker->failed, true, memory_order_release);
    }

    atomic_store_explicit(&worker->ready, true, memory_order_release);

    while (!atomic_load_explicit(&worker->stop, memory_order_acquire)) {
        if (!geo_segment_set_search_radius_count(worker->set, -23.5505, -46.6333, 5.0, &count, NULL) ||
            count != worker->expected_count) {
            atomic_store_explicit(&worker->failed, true, memory_order_release);
            break;
        }
    }

    return NULL;
}
#endif

int test_segment_set_manifest_queries_and_compaction(void)
{
#if defined(__unix__) || defined(__APPLE__)
    const size_t segment_count = 3;
    const size_t records_per_segment = 2000;
    const size_t nearest_count = 16;
    char manifest_path[128];
    char compacted_path[128];
    char segment_paths[segment_count][128];

    snprintf(manifest_path, sizeof(manifest_path), "/tmp/geobolt-manifest-%ld.bin", (long) getpid());
    snprintf(compacted_path, sizeof(compacted_path), "/tmp/geobolt-compacted-%ld.bin", (long) getpid());
    remove(manifest_path);
    remove(compacted_path);

    for (size_t segment = 0; segment < segment_count; ++segment) {
        snprintf(segment_paths[segment],
                 sizeof(segment_paths[segment]),
                 "/tmp/geobolt-segment-%ld-%zu.bin",
                 (long) getpid(),
                 segment);
        remove(segment_paths[segment]);
    }

    GeoSegmentSet *set = geo_segment_set_create(manifest_path);
    bool succeeded = set != NULL;

    for (size_t segment = 0; succeeded && segment < segment_count; ++segment) {
        GeoIndex *index = geo_index_create(records_per_segment);

        succeeded = index != NULL;

        for (size_t i = 0; succeeded && i < records_per_segment; ++i) {
            size_t global = segment * records_per_segment + i;
            double latitude = -23.5505 + ((double) ((global * 37U) % 6007U) - 3003.0) * 0.000015;
            double longitude = -46.6333 + ((double) ((global * 53U) % 6011U) - 3005.0) * 0.000015;

            succeeded = geo_index_add(index, UINT64_C(5000000) + global, latitude, longitude);
        }

        if (succeeded) {
            geo_index_build(index);
            succeeded = geo_index_save(index, segment_paths[segment]) &&
                        geo_segment_set_add_file(set, segment_paths[segment]);
        }

        geo_index_destroy(index);
    }

    size_t radius_count_before = 0;
    size_t bbox_count_before = 0;
    GeoSearchResult *radius_result = NULL;
    GeoSearchResult *nearest_before = NULL;

    if (succeeded) {
        radius_result = geo_segment_set_search_radius(set, -23.5505, -46.6333, 5.0, NULL);
        nearest_before = geo_segment_set_search_knn(set, -23.5505, -46.6333, nearest_count, 100.0, NULL);
        succeeded = geo_segment_set_count(set) == segment_count &&
                    geo_segment_set_record_count(set) == segment_count * records_per_segment &&
                    geo_segment_set_search_radius_count(set, -23.5505, -46.6333, 5.0, &radius_count_before, NULL) &&
                    geo_segment_set_search_bbox_count(set, -23.60, -23.50, -46.69, -46.58, &bbox_count_before, NULL) &&
                    radius_result &&
                    radius_result->count == radius_count_before &&
                    radius_count_before > 0 &&
                    bbox_count_before > 0 &&
                    nearest_before &&
                    nearest_before->count == nearest_count;
    }

    geo_result_destroy(radius_result);

    if (succeeded) {
        geo_segment_set_destroy(set);
        set = geo_segment_set_open(manifest_path);
        succeeded = set && geo_segment_set_count(set) == segment_count;
    }

    GeoSegmentCompactionStats compaction_stats;
    SegmentQueryWorker query_worker = {
        .set = set,
        .expected_count = radius_count_before,
    };
    pthread_t query_thread;
    bool query_thread_started = false;

    if (succeeded) {
        atomic_init(&query_worker.ready, false);
        atomic_init(&query_worker.stop, false);
        atomic_init(&query_worker.failed, false);
        query_thread_started = pthread_create(&query_thread, NULL, segment_query_during_compaction, &query_worker) == 0;
        succeeded = query_thread_started;
    }

    while (succeeded && !atomic_load_explicit(&query_worker.ready, memory_order_acquire)) {
        sched_yield();
    }

    if (succeeded) {
        succeeded = geo_segment_set_compact_with_workers(set, compacted_path, 4, &compaction_stats) &&
                    compaction_stats.records_written == segment_count * records_per_segment &&
                    compaction_stats.input_segments == segment_count &&
                    compaction_stats.worker_count >= 1U &&
                    compaction_stats.worker_count <= 4U &&
                    compaction_stats.partition_count >= compaction_stats.worker_count &&
                    geo_segment_set_count(set) == 1;
    }

    if (query_thread_started) {
        atomic_store_explicit(&query_worker.stop, true, memory_order_release);
        succeeded = pthread_join(query_thread, NULL) == 0 &&
                    !atomic_load_explicit(&query_worker.failed, memory_order_acquire) &&
                    succeeded;
    }

    for (size_t segment = 0; succeeded && segment < segment_count; ++segment) {
        succeeded = access(segment_paths[segment], F_OK) == 0;
    }

    size_t radius_count_after = 0;
    size_t bbox_count_after = 0;
    GeoSearchResult *nearest_after = NULL;

    if (succeeded) {
        nearest_after = geo_segment_set_search_knn(set, -23.5505, -46.6333, nearest_count, 100.0, NULL);
        succeeded = geo_segment_set_search_radius_count(set, -23.5505, -46.6333, 5.0, &radius_count_after, NULL) &&
                    geo_segment_set_search_bbox_count(set, -23.60, -23.50, -46.69, -46.58, &bbox_count_after, NULL) &&
                    radius_count_after == radius_count_before &&
                    bbox_count_after == bbox_count_before &&
                    nearest_after &&
                    nearest_after->count == nearest_count;
    }

    for (size_t i = 0; succeeded && i < nearest_count; ++i) {
        succeeded = nearest_before->results[i].id == nearest_after->results[i].id;
    }

    geo_result_destroy(nearest_after);
    geo_result_destroy(nearest_before);

    if (succeeded) {
        geo_segment_set_destroy(set);
        set = geo_segment_set_open(manifest_path);
        succeeded = set &&
                    geo_segment_set_count(set) == 1 &&
                    geo_segment_set_record_count(set) == segment_count * records_per_segment;
    }

    geo_segment_set_destroy(set);
    set = NULL;

    if (succeeded) {
        FILE *manifest = fopen(manifest_path, "r+b");
        int last_byte = EOF;

        if (manifest && fseek(manifest, -1, SEEK_END) == 0) {
            last_byte = fgetc(manifest);
        }

        if (last_byte != EOF && fseek(manifest, -1, SEEK_END) == 0) {
            bool wrote_byte = fputc(last_byte ^ 1, manifest) != EOF;
            bool closed_manifest = fclose(manifest) == 0;

            succeeded = wrote_byte && closed_manifest;
            manifest = NULL;
        } else {
            succeeded = false;
        }

        if (manifest) {
            fclose(manifest);
        }

        if (succeeded) {
            set = geo_segment_set_open(manifest_path);
            succeeded = set == NULL;
        }
    }

    geo_segment_set_destroy(set);
    remove(compacted_path);
    remove(manifest_path);

    for (size_t segment = 0; segment < segment_count; ++segment) {
        remove(segment_paths[segment]);
    }

    return succeeded;
#else
    return 1;
#endif
}

int test_segment_set_background_compaction(void)
{
#if defined(__unix__) || defined(__APPLE__)
    const size_t segment_count = 4;
    const size_t records_per_segment[] = { 128, 128, 128, 4096 };
    const size_t compacted_record_count = 384;
    const size_t total_record_count = 4480;
    char manifest_path[128];
    char segment_paths[segment_count][128];
    char compacted_path[512] = { 0 };

    snprintf(manifest_path, sizeof(manifest_path), "/tmp/geobolt-background-manifest-%ld.bin", (long) getpid());
    remove(manifest_path);

    for (size_t segment = 0; segment < segment_count; ++segment) {
        snprintf(segment_paths[segment],
                 sizeof(segment_paths[segment]),
                 "/tmp/geobolt-background-segment-%ld-%zu.bin",
                 (long) getpid(),
                 segment);
        remove(segment_paths[segment]);
    }

    GeoSegmentSet *set = geo_segment_set_create(manifest_path);
    bool succeeded = set && geo_segment_set_enable_background_compaction(set, "/tmp", 3);
    size_t global_offset = 0;

    for (size_t segment = 0; succeeded && segment < segment_count; ++segment) {
        GeoIndex *index = geo_index_create(records_per_segment[segment]);

        succeeded = index != NULL;

        for (size_t i = 0; succeeded && i < records_per_segment[segment]; ++i) {
            size_t global = global_offset + i;
            double latitude = -33.8688 + (double) ((global * 29U) % 1009U) * 0.00001;
            double longitude = 151.2093 + (double) ((global * 43U) % 1013U) * 0.00001;

            succeeded = geo_index_add(index, UINT64_C(7000000) + global, latitude, longitude);
        }

        if (succeeded) {
            geo_index_build(index);
            succeeded = geo_index_save(index, segment_paths[segment]) &&
                        geo_segment_set_add_file(set, segment_paths[segment]);
        }

        geo_index_destroy(index);
        global_offset += records_per_segment[segment];
    }

    GeoSegmentCompactionStats stats;
    GeoSegmentCompactionStats cumulative_stats;
    uint64_t completed_runs = 0;
    uint64_t failed_runs = 0;

    if (succeeded) {
        succeeded = geo_segment_set_wait_for_background_compaction(set, &stats) &&
                    geo_segment_set_background_compaction_totals(set,
                                                                 &cumulative_stats,
                                                                 &completed_runs,
                                                                 &failed_runs) &&
                    !geo_segment_set_background_compaction_active(set) &&
                    stats.input_segments == 3 &&
                    stats.records_written == compacted_record_count &&
                    cumulative_stats.input_segments == stats.input_segments &&
                    cumulative_stats.records_written == stats.records_written &&
                    completed_runs == 1 &&
                    failed_runs == 0 &&
                    geo_segment_set_count(set) == 2 &&
                    geo_segment_set_record_count(set) == total_record_count &&
                    geo_segment_set_copy_path(set, 1, compacted_path, sizeof(compacted_path));
    }

    for (size_t segment = 0; succeeded && segment < segment_count; ++segment) {
        succeeded = access(segment_paths[segment], F_OK) == 0;
    }

    if (succeeded) {
        geo_segment_set_destroy(set);
        set = geo_segment_set_open(manifest_path);
        succeeded = set &&
                    geo_segment_set_count(set) == 2 &&
                    geo_segment_set_record_count(set) == total_record_count;
    }

    geo_segment_set_destroy(set);

    if (compacted_path[0]) {
        remove(compacted_path);
    }

    remove(manifest_path);

    for (size_t segment = 0; segment < segment_count; ++segment) {
        remove(segment_paths[segment]);
    }

    return succeeded;
#else
    return 1;
#endif
}

int test_segment_set_background_tombstone_reclamation(void)
{
#if defined(__unix__) || defined(__APPLE__)
    const size_t record_count = 8192;
    const size_t removed_count = 4096;
    const uint64_t first_id = UINT64_C(73000000);
    char manifest_path[160];
    char segment_path[160];
    char mutation_path[192];
    char mutation_checkpoint_path[224];
    char compacted_path[512] = { 0 };
    long process_id = (long) getpid();

    snprintf(manifest_path, sizeof(manifest_path), "/tmp/geobolt-tombstone-background-manifest-%ld.bin", process_id);
    snprintf(segment_path, sizeof(segment_path), "/tmp/geobolt-tombstone-background-segment-%ld.bin", process_id);
    snprintf(mutation_path, sizeof(mutation_path), "%s.mutations", manifest_path);
    snprintf(mutation_checkpoint_path, sizeof(mutation_checkpoint_path), "%s.mutations.checkpoint", manifest_path);
    remove(manifest_path);
    remove(segment_path);
    remove(mutation_path);
    remove(mutation_checkpoint_path);

    GeoIndex *index = geo_index_create(record_count);
    uint64_t *removed_ids = malloc(removed_count * sizeof(*removed_ids));
    bool succeeded = index && removed_ids;

    for (size_t i = 0; succeeded && i < record_count; ++i) {
        double latitude = (double) (i % 1800U) * 0.1 - 90.0;
        double longitude = (double) ((i * 17U) % 3600U) * 0.1 - 180.0;

        succeeded = geo_index_add(index, first_id + i, latitude, longitude);
    }

    for (size_t i = 0; i < removed_count; ++i) {
        removed_ids[i] = first_id + i * 2U;
    }

    if (succeeded) {
        succeeded = geo_index_build(index) && geo_index_save(index, segment_path);
    }

    geo_index_destroy(index);

    GeoSegmentSet *set = succeeded ? geo_segment_set_create(manifest_path) : NULL;
    GeoSegmentCompactionPolicy policy = {
        .max_active_segments = 8,
        .minimum_mutations = removed_count,
        .maximum_mutation_rewrite_passes = 2,
        .mutation_ratio_numerator = 1,
        .mutation_ratio_denominator = 4,
    };

    if (succeeded) {
        succeeded = set &&
                    geo_segment_set_configure_background_compaction(set, "/tmp", &policy) &&
                    geo_segment_set_add_file(set, segment_path) &&
                    geo_segment_set_remove_ids(set, removed_ids, removed_count);
    }

    GeoSegmentCompactionStats compaction_stats;
    GeoSearchStats query_stats;
    size_t visible_count = 0;

    if (succeeded) {
        succeeded = geo_segment_set_wait_for_background_compaction(set, &compaction_stats) &&
                    !geo_segment_set_background_compaction_active(set) &&
                    compaction_stats.input_segments == 1 &&
                    compaction_stats.records_written == record_count - removed_count &&
                    geo_segment_set_count(set) == 1 &&
                    geo_segment_set_record_count(set) == record_count - removed_count &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &visible_count, &query_stats) &&
                    visible_count == record_count - removed_count &&
                    query_stats.records_scanned == 0 &&
                    geo_segment_set_copy_path(set, 0, compacted_path, sizeof(compacted_path)) &&
                    strcmp(compacted_path, segment_path) != 0;
    }

    geo_segment_set_destroy(set);
    set = succeeded ? geo_segment_set_open(manifest_path) : NULL;

    if (succeeded) {
        visible_count = 0;
        memset(&query_stats, 0, sizeof(query_stats));
        succeeded = set &&
                    geo_segment_set_record_count(set) == record_count - removed_count &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &visible_count, &query_stats) &&
                    visible_count == record_count - removed_count &&
                    query_stats.records_scanned == 0;
    }

    geo_segment_set_destroy(set);
    free(removed_ids);

    if (compacted_path[0]) {
        remove(compacted_path);
    }

    remove(manifest_path);
    remove(segment_path);
    remove(mutation_path);
    remove(mutation_checkpoint_path);

    return succeeded;
#else
    return 1;
#endif
}

#if defined(__unix__) || defined(__APPLE__)
typedef struct {
    GeoSegmentSet *set;
    const uint64_t *ids;
    size_t count;
    atomic_bool *start;
    atomic_bool *failed;
} ConcurrentRemoveWorker;

typedef struct {
    GeoSegmentSet *set;
    char (*paths)[128];
    size_t count;
    atomic_bool *start;
    atomic_bool *failed;
} ConcurrentInsertWorker;

typedef struct {
    GeoSegmentSet *set;
    size_t minimum_count;
    size_t maximum_count;
    atomic_bool *start;
    atomic_bool *stop;
    atomic_bool *failed;
} ConcurrentMutationReader;

typedef struct {
    GeoSegmentSet *set;
    const char *output_path;
    GeoSegmentCompactionStats stats;
    atomic_bool completed;
    bool succeeded;
} ConcurrentCompactionWorker;

static void wait_for_mutation_start(atomic_bool *start)
{
    while (!atomic_load_explicit(start, memory_order_acquire)) {
        sched_yield();
    }
}

static void *remove_ids_concurrently(void *argument)
{
    ConcurrentRemoveWorker *worker = argument;
    const size_t batch_size = 127;

    wait_for_mutation_start(worker->start);

    for (size_t position = 0; position < worker->count; position += batch_size) {
        size_t count = worker->count - position;

        if (count > batch_size) {
            count = batch_size;
        }

        if (!geo_segment_set_remove_ids(worker->set, worker->ids + position, count)) {
            atomic_store_explicit(worker->failed, true, memory_order_release);
            break;
        }
    }

    return NULL;
}

static void *insert_segments_concurrently(void *argument)
{
    ConcurrentInsertWorker *worker = argument;

    wait_for_mutation_start(worker->start);

    for (size_t i = 0; i < worker->count; ++i) {
        if (!geo_segment_set_add_file(worker->set, worker->paths[i])) {
            atomic_store_explicit(worker->failed, true, memory_order_release);
            break;
        }
    }

    return NULL;
}

static void *query_during_concurrent_mutations(void *argument)
{
    ConcurrentMutationReader *reader = argument;

    wait_for_mutation_start(reader->start);

    while (!atomic_load_explicit(reader->stop, memory_order_acquire)) {
        size_t count = 0;

        if (!geo_segment_set_search_radius_count(reader->set, 0.0, 0.0, 25000.0, &count, NULL) ||
            count < reader->minimum_count ||
            count > reader->maximum_count) {
            atomic_store_explicit(reader->failed, true, memory_order_release);
            break;
        }
    }

    return NULL;
}

static void *compact_segments_concurrently(void *argument)
{
    ConcurrentCompactionWorker *worker = argument;

    worker->succeeded = geo_segment_set_compact(worker->set, worker->output_path, &worker->stats);
    atomic_store_explicit(&worker->completed, true, memory_order_release);

    return NULL;
}

static bool write_concurrent_compaction_segment(const char *path,
                                                uint64_t first_id,
                                                size_t count,
                                                double latitude_origin,
                                                double longitude_origin)
{
    GeoIndex *index = geo_index_create(count);

    if (!index) {
        return false;
    }

    bool succeeded = true;

    for (size_t i = 0; succeeded && i < count; ++i) {
        double latitude = latitude_origin + (double) ((i * 17U) % 1009U) * 0.000001;
        double longitude = longitude_origin + (double) ((i * 31U) % 1013U) * 0.000001;

        succeeded = geo_index_add(index, first_id + i, latitude, longitude);
    }

    if (succeeded) {
        succeeded = geo_index_build(index) && geo_index_save(index, path);
    }

    geo_index_destroy(index);

    return succeeded;
}
#endif

int test_segment_set_upsert_replaces_active_locations(void)
{
#if defined(__unix__) || defined(__APPLE__)
    const uint64_t first_id = UINT64_C(910000000);
    char manifest_path[160];
    char mutation_path[192];
    char checkpoint_path[224];
    char base_path[160];
    char upsert_path[160];
    char duplicate_path[160];
    char extra_paths[3][160];
    char compacted_path[512] = { 0 };
    long process_id = (long) getpid();

    snprintf(manifest_path, sizeof(manifest_path), "/tmp/geobolt-upsert-manifest-%ld.bin", process_id);
    snprintf(mutation_path, sizeof(mutation_path), "%s.mutations", manifest_path);
    snprintf(checkpoint_path, sizeof(checkpoint_path), "%s.mutations.checkpoint", manifest_path);
    snprintf(base_path, sizeof(base_path), "/tmp/geobolt-upsert-base-%ld.bin", process_id);
    snprintf(upsert_path, sizeof(upsert_path), "/tmp/geobolt-upsert-live-%ld.bin", process_id);
    snprintf(duplicate_path, sizeof(duplicate_path), "/tmp/geobolt-upsert-duplicate-%ld.bin", process_id);
    remove(checkpoint_path);
    remove(mutation_path);
    remove(manifest_path);
    remove(duplicate_path);
    remove(upsert_path);
    remove(base_path);

    for (size_t i = 0; i < 3; ++i) {
        snprintf(extra_paths[i],
                 sizeof(extra_paths[i]),
                 "/tmp/geobolt-upsert-extra-%ld-%zu.bin",
                 process_id,
                 i);
        remove(extra_paths[i]);
    }

    bool succeeded = write_concurrent_compaction_segment(base_path, first_id, 3, 0.0, 0.0) &&
                     write_concurrent_compaction_segment(upsert_path, first_id, 2, 40.0, 40.0);
    GeoSegmentSet *set = succeeded ? geo_segment_set_create(manifest_path) : NULL;

    succeeded = set &&
                geo_segment_set_add_file(set, base_path) &&
                geo_segment_set_upsert_file(set, upsert_path);

    size_t global_count = 0;
    size_t old_location_count = 0;
    size_t new_location_count = 0;

    if (succeeded) {
        succeeded = geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &global_count, NULL) &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 10.0, &old_location_count, NULL) &&
                    geo_segment_set_search_radius_count(set, 40.0, 40.0, 10.0, &new_location_count, NULL) &&
                    global_count == 3 &&
                    old_location_count == 1 &&
                    new_location_count == 2;
    }

    GeoIndex *duplicate = succeeded ? geo_index_create(2) : NULL;

    if (succeeded) {
        succeeded = duplicate &&
                    geo_index_add(duplicate, first_id, -40.0, -40.0) &&
                    geo_index_add(duplicate, first_id, -41.0, -41.0) &&
                    geo_index_build(duplicate) &&
                    geo_index_save(duplicate, duplicate_path) &&
                    !geo_segment_set_upsert_file(set, duplicate_path);
    }

    geo_index_destroy(duplicate);

    if (succeeded) {
        geo_segment_set_destroy(set);
        set = geo_segment_set_open(manifest_path);
        global_count = 0;
        old_location_count = 0;
        new_location_count = 0;
        succeeded = set &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &global_count, NULL) &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 10.0, &old_location_count, NULL) &&
                    geo_segment_set_search_radius_count(set, 40.0, 40.0, 10.0, &new_location_count, NULL) &&
                    global_count == 3 &&
                    old_location_count == 1 &&
                    new_location_count == 2;
    }

    for (size_t i = 0; succeeded && i < 3; ++i) {
        uint64_t extra_id = first_id + UINT64_C(100) + i;

        succeeded = write_concurrent_compaction_segment(extra_paths[i], extra_id, 1, -40.0 + (double) i, -40.0) &&
                    geo_segment_set_add_file(set, extra_paths[i]);
    }

    GeoSegmentCompactionStats compaction_stats = { 0 };

    if (succeeded) {
        succeeded = geo_segment_set_enable_background_compaction(set, "/tmp", 4) &&
                    geo_segment_set_wait_for_background_compaction(set, &compaction_stats) &&
                    compaction_stats.input_segments == 4 &&
                    geo_segment_set_count(set) == 2 &&
                    geo_segment_set_record_count(set) == 8 &&
                    geo_segment_set_copy_path(set, 1, compacted_path, sizeof(compacted_path));
    }

    if (succeeded) {
        global_count = 0;
        old_location_count = 0;
        new_location_count = 0;
        succeeded = geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &global_count, NULL) &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 10.0, &old_location_count, NULL) &&
                    geo_segment_set_search_radius_count(set, 40.0, 40.0, 10.0, &new_location_count, NULL) &&
                    global_count == 6 &&
                    old_location_count == 1 &&
                    new_location_count == 2;
    }

    if (succeeded) {
        geo_segment_set_destroy(set);
        set = geo_segment_set_open(manifest_path);
        global_count = 0;
        succeeded = set &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &global_count, NULL) &&
                    global_count == 6;
    }

    geo_segment_set_destroy(set);
    remove(checkpoint_path);
    remove(mutation_path);
    remove(manifest_path);
    remove(duplicate_path);
    remove(upsert_path);
    remove(base_path);

    if (compacted_path[0]) {
        remove(compacted_path);
    }

    for (size_t i = 0; i < 3; ++i) {
        remove(extra_paths[i]);
    }

    return succeeded;
#else
    return 1;
#endif
}

int test_segment_set_concurrent_insert_remove(void)
{
#if defined(__unix__) || defined(__APPLE__)
    const size_t initial_count = 12000;
    const size_t removed_count = 5000;
    const size_t insert_segment_count = 6;
    const size_t records_per_insert = 1000;
    const uint64_t initial_id = UINT64_C(70000000);
    const uint64_t inserted_id = UINT64_C(80000000);
    char manifest_path[128];
    char initial_path[128];
    char compacted_path[128];
    char reinsert_path[128];
    char mutation_path[160];
    char mutation_checkpoint_path[192];
    char insert_paths[insert_segment_count][128];

    snprintf(manifest_path, sizeof(manifest_path), "/tmp/geobolt-mutation-manifest-%ld.bin", (long) getpid());
    snprintf(initial_path, sizeof(initial_path), "/tmp/geobolt-mutation-initial-%ld.bin", (long) getpid());
    snprintf(compacted_path, sizeof(compacted_path), "/tmp/geobolt-mutation-compact-%ld.bin", (long) getpid());
    snprintf(reinsert_path, sizeof(reinsert_path), "/tmp/geobolt-mutation-reinsert-%ld.bin", (long) getpid());
    snprintf(mutation_path, sizeof(mutation_path), "%s.mutations", manifest_path);
    snprintf(mutation_checkpoint_path, sizeof(mutation_checkpoint_path), "%s.mutations.checkpoint", manifest_path);
    remove(manifest_path);
    remove(mutation_path);
    remove(mutation_checkpoint_path);
    remove(initial_path);
    remove(compacted_path);
    remove(reinsert_path);

    for (size_t segment = 0; segment < insert_segment_count; ++segment) {
        snprintf(insert_paths[segment],
                 sizeof(insert_paths[segment]),
                 "/tmp/geobolt-mutation-insert-%ld-%zu.bin",
                 (long) getpid(),
                 segment);
        remove(insert_paths[segment]);
    }

    GeoIndex *initial = geo_index_create(initial_count);
    uint64_t *removed_ids = malloc(removed_count * sizeof(*removed_ids));
    bool succeeded = initial && removed_ids;

    for (size_t i = 0; succeeded && i < initial_count; ++i) {
        double latitude = ((double) (i % 181U) - 90.0) * 0.5;
        double longitude = ((double) ((i * 17U) % 721U) - 360.0) * 0.5;

        succeeded = geo_index_add(initial, initial_id + i, latitude, longitude);
    }

    for (size_t i = 0; i < removed_count; ++i) {
        removed_ids[i] = initial_id + i * 2U;
    }

    if (succeeded) {
        geo_index_build(initial);
        succeeded = geo_index_save(initial, initial_path);
    }

    geo_index_destroy(initial);

    for (size_t segment = 0; succeeded && segment < insert_segment_count; ++segment) {
        GeoIndex *insert = geo_index_create(records_per_insert);

        succeeded = insert != NULL;

        for (size_t i = 0; succeeded && i < records_per_insert; ++i) {
            size_t global = segment * records_per_insert + i;
            double latitude = ((double) (global % 101U) - 50.0) * 0.25;
            double longitude = ((double) ((global * 29U) % 401U) - 200.0) * 0.5;

            succeeded = geo_index_add(insert, inserted_id + global, latitude, longitude);
        }

        if (succeeded) {
            geo_index_build(insert);
            succeeded = geo_index_save(insert, insert_paths[segment]);
        }

        geo_index_destroy(insert);
    }

    GeoSegmentSet *set = succeeded ? geo_segment_set_create(manifest_path) : NULL;
    succeeded = set && geo_segment_set_add_file(set, initial_path);

    atomic_bool start;
    atomic_bool stop;
    atomic_bool failed;
    atomic_init(&start, false);
    atomic_init(&stop, false);
    atomic_init(&failed, false);

    ConcurrentRemoveWorker remove_worker = {
        .set = set,
        .ids = removed_ids,
        .count = removed_count,
        .start = &start,
        .failed = &failed,
    };
    ConcurrentInsertWorker insert_worker = {
        .set = set,
        .paths = insert_paths,
        .count = insert_segment_count,
        .start = &start,
        .failed = &failed,
    };
    ConcurrentMutationReader reader = {
        .set = set,
        .minimum_count = initial_count - removed_count,
        .maximum_count = initial_count + insert_segment_count * records_per_insert,
        .start = &start,
        .stop = &stop,
        .failed = &failed,
    };
    pthread_t remove_thread;
    pthread_t insert_thread;
    pthread_t reader_thread;
    bool remove_started = false;
    bool insert_started = false;
    bool reader_started = false;

    if (succeeded) {
        remove_started = pthread_create(&remove_thread, NULL, remove_ids_concurrently, &remove_worker) == 0;
        insert_started = pthread_create(&insert_thread, NULL, insert_segments_concurrently, &insert_worker) == 0;
        reader_started = pthread_create(&reader_thread, NULL, query_during_concurrent_mutations, &reader) == 0;
        succeeded = remove_started && insert_started && reader_started;
        atomic_store_explicit(&start, true, memory_order_release);
    }

    if (remove_started) {
        succeeded = pthread_join(remove_thread, NULL) == 0 && succeeded;
    }

    if (insert_started) {
        succeeded = pthread_join(insert_thread, NULL) == 0 && succeeded;
    }

    atomic_store_explicit(&stop, true, memory_order_release);

    if (reader_started) {
        succeeded = pthread_join(reader_thread, NULL) == 0 && succeeded;
    }

    size_t expected_count = initial_count - removed_count + insert_segment_count * records_per_insert;
    size_t radius_count = 0;
    size_t bbox_count = 0;
    GeoSearchResult *nearest = NULL;

    if (succeeded) {
        nearest = geo_segment_set_search_knn(set, 0.0, 0.0, 64, 25000.0, NULL);
        succeeded = !atomic_load_explicit(&failed, memory_order_acquire) &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &radius_count, NULL) &&
                    geo_segment_set_search_bbox_count(set, -90.0, 90.0, -180.0, 180.0, &bbox_count, NULL) &&
                    radius_count == expected_count &&
                    bbox_count == expected_count &&
                    nearest &&
                    nearest->count == 64;
    }

    geo_result_destroy(nearest);

    GeoIndex *reinsert = succeeded ? geo_index_create(1) : NULL;
    uint64_t reinserted_id = removed_ids[0];
    uint64_t deleted_after_insert_id = inserted_id;
    uint64_t reinserted_z = geo_encode(12.5, -33.25);

    if (succeeded) {
        succeeded = reinsert &&
                    geo_index_add(reinsert, reinserted_id, 12.5, -33.25);
    }

    if (succeeded) {
        geo_index_build(reinsert);
        succeeded = geo_index_save(reinsert, reinsert_path) &&
                    geo_segment_set_add_file(set, reinsert_path) &&
                    geo_segment_set_remove(set, deleted_after_insert_id);
    }

    geo_index_destroy(reinsert);

    GeoSearchResult *ordered_result = NULL;

    if (succeeded) {
        ordered_result = geo_segment_set_search_radius(set, 0.0, 0.0, 25000.0, NULL);
        size_t reinserted_occurrences = 0;
        size_t deleted_occurrences = 0;

        if (ordered_result) {
            for (size_t i = 0; i < ordered_result->count; ++i) {
                if (ordered_result->results[i].id == reinserted_id) {
                    reinserted_occurrences++;
                    succeeded = succeeded && ordered_result->results[i].z == reinserted_z;
                }

                deleted_occurrences += ordered_result->results[i].id == deleted_after_insert_id;
            }
        }

        succeeded = succeeded &&
                    ordered_result &&
                    ordered_result->count == expected_count &&
                    reinserted_occurrences == 1 &&
                    deleted_occurrences == 0;
    }

    geo_result_destroy(ordered_result);

    struct stat mutation_log_before;
    struct stat mutation_log_after;

    for (size_t repeat = 0; succeeded && repeat < 32; ++repeat) {
        succeeded = geo_segment_set_remove(set, deleted_after_insert_id);
    }

    if (succeeded) {
        succeeded = stat(mutation_path, &mutation_log_before) == 0 &&
                    geo_segment_set_checkpoint_mutations(set) &&
                    stat(mutation_checkpoint_path, &mutation_log_after) == 0 &&
                    mutation_log_after.st_size < mutation_log_before.st_size;
    }

    if (succeeded) {
        geo_segment_set_destroy(set);
        set = geo_segment_set_open(manifest_path);
        succeeded = set &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &radius_count, NULL) &&
                    radius_count == expected_count;
    }

    const size_t amplified_remove_count = 50000;
    uint64_t *amplified_remove_ids = succeeded ? malloc(amplified_remove_count * sizeof(*amplified_remove_ids)) : NULL;
    struct stat checkpoint_before_amplification;
    struct stat checkpoint_after_amplification;

    if (succeeded) {
        succeeded = amplified_remove_ids != NULL && stat(mutation_checkpoint_path, &checkpoint_before_amplification) == 0;
    }

    for (size_t i = 0; succeeded && i < amplified_remove_count; ++i) {
        amplified_remove_ids[i] = deleted_after_insert_id;
    }

    if (succeeded) {
        succeeded = geo_segment_set_remove_ids(set, amplified_remove_ids, amplified_remove_count) &&
                    stat(mutation_path, &checkpoint_after_amplification) == 0 &&
                    access(mutation_checkpoint_path, F_OK) != 0 &&
                    checkpoint_after_amplification.st_size <= checkpoint_before_amplification.st_size;
    }

    free(amplified_remove_ids);

    if (succeeded) {
        geo_segment_set_destroy(set);
        set = geo_segment_set_open(manifest_path);
        succeeded = set &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &radius_count, NULL) &&
                    radius_count == expected_count;
    }

    GeoSegmentCompactionStats stats;

    if (succeeded) {
        succeeded = geo_segment_set_compact(set, compacted_path, &stats) &&
                    stats.records_written == expected_count &&
                    geo_segment_set_record_count(set) == expected_count;
    }

    if (succeeded) {
        GeoSearchStats post_compaction_stats;
        struct stat mutation_status;

        geo_segment_set_destroy(set);
        set = geo_segment_set_open(manifest_path);
        succeeded = set &&
                    geo_segment_set_search_radius_count(set,
                                                        0.0,
                                                        0.0,
                                                        25000.0,
                                                        &radius_count,
                                                        &post_compaction_stats) &&
                    radius_count == expected_count &&
                    geo_segment_set_record_count(set) == expected_count &&
                    post_compaction_stats.records_scanned == 0 &&
                    stat(mutation_path, &mutation_status) == 0 &&
                    mutation_status.st_size == 0;
    }

    geo_segment_set_destroy(set);
    free(removed_ids);
    remove(compacted_path);
    remove(initial_path);
    remove(reinsert_path);
    remove(manifest_path);
    remove(mutation_path);
    remove(mutation_checkpoint_path);

    for (size_t segment = 0; segment < insert_segment_count; ++segment) {
        remove(insert_paths[segment]);
    }

    return succeeded;
#else
    return 1;
#endif
}

int test_segment_set_writes_during_compaction(void)
{
#if defined(__unix__) || defined(__APPLE__)
    const size_t records_per_source = 750000;
    const size_t late_record_count = 4096;
    const uint64_t first_source_id = UINT64_C(310000000);
    const uint64_t second_source_id = UINT64_C(320000000);
    const uint64_t late_id = UINT64_C(330000000);
    const uint64_t removed_id = first_source_id + records_per_source / 2U;
    char manifest_path[160];
    char source_paths[2][160];
    char late_path[160];
    char compacted_path[160];
    char mutation_path[192];
    char mutation_checkpoint_path[224];
    long process_id = (long) getpid();

    snprintf(manifest_path, sizeof(manifest_path), "/tmp/geobolt-live-compact-manifest-%ld.bin", process_id);
    snprintf(source_paths[0], sizeof(source_paths[0]), "/tmp/geobolt-live-compact-source-a-%ld.bin", process_id);
    snprintf(source_paths[1], sizeof(source_paths[1]), "/tmp/geobolt-live-compact-source-b-%ld.bin", process_id);
    snprintf(late_path, sizeof(late_path), "/tmp/geobolt-live-compact-late-%ld.bin", process_id);
    snprintf(compacted_path, sizeof(compacted_path), "/tmp/geobolt-live-compact-output-%ld.bin", process_id);
    snprintf(mutation_path, sizeof(mutation_path), "%s.mutations", manifest_path);
    snprintf(mutation_checkpoint_path, sizeof(mutation_checkpoint_path), "%s.mutations.checkpoint", manifest_path);

    remove(manifest_path);
    remove(source_paths[0]);
    remove(source_paths[1]);
    remove(late_path);
    remove(compacted_path);
    remove(mutation_path);
    remove(mutation_checkpoint_path);

    bool succeeded = write_concurrent_compaction_segment(source_paths[0], first_source_id, records_per_source, -23.55, -46.63) &&
                     write_concurrent_compaction_segment(source_paths[1], second_source_id, records_per_source, 40.71, -74.00) &&
                     write_concurrent_compaction_segment(late_path, late_id, late_record_count, 35.68, 139.69);

    GeoSegmentSet *set = succeeded ? geo_segment_set_create(manifest_path) : NULL;

    if (succeeded) {
        succeeded = set &&
                    geo_segment_set_add_file(set, source_paths[0]) &&
                    geo_segment_set_add_file(set, source_paths[1]);
    }

    ConcurrentCompactionWorker worker = {
        .set = set,
        .output_path = compacted_path,
    };
    atomic_init(&worker.completed, false);

    pthread_t compaction_thread;
    bool thread_started = false;

    if (succeeded) {
        thread_started = pthread_create(&compaction_thread, NULL, compact_segments_concurrently, &worker) == 0;
        succeeded = thread_started;
    }

    double observation_deadline = geo_get_time_ms() + 5000.0;
    bool observed_active = false;

    while (succeeded && geo_get_time_ms() < observation_deadline) {
        if (geo_segment_set_compaction_active(set)) {
            observed_active = true;
            break;
        }

        if (atomic_load_explicit(&worker.completed, memory_order_acquire)) {
            break;
        }

        sched_yield();
    }

    bool active_after_insert = false;
    bool active_after_remove = false;

    if (succeeded && observed_active) {
        succeeded = geo_segment_set_add_file(set, late_path);
        active_after_insert = succeeded && geo_segment_set_compaction_active(set);
    }

    if (succeeded && active_after_insert) {
        succeeded = geo_segment_set_remove(set, removed_id);
        active_after_remove = succeeded && geo_segment_set_compaction_active(set);
    }

    if (thread_started) {
        succeeded = pthread_join(compaction_thread, NULL) == 0 && succeeded;
    }

    size_t expected_visible_count = records_per_source * 2U + late_record_count - 1U;
    size_t visible_count = 0;

    if (succeeded) {
        succeeded = observed_active &&
                    active_after_insert &&
                    active_after_remove &&
                    worker.succeeded &&
                    worker.stats.input_segments == 2U &&
                    worker.stats.records_written == records_per_source * 2U &&
                    worker.stats.worker_count > 1U &&
                    worker.stats.partition_count > 1U &&
                    geo_segment_set_count(set) == 2U &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &visible_count, NULL) &&
                    visible_count == expected_visible_count;
    }

    if (succeeded) {
        geo_segment_set_destroy(set);
        set = geo_segment_set_open(manifest_path);
        succeeded = set &&
                    geo_segment_set_count(set) == 2U &&
                    geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &visible_count, NULL) &&
                    visible_count == expected_visible_count;
    }

    GeoSearchResult *visible = succeeded ? geo_segment_set_search_radius(set, 0.0, 0.0, 25000.0, NULL) : NULL;
    size_t removed_occurrences = 0;
    size_t late_occurrences = 0;

    if (visible) {
        for (size_t i = 0; i < visible->count; ++i) {
            removed_occurrences += visible->results[i].id == removed_id;
            late_occurrences += visible->results[i].id >= late_id && visible->results[i].id < late_id + late_record_count;
        }
    }

    succeeded = succeeded &&
                visible &&
                visible->count == expected_visible_count &&
                removed_occurrences == 0U &&
                late_occurrences == late_record_count;

    geo_result_destroy(visible);
    geo_segment_set_destroy(set);
    remove(compacted_path);
    remove(source_paths[0]);
    remove(source_paths[1]);
    remove(late_path);
    remove(manifest_path);
    remove(mutation_path);
    remove(mutation_checkpoint_path);

    return succeeded;
#else
    return 1;
#endif
}

int test_segment_batch_executor_matches_mutation_aware_queries(void)
{
#if defined(__unix__) || defined(__APPLE__)
    enum {
        segment_count = 3,
        records_per_segment = 4096,
        removed_count = 192,
        query_count = 96,
    };
    const uint64_t first_id = UINT64_C(410000000);
    char manifest_path[160];
    char segment_paths[segment_count][160];
    char reinsert_path[160];
    char mutation_path[192];
    char mutation_checkpoint_path[224];
    long process_id = (long) getpid();

    snprintf(manifest_path, sizeof(manifest_path), "/tmp/geobolt-batch-segment-manifest-%ld.bin", process_id);
    snprintf(reinsert_path, sizeof(reinsert_path), "/tmp/geobolt-batch-segment-reinsert-%ld.bin", process_id);
    snprintf(mutation_path, sizeof(mutation_path), "%s.mutations", manifest_path);
    snprintf(mutation_checkpoint_path, sizeof(mutation_checkpoint_path), "%s.mutations.checkpoint", manifest_path);

    remove(manifest_path);
    remove(reinsert_path);
    remove(mutation_path);
    remove(mutation_checkpoint_path);

    const double origins[segment_count][2] = {
        { -23.5505, -46.6333 },
        { 40.7128, -74.0060 },
        { 35.6762, 139.6503 },
    };
    bool succeeded = true;

    for (size_t segment = 0; succeeded && segment < segment_count; ++segment) {
        snprintf(segment_paths[segment],
                 sizeof(segment_paths[segment]),
                 "/tmp/geobolt-batch-segment-%ld-%zu.bin",
                 process_id,
                 segment);
        remove(segment_paths[segment]);
        succeeded = write_concurrent_compaction_segment(segment_paths[segment],
                                                        first_id + segment * records_per_segment,
                                                        records_per_segment,
                                                        origins[segment][0],
                                                        origins[segment][1]);
    }

    GeoSegmentSet *set = succeeded ? geo_segment_set_create(manifest_path) : NULL;

    for (size_t segment = 0; succeeded && segment < segment_count; ++segment) {
        succeeded = geo_segment_set_add_file(set, segment_paths[segment]);
    }

    uint64_t removed_ids[removed_count];

    for (size_t i = 0; i < removed_count; ++i) {
        removed_ids[i] = first_id + (i * 67U) % (segment_count * records_per_segment);
    }

    if (succeeded) {
        succeeded = geo_segment_set_remove_ids(set, removed_ids, removed_count) &&
                    write_concurrent_compaction_segment(reinsert_path, removed_ids[0], 1, -23.5505, -46.6333) &&
                    geo_segment_set_add_file(set, reinsert_path);
    }

    double latitudes[query_count];
    double longitudes[query_count];
    double radii_km[query_count];

    for (size_t query = 0; query < query_count; ++query) {
        size_t origin = query % segment_count;

        latitudes[query] = origins[origin][0] + (double) ((query * 7U) % 17U) * 0.00001;
        longitudes[query] = origins[origin][1] - (double) ((query * 11U) % 19U) * 0.00001;
        radii_km[query] = query % 8U == 0 ? 25000.0 : (query % 3U == 0 ? 0.05 : 5.0);
    }

    GeoQueryExecutorConfig config = {
        .thread_count = 4,
        .scheduling_chunk = 4,
        .pin_workers = false,
    };
    GeoQueryExecutor *executor = succeeded ? geo_query_executor_create_segment_set(set, &config) : NULL;
    GeoBatchIdResult *batch = succeeded ? geo_batch_id_result_create(query_count, 0) : NULL;
    GeoIdResult *exact = succeeded ? geo_id_result_create(256) : NULL;

    succeeded = succeeded && executor && batch && exact;

    for (size_t submission = 0; succeeded && submission < 2U; ++submission) {
        succeeded = geo_query_executor_search_radius_ids(executor,
                                                        latitudes,
                                                        longitudes,
                                                        radii_km,
                                                        query_count,
                                                        batch) &&
                    batch->query_count == query_count &&
                    batch->offsets[0] == 0U &&
                    batch->offsets[query_count] == batch->id_count;

        for (size_t query = 0; succeeded && query < query_count; ++query) {
            size_t begin = batch->offsets[query];
            size_t end = batch->offsets[query + 1U];
            size_t batch_count = end - begin;

            succeeded = geo_segment_set_search_radius_ids_reuse(set,
                                                                latitudes[query],
                                                                longitudes[query],
                                                                radii_km[query],
                                                                exact,
                                                                NULL) &&
                        exact->count == batch_count;

            if (succeeded && batch_count > 1U) {
                qsort(batch->ids + begin, batch_count, sizeof(*batch->ids), compare_uint64_values);
                qsort(exact->ids, exact->count, sizeof(*exact->ids), compare_uint64_values);
            }

            succeeded = succeeded &&
                        (!batch_count || memcmp(batch->ids + begin, exact->ids, batch_count * sizeof(*exact->ids)) == 0);
        }
    }

    geo_id_result_destroy(exact);
    geo_batch_id_result_destroy(batch);
    geo_query_executor_destroy(executor);
    geo_segment_set_destroy(set);
    remove(reinsert_path);
    remove(manifest_path);
    remove(mutation_path);
    remove(mutation_checkpoint_path);

    for (size_t segment = 0; segment < segment_count; ++segment) {
        remove(segment_paths[segment]);
    }

    return succeeded;
#else
    return 1;
#endif
}

#if defined(__unix__) || defined(__APPLE__)
typedef struct {
    GeoSegmentSet *set;
    uint64_t removed_id;
    atomic_bool started;
    atomic_bool completed;
    bool succeeded;
} SegmentSnapshotWriter;

typedef struct {
    const GeoSegmentSet *set;
    uint64_t *ids;
    size_t capacity;
    size_t count;
    atomic_bool completed;
    bool succeeded;
} SegmentSnapshotQuery;

static void *remove_during_segment_snapshot(void *argument)
{
    SegmentSnapshotWriter *writer = argument;

    atomic_store_explicit(&writer->started, true, memory_order_release);
    writer->succeeded = geo_segment_set_remove(writer->set, writer->removed_id);
    atomic_store_explicit(&writer->completed, true, memory_order_release);

    return NULL;
}

static void *query_with_preacquired_segment_snapshot(void *argument)
{
    SegmentSnapshotQuery *query = argument;
    size_t counted = 0;

    query->succeeded = geo_segment_set_search_radius_count_snapshot(query->set, 0.0, 0.0, 25000.0, &counted) &&
                       geo_segment_set_search_radius_ids_into_snapshot(query->set,
                                                                      0.0,
                                                                      0.0,
                                                                      25000.0,
                                                                      query->ids,
                                                                      query->capacity,
                                                                      &query->count) &&
                       query->count == counted;
    atomic_store_explicit(&query->completed, true, memory_order_release);

    return NULL;
}
#endif

int test_segment_snapshot_kernels_do_not_reenter_reader_gate(void)
{
#if defined(__unix__) || defined(__APPLE__)
    const size_t record_count = 4096;
    const uint64_t first_id = UINT64_C(520000000);
    const uint64_t removed_id = first_id + 127U;
    char manifest_path[160];
    char segment_path[160];
    char mutation_path[192];
    char mutation_checkpoint_path[224];
    long process_id = (long) getpid();

    snprintf(manifest_path, sizeof(manifest_path), "/tmp/geobolt-snapshot-manifest-%ld.bin", process_id);
    snprintf(segment_path, sizeof(segment_path), "/tmp/geobolt-snapshot-segment-%ld.bin", process_id);
    snprintf(mutation_path, sizeof(mutation_path), "%s.mutations", manifest_path);
    snprintf(mutation_checkpoint_path, sizeof(mutation_checkpoint_path), "%s.mutations.checkpoint", manifest_path);
    remove(manifest_path);
    remove(segment_path);
    remove(mutation_path);
    remove(mutation_checkpoint_path);

    bool succeeded = write_concurrent_compaction_segment(segment_path, first_id, record_count, 0.0, 0.0);
    GeoSegmentSet *set = succeeded ? geo_segment_set_create(manifest_path) : NULL;
    uint64_t *ids = malloc(record_count * sizeof(*ids));
    bool snapshot_acquired = false;

    succeeded = set && ids && geo_segment_set_add_file(set, segment_path);

    if (succeeded) {
        snapshot_acquired = geo_segment_set_read_snapshot_acquire(set);
        succeeded = snapshot_acquired;
    }

    SegmentSnapshotWriter writer = {
        .set = set,
        .removed_id = removed_id,
    };
    SegmentSnapshotQuery query = {
        .set = set,
        .ids = ids,
        .capacity = record_count,
    };
    atomic_init(&writer.started, false);
    atomic_init(&writer.completed, false);
    atomic_init(&query.completed, false);

    pthread_t writer_thread;
    pthread_t query_thread;
    bool writer_started = false;
    bool query_started = false;

    if (succeeded) {
        writer_started = pthread_create(&writer_thread, NULL, remove_during_segment_snapshot, &writer) == 0;
        succeeded = writer_started;
    }

    while (succeeded && !atomic_load_explicit(&writer.started, memory_order_acquire)) {
        sched_yield();
    }

    bool mutation_durable = false;

    for (size_t attempt = 0; succeeded && attempt < 2000 && !mutation_durable; ++attempt) {
        struct stat status;

        mutation_durable = stat(mutation_path, &status) == 0 && status.st_size >= (off_t) (sizeof(uint64_t) * 3U);

        if (!mutation_durable) {
            struct timespec delay = {
                .tv_nsec = 1000000,
            };

            nanosleep(&delay, NULL);
        }
    }

    succeeded = succeeded && mutation_durable && !atomic_load_explicit(&writer.completed, memory_order_acquire);

    if (succeeded) {
        struct timespec writer_gate_delay = {
            .tv_nsec = 20000000,
        };

        nanosleep(&writer_gate_delay, NULL);
    }

    if (succeeded) {
        query_started = pthread_create(&query_thread, NULL, query_with_preacquired_segment_snapshot, &query) == 0;
        succeeded = query_started;
    }

    bool query_completed_while_snapshot_held = false;

    for (size_t attempt = 0; query_started && attempt < 2000; ++attempt) {
        if (atomic_load_explicit(&query.completed, memory_order_acquire)) {
            query_completed_while_snapshot_held = true;
            break;
        }

        struct timespec delay = {
            .tv_nsec = 1000000,
        };

        nanosleep(&delay, NULL);
    }

    if (snapshot_acquired) {
        geo_segment_set_read_snapshot_release(set);
    }

    if (query_started) {
        succeeded = pthread_join(query_thread, NULL) == 0 && succeeded;
    }

    if (writer_started) {
        succeeded = pthread_join(writer_thread, NULL) == 0 && succeeded;
    }

    bool removed_id_was_visible = false;

    for (size_t i = 0; query.succeeded && i < query.count; ++i) {
        removed_id_was_visible |= ids[i] == removed_id;
    }

    size_t visible_after_remove = 0;

    succeeded = succeeded &&
                query_completed_while_snapshot_held &&
                query.succeeded &&
                query.count == record_count &&
                removed_id_was_visible &&
                writer.succeeded &&
                geo_segment_set_search_radius_count(set, 0.0, 0.0, 25000.0, &visible_after_remove, NULL) &&
                visible_after_remove == record_count - 1U;

    free(ids);
    geo_segment_set_destroy(set);
    remove(segment_path);
    remove(manifest_path);
    remove(mutation_path);
    remove(mutation_checkpoint_path);

    return succeeded;
#else
    return 1;
#endif
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

    size_t count_only = 0;

    ASSERT_TRUE(geo_search_radius_count(index, -23.5505, -46.6333, 500.0, &count_only, NULL),
                "count-only radius search should succeed");
    ASSERT_EQ(count_only, result->count, "count-only radius search should match materialized search");
    
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

    size_t count_only = 0;

    ASSERT_TRUE(geo_search_bbox_count(index, -35.0, -20.0, -60.0, -40.0, &count_only, NULL),
                "count-only bbox search should succeed");
    ASSERT_EQ(count_only, result->count, "count-only bbox search should match materialized search");

    printf("    Found %zu results in South America bbox\n", result->count);
    
    geo_result_destroy(result);
    geo_index_destroy(index);
    return 1;
}

int test_search_knn_basic(void)
{
    GeoIndex *index = geo_index_create(100);

    for (int i = 0; i < NUM_KNOWN_CITIES; i++) {
        geo_index_add(index, i, KNOWN_CITIES[i].lat, KNOWN_CITIES[i].lng);
    }

    geo_index_build(index);

    GeoSearchResult *result = geo_search_knn(index, -23.5505, -46.6333, 3, 10000.0, NULL);

    ASSERT_TRUE(result != NULL, "search should return result");
    ASSERT_EQ(result->count, 3, "should return exactly k results when enough candidates exist");
    ASSERT_EQ(result->results[0].id, 0, "the nearest result should be Sao Paulo");

    printf("    Found %zu nearest neighbors\n", result->count);

    geo_result_destroy(result);

    GeoSearchResult *limited = geo_search_knn(index, -23.5505, -46.6333, 3, 1.0, NULL);

    ASSERT_TRUE(limited != NULL, "a radius smaller than k coverage should still return a result object");
    ASSERT_EQ(limited->count, 1, "the radius limit should return only the center point");

    geo_result_destroy(limited);
    geo_index_destroy(index);

    return 1;
}

int test_search_knn_matches_brute_force(void)
{
    const size_t point_count = 10000;
    const size_t query_count = 24;
    GeoIndex *index = geo_index_create(point_count);

    ASSERT_TRUE(index != NULL, "kNN brute-force test index allocation should succeed");

    for (size_t i = 0; i < point_count; ++i) {
        double latitude = (double) ((i * 104729U) % 180000U) / 1000.0 - 90.0;
        double longitude = (double) ((i * 130363U) % 360000U) / 1000.0 - 180.0;

        ASSERT_TRUE(geo_index_add(index, i, latitude, longitude), "kNN brute-force point insertion should succeed");
    }

    geo_index_build(index);

    for (size_t query = 0; query < query_count; ++query) {
        double latitude = (double) ((query * 7919U) % 178000U) / 1000.0 - 89.0;
        double longitude = (double) ((query * 104729U) % 360000U) / 1000.0 - 180.0;
        double max_radius = 100.0 + (double) ((query * 1543U) % 12000U);
        size_t k = query % 31 + 1;

        if (query == 0) {
            latitude = 89.0;
            longitude = 179.9;
        } else if (query == 1) {
            latitude = -89.0;
            longitude = -179.9;
        }

        GeoSearchResult *expected = geo_result_create(64);

        ASSERT_TRUE(expected != NULL, "kNN brute-force result allocation should succeed");

        for (size_t i = 0; i < index->count; ++i) {
            GeoPoint point = geo_decode(index->records[i].z);

            if (geo_haversine_km(latitude, longitude, point.lat, point.lng) <= max_radius) {
                ASSERT_TRUE(geo_result_add(expected, index->records + i), "kNN brute-force result growth should succeed");
            }
        }

        geo_result_sort_by_distance(expected, latitude, longitude);

        GeoSearchResult *actual = geo_search_knn(index, latitude, longitude, k, max_radius, NULL);
        size_t expected_count = expected->count < k ? expected->count : k;

        ASSERT_TRUE(actual != NULL, "best-first kNN query should succeed");
        ASSERT_EQ(actual->count, expected_count, "best-first kNN count should match brute force");

        for (size_t i = 0; i < expected_count; ++i) {
            ASSERT_EQ(actual->results[i].id, expected->results[i].id, "best-first kNN ordering should match brute force");
        }

        geo_result_destroy(actual);
        geo_result_destroy(expected);
    }

    geo_index_destroy(index);

    return 1;
}

int test_search_matches_brute_force(void)
{
    const size_t point_count = 10000;
    const size_t query_count = 100;
    GeoIndex *index = geo_index_create(point_count);

    if (!index) {
        return 0;
    }

    for (size_t i = 0; i < point_count; ++i) {
        double latitude = (double)((i * 104729U) % 180000U) / 1000.0 - 90.0;
        double longitude = (double)((i * 130363U) % 360000U) / 1000.0 - 180.0;

        if (!geo_index_add(index, i, latitude, longitude)) {
            geo_index_destroy(index);

            return 0;
        }
    }

    geo_index_build(index);

    for (size_t query = 0; query < query_count; ++query) {
        double center_latitude = (double)((query * 7919U) % 170000U) / 1000.0 - 85.0;
        double center_longitude = (double)((query * 104729U) % 360000U) / 1000.0 - 180.0;
        double radius = 1.0 + (double)((query * 97U) % 2000U);
        size_t brute_force_count = 0;

        for (size_t i = 0; i < index->count; ++i) {
            GeoPoint point = geo_decode(index->records[i].z);

            if (geo_haversine_km(center_latitude, center_longitude, point.lat, point.lng) <= radius) {
                brute_force_count++;
            }
        }

        GeoSearchResult *result = geo_search_radius(index,
                                                    center_latitude,
                                                    center_longitude,
                                                    radius,
                                                    NULL);

        bool matches = result != NULL && result->count == brute_force_count;

        geo_result_destroy(result);

        if (!matches) {
            geo_index_destroy(index);

            return 0;
        }

        double min_latitude = geo_clamp_lat(center_latitude - 10.0);
        double max_latitude = geo_clamp_lat(center_latitude + 10.0);
        double min_longitude = geo_wrap_lng(center_longitude - 20.0);
        double max_longitude = geo_wrap_lng(center_longitude + 20.0);
        bool wraps_antimeridian = min_longitude > max_longitude;

        brute_force_count = 0;

        for (size_t i = 0; i < index->count; ++i) {
            GeoPoint point = geo_decode(index->records[i].z);
            bool longitude_matches = wraps_antimeridian ? point.lng >= min_longitude || point.lng <= max_longitude
                                                        : point.lng >= min_longitude && point.lng <= max_longitude;

            if (point.lat >= min_latitude && point.lat <= max_latitude && longitude_matches) {
                brute_force_count++;
            }
        }

        result = geo_search_bbox(index,
                                 min_latitude,
                                 max_latitude,
                                 min_longitude,
                                 max_longitude,
                                 NULL);

        matches = result != NULL && result->count == brute_force_count;

        geo_result_destroy(result);

        if (!matches) {
            geo_index_destroy(index);

            return 0;
        }
    }

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
    (void)sum;
    
    start = geo_get_time_ms();
    volatile double lat_sum = 0;
    
    for (int i = 0; i < iterations; i++) {
        GeoPoint p = geo_decode((uint64_t)i * 12345);
        lat_sum += p.lat;
    }
    
    double decode_time = geo_get_time_ms() - start;
    (void)lat_sum;
    
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
    (void)sum;
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
    
    // Allocate arrays for batch operations
    double *lats = (double*)malloc(n * sizeof(double));
    double *lngs = (double*)malloc(n * sizeof(double));
    uint64_t *ids = (uint64_t*)malloc(n * sizeof(uint64_t));

    if (!lats || !lngs || !ids) {
        free(ids);
        free(lngs);
        free(lats);

        return 0;
    }
    
    // Generate random data
    srand(42);
    for (int i = 0; i < n; i++) {
        lats[i] = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        lngs[i] = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
        ids[i] = i;
    }
    
    // === SCALAR INSERT ===
    GeoIndex *index = geo_index_create(n);
    ASSERT_TRUE(index != NULL, "failed to create large index");
    
    double start = geo_get_time_ms();
    for (int i = 0; i < n; i++) {
        geo_index_add(index, ids[i], lats[i], lngs[i]);
    }
    double scalar_insert_time = geo_get_time_ms() - start;
    
    start = geo_get_time_ms();
    geo_index_build(index);
    double build_time = geo_get_time_ms() - start;
    
    printf("    ┌─────────────────────────────────────────────────────────┐\n");
    printf("    │ LARGE DATASET TEST (%d records)                    │\n", n);
    printf("    ├─────────────────────────────────────────────────────────┤\n");
    printf("    │ Scalar Insert: %8.2f ms (%6.2f M ops/sec)           │\n",
           scalar_insert_time, n / scalar_insert_time / 1000.0);
    
#if GEO_SIMD_ENABLED
    // === SIMD ENCODE ===
    uint64_t *z_codes = (uint64_t*)malloc(n * sizeof(uint64_t));
    ASSERT_TRUE(z_codes != NULL, "failed to allocate z_codes");
    
    start = geo_get_time_ms();
    geo_simd_encode_batch(lats, lngs, z_codes, n);
    double simd_encode_time = geo_get_time_ms() - start;
    
    printf("    │ SIMD Encode:   %8.2f ms (%6.2f M ops/sec)           │\n",
           simd_encode_time, n / simd_encode_time / 1000.0);
    printf("    │ Encode Speedup: %6.2fx                                 │\n",
           scalar_insert_time / simd_encode_time);
    
    free(z_codes);
#endif
    
    printf("    │ Build time:    %8.2f ms                              │\n", build_time);
    printf("    ├─────────────────────────────────────────────────────────┤\n");
    
    // === SEARCH TESTS ===
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
    
    printf("    │ Avg search (50km): %.3f ms                            │\n", 
           total_search_time / num_searches);
    printf("    └─────────────────────────────────────────────────────────┘\n");
    
    geo_index_destroy(index);
    free(lats);
    free(lngs);
    free(ids);
    return 1;
}

int test_stress_many_searches(void)
{
    const size_t point_count = 100000;
    const size_t search_count = 10000;

    double *latitudes = malloc(point_count * sizeof(*latitudes));
    double *longitudes = malloc(point_count * sizeof(*longitudes));
    double *search_latitudes = malloc(search_count * sizeof(*search_latitudes));
    double *search_longitudes = malloc(search_count * sizeof(*search_longitudes));
    GeoIndex *index = geo_index_create(point_count);
    GeoSearchResult *reusable_result = geo_result_create(256);

    if (!latitudes || !longitudes || !search_latitudes || !search_longitudes || !index || !reusable_result) {
        free(latitudes);
        free(longitudes);
        free(search_latitudes);
        free(search_longitudes);
        geo_index_destroy(index);
        geo_result_destroy(reusable_result);

        return 0;
    }

    srand(42);

    for (size_t i = 0; i < point_count; ++i) {
        latitudes[i] = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        longitudes[i] = ((double)rand() / RAND_MAX) * 360.0 - 180.0;

        if (!geo_index_add(index, i, latitudes[i], longitudes[i])) {
            geo_result_destroy(reusable_result);
            geo_index_destroy(index);
            free(latitudes);
            free(longitudes);
            free(search_latitudes);
            free(search_longitudes);

            return 0;
        }
    }

    geo_index_build(index);
    srand(123);

    for (size_t i = 0; i < search_count; ++i) {
        search_latitudes[i] = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        search_longitudes[i] = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
    }

    size_t allocated_checksum = 0;
    double start = geo_get_time_ms();

    for (size_t i = 0; i < search_count; ++i) {
        GeoSearchResult *result = geo_search_radius(index, search_latitudes[i], search_longitudes[i], 100.0, NULL);

        if (!result) {
            geo_result_destroy(reusable_result);
            geo_index_destroy(index);
            free(latitudes);
            free(longitudes);
            free(search_latitudes);
            free(search_longitudes);

            return 0;
        }

        allocated_checksum += result->count;
        geo_result_destroy(result);
    }

    double allocated_time = geo_get_time_ms() - start;
    size_t reused_checksum = 0;

    start = geo_get_time_ms();

    for (size_t i = 0; i < search_count; ++i) {
        bool succeeded = geo_search_radius_reuse(index,
                                                 search_latitudes[i],
                                                 search_longitudes[i],
                                                 100.0,
                                                 reusable_result,
                                                 NULL);

        if (!succeeded) {
            geo_result_destroy(reusable_result);
            geo_index_destroy(index);
            free(latitudes);
            free(longitudes);
            free(search_latitudes);
            free(search_longitudes);

            return 0;
        }

        reused_checksum += reusable_result->count;
    }

    double reused_time = geo_get_time_ms() - start;

    printf("    ┌─────────────────────────────────────────────────────────┐\n");
    printf("    │ MANY SEARCHES TEST (%zu searches)                    │\n", search_count);
    printf("    ├─────────────────────────────────────────────────────────┤\n");
    printf("    │ Allocating result: %8.2f ms (%8.0f searches/sec)  │\n",
           allocated_time,
           search_count / allocated_time * 1000.0);
    printf("    │ Reusing result:    %8.2f ms (%8.0f searches/sec)  │\n",
           reused_time,
           search_count / reused_time * 1000.0);
    printf("    │ Reuse speedup:        %5.2fx                          │\n", allocated_time / reused_time);
    printf("    └─────────────────────────────────────────────────────────┘\n");

    geo_result_destroy(reusable_result);
    geo_index_destroy(index);
    free(latitudes);
    free(longitudes);
    free(search_latitudes);
    free(search_longitudes);

    return allocated_checksum == reused_checksum;
}

int test_stress_dense_area(void) {
    int n = 100000;
    
    double center_lat = -23.5505;
    double center_lng = -46.6333;
    double radius_km = 5.0;
    
    // Allocate arrays
    double *lats = (double*)malloc(n * sizeof(double));
    double *lngs = (double*)malloc(n * sizeof(double));
    
    srand(42);
    for (int i = 0; i < n; i++) {
        lats[i] = center_lat + ((double)rand() / RAND_MAX - 0.5) * 0.1;
        lngs[i] = center_lng + ((double)rand() / RAND_MAX - 0.5) * 0.1;
    }
    
    GeoIndex *index = geo_index_create(n);
    for (int i = 0; i < n; i++) {
        geo_index_add(index, i, lats[i], lngs[i]);
    }
    geo_index_build(index);
    
    // === SCALAR SEARCH ===
    GeoSearchStats stats;
    double start = geo_get_time_ms();
    GeoSearchResult *result = geo_search_radius(index, center_lat, center_lng, radius_km, &stats);
    double scalar_time = geo_get_time_ms() - start;
    size_t scalar_count = result->count;
    geo_result_destroy(result);
    
    printf("    ┌─────────────────────────────────────────────────────────┐\n");
    printf("    │ DENSE AREA TEST (%d points, %.1fkm radius)           │\n", n, radius_km);
    printf("    ├─────────────────────────────────────────────────────────┤\n");
    printf("    │ Scalar Search: %8.3f ms, found %zu points           │\n",
           scalar_time, scalar_count);
    
#if GEO_SIMD_ENABLED
    // === SIMD FILTER RADIUS ===
    uint8_t *mask = (uint8_t*)malloc(n * sizeof(uint8_t));
    
    start = geo_get_time_ms();
    size_t simd_count = geo_simd_filter_radius(lats, lngs, n, center_lat, center_lng, radius_km, mask);
    double simd_time = geo_get_time_ms() - start;
    
    printf("    │ SIMD Filter:   %8.3f ms, found %zu points           │\n",
           simd_time, simd_count);
    printf("    │ Filter Speedup: %6.2fx                                 │\n",
           scalar_time / simd_time);
    
    // === SIMD HAVERSINE BATCH ===
    double *dists = (double*)malloc(n * sizeof(double));
    
    start = geo_get_time_ms();
    geo_simd_haversine_batch(center_lat, center_lng, lats, lngs, dists, n);
    double haversine_time = geo_get_time_ms() - start;
    
    // Count matches
    size_t haversine_count = 0;
    for (int i = 0; i < n; i++) {
        if (dists[i] <= radius_km) haversine_count++;
    }
    
    printf("    │ SIMD Haversine: %7.3f ms, found %zu points           │\n",
           haversine_time, haversine_count);
    printf("    │ Haversine Speedup: %5.2fx                              │\n",
           scalar_time / haversine_time);
    
    free(mask);
    free(dists);
#endif
    
    printf("    └─────────────────────────────────────────────────────────┘\n");
    
    geo_index_destroy(index);
    free(lats);
    free(lngs);
    return 1;
}

// =========================================================
// Edge Case Tests
// =========================================================

int test_edge_antimeridian(void)
{
    GeoIndex *index = geo_index_create(10);

    geo_index_add(index, 1, 0.0, 179.9);
    geo_index_add(index, 2, 0.0, -179.9);
    geo_index_add(index, 3, 0.0, 180.0);
    geo_index_add(index, 4, 0.0, -180.0);
    geo_index_build(index);

    GeoSearchResult *result = geo_search_radius(index, 0.0, 180.0, 100.0, NULL);

    ASSERT_TRUE(result != NULL, "antimeridian search should succeed");
    ASSERT_EQ(result->count, 4, "antimeridian search should include points on both sides");

    printf("    Points near antimeridian: %zu\n", result->count);

    geo_result_destroy(result);
    geo_index_destroy(index);

    return 1;
}

int test_edge_poles(void)
{
    GeoIndex *index = geo_index_create(10);

    geo_index_add(index, 1, 89.9, 0.0);
    geo_index_add(index, 2, 89.9, 90.0);
    geo_index_add(index, 3, 89.9, 180.0);
    geo_index_add(index, 4, 89.9, -90.0);
    geo_index_build(index);

    GeoSearchResult *result = geo_search_radius(index, 90.0, 0.0, 100.0, NULL);

    ASSERT_TRUE(result != NULL, "polar search should succeed");
    ASSERT_EQ(result->count, 4, "polar search should span every longitude");

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
// SIMD Benchmark Tests
// =========================================================

#if GEO_SIMD_ENABLED

int test_simd_available(void) {
    bool available = geo_simd_available();
    const char *name = geo_simd_get_name();
    size_t batch_size = geo_simd_optimal_batch_size();
    
    printf("    SIMD Available: %s\n", available ? "YES" : "NO");
    printf("    SIMD Implementation: %s\n", name);
    printf("    Optimal Batch Size: %zu\n", batch_size);
    
    return 1;
}

int test_simd_encode_correctness(void)
{
    const size_t count = 1003;
    double *latitudes = malloc(count * sizeof(*latitudes));
    double *longitudes = malloc(count * sizeof(*longitudes));
    uint64_t *scalar_codes = malloc(count * sizeof(*scalar_codes));
    uint64_t *simd_codes = malloc(count * sizeof(*simd_codes));

    if (!latitudes || !longitudes || !scalar_codes || !simd_codes) {
        free(latitudes);
        free(longitudes);
        free(scalar_codes);
        free(simd_codes);

        return 0;
    }

    for (size_t i = 0; i < count; ++i) {
        latitudes[i] = (double)((i * 104729U) % 180000U) / 1000.0 - 90.0;
        longitudes[i] = (double)((i * 130363U) % 360000U) / 1000.0 - 180.0;
    }

    latitudes[0] = -90.0;
    longitudes[0] = -180.0;
    latitudes[1] = 90.0;
    longitudes[1] = 180.0;
    latitudes[2] = 0.0;
    longitudes[2] = 0.0;

    for (size_t i = 0; i < count; ++i) {
        scalar_codes[i] = geo_encode(latitudes[i], longitudes[i]);
    }

    geo_simd_encode_batch(latitudes, longitudes, simd_codes, count);

    int mismatches = 0;

    for (size_t i = 0; i < count; ++i) {
        if (scalar_codes[i] != simd_codes[i]) {
            mismatches++;

            if (mismatches <= 3) {
                printf("    Mismatch at %zu: scalar=%" PRIu64 " simd=%" PRIu64 "\n",
                       i,
                       scalar_codes[i],
                       simd_codes[i]);
            }
        }
    }

    printf("    Compared %zu encodings, %d mismatches\n", count, mismatches);

    free(latitudes);
    free(longitudes);
    free(scalar_codes);
    free(simd_codes);

    return mismatches == 0;
}

int test_simd_decode_correctness(void) {
    size_t count = 1000;
    uint64_t *z_codes = (uint64_t*)malloc(count * sizeof(uint64_t));
    double *lats_scalar = (double*)malloc(count * sizeof(double));
    double *lngs_scalar = (double*)malloc(count * sizeof(double));
    double *lats_simd = (double*)malloc(count * sizeof(double));
    double *lngs_simd = (double*)malloc(count * sizeof(double));
    
    if (!z_codes || !lats_scalar || !lngs_scalar || !lats_simd || !lngs_simd) {
        free(z_codes); free(lats_scalar); free(lngs_scalar);
        free(lats_simd); free(lngs_simd);
        return 0;
    }
    
    // Generate test data
    for (size_t i = 0; i < count; i++) {
        double lat = ((double)(i % 18000) / 100.0) - 90.0;
        double lng = ((double)(i % 36000) / 100.0) - 180.0;
        z_codes[i] = geo_encode(lat, lng);
    }
    
    // Scalar decoding
    for (size_t i = 0; i < count; i++) {
        GeoPoint p = geo_decode(z_codes[i]);
        lats_scalar[i] = p.lat;
        lngs_scalar[i] = p.lng;
    }
    
    // SIMD decoding
    geo_simd_decode_batch(z_codes, lats_simd, lngs_simd, count);
    
    // Compare results
    double max_lat_diff = 0, max_lng_diff = 0;
    for (size_t i = 0; i < count; i++) {
        double lat_diff = fabs(lats_scalar[i] - lats_simd[i]);
        double lng_diff = fabs(lngs_scalar[i] - lngs_simd[i]);
        if (lat_diff > max_lat_diff) max_lat_diff = lat_diff;
        if (lng_diff > max_lng_diff) max_lng_diff = lng_diff;
    }
    
    printf("    Max lat diff: %.10f, max lng diff: %.10f\n", max_lat_diff, max_lng_diff);
    
    free(z_codes); free(lats_scalar); free(lngs_scalar);
    free(lats_simd); free(lngs_simd);
    
    // Allow small precision differences due to SIMD approximations
    return max_lat_diff < 0.0001 && max_lng_diff < 0.0001;
}

int test_simd_haversine_accuracy(void)
{
    const size_t count = 1000;
    double *latitudes = malloc(count * sizeof(*latitudes));
    double *longitudes = malloc(count * sizeof(*longitudes));
    double *scalar_distances = malloc(count * sizeof(*scalar_distances));
    double *simd_distances = malloc(count * sizeof(*simd_distances));

    if (!latitudes || !longitudes || !scalar_distances || !simd_distances) {
        free(latitudes);
        free(longitudes);
        free(scalar_distances);
        free(simd_distances);

        return 0;
    }

    const double center_latitude = -23.5505;
    const double center_longitude = -46.6333;

    for (size_t i = 0; i < count; ++i) {
        latitudes[i] = (double)((i * 104729U) % 180000U) / 1000.0 - 90.0;
        longitudes[i] = (double)((i * 130363U) % 360000U) / 1000.0 - 180.0;
    }

    latitudes[0] = center_latitude;
    longitudes[0] = center_longitude;
    latitudes[1] = -center_latitude;
    longitudes[1] = geo_wrap_lng(center_longitude + 180.0);
    latitudes[2] = 90.0;
    longitudes[2] = 180.0;
    latitudes[3] = -90.0;
    longitudes[3] = -180.0;

    for (size_t i = 0; i < count; ++i) {
        scalar_distances[i] = geo_haversine_km(center_latitude,
                                               center_longitude,
                                               latitudes[i],
                                               longitudes[i]);
    }

    geo_simd_haversine_batch(center_latitude,
                             center_longitude,
                             latitudes,
                             longitudes,
                             simd_distances,
                             count);

    double max_absolute_error = 0.0;
    double max_relative_error = 0.0;

    for (size_t i = 0; i < count; ++i) {
        double absolute_error = fabs(scalar_distances[i] - simd_distances[i]);
        double relative_error = scalar_distances[i] > 1.0 ? absolute_error / scalar_distances[i] : 0.0;

        if (absolute_error > max_absolute_error) {
            max_absolute_error = absolute_error;
        }

        if (relative_error > max_relative_error) {
            max_relative_error = relative_error;
        }
    }

    printf("    Global max absolute diff: %.9f km\n", max_absolute_error);
    printf("    Global max relative error: %.9f%%\n", max_relative_error * 100.0);

    free(latitudes);
    free(longitudes);
    free(scalar_distances);
    free(simd_distances);

    return max_absolute_error < 0.00001 && max_relative_error < 0.00000001;
}

int test_simd_filter_bitset_correctness(void)
{
    const size_t count = 1003;
    double *latitudes = malloc(count * sizeof(*latitudes));
    double *longitudes = malloc(count * sizeof(*longitudes));
    uint64_t *codes = malloc(count * sizeof(*codes));
    uint8_t *byte_mask = malloc(count);
    uint64_t *bit_mask = calloc((count + 63U) / 64U, sizeof(*bit_mask));

    if (!latitudes || !longitudes || !codes || !byte_mask || !bit_mask) {
        free(bit_mask);
        free(byte_mask);
        free(codes);
        free(longitudes);
        free(latitudes);

        return 0;
    }

    for (size_t i = 0; i < count; ++i) {
        latitudes[i] = (double) ((i * 104729U) % 180000U) / 1000.0 - 90.0;
        longitudes[i] = (double) ((i * 130363U) % 360000U) / 1000.0 - 180.0;
        codes[i] = geo_encode(latitudes[i], longitudes[i]);
    }

    size_t byte_count = geo_simd_filter_bbox_codes(codes, count, -35.0, 18.0, 140.0, -155.0, byte_mask);
    size_t bit_count = geo_simd_filter_bbox_codes_bits(codes, count, -35.0, 18.0, 140.0, -155.0, bit_mask);
    bool succeeded = byte_count == bit_count;

    for (size_t i = 0; succeeded && i < count; ++i) {
        succeeded = byte_mask[i] == ((bit_mask[i >> 6] >> (i & 63U)) & 1U);
    }

    memset(bit_mask, 0, ((count + 63U) / 64U) * sizeof(*bit_mask));
    byte_count = geo_simd_filter_radius(latitudes, longitudes, count, -23.5505, -46.6333, 5000.0, byte_mask);
    bit_count = geo_simd_filter_radius_bits(latitudes, longitudes, count, -23.5505, -46.6333, 5000.0, bit_mask);
    succeeded = succeeded && byte_count == bit_count;

    for (size_t i = 0; succeeded && i < count; ++i) {
        succeeded = byte_mask[i] == ((bit_mask[i >> 6] >> (i & 63U)) & 1U);
    }

    free(bit_mask);
    free(byte_mask);
    free(codes);
    free(longitudes);
    free(latitudes);

    return succeeded;
}

int test_simd_exclude_id_bitset_tails(void)
{
    GeoRecord records[130];
    uint64_t actual_bits[3];
    uint64_t expected_bits[3];

    for (size_t i = 0; i < sizeof(records) / sizeof(records[0]); ++i) {
        records[i] = (GeoRecord) {
            .id = i % 11U,
            .z = i,
        };
    }

    for (size_t count = 0; count <= sizeof(records) / sizeof(records[0]); ++count) {
        memset(actual_bits, 0, sizeof(actual_bits));

        for (size_t i = 0; i < count; ++i) {
            if ((i % 5U) != 2U) {
                actual_bits[i >> 6U] |= UINT64_C(1) << (i & 63U);
            }
        }

        memcpy(expected_bits, actual_bits, sizeof(actual_bits));
        size_t expected_count = 0;

        for (size_t i = 0; i < count; ++i) {
            uint64_t bit = UINT64_C(1) << (i & 63U);

            if (records[i].id == 7U) {
                expected_bits[i >> 6U] &= ~bit;
            }

            expected_count += (expected_bits[i >> 6U] & bit) != 0;
        }

        size_t actual_count = geo_simd_exclude_id_bits(records, count, 7U, actual_bits);

        if (actual_count != expected_count || memcmp(actual_bits, expected_bits, sizeof(actual_bits)) != 0) {
            return 0;
        }
    }

    return 1;
}

int test_simd_benchmark_encode(void) {
    printf("    Running SIMD encode benchmark...\n");
    
    GeoSimdBenchmark result = geo_simd_benchmark_encode(100000, 10);
    
    printf("    ┌─────────────────────────────────────────────────────────┐\n");
    printf("    │ ENCODE BENCHMARK (100K points x 10 iterations)          │\n");
    printf("    ├─────────────────────────────────────────────────────────┤\n");
    printf("    │ Scalar:  %8.2f ms  (%6.2f M ops/sec)                 │\n",
           result.scalar_time_ms, 
           (result.operations / result.scalar_time_ms) / 1000.0);
    printf("    │ SIMD:    %8.2f ms  (%6.2f M ops/sec)                 │\n",
           result.simd_time_ms,
           (result.operations / result.simd_time_ms) / 1000.0);
    printf("    │ Speedup: %8.2fx                                      │\n", result.speedup);
    printf("    └─────────────────────────────────────────────────────────┘\n");
    
    return 1;
}

int test_simd_benchmark_decode(void) {
    printf("    Running SIMD decode benchmark...\n");
    
    GeoSimdBenchmark result = geo_simd_benchmark_decode(100000, 10);
    
    printf("    ┌─────────────────────────────────────────────────────────┐\n");
    printf("    │ DECODE BENCHMARK (100K codes x 10 iterations)           │\n");
    printf("    ├─────────────────────────────────────────────────────────┤\n");
    printf("    │ Scalar:  %8.2f ms  (%6.2f M ops/sec)                 │\n",
           result.scalar_time_ms, 
           (result.operations / result.scalar_time_ms) / 1000.0);
    printf("    │ SIMD:    %8.2f ms  (%6.2f M ops/sec)                 │\n",
           result.simd_time_ms,
           (result.operations / result.simd_time_ms) / 1000.0);
    printf("    │ Speedup: %8.2fx                                      │\n", result.speedup);
    printf("    └─────────────────────────────────────────────────────────┘\n");
    
    return 1;
}

int test_simd_benchmark_haversine(void) {
    printf("    Running SIMD haversine benchmark...\n");
    
    GeoSimdBenchmark result = geo_simd_benchmark_haversine(100000, 10);
    
    printf("    ┌─────────────────────────────────────────────────────────┐\n");
    printf("    │ HAVERSINE BENCHMARK (100K distances x 10 iterations)    │\n");
    printf("    ├─────────────────────────────────────────────────────────┤\n");
    printf("    │ Scalar:  %8.2f ms  (%6.2f M ops/sec)                 │\n",
           result.scalar_time_ms, 
           (result.operations / result.scalar_time_ms) / 1000.0);
    printf("    │ SIMD:    %8.2f ms  (%6.2f M ops/sec)                 │\n",
           result.simd_time_ms,
           (result.operations / result.simd_time_ms) / 1000.0);
    printf("    │ Speedup: %8.2fx                                      │\n", result.speedup);
    printf("    └─────────────────────────────────────────────────────────┘\n");
    
    return 1;
}

int test_simd_benchmark_filter_radius(void) {
    printf("    Running SIMD filter_radius benchmark...\n");
    
    GeoSimdBenchmark result = geo_simd_benchmark_filter_radius(100000, 5);
    
    printf("    ┌─────────────────────────────────────────────────────────┐\n");
    printf("    │ FILTER RADIUS BENCHMARK (100K points x 5 iterations)    │\n");
    printf("    ├─────────────────────────────────────────────────────────┤\n");
    printf("    │ Scalar:  %8.2f ms  (%6.2f M ops/sec)                 │\n",
           result.scalar_time_ms, 
           (result.operations / result.scalar_time_ms) / 1000.0);
    printf("    │ SIMD:    %8.2f ms  (%6.2f M ops/sec)                 │\n",
           result.simd_time_ms,
           (result.operations / result.simd_time_ms) / 1000.0);
    printf("    │ Speedup: %8.2fx                                      │\n", result.speedup);
    printf("    └─────────────────────────────────────────────────────────┘\n");
    
    return 1;
}

int test_simd_comprehensive_benchmark(void) {
    printf("\n");
    printf("    ╔═════════════════════════════════════════════════════════╗\n");
    printf("    ║     COMPREHENSIVE SIMD vs SCALAR BENCHMARK SUMMARY      ║\n");
    printf("    ╠═════════════════════════════════════════════════════════╣\n");
    printf("    ║ Architecture: %-42s ║\n", geo_simd_get_name());
    printf("    ╠═══════════════╦═══════════════╦═══════════════╦═════════╣\n");
    printf("    ║   Operation   ║  Scalar (ms)  ║   SIMD (ms)   ║ Speedup ║\n");
    printf("    ╠═══════════════╬═══════════════╬═══════════════╬═════════╣\n");
    
    GeoSimdBenchmark enc = geo_simd_benchmark_encode(50000, 20);
    GeoSimdBenchmark dec = geo_simd_benchmark_decode(50000, 20);
    GeoSimdBenchmark hav = geo_simd_benchmark_haversine(50000, 20);
    GeoSimdBenchmark flt = geo_simd_benchmark_filter_radius(50000, 10);
    
    printf("    ║ Encode        ║ %13.2f ║ %13.2f ║ %6.2fx ║\n",
           enc.scalar_time_ms, enc.simd_time_ms, enc.speedup);
    printf("    ║ Decode        ║ %13.2f ║ %13.2f ║ %6.2fx ║\n",
           dec.scalar_time_ms, dec.simd_time_ms, dec.speedup);
    printf("    ║ Haversine     ║ %13.2f ║ %13.2f ║ %6.2fx ║\n",
           hav.scalar_time_ms, hav.simd_time_ms, hav.speedup);
    printf("    ║ Filter Radius ║ %13.2f ║ %13.2f ║ %6.2fx ║\n",
           flt.scalar_time_ms, flt.simd_time_ms, flt.speedup);
    printf("    ╚═══════════════╩═══════════════╩═══════════════╩═════════╝\n");
    
    double avg_speedup = (enc.speedup + dec.speedup + hav.speedup + flt.speedup) / 4.0;
    printf("\n    Average Speedup: %.2fx\n", avg_speedup);
    
    return 1;
}

#endif // GEO_SIMD_ENABLED

// =========================================================
// Main
// =========================================================

int main(int argc, char **argv)
{
    g_test_filter = argc > 1 ? argv[1] : getenv("GEO_TEST_FILTER");

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
    RUN_TEST(test_index_persistence_roundtrip);
    RUN_TEST(test_persisted_checksum_is_independent_of_write_blocks);
    RUN_TEST(test_index_persistence_rejects_corruption);
    RUN_TEST(test_density_metadata_mmap_roundtrip);
    RUN_TEST(test_batch_executor_ids_match_exact_queries);
    RUN_TEST(test_stream_builder_multi_run_roundtrip);
    RUN_TEST(test_stream_builder_single_run_direct_copy_roundtrip);
    RUN_TEST(test_parallel_radix_build_preserves_records);
    RUN_TEST(test_parallel_radix_sorter_reuses_workers_and_scratch);
    RUN_TEST(test_segment_set_manifest_queries_and_compaction);
    RUN_TEST(test_segment_set_background_compaction);
    RUN_TEST(test_segment_set_background_tombstone_reclamation);
    RUN_TEST(test_segment_set_upsert_replaces_active_locations);
    RUN_TEST(test_segment_set_concurrent_insert_remove);
    RUN_TEST(test_segment_set_writes_during_compaction);
    RUN_TEST(test_segment_batch_executor_matches_mutation_aware_queries);
    RUN_TEST(test_segment_snapshot_kernels_do_not_reenter_reader_gate);
    
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
    RUN_TEST(test_search_knn_matches_brute_force);
    RUN_TEST(test_search_matches_brute_force);
    
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
    
#if GEO_SIMD_ENABLED
    // SIMD Tests
    printf("\n--- SIMD TESTS ---\n");
    printf("SIMD Implementation: %s\n", geo_simd_get_name());
    RUN_TEST(test_simd_available);
    RUN_TEST(test_simd_encode_correctness);
    RUN_TEST(test_simd_decode_correctness);
    RUN_TEST(test_simd_haversine_accuracy);
    RUN_TEST(test_simd_filter_bitset_correctness);
    RUN_TEST(test_simd_exclude_id_bitset_tails);
    
    // SIMD Benchmarks
    printf("\n--- SIMD BENCHMARKS (Scalar vs SIMD) ---\n");
    RUN_TEST(test_simd_benchmark_encode);
    RUN_TEST(test_simd_benchmark_decode);
    RUN_TEST(test_simd_benchmark_haversine);
    RUN_TEST(test_simd_benchmark_filter_radius);
    RUN_TEST(test_simd_comprehensive_benchmark);
#else
    printf("\n--- SIMD TESTS ---\n");
    printf("SIMD not available on this platform\n");
#endif
    
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
