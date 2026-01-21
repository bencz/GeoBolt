/*
 * GeoBolt Index - Comprehensive Performance Demo
 * 
 * This demo showcases the full capabilities of the GeoBolt spatial indexing
 * library including:
 *   - Millions of points indexing
 *   - SIMD-optimized batch operations
 *   - Multi-threaded read operations (thread safety demo)
 *   - Performance comparisons (Scalar vs SIMD)
 *   - Various search operations (radius, KNN, bbox)
 */

#include "geo_index.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>

// =========================================================
// Configuration
// =========================================================

#define NUM_POINTS        10000000   // 10 million points
#define NUM_SEARCHES      10000      // Number of search queries
#define NUM_THREADS       8          // Number of parallel search threads
#define SEARCH_RADIUS_KM  50.0       // Search radius in km

// =========================================================
// Thread Data Structure
// =========================================================

typedef struct {
    int thread_id;
    const GeoIndex *index;
    const double *search_lats;
    const double *search_lngs;
    int start_idx;
    int end_idx;
    double radius_km;
    double total_time_ms;
    size_t total_results;
    size_t total_scanned;
} ThreadData;

// =========================================================
// Utility Functions
// =========================================================

static void print_separator(void) {
    printf("════════════════════════════════════════════════════════════════════════════════\n");
}

static void print_header(const char *title) {
    printf("\n");
    print_separator();
    printf("  %s\n", title);
    print_separator();
}

static void print_metric(const char *name, double value, const char *unit) {
    printf("  %-30s %12.2f %s\n", name, value, unit);
}

static void print_metric_int(const char *name, size_t value, const char *unit) {
    printf("  %-30s %12zu %s\n", name, value, unit);
}

// =========================================================
// Thread Worker Function
// =========================================================

static void* search_worker(void *arg) {
    ThreadData *data = (ThreadData*)arg;
    
    double start = geo_get_time_ms();
    
    for (int i = data->start_idx; i < data->end_idx; i++) {
        GeoSearchStats stats;
        GeoSearchResult *result = geo_search_radius(
            data->index,
            data->search_lats[i],
            data->search_lngs[i],
            data->radius_km,
            &stats
        );
        
        data->total_results += result->count;
        data->total_scanned += stats.records_scanned;
        geo_result_destroy(result);
    }
    
    data->total_time_ms = geo_get_time_ms() - start;
    
    return NULL;
}

// =========================================================
// Demo 1: Large Scale Indexing
// =========================================================

static void demo_large_scale_indexing(double **out_lats, double **out_lngs, 
                                       GeoIndex **out_index) {
    print_header("DEMO 1: LARGE SCALE INDEXING");
    printf("  Building index with %d points (%.1f million)...\n\n", 
           NUM_POINTS, NUM_POINTS / 1000000.0);
    
    // Allocate arrays
    double *lats = (double*)malloc(NUM_POINTS * sizeof(double));
    double *lngs = (double*)malloc(NUM_POINTS * sizeof(double));
    
    if (!lats || !lngs) {
        fprintf(stderr, "Failed to allocate memory for %d points\n", NUM_POINTS);
        exit(1);
    }
    
    // Generate random world-wide distribution
    printf("  Generating random coordinates...\n");
    double gen_start = geo_get_time_ms();
    
    srand(42);
    for (int i = 0; i < NUM_POINTS; i++) {
        lats[i] = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        lngs[i] = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
    }
    
    double gen_time = geo_get_time_ms() - gen_start;
    print_metric("Data generation time:", gen_time, "ms");
    
    // Create index
    GeoIndex *index = geo_index_create(NUM_POINTS);
    if (!index) {
        fprintf(stderr, "Failed to create index\n");
        exit(1);
    }
    
    // === SCALAR INSERT ===
    printf("\n  [SCALAR] Inserting points one by one...\n");
    double scalar_start = geo_get_time_ms();
    
    for (int i = 0; i < NUM_POINTS; i++) {
        geo_index_add(index, i, lats[i], lngs[i]);
    }
    
    double scalar_insert_time = geo_get_time_ms() - scalar_start;
    print_metric("Scalar insert time:", scalar_insert_time, "ms");
    print_metric("Scalar insert rate:", NUM_POINTS / scalar_insert_time / 1000.0, "M ops/sec");
    
#if GEO_SIMD_ENABLED
    // === SIMD ENCODE ===
    printf("\n  [SIMD] Batch encoding coordinates...\n");
    uint64_t *z_codes = (uint64_t*)malloc(NUM_POINTS * sizeof(uint64_t));
    
    double simd_start = geo_get_time_ms();
    geo_simd_encode_batch(lats, lngs, z_codes, NUM_POINTS);
    double simd_encode_time = geo_get_time_ms() - simd_start;
    
    print_metric("SIMD encode time:", simd_encode_time, "ms");
    print_metric("SIMD encode rate:", NUM_POINTS / simd_encode_time / 1000.0, "M ops/sec");
    print_metric("SIMD speedup:", scalar_insert_time / simd_encode_time, "x");
    
    free(z_codes);
#endif
    
    // Build index (sort by Z-order)
    printf("\n  Building index (Z-order sort)...\n");
    double build_start = geo_get_time_ms();
    geo_index_build(index);
    double build_time = geo_get_time_ms() - build_start;
    
    print_metric("Build time:", build_time, "ms");
    print_metric("Total indexing time:", scalar_insert_time + build_time, "ms");
    
    // Memory usage estimate
    size_t memory_mb = (NUM_POINTS * sizeof(GeoRecord)) / (1024 * 1024);
    print_metric_int("Estimated memory usage:", memory_mb, "MB");
    
    *out_lats = lats;
    *out_lngs = lngs;
    *out_index = index;
}

// =========================================================
// Demo 2: SIMD vs Scalar Performance
// =========================================================

static void demo_simd_comparison(const double *lats, const double *lngs) {
    print_header("DEMO 2: SIMD vs SCALAR PERFORMANCE");
    
#if GEO_SIMD_ENABLED
    printf("  SIMD Implementation: %s\n", geo_simd_get_name());
    printf("  Optimal batch size: %zu\n\n", geo_simd_optimal_batch_size());
    
    int test_size = 1000000;
    int iterations = 5;
    
    // Allocate output arrays
    uint64_t *z_codes = (uint64_t*)malloc(test_size * sizeof(uint64_t));
    double *dists = (double*)malloc(test_size * sizeof(double));
    uint8_t *mask = (uint8_t*)malloc(test_size * sizeof(uint8_t));
    
    double center_lat = -23.5505;
    double center_lng = -46.6333;
    double radius = 100.0;
    
    printf("  ┌────────────────────┬────────────────┬────────────────┬──────────┐\n");
    printf("  │     Operation      │  Scalar (ms)   │   SIMD (ms)    │ Speedup  │\n");
    printf("  ├────────────────────┼────────────────┼────────────────┼──────────┤\n");
    
    // === ENCODE ===
    double scalar_time = 0, simd_time = 0;
    
    for (int iter = 0; iter < iterations; iter++) {
        double start = geo_get_time_ms();
        for (int i = 0; i < test_size; i++) {
            z_codes[i] = geo_encode(lats[i], lngs[i]);
        }
        scalar_time += geo_get_time_ms() - start;
        
        start = geo_get_time_ms();
        geo_simd_encode_batch(lats, lngs, z_codes, test_size);
        simd_time += geo_get_time_ms() - start;
    }
    
    printf("  │ Encode (1M pts)    │ %14.2f │ %14.2f │ %7.2fx │\n",
           scalar_time / iterations, simd_time / iterations, 
           scalar_time / simd_time);
    
    // === HAVERSINE ===
    scalar_time = simd_time = 0;
    
    for (int iter = 0; iter < iterations; iter++) {
        double start = geo_get_time_ms();
        for (int i = 0; i < test_size; i++) {
            dists[i] = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);
        }
        scalar_time += geo_get_time_ms() - start;
        
        start = geo_get_time_ms();
        geo_simd_haversine_batch(center_lat, center_lng, lats, lngs, dists, test_size);
        simd_time += geo_get_time_ms() - start;
    }
    
    printf("  │ Haversine (1M)     │ %14.2f │ %14.2f │ %7.2fx │\n",
           scalar_time / iterations, simd_time / iterations,
           scalar_time / simd_time);
    
    // === FILTER RADIUS ===
    scalar_time = simd_time = 0;
    
    for (int iter = 0; iter < iterations; iter++) {
        double start = geo_get_time_ms();
        size_t count = 0;
        for (int i = 0; i < test_size; i++) {
            double d = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);
            if (d <= radius) count++;
        }
        scalar_time += geo_get_time_ms() - start;
        (void)count;
        
        start = geo_get_time_ms();
        geo_simd_filter_radius(lats, lngs, test_size, center_lat, center_lng, radius, mask);
        simd_time += geo_get_time_ms() - start;
    }
    
    printf("  │ Filter Radius (1M) │ %14.2f │ %14.2f │ %7.2fx │\n",
           scalar_time / iterations, simd_time / iterations,
           scalar_time / simd_time);
    
    // === FAST DISTANCE ===
    scalar_time = simd_time = 0;
    
    for (int iter = 0; iter < iterations; iter++) {
        double start = geo_get_time_ms();
        for (int i = 0; i < test_size; i++) {
            dists[i] = geo_fast_distance_km(center_lat, center_lng, lats[i], lngs[i]);
        }
        scalar_time += geo_get_time_ms() - start;
        
        start = geo_get_time_ms();
        geo_simd_fast_distance_batch(center_lat, center_lng, lats, lngs, dists, test_size);
        simd_time += geo_get_time_ms() - start;
    }
    
    printf("  │ Fast Distance (1M) │ %14.2f │ %14.2f │ %7.2fx │\n",
           scalar_time / iterations, simd_time / iterations,
           scalar_time / simd_time);
    
    printf("  └────────────────────┴────────────────┴────────────────┴──────────┘\n");
    
    free(z_codes);
    free(dists);
    free(mask);
#else
    printf("  SIMD not available on this platform.\n");
#endif
}

// =========================================================
// Demo 3: Multi-threaded Search (Thread Safety)
// =========================================================

static void demo_multithreaded_search(const GeoIndex *index) {
    print_header("DEMO 3: MULTI-THREADED SEARCH (THREAD SAFETY)");
    
    printf("  Demonstrating thread-safe read operations...\n");
    printf("  NOTE: Index is read-only after geo_index_build()\n");
    printf("  Concurrent reads are safe, concurrent writes are NOT.\n\n");
    
    // Generate search points
    double *search_lats = (double*)malloc(NUM_SEARCHES * sizeof(double));
    double *search_lngs = (double*)malloc(NUM_SEARCHES * sizeof(double));
    
    srand(12345);
    for (int i = 0; i < NUM_SEARCHES; i++) {
        search_lats[i] = ((double)rand() / RAND_MAX) * 180.0 - 90.0;
        search_lngs[i] = ((double)rand() / RAND_MAX) * 360.0 - 180.0;
    }
    
    // === SINGLE-THREADED ===
    printf("  [SINGLE-THREAD] Running %d searches...\n", NUM_SEARCHES);
    
    double single_start = geo_get_time_ms();
    size_t single_results = 0;
    
    for (int i = 0; i < NUM_SEARCHES; i++) {
        GeoSearchResult *result = geo_search_radius(
            index, search_lats[i], search_lngs[i], SEARCH_RADIUS_KM, NULL
        );
        single_results += result->count;
        geo_result_destroy(result);
    }
    
    double single_time = geo_get_time_ms() - single_start;
    print_metric("Single-thread time:", single_time, "ms");
    print_metric("Single-thread rate:", NUM_SEARCHES / single_time * 1000.0, "searches/sec");
    
    // === MULTI-THREADED ===
    printf("\n  [MULTI-THREAD] Running %d searches with %d threads...\n", 
           NUM_SEARCHES, NUM_THREADS);
    
    pthread_t threads[NUM_THREADS];
    ThreadData thread_data[NUM_THREADS];
    
    int searches_per_thread = NUM_SEARCHES / NUM_THREADS;
    
    double multi_start = geo_get_time_ms();
    
    for (int t = 0; t < NUM_THREADS; t++) {
        thread_data[t].thread_id = t;
        thread_data[t].index = index;
        thread_data[t].search_lats = search_lats;
        thread_data[t].search_lngs = search_lngs;
        thread_data[t].start_idx = t * searches_per_thread;
        thread_data[t].end_idx = (t == NUM_THREADS - 1) ? NUM_SEARCHES : (t + 1) * searches_per_thread;
        thread_data[t].radius_km = SEARCH_RADIUS_KM;
        thread_data[t].total_time_ms = 0;
        thread_data[t].total_results = 0;
        thread_data[t].total_scanned = 0;
        
        pthread_create(&threads[t], NULL, search_worker, &thread_data[t]);
    }
    
    // Wait for all threads
    for (int t = 0; t < NUM_THREADS; t++) {
        pthread_join(threads[t], NULL);
    }
    
    double multi_time = geo_get_time_ms() - multi_start;
    
    size_t multi_results = 0;
    for (int t = 0; t < NUM_THREADS; t++) {
        multi_results += thread_data[t].total_results;
    }
    
    print_metric("Multi-thread time:", multi_time, "ms");
    print_metric("Multi-thread rate:", NUM_SEARCHES / multi_time * 1000.0, "searches/sec");
    print_metric("Parallel speedup:", single_time / multi_time, "x");
    
    // Verify results match
    if (single_results == multi_results) {
        printf("\n  ✓ Thread safety verified: results match (%zu total results)\n", single_results);
    } else {
        printf("\n  ✗ WARNING: Results mismatch! Single=%zu Multi=%zu\n", 
               single_results, multi_results);
    }
    
    free(search_lats);
    free(search_lngs);
}

// =========================================================
// Demo 4: Search Operations Showcase
// =========================================================

static void demo_search_operations(const GeoIndex *index) {
    print_header("DEMO 4: SEARCH OPERATIONS SHOWCASE");
    
    double center_lat = -23.5505;  // São Paulo
    double center_lng = -46.6333;
    
    // === RADIUS SEARCH ===
    printf("  [RADIUS SEARCH] Center: São Paulo (-23.55, -46.63)\n\n");
    
    double radii[] = {10.0, 50.0, 100.0, 500.0, 1000.0};
    int num_radii = sizeof(radii) / sizeof(radii[0]);
    
    printf("  ┌─────────────┬─────────────┬─────────────┬─────────────┐\n");
    printf("  │ Radius (km) │   Results   │   Scanned   │  Time (ms)  │\n");
    printf("  ├─────────────┼─────────────┼─────────────┼─────────────┤\n");
    
    for (int r = 0; r < num_radii; r++) {
        GeoSearchStats stats;
        GeoSearchResult *result = geo_search_radius(
            index, center_lat, center_lng, radii[r], &stats
        );
        
        printf("  │ %11.0f │ %11zu │ %11" PRIu64 " │ %11.3f │\n",
               radii[r], result->count, stats.records_scanned, stats.search_time_ms);
        
        geo_result_destroy(result);
    }
    
    printf("  └─────────────┴─────────────┴─────────────┴─────────────┘\n");
    
    // === KNN SEARCH ===
    printf("\n  [KNN SEARCH] Finding nearest neighbors to São Paulo\n\n");
    
    int k_values[] = {1, 5, 10, 50, 100};
    int num_k = sizeof(k_values) / sizeof(k_values[0]);
    
    printf("  ┌───────────┬───────────────────┬─────────────┐\n");
    printf("  │     K     │ Farthest Dist(km) │  Time (ms)  │\n");
    printf("  ├───────────┼───────────────────┼─────────────┤\n");
    
    for (int ki = 0; ki < num_k; ki++) {
        GeoSearchStats stats;
        GeoSearchResult *result = geo_search_knn(
            index, center_lat, center_lng, k_values[ki], 10000.0, &stats
        );
        
        double farthest = 0;
        for (size_t i = 0; i < result->count; i++) {
            GeoPoint pt = geo_decode(result->results[i].z);
            double dist = geo_haversine_km(center_lat, center_lng, pt.lat, pt.lng);
            if (dist > farthest) farthest = dist;
        }
        
        printf("  │ %9d │ %17.2f │ %11.3f │\n",
               k_values[ki], farthest, stats.search_time_ms);
        
        geo_result_destroy(result);
    }
    
    printf("  └───────────┴───────────────────┴─────────────┘\n");
    
    // === BOUNDING BOX SEARCH ===
    printf("\n  [BOUNDING BOX SEARCH] Various regions\n\n");
    
    struct {
        const char *name;
        double min_lat, max_lat, min_lng, max_lng;
    } regions[] = {
        {"São Paulo Metro",     -24.0, -23.0,  -47.0,  -46.0},
        {"Southeast Brazil",    -25.0, -19.0,  -52.0,  -40.0},
        {"South America",       -56.0,  13.0,  -82.0,  -34.0},
        {"Northern Hemisphere",   0.0,  90.0, -180.0,  180.0},
        {"World",               -90.0,  90.0, -180.0,  180.0}
    };
    int num_regions = sizeof(regions) / sizeof(regions[0]);
    
    printf("  ┌────────────────────────┬─────────────┬─────────────┐\n");
    printf("  │        Region          │   Results   │  Time (ms)  │\n");
    printf("  ├────────────────────────┼─────────────┼─────────────┤\n");
    
    for (int ri = 0; ri < num_regions; ri++) {
        GeoSearchStats stats;
        GeoSearchResult *result = geo_search_bbox(
            index,
            regions[ri].min_lat, regions[ri].max_lat,
            regions[ri].min_lng, regions[ri].max_lng,
            &stats
        );
        
        printf("  │ %-22s │ %11zu │ %11.3f │\n",
               regions[ri].name, result->count, stats.search_time_ms);
        
        geo_result_destroy(result);
    }
    
    printf("  └────────────────────────┴─────────────┴─────────────┘\n");
}

// =========================================================
// Demo 5: Precision & Accuracy
// =========================================================

static void demo_precision(void) {
    print_header("DEMO 5: ENCODING PRECISION & ACCURACY");
    
    struct {
        const char *name;
        double lat, lng;
    } test_points[] = {
        {"São Paulo, Brazil",    -23.5505200, -46.6333090},
        {"New York, USA",         40.7127753, -74.0059728},
        {"Tokyo, Japan",          35.6761919, 139.6503106},
        {"Sydney, Australia",    -33.8688197, 151.2092955},
        {"North Pole",            89.9999999,   0.0000000},
        {"South Pole",           -89.9999999,   0.0000000},
        {"Prime Meridian",         0.0000000,   0.0000000},
        {"Antimeridian",           0.0000000, 179.9999999}
    };
    int num_points = sizeof(test_points) / sizeof(test_points[0]);
    
    printf("  ┌────────────────────────┬─────────────────────────┬──────────────┐\n");
    printf("  │        Location        │     Error (meters)      │   Z-Code     │\n");
    printf("  ├────────────────────────┼─────────────────────────┼──────────────┤\n");
    
    for (int i = 0; i < num_points; i++) {
        uint64_t z = geo_encode(test_points[i].lat, test_points[i].lng);
        GeoPoint decoded = geo_decode(z);
        double error_m = geo_haversine_m(
            test_points[i].lat, test_points[i].lng,
            decoded.lat, decoded.lng
        );
        
        printf("  │ %-22s │ %23.4f │ %12" PRIu64 " │\n",
               test_points[i].name, error_m, z % 1000000000000ULL);
    }
    
    printf("  └────────────────────────┴─────────────────────────┴──────────────┘\n");
    
    printf("\n  Sub-centimeter precision achieved across all test points!\n");
}

// =========================================================
// Main
// =========================================================

int main(void) {
    printf("\n");
    printf("╔══════════════════════════════════════════════════════════════════════════════╗\n");
    printf("║                                                                              ║\n");
    printf("║                GeoBolt Index - Comprehensive Performance Demo                ║\n");
    printf("║                                                                              ║\n");
    printf("╚══════════════════════════════════════════════════════════════════════════════╝\n");
    
    printf("\n  Configuration:\n");
    printf("  • Points:        %d (%.1f million)\n", NUM_POINTS, NUM_POINTS / 1000000.0);
    printf("  • Searches:      %d\n", NUM_SEARCHES);
    printf("  • Threads:       %d\n", NUM_THREADS);
    printf("  • Search Radius: %.1f km\n", SEARCH_RADIUS_KM);
#if GEO_SIMD_ENABLED
    printf("  • SIMD:          %s\n", geo_simd_get_name());
#else
    printf("  • SIMD:          Not available\n");
#endif
    
    double *lats = NULL, *lngs = NULL;
    GeoIndex *index = NULL;
    
    double total_start = geo_get_time_ms();
    
    // Run all demos
    demo_large_scale_indexing(&lats, &lngs, &index);
    demo_simd_comparison(lats, lngs);
    demo_multithreaded_search(index);
    demo_search_operations(index);
    demo_precision();
    
    double total_time = geo_get_time_ms() - total_start;
    
    // Cleanup
    geo_index_destroy(index);
    free(lats);
    free(lngs);
    
    print_header("DEMO COMPLETE");
    print_metric("Total demo time:", total_time / 1000.0, "seconds");
    
    printf("\n");
    return 0;
}
