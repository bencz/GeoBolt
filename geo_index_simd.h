/*
 * GeoIndex SIMD Optimizations Header
 * 
 * This header provides architecture detection and declares SIMD-optimized
 * functions for geo-indexing operations. The actual implementations are
 * in separate files for each architecture:
 *   - geo_index_simd_arm64.c  (ARM64 NEON)
 *   - geo_index_simd_x86.c    (x86-64 AVX2/SSE4)
 */

#ifndef GEO_INDEX_SIMD_H
#define GEO_INDEX_SIMD_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// =========================================================
// Architecture Detection
// =========================================================

#if defined(__aarch64__) || defined(_M_ARM64)
    #define GEO_SIMD_ARM64 1
    #define GEO_SIMD_ENABLED 1
    #define GEO_SIMD_NAME "ARM64 NEON"
#elif defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    #if defined(__AVX2__)
        #define GEO_SIMD_X86_AVX2 1
        #define GEO_SIMD_ENABLED 1
        #define GEO_SIMD_NAME "x86-64 AVX2"
    #elif defined(__SSE4_1__) || defined(__SSE4_2__)
        #define GEO_SIMD_X86_SSE4 1
        #define GEO_SIMD_ENABLED 1
        #define GEO_SIMD_NAME "x86-64 SSE4"
    #elif defined(__SSE2__)
        #define GEO_SIMD_X86_SSE2 1
        #define GEO_SIMD_ENABLED 1
        #define GEO_SIMD_NAME "x86-64 SSE2"
    #endif
#endif

#ifndef GEO_SIMD_ENABLED
    #define GEO_SIMD_ENABLED 0
    #define GEO_SIMD_NAME "Scalar (no SIMD)"
#endif

// =========================================================
// SIMD Configuration
// =========================================================

#define GEO_SIMD_ALIGN_ARM64 16
#define GEO_SIMD_ALIGN_X86   32

#if defined(GEO_SIMD_ARM64)
    #define GEO_SIMD_ALIGN GEO_SIMD_ALIGN_ARM64
    #define GEO_SIMD_VECTOR_SIZE 4
#elif defined(GEO_SIMD_X86_AVX2)
    #define GEO_SIMD_ALIGN GEO_SIMD_ALIGN_X86
    #define GEO_SIMD_VECTOR_SIZE 8
#else
    #define GEO_SIMD_ALIGN 16
    #define GEO_SIMD_VECTOR_SIZE 4
#endif

// Alignment macro
#define GEO_SIMD_ALIGNED __attribute__((aligned(GEO_SIMD_ALIGN)))

// =========================================================
// SIMD Function Declarations - Batch Encoding/Decoding
// =========================================================

/*
 * Batch encode multiple lat/lng pairs to Morton codes using SIMD.
 * 
 * @param lats      Array of latitudes (count elements)
 * @param lngs      Array of longitudes (count elements)
 * @param out_z     Output array for Morton codes (count elements)
 * @param count     Number of points to encode
 */
void geo_simd_encode_batch(const double *lats, const double *lngs,
                           uint64_t *out_z, size_t count);

/*
 * Batch decode Morton codes to lat/lng pairs using SIMD.
 * 
 * @param z_codes   Array of Morton codes (count elements)
 * @param out_lats  Output array for latitudes (count elements)
 * @param out_lngs  Output array for longitudes (count elements)
 * @param count     Number of codes to decode
 */
void geo_simd_decode_batch(const uint64_t *z_codes,
                           double *out_lats, double *out_lngs, size_t count);

// =========================================================
// SIMD Function Declarations - Distance Calculations
// =========================================================

/*
 * Calculate Haversine distances from one point to multiple points using SIMD.
 * 
 * @param lat1      Latitude of reference point
 * @param lng1      Longitude of reference point
 * @param lats      Array of latitudes (count elements)
 * @param lngs      Array of longitudes (count elements)
 * @param out_dist  Output array for distances in km (count elements)
 * @param count     Number of distances to calculate
 */
void geo_simd_haversine_batch(double lat1, double lng1,
                              const double *lats, const double *lngs,
                              double *out_dist, size_t count);

/*
 * Fast approximate distances (equirectangular) using SIMD.
 * 
 * @param lat1      Latitude of reference point
 * @param lng1      Longitude of reference point
 * @param lats      Array of latitudes (count elements)
 * @param lngs      Array of longitudes (count elements)
 * @param out_dist  Output array for distances in km (count elements)
 * @param count     Number of distances to calculate
 */
void geo_simd_fast_distance_batch(double lat1, double lng1,
                                  const double *lats, const double *lngs,
                                  double *out_dist, size_t count);

// =========================================================
// SIMD Function Declarations - Range Filtering
// =========================================================

/*
 * Filter Morton codes within a range using SIMD.
 * Returns indices of codes within [z_min, z_max].
 * 
 * @param z_codes    Array of Morton codes (count elements)
 * @param count      Number of codes
 * @param z_min      Minimum Z value (inclusive)
 * @param z_max      Maximum Z value (inclusive)
 * @param out_mask   Output bitmask or indices
 * @return           Number of matches
 */
size_t geo_simd_filter_range(const uint64_t *z_codes, size_t count,
                             uint64_t z_min, uint64_t z_max,
                             uint8_t *out_mask);

/*
 * Filter points within bounding box using SIMD.
 * 
 * @param lats       Array of latitudes (count elements)
 * @param lngs       Array of longitudes (count elements)
 * @param count      Number of points
 * @param min_lat    Minimum latitude
 * @param max_lat    Maximum latitude
 * @param min_lng    Minimum longitude
 * @param max_lng    Maximum longitude
 * @param out_mask   Output bitmask (1 = inside, 0 = outside)
 * @return           Number of points inside
 */
size_t geo_simd_filter_bbox(const double *lats, const double *lngs, size_t count,
                            double min_lat, double max_lat,
                            double min_lng, double max_lng,
                            uint8_t *out_mask);

/*
 * Filter points within radius using SIMD.
 * Combines bounding box pre-filter with precise distance check.
 * 
 * @param lats       Array of latitudes (count elements)
 * @param lngs       Array of longitudes (count elements)
 * @param count      Number of points
 * @param center_lat Center latitude
 * @param center_lng Center longitude
 * @param radius_km  Radius in kilometers
 * @param out_mask   Output bitmask (1 = inside, 0 = outside)
 * @return           Number of points inside radius
 */
size_t geo_simd_filter_radius(const double *lats, const double *lngs, size_t count,
                              double center_lat, double center_lng, double radius_km,
                              uint8_t *out_mask);

// =========================================================
// SIMD Function Declarations - Morton Code Utilities
// =========================================================

/*
 * Batch spread bits for Morton encoding using SIMD.
 * Interleaves bits of 32-bit values into 64-bit values.
 * 
 * @param values     Array of 32-bit values (count elements)
 * @param out        Output array of 64-bit spread values (count elements)
 * @param count      Number of values
 */
void geo_simd_spread_bits_batch(const uint32_t *values, uint64_t *out, size_t count);

/*
 * Batch compact bits for Morton decoding using SIMD.
 * Extracts bits from 64-bit values to 32-bit values.
 * 
 * @param values     Array of 64-bit values (count elements)
 * @param out        Output array of 32-bit compact values (count elements)
 * @param count      Number of values
 */
void geo_simd_compact_bits_batch(const uint64_t *values, uint32_t *out, size_t count);

// =========================================================
// SIMD Utility Functions
// =========================================================

/*
 * Check if SIMD is available at runtime.
 */
bool geo_simd_available(void);

/*
 * Get SIMD implementation name.
 */
const char* geo_simd_get_name(void);

/*
 * Get optimal batch size for current SIMD implementation.
 */
size_t geo_simd_optimal_batch_size(void);

// =========================================================
// Benchmark Utilities
// =========================================================

typedef struct {
    double scalar_time_ms;
    double simd_time_ms;
    double speedup;
    size_t operations;
    const char *operation_name;
} GeoSimdBenchmark;

/*
 * Run a comparative benchmark between scalar and SIMD implementations.
 */
GeoSimdBenchmark geo_simd_benchmark_encode(size_t count, int iterations);
GeoSimdBenchmark geo_simd_benchmark_decode(size_t count, int iterations);
GeoSimdBenchmark geo_simd_benchmark_haversine(size_t count, int iterations);
GeoSimdBenchmark geo_simd_benchmark_filter_radius(size_t count, int iterations);

#ifdef __cplusplus
}
#endif

#endif // GEO_INDEX_SIMD_H
