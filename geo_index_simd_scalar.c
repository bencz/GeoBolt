/*
 * GeoIndex SIMD Optimizations - Scalar Fallback Implementation
 * 
 * This file provides scalar implementations of the SIMD functions
 * for platforms where SIMD is not available or not detected.
 */

#include "geo_index_simd.h"
#include "geo_index.h"
#include "geo_index_internal.h"

#if !defined(GEO_SIMD_ARM64) && !defined(GEO_SIMD_X86_AVX2) && !defined(GEO_SIMD_X86_SSE4) && !defined(GEO_SIMD_X86_SSE2)

#include <stdlib.h>
#include <string.h>

// =========================================================
// Scalar Batch Encoding
// =========================================================

void geo_simd_encode_batch(const double *lats, const double *lngs,
                           uint64_t *out_z, size_t count) {
    for (size_t i = 0; i < count; i++) {
        out_z[i] = geo_encode(lats[i], lngs[i]);
    }
}

// =========================================================
// Scalar Batch Decoding
// =========================================================

void geo_simd_decode_batch(const uint64_t *z_codes,
                           double *out_lats, double *out_lngs, size_t count) {
    for (size_t i = 0; i < count; i++) {
        GeoPoint p = geo_decode(z_codes[i]);
        out_lats[i] = p.lat;
        out_lngs[i] = p.lng;
    }
}

// =========================================================
// Scalar Haversine Distance
// =========================================================

void geo_simd_haversine_batch(double lat1, double lng1,
                              const double *lats, const double *lngs,
                              double *out_dist, size_t count) {
    for (size_t i = 0; i < count; i++) {
        out_dist[i] = geo_haversine_km(lat1, lng1, lats[i], lngs[i]);
    }
}

// =========================================================
// Scalar Fast Distance
// =========================================================

void geo_simd_fast_distance_batch(double lat1, double lng1,
                                  const double *lats, const double *lngs,
                                  double *out_dist, size_t count) {
    for (size_t i = 0; i < count; i++) {
        out_dist[i] = geo_fast_distance_km(lat1, lng1, lats[i], lngs[i]);
    }
}

// =========================================================
// Scalar Range Filtering
// =========================================================

size_t geo_simd_filter_range(const uint64_t *z_codes, size_t count,
                             uint64_t z_min, uint64_t z_max,
                             uint8_t *out_mask) {
    size_t matches = 0;
    for (size_t i = 0; i < count; i++) {
        uint8_t in = (z_codes[i] >= z_min && z_codes[i] <= z_max) ? 1 : 0;
        out_mask[i] = in;
        matches += in;
    }
    return matches;
}

// =========================================================
// Scalar Bounding Box Filtering
// =========================================================

size_t geo_simd_filter_bbox(const double *lats, const double *lngs, size_t count,
                            double min_lat, double max_lat,
                            double min_lng, double max_lng,
                            uint8_t *out_mask) {
    size_t matches = 0;
    for (size_t i = 0; i < count; i++) {
        uint8_t in = (lats[i] >= min_lat && lats[i] <= max_lat &&
                      lngs[i] >= min_lng && lngs[i] <= max_lng) ? 1 : 0;
        out_mask[i] = in;
        matches += in;
    }
    return matches;
}

// =========================================================
// Scalar Radius Filtering
// =========================================================

size_t geo_simd_filter_radius(const double *lats, const double *lngs, size_t count,
                              double center_lat, double center_lng, double radius_km,
                              uint8_t *out_mask) {
    size_t matches = 0;
    for (size_t i = 0; i < count; i++) {
        double dist = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);
        uint8_t in = (dist <= radius_km) ? 1 : 0;
        out_mask[i] = in;
        matches += in;
    }
    return matches;
}

// =========================================================
// Scalar Bit Manipulation
// =========================================================

void geo_simd_spread_bits_batch(const uint32_t *values, uint64_t *out, size_t count) {
    for (size_t i = 0; i < count; i++) {
        out[i] = geo_internal_spread_bits(values[i]);
    }
}

void geo_simd_compact_bits_batch(const uint64_t *values, uint32_t *out, size_t count) {
    for (size_t i = 0; i < count; i++) {
        out[i] = geo_internal_compact_bits(values[i]);
    }
}

// =========================================================
// Utility Functions - Scalar
// =========================================================

bool geo_simd_available(void) {
    return false;
}

const char* geo_simd_get_name(void) {
    return "Scalar (no SIMD)";
}

size_t geo_simd_optimal_batch_size(void) {
    return 64;
}

// =========================================================
// Benchmark Functions - Scalar
// =========================================================

GeoSimdBenchmark geo_simd_benchmark_encode(size_t count, int iterations) {
    GeoSimdBenchmark result = {0};
    result.operation_name = "encode";
    result.operations = count * iterations;
    
    double *lats = (double*)malloc(count * sizeof(double));
    double *lngs = (double*)malloc(count * sizeof(double));
    uint64_t *z_codes = (uint64_t*)malloc(count * sizeof(uint64_t));
    
    if (!lats || !lngs || !z_codes) {
        free(lats); free(lngs); free(z_codes);
        return result;
    }
    
    for (size_t i = 0; i < count; i++) {
        lats[i] = ((double)(i % 18000) / 100.0) - 90.0;
        lngs[i] = ((double)(i % 36000) / 100.0) - 180.0;
    }
    
    double start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            z_codes[i] = geo_encode(lats[i], lngs[i]);
        }
    }
    result.scalar_time_ms = geo_get_time_ms() - start;
    result.simd_time_ms = result.scalar_time_ms; // Same as scalar
    result.speedup = 1.0;
    
    free(lats); free(lngs); free(z_codes);
    return result;
}

GeoSimdBenchmark geo_simd_benchmark_decode(size_t count, int iterations) {
    GeoSimdBenchmark result = {0};
    result.operation_name = "decode";
    result.operations = count * iterations;
    
    uint64_t *z_codes = (uint64_t*)malloc(count * sizeof(uint64_t));
    double *lats = (double*)malloc(count * sizeof(double));
    double *lngs = (double*)malloc(count * sizeof(double));
    
    if (!z_codes || !lats || !lngs) {
        free(z_codes); free(lats); free(lngs);
        return result;
    }
    
    for (size_t i = 0; i < count; i++) {
        z_codes[i] = (uint64_t)i * 12345ULL;
    }
    
    double start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            GeoPoint p = geo_decode(z_codes[i]);
            lats[i] = p.lat;
            lngs[i] = p.lng;
        }
    }
    result.scalar_time_ms = geo_get_time_ms() - start;
    result.simd_time_ms = result.scalar_time_ms;
    result.speedup = 1.0;
    
    free(z_codes); free(lats); free(lngs);
    return result;
}

GeoSimdBenchmark geo_simd_benchmark_haversine(size_t count, int iterations) {
    GeoSimdBenchmark result = {0};
    result.operation_name = "haversine";
    result.operations = count * iterations;
    
    double *lats = (double*)malloc(count * sizeof(double));
    double *lngs = (double*)malloc(count * sizeof(double));
    double *dists = (double*)malloc(count * sizeof(double));
    
    if (!lats || !lngs || !dists) {
        free(lats); free(lngs); free(dists);
        return result;
    }
    
    double center_lat = -23.5505;
    double center_lng = -46.6333;
    
    for (size_t i = 0; i < count; i++) {
        lats[i] = center_lat + ((double)(i % 1000) / 1000.0 - 0.5) * 10.0;
        lngs[i] = center_lng + ((double)(i % 1000) / 1000.0 - 0.5) * 10.0;
    }
    
    double start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            dists[i] = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);
        }
    }
    result.scalar_time_ms = geo_get_time_ms() - start;
    result.simd_time_ms = result.scalar_time_ms;
    result.speedup = 1.0;
    
    free(lats); free(lngs); free(dists);
    return result;
}

GeoSimdBenchmark geo_simd_benchmark_filter_radius(size_t count, int iterations) {
    GeoSimdBenchmark result = {0};
    result.operation_name = "filter_radius";
    result.operations = count * iterations;
    
    double *lats = (double*)malloc(count * sizeof(double));
    double *lngs = (double*)malloc(count * sizeof(double));
    uint8_t *mask = (uint8_t*)malloc(count * sizeof(uint8_t));
    
    if (!lats || !lngs || !mask) {
        free(lats); free(lngs); free(mask);
        return result;
    }
    
    double center_lat = -23.5505;
    double center_lng = -46.6333;
    double radius_km = 50.0;
    
    for (size_t i = 0; i < count; i++) {
        lats[i] = ((double)(i % 18000) / 100.0) - 90.0;
        lngs[i] = ((double)(i % 36000) / 100.0) - 180.0;
    }
    
    double start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            double dist = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);
            mask[i] = (dist <= radius_km) ? 1 : 0;
        }
    }
    result.scalar_time_ms = geo_get_time_ms() - start;
    result.simd_time_ms = result.scalar_time_ms;
    result.speedup = 1.0;
    
    free(lats); free(lngs); free(mask);
    return result;
}

#endif // No SIMD
