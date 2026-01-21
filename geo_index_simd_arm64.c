/*
 * GeoIndex SIMD Optimizations - ARM64 NEON Implementation
 * 
 * Optimized implementations using ARM NEON intrinsics for:
 *   - Apple Silicon (M1/M2/M3)
 *   - ARM Cortex-A series
 *   - AWS Graviton
 */

#include "geo_index_simd.h"
#include "geo_index.h"
#include "geo_index_internal.h"

#if defined(GEO_SIMD_ARM64)

#include <arm_neon.h>
#include <stdlib.h>
#include <string.h>

// =========================================================
// SIMD Batch Encoding - ARM64 NEON
// =========================================================

void geo_simd_encode_batch(const double *lats, const double *lngs,
                           uint64_t *out_z, size_t count) {
    if (count == 0) return;
    
    float64x2_t lat_scale_v = vdupq_n_f64(GEO_INTERNAL_LAT_SCALE);
    float64x2_t lng_scale_v = vdupq_n_f64(GEO_INTERNAL_LNG_SCALE);
    float64x2_t lat_offset_v = vdupq_n_f64(GEO_INTERNAL_LAT_OFFSET);
    float64x2_t lng_offset_v = vdupq_n_f64(GEO_INTERNAL_LNG_OFFSET);
    float64x2_t zero_v = vdupq_n_f64(0.0);
    float64x2_t max_v = vdupq_n_f64(GEO_INTERNAL_COORD_MAX);
    
    size_t i = 0;
    
    // Process 2 points at a time (NEON f64 lanes)
    for (; i + 1 < count; i += 2) {
        // Load 2 lat/lng pairs
        float64x2_t lat_v = vld1q_f64(&lats[i]);
        float64x2_t lng_v = vld1q_f64(&lngs[i]);
        
        // Normalize: (coord + offset) * scale
        float64x2_t lat_norm = vmulq_f64(vaddq_f64(lat_v, lat_offset_v), lat_scale_v);
        float64x2_t lng_norm = vmulq_f64(vaddq_f64(lng_v, lng_offset_v), lng_scale_v);
        
        // Clamp to [0, 4294967295]
        lat_norm = vmaxq_f64(zero_v, vminq_f64(lat_norm, max_v));
        lng_norm = vmaxq_f64(zero_v, vminq_f64(lng_norm, max_v));
        
        // Convert to uint32 and spread bits
        uint32_t lat0 = (uint32_t)vgetq_lane_f64(lat_norm, 0);
        uint32_t lat1 = (uint32_t)vgetq_lane_f64(lat_norm, 1);
        uint32_t lng0 = (uint32_t)vgetq_lane_f64(lng_norm, 0);
        uint32_t lng1 = (uint32_t)vgetq_lane_f64(lng_norm, 1);
        
        out_z[i]     = geo_internal_spread_bits(lat0) | (geo_internal_spread_bits(lng0) << 1);
        out_z[i + 1] = geo_internal_spread_bits(lat1) | (geo_internal_spread_bits(lng1) << 1);
    }
    
    // Handle remaining point
    for (; i < count; i++) {
        out_z[i] = geo_encode(lats[i], lngs[i]);
    }
}

// =========================================================
// SIMD Batch Decoding - ARM64 NEON
// =========================================================

void geo_simd_decode_batch(const uint64_t *z_codes,
                           double *out_lats, double *out_lngs, size_t count) {
    if (count == 0) return;
    
    float64x2_t lat_scale_v = vdupq_n_f64(GEO_INTERNAL_LAT_DENORM_SCALE);
    float64x2_t lng_scale_v = vdupq_n_f64(GEO_INTERNAL_LNG_DENORM_SCALE);
    float64x2_t lat_offset_v = vdupq_n_f64(-90.0);
    float64x2_t lng_offset_v = vdupq_n_f64(-180.0);
    
    size_t i = 0;
    
    // Process 2 codes at a time
    for (; i + 1 < count; i += 2) {
        // Compact bits (scalar, then vectorize conversion)
        uint32_t lat0 = geo_internal_compact_bits(z_codes[i]);
        uint32_t lng0 = geo_internal_compact_bits(z_codes[i] >> 1);
        uint32_t lat1 = geo_internal_compact_bits(z_codes[i + 1]);
        uint32_t lng1 = geo_internal_compact_bits(z_codes[i + 1] >> 1);
        
        // Convert to double
        float64x2_t lat_v = {(double)lat0, (double)lat1};
        float64x2_t lng_v = {(double)lng0, (double)lng1};
        
        // Denormalize
        lat_v = vaddq_f64(vmulq_f64(lat_v, lat_scale_v), lat_offset_v);
        lng_v = vaddq_f64(vmulq_f64(lng_v, lng_scale_v), lng_offset_v);
        
        // Store
        vst1q_f64(&out_lats[i], lat_v);
        vst1q_f64(&out_lngs[i], lng_v);
    }
    
    // Handle remaining
    for (; i < count; i++) {
        GeoPoint p = geo_decode(z_codes[i]);
        out_lats[i] = p.lat;
        out_lngs[i] = p.lng;
    }
}

// =========================================================
// SIMD Haversine Distance - ARM64 NEON
// =========================================================

// Fast sine approximation for NEON (Bhaskara I's approximation enhanced)
static inline float64x2_t neon_sin_approx(float64x2_t x) {
    // Reduce x to [-π, π] range
    float64x2_t pi = vdupq_n_f64(M_PI);
    float64x2_t two_pi = vdupq_n_f64(2.0 * M_PI);
    
    // Normalize to [-π, π]
    x = vsubq_f64(x, vmulq_f64(vrndnq_f64(vmulq_f64(x, vdupq_n_f64(1.0 / (2.0 * M_PI)))), two_pi));
    
    // Parabolic approximation: sin(x) ≈ 16x(π - x) / (5π² - 4x(π - x))
    float64x2_t x_sub_pi = vsubq_f64(pi, vabsq_f64(x));
    float64x2_t numerator = vmulq_f64(vmulq_f64(vdupq_n_f64(16.0), vabsq_f64(x)), x_sub_pi);
    float64x2_t inner = vmulq_f64(vmulq_f64(vdupq_n_f64(4.0), vabsq_f64(x)), x_sub_pi);
    float64x2_t denominator = vsubq_f64(vdupq_n_f64(5.0 * M_PI * M_PI), inner);
    
    float64x2_t result = vdivq_f64(numerator, denominator);
    
    // Apply sign
    uint64x2_t sign_mask = vcltq_f64(x, vdupq_n_f64(0.0));
    return vbslq_f64(sign_mask, vnegq_f64(result), result);
}

// Fast cosine via sin(x + π/2)
static inline float64x2_t neon_cos_approx(float64x2_t x) {
    return neon_sin_approx(vaddq_f64(x, vdupq_n_f64(M_PI / 2.0)));
}

void geo_simd_haversine_batch(double lat1, double lng1,
                              const double *lats, const double *lngs,
                              double *out_dist, size_t count) {
    if (count == 0) return;
    
    float64x2_t deg_to_rad = vdupq_n_f64(GEO_INTERNAL_DEG_TO_RAD);
    float64x2_t earth_radius = vdupq_n_f64(GEO_INTERNAL_EARTH_RADIUS_KM);
    float64x2_t lat1_rad_v = vdupq_n_f64(lat1 * GEO_INTERNAL_DEG_TO_RAD);
    float64x2_t lng1_rad_v = vdupq_n_f64(lng1 * GEO_INTERNAL_DEG_TO_RAD);
    float64x2_t half = vdupq_n_f64(0.5);
    float64x2_t two = vdupq_n_f64(2.0);
    float64x2_t one = vdupq_n_f64(1.0);
    
    // Precompute cos(lat1)
    float64x2_t cos_lat1 = neon_cos_approx(lat1_rad_v);
    
    size_t i = 0;
    
    // Process 2 points at a time
    for (; i + 1 < count; i += 2) {
        float64x2_t lat2_v = vld1q_f64(&lats[i]);
        float64x2_t lng2_v = vld1q_f64(&lngs[i]);
        
        // Convert to radians
        float64x2_t lat2_rad = vmulq_f64(lat2_v, deg_to_rad);
        float64x2_t lng2_rad = vmulq_f64(lng2_v, deg_to_rad);
        
        // Delta lat/lng
        float64x2_t dlat = vsubq_f64(lat2_rad, lat1_rad_v);
        float64x2_t dlng = vsubq_f64(lng2_rad, lng1_rad_v);
        
        // sin(dlat/2)^2
        float64x2_t sin_dlat_2 = neon_sin_approx(vmulq_f64(dlat, half));
        float64x2_t sin_dlat_2_sq = vmulq_f64(sin_dlat_2, sin_dlat_2);
        
        // sin(dlng/2)^2
        float64x2_t sin_dlng_2 = neon_sin_approx(vmulq_f64(dlng, half));
        float64x2_t sin_dlng_2_sq = vmulq_f64(sin_dlng_2, sin_dlng_2);
        
        // cos(lat2)
        float64x2_t cos_lat2 = neon_cos_approx(lat2_rad);
        
        // a = sin²(Δlat/2) + cos(lat1)·cos(lat2)·sin²(Δlng/2)
        float64x2_t a = vaddq_f64(sin_dlat_2_sq,
                                   vmulq_f64(vmulq_f64(cos_lat1, cos_lat2), sin_dlng_2_sq));
        
        // c = 2·atan2(√a, √(1-a)) ≈ 2·asin(√a) for small angles
        // Using approximation: asin(x) ≈ x + x³/6 for small x
        float64x2_t sqrt_a = vsqrtq_f64(a);
        float64x2_t sqrt_1_minus_a = vsqrtq_f64(vsubq_f64(one, a));
        
        // atan2 approximation
        float64x2_t c = vmulq_f64(two, vmulq_f64(sqrt_a, 
                        vdivq_f64(one, vaddq_f64(sqrt_1_minus_a, vdupq_n_f64(0.00001)))));
        
        // More accurate: c = 2 * asin(sqrt(a))
        // asin(x) ≈ x(1 + x²(1/6 + x²(3/40)))
        float64x2_t x2 = vmulq_f64(sqrt_a, sqrt_a);
        float64x2_t asin_approx = vmulq_f64(sqrt_a,
                                   vaddq_f64(one,
                                   vmulq_f64(x2,
                                   vaddq_f64(vdupq_n_f64(1.0/6.0),
                                   vmulq_f64(x2, vdupq_n_f64(3.0/40.0))))));
        c = vmulq_f64(two, asin_approx);
        
        // Distance = R * c
        float64x2_t dist = vmulq_f64(earth_radius, c);
        
        vst1q_f64(&out_dist[i], dist);
    }
    
    // Handle remaining with scalar
    for (; i < count; i++) {
        out_dist[i] = geo_haversine_km(lat1, lng1, lats[i], lngs[i]);
    }
}

// =========================================================
// SIMD Fast Distance (Equirectangular) - ARM64 NEON
// =========================================================

void geo_simd_fast_distance_batch(double lat1, double lng1,
                                  const double *lats, const double *lngs,
                                  double *out_dist, size_t count) {
    if (count == 0) return;
    
    float64x2_t km_per_deg = vdupq_n_f64(GEO_INTERNAL_KM_PER_DEG);
    float64x2_t lat1_v = vdupq_n_f64(lat1);
    float64x2_t lng1_v = vdupq_n_f64(lng1);
    float64x2_t deg_to_rad = vdupq_n_f64(GEO_INTERNAL_DEG_TO_RAD);
    float64x2_t half = vdupq_n_f64(0.5);
    
    size_t i = 0;
    
    // Process 2 points at a time
    for (; i + 1 < count; i += 2) {
        float64x2_t lat2_v = vld1q_f64(&lats[i]);
        float64x2_t lng2_v = vld1q_f64(&lngs[i]);
        
        // Mid latitude for cos adjustment
        float64x2_t lat_mid = vmulq_f64(vaddq_f64(lat1_v, lat2_v), half);
        float64x2_t lat_mid_rad = vmulq_f64(lat_mid, deg_to_rad);
        float64x2_t cos_lat = neon_cos_approx(lat_mid_rad);
        
        // Delta calculations
        float64x2_t dlat = vmulq_f64(vsubq_f64(lat2_v, lat1_v), km_per_deg);
        float64x2_t dlng = vmulq_f64(vmulq_f64(vsubq_f64(lng2_v, lng1_v), km_per_deg), cos_lat);
        
        // Distance = sqrt(dlat² + dlng²)
        float64x2_t dist_sq = vaddq_f64(vmulq_f64(dlat, dlat), vmulq_f64(dlng, dlng));
        float64x2_t dist = vsqrtq_f64(dist_sq);
        
        vst1q_f64(&out_dist[i], dist);
    }
    
    // Handle remaining
    for (; i < count; i++) {
        out_dist[i] = geo_fast_distance_km(lat1, lng1, lats[i], lngs[i]);
    }
}

// =========================================================
// SIMD Range Filtering - ARM64 NEON
// =========================================================

size_t geo_simd_filter_range(const uint64_t *z_codes, size_t count,
                             uint64_t z_min, uint64_t z_max,
                             uint8_t *out_mask) {
    if (count == 0) return 0;
    
    uint64x2_t min_v = vdupq_n_u64(z_min);
    uint64x2_t max_v = vdupq_n_u64(z_max);
    
    size_t matches = 0;
    size_t i = 0;
    
    // Process 2 codes at a time
    for (; i + 1 < count; i += 2) {
        uint64x2_t z_v = vld1q_u64(&z_codes[i]);
        
        // Check z >= z_min AND z <= z_max
        uint64x2_t ge_min = vcgeq_u64(z_v, min_v);
        uint64x2_t le_max = vcleq_u64(z_v, max_v);
        uint64x2_t in_range = vandq_u64(ge_min, le_max);
        
        // Extract results
        uint64_t mask0 = vgetq_lane_u64(in_range, 0);
        uint64_t mask1 = vgetq_lane_u64(in_range, 1);
        
        out_mask[i]     = mask0 ? 1 : 0;
        out_mask[i + 1] = mask1 ? 1 : 0;
        
        matches += (mask0 ? 1 : 0) + (mask1 ? 1 : 0);
    }
    
    // Handle remaining
    for (; i < count; i++) {
        uint8_t in = (z_codes[i] >= z_min && z_codes[i] <= z_max) ? 1 : 0;
        out_mask[i] = in;
        matches += in;
    }
    
    return matches;
}

// =========================================================
// SIMD Bounding Box Filtering - ARM64 NEON
// =========================================================

size_t geo_simd_filter_bbox(const double *lats, const double *lngs, size_t count,
                            double min_lat, double max_lat,
                            double min_lng, double max_lng,
                            uint8_t *out_mask) {
    if (count == 0) return 0;
    
    float64x2_t min_lat_v = vdupq_n_f64(min_lat);
    float64x2_t max_lat_v = vdupq_n_f64(max_lat);
    float64x2_t min_lng_v = vdupq_n_f64(min_lng);
    float64x2_t max_lng_v = vdupq_n_f64(max_lng);
    
    size_t matches = 0;
    size_t i = 0;
    
    // Process 2 points at a time
    for (; i + 1 < count; i += 2) {
        float64x2_t lat_v = vld1q_f64(&lats[i]);
        float64x2_t lng_v = vld1q_f64(&lngs[i]);
        
        // Check lat in range
        uint64x2_t lat_ge_min = vcgeq_f64(lat_v, min_lat_v);
        uint64x2_t lat_le_max = vcleq_f64(lat_v, max_lat_v);
        uint64x2_t lat_ok = vandq_u64(lat_ge_min, lat_le_max);
        
        // Check lng in range
        uint64x2_t lng_ge_min = vcgeq_f64(lng_v, min_lng_v);
        uint64x2_t lng_le_max = vcleq_f64(lng_v, max_lng_v);
        uint64x2_t lng_ok = vandq_u64(lng_ge_min, lng_le_max);
        
        // Combined
        uint64x2_t in_bbox = vandq_u64(lat_ok, lng_ok);
        
        uint64_t mask0 = vgetq_lane_u64(in_bbox, 0);
        uint64_t mask1 = vgetq_lane_u64(in_bbox, 1);
        
        out_mask[i]     = mask0 ? 1 : 0;
        out_mask[i + 1] = mask1 ? 1 : 0;
        
        matches += (mask0 ? 1 : 0) + (mask1 ? 1 : 0);
    }
    
    // Handle remaining
    for (; i < count; i++) {
        uint8_t in = (lats[i] >= min_lat && lats[i] <= max_lat &&
                      lngs[i] >= min_lng && lngs[i] <= max_lng) ? 1 : 0;
        out_mask[i] = in;
        matches += in;
    }
    
    return matches;
}

// =========================================================
// SIMD Radius Filtering - ARM64 NEON
// =========================================================

size_t geo_simd_filter_radius(const double *lats, const double *lngs, size_t count,
                              double center_lat, double center_lng, double radius_km,
                              uint8_t *out_mask) {
    if (count == 0) return 0;
    
    // Pre-filter with bounding box
    double min_lat, max_lat, min_lng, max_lng;
    geo_bounding_box(center_lat, center_lng, radius_km,
                     &min_lat, &max_lat, &min_lng, &max_lng);
    
    // Allocate temporary arrays for candidates
    double *candidate_lats = (double*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(double));
    double *candidate_lngs = (double*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(double));
    size_t *candidate_idx = (size_t*)malloc(count * sizeof(size_t));
    
    if (!candidate_lats || !candidate_lngs || !candidate_idx) {
        free(candidate_lats);
        free(candidate_lngs);
        free(candidate_idx);
        // Fallback to scalar
        size_t matches = 0;
        for (size_t i = 0; i < count; i++) {
            double dist = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);
            out_mask[i] = (dist <= radius_km) ? 1 : 0;
            matches += out_mask[i];
        }
        return matches;
    }
    
    // Initialize mask to 0
    memset(out_mask, 0, count);
    
    // First pass: bounding box filter
    size_t num_candidates = 0;
    size_t i = 0;
    
    float64x2_t min_lat_v = vdupq_n_f64(min_lat);
    float64x2_t max_lat_v = vdupq_n_f64(max_lat);
    float64x2_t min_lng_v = vdupq_n_f64(min_lng);
    float64x2_t max_lng_v = vdupq_n_f64(max_lng);
    
    for (; i + 1 < count; i += 2) {
        float64x2_t lat_v = vld1q_f64(&lats[i]);
        float64x2_t lng_v = vld1q_f64(&lngs[i]);
        
        uint64x2_t lat_ok = vandq_u64(vcgeq_f64(lat_v, min_lat_v), vcleq_f64(lat_v, max_lat_v));
        uint64x2_t lng_ok = vandq_u64(vcgeq_f64(lng_v, min_lng_v), vcleq_f64(lng_v, max_lng_v));
        uint64x2_t in_bbox = vandq_u64(lat_ok, lng_ok);
        
        if (vgetq_lane_u64(in_bbox, 0)) {
            candidate_lats[num_candidates] = lats[i];
            candidate_lngs[num_candidates] = lngs[i];
            candidate_idx[num_candidates] = i;
            num_candidates++;
        }
        if (vgetq_lane_u64(in_bbox, 1)) {
            candidate_lats[num_candidates] = lats[i + 1];
            candidate_lngs[num_candidates] = lngs[i + 1];
            candidate_idx[num_candidates] = i + 1;
            num_candidates++;
        }
    }
    
    for (; i < count; i++) {
        if (lats[i] >= min_lat && lats[i] <= max_lat &&
            lngs[i] >= min_lng && lngs[i] <= max_lng) {
            candidate_lats[num_candidates] = lats[i];
            candidate_lngs[num_candidates] = lngs[i];
            candidate_idx[num_candidates] = i;
            num_candidates++;
        }
    }
    
    // Second pass: precise distance check on candidates
    double *distances = (double*)aligned_alloc(GEO_SIMD_ALIGN, num_candidates * sizeof(double));
    if (distances) {
        geo_simd_haversine_batch(center_lat, center_lng,
                                  candidate_lats, candidate_lngs,
                                  distances, num_candidates);
        
        float64x2_t radius_v = vdupq_n_f64(radius_km);
        size_t j = 0;
        
        for (; j + 1 < num_candidates; j += 2) {
            float64x2_t dist_v = vld1q_f64(&distances[j]);
            uint64x2_t in_radius = vcleq_f64(dist_v, radius_v);
            
            if (vgetq_lane_u64(in_radius, 0)) out_mask[candidate_idx[j]] = 1;
            if (vgetq_lane_u64(in_radius, 1)) out_mask[candidate_idx[j + 1]] = 1;
        }
        
        for (; j < num_candidates; j++) {
            if (distances[j] <= radius_km) {
                out_mask[candidate_idx[j]] = 1;
            }
        }
        
        free(distances);
    }
    
    free(candidate_lats);
    free(candidate_lngs);
    free(candidate_idx);
    
    // Count matches
    size_t matches = 0;
    for (i = 0; i < count; i++) {
        matches += out_mask[i];
    }
    
    return matches;
}

// =========================================================
// SIMD Bit Manipulation - ARM64 NEON
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
// Utility Functions - ARM64
// =========================================================

bool geo_simd_available(void) {
    return true; // NEON is always available on ARM64
}

const char* geo_simd_get_name(void) {
    return GEO_SIMD_NAME;
}

size_t geo_simd_optimal_batch_size(void) {
    return 256; // Good balance for NEON
}

// =========================================================
// Benchmark Functions - ARM64
// =========================================================

GeoSimdBenchmark geo_simd_benchmark_encode(size_t count, int iterations) {
    GeoSimdBenchmark result = {0};
    result.operation_name = "encode";
    result.operations = count * iterations;
    
    double *lats = (double*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(double));
    double *lngs = (double*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(double));
    uint64_t *z_codes = (uint64_t*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(uint64_t));
    
    if (!lats || !lngs || !z_codes) {
        free(lats); free(lngs); free(z_codes);
        return result;
    }
    
    // Initialize test data
    for (size_t i = 0; i < count; i++) {
        lats[i] = ((double)(i % 18000) / 100.0) - 90.0;
        lngs[i] = ((double)(i % 36000) / 100.0) - 180.0;
    }
    
    // Scalar benchmark
    double start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            z_codes[i] = geo_encode(lats[i], lngs[i]);
        }
    }
    result.scalar_time_ms = geo_get_time_ms() - start;
    
    // SIMD benchmark
    start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        geo_simd_encode_batch(lats, lngs, z_codes, count);
    }
    result.simd_time_ms = geo_get_time_ms() - start;
    
    result.speedup = result.scalar_time_ms / result.simd_time_ms;
    
    free(lats); free(lngs); free(z_codes);
    return result;
}

GeoSimdBenchmark geo_simd_benchmark_decode(size_t count, int iterations) {
    GeoSimdBenchmark result = {0};
    result.operation_name = "decode";
    result.operations = count * iterations;
    
    uint64_t *z_codes = (uint64_t*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(uint64_t));
    double *lats = (double*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(double));
    double *lngs = (double*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(double));
    
    if (!z_codes || !lats || !lngs) {
        free(z_codes); free(lats); free(lngs);
        return result;
    }
    
    // Initialize test data
    for (size_t i = 0; i < count; i++) {
        z_codes[i] = (uint64_t)i * 12345ULL;
    }
    
    // Scalar benchmark
    double start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            GeoPoint p = geo_decode(z_codes[i]);
            lats[i] = p.lat;
            lngs[i] = p.lng;
        }
    }
    result.scalar_time_ms = geo_get_time_ms() - start;
    
    // SIMD benchmark
    start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        geo_simd_decode_batch(z_codes, lats, lngs, count);
    }
    result.simd_time_ms = geo_get_time_ms() - start;
    
    result.speedup = result.scalar_time_ms / result.simd_time_ms;
    
    free(z_codes); free(lats); free(lngs);
    return result;
}

GeoSimdBenchmark geo_simd_benchmark_haversine(size_t count, int iterations) {
    GeoSimdBenchmark result = {0};
    result.operation_name = "haversine";
    result.operations = count * iterations;
    
    double *lats = (double*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(double));
    double *lngs = (double*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(double));
    double *dists = (double*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(double));
    
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
    
    // Scalar benchmark
    double start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            dists[i] = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);
        }
    }
    result.scalar_time_ms = geo_get_time_ms() - start;
    
    // SIMD benchmark
    start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        geo_simd_haversine_batch(center_lat, center_lng, lats, lngs, dists, count);
    }
    result.simd_time_ms = geo_get_time_ms() - start;
    
    result.speedup = result.scalar_time_ms / result.simd_time_ms;
    
    free(lats); free(lngs); free(dists);
    return result;
}

GeoSimdBenchmark geo_simd_benchmark_filter_radius(size_t count, int iterations) {
    GeoSimdBenchmark result = {0};
    result.operation_name = "filter_radius";
    result.operations = count * iterations;
    
    double *lats = (double*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(double));
    double *lngs = (double*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(double));
    uint8_t *mask = (uint8_t*)aligned_alloc(GEO_SIMD_ALIGN, count * sizeof(uint8_t));
    
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
    
    // Scalar benchmark
    double start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            double dist = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);
            mask[i] = (dist <= radius_km) ? 1 : 0;
        }
    }
    result.scalar_time_ms = geo_get_time_ms() - start;
    
    // SIMD benchmark
    start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        geo_simd_filter_radius(lats, lngs, count, center_lat, center_lng, radius_km, mask);
    }
    result.simd_time_ms = geo_get_time_ms() - start;
    
    result.speedup = result.scalar_time_ms / result.simd_time_ms;
    
    free(lats); free(lngs); free(mask);
    return result;
}

#endif // GEO_SIMD_ARM64
