/*
 * GeoIndex SIMD Optimizations - x86-64 AVX2/SSE Implementation
 * 
 * Optimized implementations using x86 SIMD intrinsics for:
 *   - Intel Core series (Haswell and newer for AVX2)
 *   - AMD Ryzen series
 *   - SSE4 fallback for older processors
 */

#include "geo_index_simd.h"
#include "geo_index.h"
#include "geo_index_internal.h"

#if defined(GEO_SIMD_X86_AVX2) || defined(GEO_SIMD_X86_SSE4) || defined(GEO_SIMD_X86_SSE2)

#include <immintrin.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#include <intrin.h>
#else
#include <cpuid.h>
#endif

// =========================================================
// Runtime CPU Feature Detection
// =========================================================

static int cpu_has_avx2 = -1;
static int cpu_has_avx = -1;
static int cpu_has_sse4 = -1;

static void detect_cpu_features(void) {
    if (cpu_has_avx2 >= 0) return; // Already detected
    
#ifdef _MSC_VER
    int cpuinfo[4];
    __cpuid(cpuinfo, 1);
    cpu_has_sse4 = (cpuinfo[2] & (1 << 19)) != 0; // SSE4.1
    cpu_has_avx = (cpuinfo[2] & (1 << 28)) != 0;  // AVX
    
    if (cpu_has_avx) {
        __cpuidex(cpuinfo, 7, 0);
        cpu_has_avx2 = (cpuinfo[1] & (1 << 5)) != 0; // AVX2
    } else {
        cpu_has_avx2 = 0;
    }
#else
    unsigned int eax, ebx, ecx, edx;
    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        cpu_has_sse4 = (ecx & (1 << 19)) != 0;
        cpu_has_avx = (ecx & (1 << 28)) != 0;
    } else {
        cpu_has_sse4 = 0;
        cpu_has_avx = 0;
    }
    
    if (cpu_has_avx && __get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
        cpu_has_avx2 = (ebx & (1 << 5)) != 0;
    } else {
        cpu_has_avx2 = 0;
    }
#endif
}

// =========================================================
// SIMD Batch Encoding - x86 AVX2
// =========================================================

#if defined(GEO_SIMD_X86_AVX2)

void geo_simd_encode_batch(const double *lats, const double *lngs,
                           uint64_t *out_z, size_t count) {
    if (count == 0) return;
    
    __m256d lat_scale_v = _mm256_set1_pd(GEO_INTERNAL_LAT_SCALE);
    __m256d lng_scale_v = _mm256_set1_pd(GEO_INTERNAL_LNG_SCALE);
    __m256d lat_offset_v = _mm256_set1_pd(GEO_INTERNAL_LAT_OFFSET);
    __m256d lng_offset_v = _mm256_set1_pd(GEO_INTERNAL_LNG_OFFSET);
    __m256d zero_v = _mm256_setzero_pd();
    __m256d max_v = _mm256_set1_pd(GEO_INTERNAL_COORD_MAX);
    
    size_t i = 0;
    
    // Process 4 points at a time (AVX2 4x double)
    for (; i + 3 < count; i += 4) {
        // Load 4 lat/lng pairs
        __m256d lat_v = _mm256_loadu_pd(&lats[i]);
        __m256d lng_v = _mm256_loadu_pd(&lngs[i]);
        
        // Normalize: (coord + offset) * scale
        __m256d lat_norm = _mm256_mul_pd(_mm256_add_pd(lat_v, lat_offset_v), lat_scale_v);
        __m256d lng_norm = _mm256_mul_pd(_mm256_add_pd(lng_v, lng_offset_v), lng_scale_v);
        
        // Clamp to [0, 4294967295]
        lat_norm = _mm256_max_pd(zero_v, _mm256_min_pd(lat_norm, max_v));
        lng_norm = _mm256_max_pd(zero_v, _mm256_min_pd(lng_norm, max_v));
        
        // Extract to scalar and process Morton encoding
        double lat_arr[4], lng_arr[4];
        _mm256_storeu_pd(lat_arr, lat_norm);
        _mm256_storeu_pd(lng_arr, lng_norm);
        
        for (int j = 0; j < 4; j++) {
            uint32_t lat_u = (uint32_t)lat_arr[j];
            uint32_t lng_u = (uint32_t)lng_arr[j];
            out_z[i + j] = geo_internal_spread_bits(lat_u) | (geo_internal_spread_bits(lng_u) << 1);
        }
    }
    
    // Handle remaining points
    for (; i < count; i++) {
        out_z[i] = geo_encode(lats[i], lngs[i]);
    }
}

#elif defined(GEO_SIMD_X86_SSE4) || defined(GEO_SIMD_X86_SSE2)

void geo_simd_encode_batch(const double *lats, const double *lngs,
                           uint64_t *out_z, size_t count) {
    if (count == 0) return;
    
    __m128d lat_scale_v = _mm_set1_pd(GEO_INTERNAL_LAT_SCALE);
    __m128d lng_scale_v = _mm_set1_pd(GEO_INTERNAL_LNG_SCALE);
    __m128d lat_offset_v = _mm_set1_pd(GEO_INTERNAL_LAT_OFFSET);
    __m128d lng_offset_v = _mm_set1_pd(GEO_INTERNAL_LNG_OFFSET);
    __m128d zero_v = _mm_setzero_pd();
    __m128d max_v = _mm_set1_pd(GEO_INTERNAL_COORD_MAX);
    
    size_t i = 0;
    
    // Process 2 points at a time (SSE 2x double)
    for (; i + 1 < count; i += 2) {
        __m128d lat_v = _mm_loadu_pd(&lats[i]);
        __m128d lng_v = _mm_loadu_pd(&lngs[i]);
        
        __m128d lat_norm = _mm_mul_pd(_mm_add_pd(lat_v, lat_offset_v), lat_scale_v);
        __m128d lng_norm = _mm_mul_pd(_mm_add_pd(lng_v, lng_offset_v), lng_scale_v);
        
        lat_norm = _mm_max_pd(zero_v, _mm_min_pd(lat_norm, max_v));
        lng_norm = _mm_max_pd(zero_v, _mm_min_pd(lng_norm, max_v));
        
        double lat_arr[2], lng_arr[2];
        _mm_storeu_pd(lat_arr, lat_norm);
        _mm_storeu_pd(lng_arr, lng_norm);
        
        for (int j = 0; j < 2; j++) {
            uint32_t lat_u = (uint32_t)lat_arr[j];
            uint32_t lng_u = (uint32_t)lng_arr[j];
            out_z[i + j] = geo_internal_spread_bits(lat_u) | (geo_internal_spread_bits(lng_u) << 1);
        }
    }
    
    for (; i < count; i++) {
        out_z[i] = geo_encode(lats[i], lngs[i]);
    }
}

#endif

// =========================================================
// SIMD Batch Decoding - x86
// =========================================================

#if defined(GEO_SIMD_X86_AVX2)

void geo_simd_decode_batch(const uint64_t *z_codes,
                           double *out_lats, double *out_lngs, size_t count) {
    if (count == 0) return;
    
    __m256d lat_scale_v = _mm256_set1_pd(GEO_INTERNAL_LAT_DENORM_SCALE);
    __m256d lng_scale_v = _mm256_set1_pd(GEO_INTERNAL_LNG_DENORM_SCALE);
    __m256d lat_offset_v = _mm256_set1_pd(-90.0);
    __m256d lng_offset_v = _mm256_set1_pd(-180.0);
    
    size_t i = 0;
    
    // Process 4 codes at a time
    for (; i + 3 < count; i += 4) {
        double lat_arr[4], lng_arr[4];
        
        for (int j = 0; j < 4; j++) {
            lat_arr[j] = (double)geo_internal_compact_bits(z_codes[i + j]);
            lng_arr[j] = (double)geo_internal_compact_bits(z_codes[i + j] >> 1);
        }
        
        __m256d lat_v = _mm256_loadu_pd(lat_arr);
        __m256d lng_v = _mm256_loadu_pd(lng_arr);
        
        lat_v = _mm256_add_pd(_mm256_mul_pd(lat_v, lat_scale_v), lat_offset_v);
        lng_v = _mm256_add_pd(_mm256_mul_pd(lng_v, lng_scale_v), lng_offset_v);
        
        _mm256_storeu_pd(&out_lats[i], lat_v);
        _mm256_storeu_pd(&out_lngs[i], lng_v);
    }
    
    for (; i < count; i++) {
        GeoPoint p = geo_decode(z_codes[i]);
        out_lats[i] = p.lat;
        out_lngs[i] = p.lng;
    }
}

#elif defined(GEO_SIMD_X86_SSE4) || defined(GEO_SIMD_X86_SSE2)

void geo_simd_decode_batch(const uint64_t *z_codes,
                           double *out_lats, double *out_lngs, size_t count) {
    if (count == 0) return;
    
    __m128d lat_scale_v = _mm_set1_pd(GEO_INTERNAL_LAT_DENORM_SCALE);
    __m128d lng_scale_v = _mm_set1_pd(GEO_INTERNAL_LNG_DENORM_SCALE);
    __m128d lat_offset_v = _mm_set1_pd(-90.0);
    __m128d lng_offset_v = _mm_set1_pd(-180.0);
    
    size_t i = 0;
    
    for (; i + 1 < count; i += 2) {
        double lat_arr[2], lng_arr[2];
        
        for (int j = 0; j < 2; j++) {
            lat_arr[j] = (double)geo_internal_compact_bits(z_codes[i + j]);
            lng_arr[j] = (double)geo_internal_compact_bits(z_codes[i + j] >> 1);
        }
        
        __m128d lat_v = _mm_loadu_pd(lat_arr);
        __m128d lng_v = _mm_loadu_pd(lng_arr);
        
        lat_v = _mm_add_pd(_mm_mul_pd(lat_v, lat_scale_v), lat_offset_v);
        lng_v = _mm_add_pd(_mm_mul_pd(lng_v, lng_scale_v), lng_offset_v);
        
        _mm_storeu_pd(&out_lats[i], lat_v);
        _mm_storeu_pd(&out_lngs[i], lng_v);
    }
    
    for (; i < count; i++) {
        GeoPoint p = geo_decode(z_codes[i]);
        out_lats[i] = p.lat;
        out_lngs[i] = p.lng;
    }
}

#endif

// =========================================================
// SIMD Haversine Distance - x86 AVX2
// =========================================================

#if defined(GEO_SIMD_X86_AVX2)

// Polynomial approximation for sin (optimized for range reduction)
static inline __m256d avx2_sin_approx(__m256d x) {
    // Constants for sin approximation
    __m256d pi = _mm256_set1_pd(M_PI);
    __m256d two_pi = _mm256_set1_pd(2.0 * M_PI);
    __m256d inv_two_pi = _mm256_set1_pd(1.0 / (2.0 * M_PI));
    
    // Range reduction to [-π, π]
    __m256d n = _mm256_round_pd(_mm256_mul_pd(x, inv_two_pi), _MM_FROUND_TO_NEAREST_INT);
    x = _mm256_sub_pd(x, _mm256_mul_pd(n, two_pi));
    
    // Taylor series coefficients
    __m256d c1 = _mm256_set1_pd(-1.0 / 6.0);
    __m256d c2 = _mm256_set1_pd(1.0 / 120.0);
    __m256d c3 = _mm256_set1_pd(-1.0 / 5040.0);
    
    __m256d x2 = _mm256_mul_pd(x, x);
    __m256d x3 = _mm256_mul_pd(x2, x);
    __m256d x5 = _mm256_mul_pd(x3, x2);
    __m256d x7 = _mm256_mul_pd(x5, x2);
    
    // sin(x) ≈ x - x³/6 + x⁵/120 - x⁷/5040
    __m256d result = _mm256_add_pd(x,
                      _mm256_add_pd(_mm256_mul_pd(c1, x3),
                      _mm256_add_pd(_mm256_mul_pd(c2, x5),
                                    _mm256_mul_pd(c3, x7))));
    
    return result;
}

static inline __m256d avx2_cos_approx(__m256d x) {
    __m256d half_pi = _mm256_set1_pd(M_PI / 2.0);
    return avx2_sin_approx(_mm256_add_pd(x, half_pi));
}

static inline __m256d avx2_sqrt(__m256d x) {
    return _mm256_sqrt_pd(x);
}

void geo_simd_haversine_batch(double lat1, double lng1,
                              const double *lats, const double *lngs,
                              double *out_dist, size_t count) {
    if (count == 0) return;
    
    __m256d deg_to_rad = _mm256_set1_pd(GEO_INTERNAL_DEG_TO_RAD);
    __m256d earth_radius = _mm256_set1_pd(GEO_INTERNAL_EARTH_RADIUS_KM);
    __m256d lat1_rad_v = _mm256_set1_pd(lat1 * GEO_INTERNAL_DEG_TO_RAD);
    __m256d lng1_rad_v = _mm256_set1_pd(lng1 * GEO_INTERNAL_DEG_TO_RAD);
    __m256d half = _mm256_set1_pd(0.5);
    __m256d two = _mm256_set1_pd(2.0);
    __m256d one = _mm256_set1_pd(1.0);
    
    __m256d cos_lat1 = avx2_cos_approx(lat1_rad_v);
    
    size_t i = 0;
    
    // Process 4 points at a time
    for (; i + 3 < count; i += 4) {
        __m256d lat2_v = _mm256_loadu_pd(&lats[i]);
        __m256d lng2_v = _mm256_loadu_pd(&lngs[i]);
        
        __m256d lat2_rad = _mm256_mul_pd(lat2_v, deg_to_rad);
        __m256d lng2_rad = _mm256_mul_pd(lng2_v, deg_to_rad);
        
        __m256d dlat = _mm256_sub_pd(lat2_rad, lat1_rad_v);
        __m256d dlng = _mm256_sub_pd(lng2_rad, lng1_rad_v);
        
        __m256d sin_dlat_2 = avx2_sin_approx(_mm256_mul_pd(dlat, half));
        __m256d sin_dlng_2 = avx2_sin_approx(_mm256_mul_pd(dlng, half));
        
        __m256d sin_dlat_2_sq = _mm256_mul_pd(sin_dlat_2, sin_dlat_2);
        __m256d sin_dlng_2_sq = _mm256_mul_pd(sin_dlng_2, sin_dlng_2);
        
        __m256d cos_lat2 = avx2_cos_approx(lat2_rad);
        
        __m256d a = _mm256_add_pd(sin_dlat_2_sq,
                     _mm256_mul_pd(_mm256_mul_pd(cos_lat1, cos_lat2), sin_dlng_2_sq));
        
        __m256d sqrt_a = avx2_sqrt(a);
        
        // asin approximation for small angles
        __m256d x2 = _mm256_mul_pd(sqrt_a, sqrt_a);
        __m256d c1 = _mm256_set1_pd(1.0 / 6.0);
        __m256d c2 = _mm256_set1_pd(3.0 / 40.0);
        __m256d asin_approx = _mm256_mul_pd(sqrt_a,
                               _mm256_add_pd(one,
                               _mm256_mul_pd(x2,
                               _mm256_add_pd(c1,
                               _mm256_mul_pd(x2, c2)))));
        
        __m256d c = _mm256_mul_pd(two, asin_approx);
        __m256d dist = _mm256_mul_pd(earth_radius, c);
        
        _mm256_storeu_pd(&out_dist[i], dist);
    }
    
    for (; i < count; i++) {
        out_dist[i] = geo_haversine_km(lat1, lng1, lats[i], lngs[i]);
    }
}

#elif defined(GEO_SIMD_X86_SSE4) || defined(GEO_SIMD_X86_SSE2)

// SSE versions of trig approximations
static inline __m128d sse_sin_approx(__m128d x) {
    __m128d two_pi = _mm_set1_pd(2.0 * M_PI);
    __m128d inv_two_pi = _mm_set1_pd(1.0 / (2.0 * M_PI));
    
    __m128d n = _mm_round_pd(_mm_mul_pd(x, inv_two_pi), _MM_FROUND_TO_NEAREST_INT);
    x = _mm_sub_pd(x, _mm_mul_pd(n, two_pi));
    
    __m128d c1 = _mm_set1_pd(-1.0 / 6.0);
    __m128d c2 = _mm_set1_pd(1.0 / 120.0);
    
    __m128d x2 = _mm_mul_pd(x, x);
    __m128d x3 = _mm_mul_pd(x2, x);
    __m128d x5 = _mm_mul_pd(x3, x2);
    
    return _mm_add_pd(x, _mm_add_pd(_mm_mul_pd(c1, x3), _mm_mul_pd(c2, x5)));
}

static inline __m128d sse_cos_approx(__m128d x) {
    return sse_sin_approx(_mm_add_pd(x, _mm_set1_pd(M_PI / 2.0)));
}

void geo_simd_haversine_batch(double lat1, double lng1,
                              const double *lats, const double *lngs,
                              double *out_dist, size_t count) {
    if (count == 0) return;
    
    __m128d deg_to_rad = _mm_set1_pd(GEO_INTERNAL_DEG_TO_RAD);
    __m128d earth_radius = _mm_set1_pd(GEO_INTERNAL_EARTH_RADIUS_KM);
    __m128d lat1_rad_v = _mm_set1_pd(lat1 * GEO_INTERNAL_DEG_TO_RAD);
    __m128d lng1_rad_v = _mm_set1_pd(lng1 * GEO_INTERNAL_DEG_TO_RAD);
    __m128d half = _mm_set1_pd(0.5);
    __m128d two = _mm_set1_pd(2.0);
    __m128d one = _mm_set1_pd(1.0);
    
    __m128d cos_lat1 = sse_cos_approx(lat1_rad_v);
    
    size_t i = 0;
    
    for (; i + 1 < count; i += 2) {
        __m128d lat2_v = _mm_loadu_pd(&lats[i]);
        __m128d lng2_v = _mm_loadu_pd(&lngs[i]);
        
        __m128d lat2_rad = _mm_mul_pd(lat2_v, deg_to_rad);
        __m128d lng2_rad = _mm_mul_pd(lng2_v, deg_to_rad);
        
        __m128d dlat = _mm_sub_pd(lat2_rad, lat1_rad_v);
        __m128d dlng = _mm_sub_pd(lng2_rad, lng1_rad_v);
        
        __m128d sin_dlat_2 = sse_sin_approx(_mm_mul_pd(dlat, half));
        __m128d sin_dlng_2 = sse_sin_approx(_mm_mul_pd(dlng, half));
        
        __m128d sin_dlat_2_sq = _mm_mul_pd(sin_dlat_2, sin_dlat_2);
        __m128d sin_dlng_2_sq = _mm_mul_pd(sin_dlng_2, sin_dlng_2);
        
        __m128d cos_lat2 = sse_cos_approx(lat2_rad);
        
        __m128d a = _mm_add_pd(sin_dlat_2_sq,
                    _mm_mul_pd(_mm_mul_pd(cos_lat1, cos_lat2), sin_dlng_2_sq));
        
        __m128d sqrt_a = _mm_sqrt_pd(a);
        
        __m128d x2 = _mm_mul_pd(sqrt_a, sqrt_a);
        __m128d asin_approx = _mm_mul_pd(sqrt_a,
                              _mm_add_pd(one,
                              _mm_mul_pd(x2, _mm_set1_pd(1.0 / 6.0))));
        
        __m128d c = _mm_mul_pd(two, asin_approx);
        __m128d dist = _mm_mul_pd(earth_radius, c);
        
        _mm_storeu_pd(&out_dist[i], dist);
    }
    
    for (; i < count; i++) {
        out_dist[i] = geo_haversine_km(lat1, lng1, lats[i], lngs[i]);
    }
}

#endif

// =========================================================
// SIMD Fast Distance - x86
// =========================================================

#if defined(GEO_SIMD_X86_AVX2)

void geo_simd_fast_distance_batch(double lat1, double lng1,
                                  const double *lats, const double *lngs,
                                  double *out_dist, size_t count) {
    if (count == 0) return;
    
    __m256d km_per_deg = _mm256_set1_pd(GEO_INTERNAL_KM_PER_DEG);
    __m256d lat1_v = _mm256_set1_pd(lat1);
    __m256d lng1_v = _mm256_set1_pd(lng1);
    __m256d deg_to_rad = _mm256_set1_pd(GEO_INTERNAL_DEG_TO_RAD);
    __m256d half = _mm256_set1_pd(0.5);
    
    size_t i = 0;
    
    for (; i + 3 < count; i += 4) {
        __m256d lat2_v = _mm256_loadu_pd(&lats[i]);
        __m256d lng2_v = _mm256_loadu_pd(&lngs[i]);
        
        __m256d lat_mid = _mm256_mul_pd(_mm256_add_pd(lat1_v, lat2_v), half);
        __m256d lat_mid_rad = _mm256_mul_pd(lat_mid, deg_to_rad);
        __m256d cos_lat = avx2_cos_approx(lat_mid_rad);
        
        __m256d dlat = _mm256_mul_pd(_mm256_sub_pd(lat2_v, lat1_v), km_per_deg);
        __m256d dlng = _mm256_mul_pd(_mm256_mul_pd(_mm256_sub_pd(lng2_v, lng1_v), km_per_deg), cos_lat);
        
        __m256d dist_sq = _mm256_add_pd(_mm256_mul_pd(dlat, dlat), _mm256_mul_pd(dlng, dlng));
        __m256d dist = _mm256_sqrt_pd(dist_sq);
        
        _mm256_storeu_pd(&out_dist[i], dist);
    }
    
    for (; i < count; i++) {
        out_dist[i] = geo_fast_distance_km(lat1, lng1, lats[i], lngs[i]);
    }
}

#elif defined(GEO_SIMD_X86_SSE4) || defined(GEO_SIMD_X86_SSE2)

void geo_simd_fast_distance_batch(double lat1, double lng1,
                                  const double *lats, const double *lngs,
                                  double *out_dist, size_t count) {
    if (count == 0) return;
    
    __m128d km_per_deg = _mm_set1_pd(GEO_INTERNAL_KM_PER_DEG);
    __m128d lat1_v = _mm_set1_pd(lat1);
    __m128d lng1_v = _mm_set1_pd(lng1);
    __m128d deg_to_rad = _mm_set1_pd(GEO_INTERNAL_DEG_TO_RAD);
    __m128d half = _mm_set1_pd(0.5);
    
    size_t i = 0;
    
    for (; i + 1 < count; i += 2) {
        __m128d lat2_v = _mm_loadu_pd(&lats[i]);
        __m128d lng2_v = _mm_loadu_pd(&lngs[i]);
        
        __m128d lat_mid = _mm_mul_pd(_mm_add_pd(lat1_v, lat2_v), half);
        __m128d lat_mid_rad = _mm_mul_pd(lat_mid, deg_to_rad);
        __m128d cos_lat = sse_cos_approx(lat_mid_rad);
        
        __m128d dlat = _mm_mul_pd(_mm_sub_pd(lat2_v, lat1_v), km_per_deg);
        __m128d dlng = _mm_mul_pd(_mm_mul_pd(_mm_sub_pd(lng2_v, lng1_v), km_per_deg), cos_lat);
        
        __m128d dist_sq = _mm_add_pd(_mm_mul_pd(dlat, dlat), _mm_mul_pd(dlng, dlng));
        __m128d dist = _mm_sqrt_pd(dist_sq);
        
        _mm_storeu_pd(&out_dist[i], dist);
    }
    
    for (; i < count; i++) {
        out_dist[i] = geo_fast_distance_km(lat1, lng1, lats[i], lngs[i]);
    }
}

#endif

// =========================================================
// SIMD Range Filtering - x86
// =========================================================

#if defined(GEO_SIMD_X86_AVX2)

size_t geo_simd_filter_range(const uint64_t *z_codes, size_t count,
                             uint64_t z_min, uint64_t z_max,
                             uint8_t *out_mask) {
    if (count == 0) return 0;
    
    __m256i min_v = _mm256_set1_epi64x((long long)z_min);
    __m256i max_v = _mm256_set1_epi64x((long long)z_max);
    
    size_t matches = 0;
    size_t i = 0;
    
    // Process 4 codes at a time
    for (; i + 3 < count; i += 4) {
        __m256i z_v = _mm256_loadu_si256((const __m256i*)&z_codes[i]);
        
        // For unsigned comparison, we use signed comparison with offset
        // or use AVX-512 if available. Here using a workaround:
        // Check z >= min: (z - min) has no underflow
        // Check z <= max: (max - z) has no underflow
        
        // Simple approach: extract and compare
        uint64_t z_arr[4];
        _mm256_storeu_si256((__m256i*)z_arr, z_v);
        
        for (int j = 0; j < 4; j++) {
            uint8_t in = (z_arr[j] >= z_min && z_arr[j] <= z_max) ? 1 : 0;
            out_mask[i + j] = in;
            matches += in;
        }
    }
    
    for (; i < count; i++) {
        uint8_t in = (z_codes[i] >= z_min && z_codes[i] <= z_max) ? 1 : 0;
        out_mask[i] = in;
        matches += in;
    }
    
    return matches;
}

#elif defined(GEO_SIMD_X86_SSE4) || defined(GEO_SIMD_X86_SSE2)

size_t geo_simd_filter_range(const uint64_t *z_codes, size_t count,
                             uint64_t z_min, uint64_t z_max,
                             uint8_t *out_mask) {
    if (count == 0) return 0;
    
    size_t matches = 0;
    
    // SSE2 doesn't have good 64-bit comparison, fall back to scalar
    for (size_t i = 0; i < count; i++) {
        uint8_t in = (z_codes[i] >= z_min && z_codes[i] <= z_max) ? 1 : 0;
        out_mask[i] = in;
        matches += in;
    }
    
    return matches;
}

#endif

// =========================================================
// SIMD Bounding Box Filtering - x86
// =========================================================

#if defined(GEO_SIMD_X86_AVX2)

size_t geo_simd_filter_bbox(const double *lats, const double *lngs, size_t count,
                            double min_lat, double max_lat,
                            double min_lng, double max_lng,
                            uint8_t *out_mask) {
    if (count == 0) return 0;
    
    __m256d min_lat_v = _mm256_set1_pd(min_lat);
    __m256d max_lat_v = _mm256_set1_pd(max_lat);
    __m256d min_lng_v = _mm256_set1_pd(min_lng);
    __m256d max_lng_v = _mm256_set1_pd(max_lng);
    
    size_t matches = 0;
    size_t i = 0;
    
    for (; i + 3 < count; i += 4) {
        __m256d lat_v = _mm256_loadu_pd(&lats[i]);
        __m256d lng_v = _mm256_loadu_pd(&lngs[i]);
        
        __m256d lat_ge_min = _mm256_cmp_pd(lat_v, min_lat_v, _CMP_GE_OQ);
        __m256d lat_le_max = _mm256_cmp_pd(lat_v, max_lat_v, _CMP_LE_OQ);
        __m256d lat_ok = _mm256_and_pd(lat_ge_min, lat_le_max);
        
        __m256d lng_ge_min = _mm256_cmp_pd(lng_v, min_lng_v, _CMP_GE_OQ);
        __m256d lng_le_max = _mm256_cmp_pd(lng_v, max_lng_v, _CMP_LE_OQ);
        __m256d lng_ok = _mm256_and_pd(lng_ge_min, lng_le_max);
        
        __m256d in_bbox = _mm256_and_pd(lat_ok, lng_ok);
        
        int mask = _mm256_movemask_pd(in_bbox);
        
        out_mask[i]     = (mask & 1) ? 1 : 0;
        out_mask[i + 1] = (mask & 2) ? 1 : 0;
        out_mask[i + 2] = (mask & 4) ? 1 : 0;
        out_mask[i + 3] = (mask & 8) ? 1 : 0;
        
        matches += __builtin_popcount(mask);
    }
    
    for (; i < count; i++) {
        uint8_t in = (lats[i] >= min_lat && lats[i] <= max_lat &&
                      lngs[i] >= min_lng && lngs[i] <= max_lng) ? 1 : 0;
        out_mask[i] = in;
        matches += in;
    }
    
    return matches;
}

#elif defined(GEO_SIMD_X86_SSE4) || defined(GEO_SIMD_X86_SSE2)

size_t geo_simd_filter_bbox(const double *lats, const double *lngs, size_t count,
                            double min_lat, double max_lat,
                            double min_lng, double max_lng,
                            uint8_t *out_mask) {
    if (count == 0) return 0;
    
    __m128d min_lat_v = _mm_set1_pd(min_lat);
    __m128d max_lat_v = _mm_set1_pd(max_lat);
    __m128d min_lng_v = _mm_set1_pd(min_lng);
    __m128d max_lng_v = _mm_set1_pd(max_lng);
    
    size_t matches = 0;
    size_t i = 0;
    
    for (; i + 1 < count; i += 2) {
        __m128d lat_v = _mm_loadu_pd(&lats[i]);
        __m128d lng_v = _mm_loadu_pd(&lngs[i]);
        
        __m128d lat_ge_min = _mm_cmpge_pd(lat_v, min_lat_v);
        __m128d lat_le_max = _mm_cmple_pd(lat_v, max_lat_v);
        __m128d lat_ok = _mm_and_pd(lat_ge_min, lat_le_max);
        
        __m128d lng_ge_min = _mm_cmpge_pd(lng_v, min_lng_v);
        __m128d lng_le_max = _mm_cmple_pd(lng_v, max_lng_v);
        __m128d lng_ok = _mm_and_pd(lng_ge_min, lng_le_max);
        
        __m128d in_bbox = _mm_and_pd(lat_ok, lng_ok);
        
        int mask = _mm_movemask_pd(in_bbox);
        
        out_mask[i]     = (mask & 1) ? 1 : 0;
        out_mask[i + 1] = (mask & 2) ? 1 : 0;
        
        matches += __builtin_popcount(mask);
    }
    
    for (; i < count; i++) {
        uint8_t in = (lats[i] >= min_lat && lats[i] <= max_lat &&
                      lngs[i] >= min_lng && lngs[i] <= max_lng) ? 1 : 0;
        out_mask[i] = in;
        matches += in;
    }
    
    return matches;
}

#endif

// =========================================================
// SIMD Radius Filtering - x86
// =========================================================

size_t geo_simd_filter_radius(const double *lats, const double *lngs, size_t count,
                              double center_lat, double center_lng, double radius_km,
                              uint8_t *out_mask) {
    if (count == 0) return 0;
    
    double min_lat, max_lat, min_lng, max_lng;
    geo_bounding_box(center_lat, center_lng, radius_km,
                     &min_lat, &max_lat, &min_lng, &max_lng);
    
    // Allocate temporary arrays
    double *candidate_lats = (double*)aligned_alloc(32, count * sizeof(double));
    double *candidate_lngs = (double*)aligned_alloc(32, count * sizeof(double));
    size_t *candidate_idx = (size_t*)malloc(count * sizeof(size_t));
    
    if (!candidate_lats || !candidate_lngs || !candidate_idx) {
        free(candidate_lats);
        free(candidate_lngs);
        free(candidate_idx);
        
        size_t matches = 0;
        for (size_t i = 0; i < count; i++) {
            double dist = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);
            out_mask[i] = (dist <= radius_km) ? 1 : 0;
            matches += out_mask[i];
        }
        return matches;
    }
    
    memset(out_mask, 0, count);
    
    // First pass: bounding box filter with SIMD
    uint8_t *bbox_mask = (uint8_t*)aligned_alloc(32, count * sizeof(uint8_t));
    geo_simd_filter_bbox(lats, lngs, count, min_lat, max_lat, min_lng, max_lng, bbox_mask);
    
    size_t num_candidates = 0;
    for (size_t i = 0; i < count; i++) {
        if (bbox_mask[i]) {
            candidate_lats[num_candidates] = lats[i];
            candidate_lngs[num_candidates] = lngs[i];
            candidate_idx[num_candidates] = i;
            num_candidates++;
        }
    }
    free(bbox_mask);
    
    // Second pass: precise distance check
    if (num_candidates > 0) {
        double *distances = (double*)aligned_alloc(32, num_candidates * sizeof(double));
        if (distances) {
            geo_simd_haversine_batch(center_lat, center_lng,
                                      candidate_lats, candidate_lngs,
                                      distances, num_candidates);
            
#if defined(GEO_SIMD_X86_AVX2)
            __m256d radius_v = _mm256_set1_pd(radius_km);
            size_t j = 0;
            
            for (; j + 3 < num_candidates; j += 4) {
                __m256d dist_v = _mm256_loadu_pd(&distances[j]);
                __m256d in_radius = _mm256_cmp_pd(dist_v, radius_v, _CMP_LE_OQ);
                int mask = _mm256_movemask_pd(in_radius);
                
                if (mask & 1) out_mask[candidate_idx[j]] = 1;
                if (mask & 2) out_mask[candidate_idx[j + 1]] = 1;
                if (mask & 4) out_mask[candidate_idx[j + 2]] = 1;
                if (mask & 8) out_mask[candidate_idx[j + 3]] = 1;
            }
            
            for (; j < num_candidates; j++) {
                if (distances[j] <= radius_km) {
                    out_mask[candidate_idx[j]] = 1;
                }
            }
#else
            for (size_t j = 0; j < num_candidates; j++) {
                if (distances[j] <= radius_km) {
                    out_mask[candidate_idx[j]] = 1;
                }
            }
#endif
            free(distances);
        }
    }
    
    free(candidate_lats);
    free(candidate_lngs);
    free(candidate_idx);
    
    size_t matches = 0;
    for (size_t i = 0; i < count; i++) {
        matches += out_mask[i];
    }
    
    return matches;
}

// =========================================================
// SIMD Bit Manipulation - x86
// =========================================================

void geo_simd_spread_bits_batch(const uint32_t *values, uint64_t *out, size_t count) {
    // Morton encoding bit manipulation is tricky to vectorize effectively
    // Using scalar with loop unrolling
    size_t i = 0;
    for (; i + 3 < count; i += 4) {
        out[i]     = geo_internal_spread_bits(values[i]);
        out[i + 1] = geo_internal_spread_bits(values[i + 1]);
        out[i + 2] = geo_internal_spread_bits(values[i + 2]);
        out[i + 3] = geo_internal_spread_bits(values[i + 3]);
    }
    for (; i < count; i++) {
        out[i] = geo_internal_spread_bits(values[i]);
    }
}

void geo_simd_compact_bits_batch(const uint64_t *values, uint32_t *out, size_t count) {
    size_t i = 0;
    for (; i + 3 < count; i += 4) {
        out[i]     = geo_internal_compact_bits(values[i]);
        out[i + 1] = geo_internal_compact_bits(values[i + 1]);
        out[i + 2] = geo_internal_compact_bits(values[i + 2]);
        out[i + 3] = geo_internal_compact_bits(values[i + 3]);
    }
    for (; i < count; i++) {
        out[i] = geo_internal_compact_bits(values[i]);
    }
}

// =========================================================
// Utility Functions - x86
// =========================================================

bool geo_simd_available(void) {
    detect_cpu_features();
#if defined(GEO_SIMD_X86_AVX2)
    return cpu_has_avx2 != 0;
#elif defined(GEO_SIMD_X86_SSE4)
    return cpu_has_sse4 != 0;
#else
    return true; // SSE2 is baseline for x86-64
#endif
}

const char* geo_simd_get_name(void) {
    detect_cpu_features();
    if (cpu_has_avx2) return "x86-64 AVX2";
    if (cpu_has_avx) return "x86-64 AVX";
    if (cpu_has_sse4) return "x86-64 SSE4";
    return "x86-64 SSE2";
}

size_t geo_simd_optimal_batch_size(void) {
    detect_cpu_features();
    if (cpu_has_avx2) return 512;
    if (cpu_has_sse4) return 256;
    return 128;
}

// =========================================================
// Benchmark Functions - x86
// =========================================================

GeoSimdBenchmark geo_simd_benchmark_encode(size_t count, int iterations) {
    GeoSimdBenchmark result = {0};
    result.operation_name = "encode";
    result.operations = count * iterations;
    
    double *lats = (double*)aligned_alloc(32, count * sizeof(double));
    double *lngs = (double*)aligned_alloc(32, count * sizeof(double));
    uint64_t *z_codes = (uint64_t*)aligned_alloc(32, count * sizeof(uint64_t));
    
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
    
    uint64_t *z_codes = (uint64_t*)aligned_alloc(32, count * sizeof(uint64_t));
    double *lats = (double*)aligned_alloc(32, count * sizeof(double));
    double *lngs = (double*)aligned_alloc(32, count * sizeof(double));
    
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
    
    double *lats = (double*)aligned_alloc(32, count * sizeof(double));
    double *lngs = (double*)aligned_alloc(32, count * sizeof(double));
    double *dists = (double*)aligned_alloc(32, count * sizeof(double));
    
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
    
    double *lats = (double*)aligned_alloc(32, count * sizeof(double));
    double *lngs = (double*)aligned_alloc(32, count * sizeof(double));
    uint8_t *mask = (uint8_t*)aligned_alloc(32, count * sizeof(uint8_t));
    
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
    
    start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        geo_simd_filter_radius(lats, lngs, count, center_lat, center_lng, radius_km, mask);
    }
    result.simd_time_ms = geo_get_time_ms() - start;
    
    result.speedup = result.scalar_time_ms / result.simd_time_ms;
    
    free(lats); free(lngs); free(mask);
    return result;
}

#endif // GEO_SIMD_X86_*
