/*
 * GeoIndex SIMD Optimizations - ARM64 NEON Implementation
 * 
 * Optimized implementations using ARM NEON intrinsics for:
 *   - Apple Silicon (M1/M2/M3)
 *   - ARM Cortex-A series
 *   - AWS Graviton
 */

#include "geobolt/geo_index_simd.h"
#include "geobolt/geo_index.h"
#include "geo_index_internal.h"
#include "geo_index_simd_scalar_kernels.h"

#if defined(GEO_SIMD_ARM64)

#include <arm_neon.h>
#include <stdlib.h>
#include <string.h>

static volatile uint64_t g_benchmark_sink_u64;
static volatile double g_benchmark_sink_double;

// =========================================================
// Morton bit manipulation - ARM64 NEON
// =========================================================

static inline uint64x2_t neon_spread_bits(uint64x2_t values)
{
    values = vandq_u64(vorrq_u64(values, vshlq_n_u64(values, 16)), vdupq_n_u64(UINT64_C(0x0000FFFF0000FFFF)));
    values = vandq_u64(vorrq_u64(values, vshlq_n_u64(values, 8)), vdupq_n_u64(UINT64_C(0x00FF00FF00FF00FF)));
    values = vandq_u64(vorrq_u64(values, vshlq_n_u64(values, 4)), vdupq_n_u64(UINT64_C(0x0F0F0F0F0F0F0F0F)));
    values = vandq_u64(vorrq_u64(values, vshlq_n_u64(values, 2)), vdupq_n_u64(UINT64_C(0x3333333333333333)));

    return vandq_u64(vorrq_u64(values, vshlq_n_u64(values, 1)), vdupq_n_u64(UINT64_C(0x5555555555555555)));
}

static inline uint64x2_t neon_compact_bits(uint64x2_t values)
{
    values = vandq_u64(values, vdupq_n_u64(UINT64_C(0x5555555555555555)));
    values = vandq_u64(vorrq_u64(values, vshrq_n_u64(values, 1)), vdupq_n_u64(UINT64_C(0x3333333333333333)));
    values = vandq_u64(vorrq_u64(values, vshrq_n_u64(values, 2)), vdupq_n_u64(UINT64_C(0x0F0F0F0F0F0F0F0F)));
    values = vandq_u64(vorrq_u64(values, vshrq_n_u64(values, 4)), vdupq_n_u64(UINT64_C(0x00FF00FF00FF00FF)));
    values = vandq_u64(vorrq_u64(values, vshrq_n_u64(values, 8)), vdupq_n_u64(UINT64_C(0x0000FFFF0000FFFF)));

    return vandq_u64(vorrq_u64(values, vshrq_n_u64(values, 16)), vdupq_n_u64(UINT64_C(0x00000000FFFFFFFF)));
}

static inline void neon_decode2(uint64x2_t morton_codes, float64x2_t *latitudes, float64x2_t *longitudes)
{
    uint64x2_t normalized_latitudes = neon_compact_bits(morton_codes);
    uint64x2_t normalized_longitudes = neon_compact_bits(vshrq_n_u64(morton_codes, 1));

    *latitudes = vfmaq_f64(vdupq_n_f64(-90.0),
                           vcvtq_f64_u64(normalized_latitudes),
                           vdupq_n_f64(GEO_INTERNAL_LAT_DENORM_SCALE));

    *longitudes = vfmaq_f64(vdupq_n_f64(-180.0),
                            vcvtq_f64_u64(normalized_longitudes),
                            vdupq_n_f64(GEO_INTERNAL_LNG_DENORM_SCALE));
}

// =========================================================
// SIMD Batch Encoding - ARM64 NEON
// =========================================================

static inline uint64x2_t neon_encode2(float64x2_t latitudes, float64x2_t longitudes)
{
    float64x2_t normalized_latitudes = vmulq_f64(vaddq_f64(latitudes, vdupq_n_f64(GEO_INTERNAL_LAT_OFFSET)),
                                                  vdupq_n_f64(GEO_INTERNAL_LAT_SCALE));

    float64x2_t normalized_longitudes = vmulq_f64(vaddq_f64(longitudes, vdupq_n_f64(GEO_INTERNAL_LNG_OFFSET)),
                                                   vdupq_n_f64(GEO_INTERNAL_LNG_SCALE));

    normalized_latitudes = vmaxq_f64(vdupq_n_f64(0.0),
                                     vminq_f64(normalized_latitudes, vdupq_n_f64(GEO_INTERNAL_COORD_MAX)));

    normalized_longitudes = vmaxq_f64(vdupq_n_f64(0.0),
                                      vminq_f64(normalized_longitudes, vdupq_n_f64(GEO_INTERNAL_COORD_MAX)));

    uint64x2_t latitude_bits = neon_spread_bits(vcvtq_u64_f64(normalized_latitudes));
    uint64x2_t longitude_bits = neon_spread_bits(vcvtq_u64_f64(normalized_longitudes));

    return vorrq_u64(latitude_bits, vshlq_n_u64(longitude_bits, 1));
}

void geo_simd_encode_batch(const double *lats, const double *lngs, uint64_t *out_z, size_t count)
{
    if (count == 0) {
        return;
    }

    size_t i = 0;

    for (; i + 1 < count; i += 2) {
        vst1q_u64(out_z + i, neon_encode2(vld1q_f64(lats + i), vld1q_f64(lngs + i)));
    }

    geo_scalar_encode_batch(lats + i, lngs + i, out_z + i, count - i);
}

void geo_simd_encode_records(const uint64_t *ids,
                             const double *lats,
                             const double *lngs,
                             GeoRecord *records,
                             size_t count)
{
    size_t i = 0;

    for (; i + 1 < count; i += 2) {
        uint64x2x2_t record_fields = {
            .val = {
                vld1q_u64(ids + i),
                neon_encode2(vld1q_f64(lats + i), vld1q_f64(lngs + i)),
            },
        };

        vst2q_u64((uint64_t *) (records + i), record_fields);
    }

    geo_scalar_encode_records(ids + i, lats + i, lngs + i, records + i, count - i);
}

bool geo_simd_validate_points(const double *lats, const double *lngs, size_t count)
{
    float64x2_t minimum_latitude = vdupq_n_f64(GEO_MIN_LAT);
    float64x2_t maximum_latitude = vdupq_n_f64(GEO_MAX_LAT);
    float64x2_t minimum_longitude = vdupq_n_f64(GEO_MIN_LNG);
    float64x2_t maximum_longitude = vdupq_n_f64(GEO_MAX_LNG);
    size_t i = 0;

    for (; i + 1 < count; i += 2) {
        float64x2_t latitudes = vld1q_f64(lats + i);
        float64x2_t longitudes = vld1q_f64(lngs + i);
        uint64x2_t valid_latitudes = vandq_u64(vcgeq_f64(latitudes, minimum_latitude),
                                               vcleq_f64(latitudes, maximum_latitude));

        uint64x2_t valid_longitudes = vandq_u64(vcgeq_f64(longitudes, minimum_longitude),
                                                vcleq_f64(longitudes, maximum_longitude));

        uint64x2_t valid = vandq_u64(valid_latitudes, valid_longitudes);

        if (vgetq_lane_u64(valid, 0) != UINT64_MAX || vgetq_lane_u64(valid, 1) != UINT64_MAX) {
            return false;
        }
    }

    return geo_scalar_validate_points(lats + i, lngs + i, count - i);
}

// =========================================================
// SIMD Batch Decoding - ARM64 NEON
// =========================================================

void geo_simd_decode_batch(const uint64_t *z_codes, double *out_lats, double *out_lngs, size_t count)
{
    if (count == 0) {
        return;
    }

    size_t i = 0;

    // Process 2 codes at a time
    for (; i + 1 < count; i += 2) {
        uint64x2_t morton_codes = vld1q_u64(&z_codes[i]);
        float64x2_t lat_v;
        float64x2_t lng_v;

        neon_decode2(morton_codes, &lat_v, &lng_v);

        // Store
        vst1q_f64(&out_lats[i], lat_v);
        vst1q_f64(&out_lngs[i], lng_v);
    }

    geo_scalar_decode_batch(z_codes + i, out_lats + i, out_lngs + i, count - i);
}

void geo_simd_extract_interleaved_codes(const GeoRecord *records, uint64_t *out_codes, size_t count)
{
    size_t i = 0;

    for (; i + 1 < count; i += 2) {
        uint64x2x2_t pair_vectors = vld2q_u64((const uint64_t *) (records + i));

        vst1q_u64(out_codes + i, pair_vectors.val[1]);
    }

    geo_scalar_extract_interleaved_codes(records + i, out_codes + i, count - i);
}

void geo_simd_decode_batch_narrow(const uint64_t *z_codes, double *out_lats, double *out_lngs, size_t count)
{
    geo_simd_decode_batch(z_codes, out_lats, out_lngs, count);
}

void geo_simd_extract_interleaved_codes_narrow(const GeoRecord *records, uint64_t *out_codes, size_t count)
{
    geo_simd_extract_interleaved_codes(records, out_codes, count);
}

void geo_simd_decode_interleaved_records_narrow(const GeoRecord *records,
                                                double *out_lats,
                                                double *out_lngs,
                                                size_t count)
{
    size_t i = 0;

    for (; i + 1 < count; i += 2) {
        uint64x2x2_t record_fields = vld2q_u64((const uint64_t *) (records + i));
        float64x2_t latitudes;
        float64x2_t longitudes;

        neon_decode2(record_fields.val[1], &latitudes, &longitudes);
        vst1q_f64(out_lats + i, latitudes);
        vst1q_f64(out_lngs + i, longitudes);
    }

    geo_scalar_decode_interleaved_records(records + i, out_lats + i, out_lngs + i, count - i);
}

// =========================================================
// SIMD Haversine Distance - ARM64 NEON
// =========================================================

// The inputs used by Haversine are reduced to [-pi/2, pi/2]. The degree-17
// sine and degree-16 cosine polynomials keep global distance calculations
// accurate while avoiding scalar sin/cos calls in the vector loop.
static inline float64x2_t neon_sin_reduced(float64x2_t value)
{
    const float64x2_t pi = vdupq_n_f64(M_PI);
    const float64x2_t negative_pi = vdupq_n_f64(-M_PI);
    const float64x2_t half_pi = vdupq_n_f64(M_PI * 0.5);
    const float64x2_t negative_half_pi = vdupq_n_f64(-M_PI * 0.5);

    value = vbslq_f64(vcgtq_f64(value, half_pi), vsubq_f64(pi, value), value);
    value = vbslq_f64(vcltq_f64(value, negative_half_pi), vsubq_f64(negative_pi, value), value);

    float64x2_t squared = vmulq_f64(value, value);
    float64x2_t polynomial = vdupq_n_f64(1.0 / 355687428096000.0);

    polynomial = vfmaq_f64(vdupq_n_f64(-1.0 / 1307674368000.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(1.0 / 6227020800.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(-1.0 / 39916800.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(1.0 / 362880.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(-1.0 / 5040.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(1.0 / 120.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(-1.0 / 6.0), polynomial, squared);

    return vfmaq_f64(value, vmulq_f64(value, squared), polynomial);
}

static inline float64x2_t neon_cos_reduced(float64x2_t value)
{
    float64x2_t squared = vmulq_f64(value, value);
    float64x2_t polynomial = vdupq_n_f64(1.0 / 20922789888000.0);

    polynomial = vfmaq_f64(vdupq_n_f64(-1.0 / 87178291200.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(1.0 / 479001600.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(-1.0 / 3628800.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(1.0 / 40320.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(-1.0 / 720.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(1.0 / 24.0), polynomial, squared);
    polynomial = vfmaq_f64(vdupq_n_f64(-0.5), polynomial, squared);

    return vfmaq_f64(vdupq_n_f64(1.0), squared, polynomial);
}

static inline float64x2_t neon_haversine_a(const GeoSimdRadiusQuery *query,
                                           float64x2_t latitudes,
                                           float64x2_t longitudes)
{
    const float64x2_t degrees_to_radians = vdupq_n_f64(GEO_INTERNAL_DEG_TO_RAD);
    const float64x2_t half = vdupq_n_f64(0.5);
    const float64x2_t zero = vdupq_n_f64(0.0);
    const float64x2_t one = vdupq_n_f64(1.0);
    const float64x2_t center_latitude = vdupq_n_f64(query->center_latitude_radians);
    const float64x2_t center_longitude = vdupq_n_f64(query->center_longitude_radians);
    const float64x2_t cosine_center_latitude = vdupq_n_f64(query->cosine_center_latitude);

    float64x2_t latitude = vmulq_f64(latitudes, degrees_to_radians);
    float64x2_t longitude = vmulq_f64(longitudes, degrees_to_radians);
    float64x2_t sin_delta_latitude = neon_sin_reduced(vmulq_f64(vsubq_f64(latitude, center_latitude), half));
    float64x2_t sin_delta_longitude = neon_sin_reduced(vmulq_f64(vsubq_f64(longitude, center_longitude), half));
    float64x2_t cosine_product = vmulq_f64(cosine_center_latitude, neon_cos_reduced(latitude));
    float64x2_t haversine = vfmaq_f64(vmulq_f64(sin_delta_latitude, sin_delta_latitude),
                                      cosine_product,
                                      vmulq_f64(sin_delta_longitude, sin_delta_longitude));

    return vmaxq_f64(zero, vminq_f64(one, haversine));
}

void geo_simd_haversine_batch(double lat1, double lng1, const double *lats, const double *lngs, double *out_dist, size_t count)
{
    if (count == 0) {
        return;
    }

    GeoSimdRadiusQuery query = {
        .center_latitude_radians = lat1 * GEO_INTERNAL_DEG_TO_RAD,
        .center_longitude_radians = lng1 * GEO_INTERNAL_DEG_TO_RAD,
        .cosine_center_latitude = cos(lat1 * GEO_INTERNAL_DEG_TO_RAD),
        .sine_angular_radius = 0.0,
        .haversine_limit = 0.0,
    };

    size_t i = 0;

    for (; i + 1 < count; i += 2) {
        double haversine[2];

        vst1q_f64(haversine, neon_haversine_a(&query, vld1q_f64(&lats[i]), vld1q_f64(&lngs[i])));

        for (size_t lane = 0; lane < 2; ++lane) {
            out_dist[i + lane] = 2.0 * GEO_INTERNAL_EARTH_RADIUS_KM * atan2(sqrt(haversine[lane]), sqrt(1.0 - haversine[lane]));
        }
    }

    geo_scalar_haversine_batch(lat1, lng1, lats + i, lngs + i, out_dist + i, count - i);
}

// =========================================================
// SIMD Fast Distance (Equirectangular) - ARM64 NEON
// =========================================================

void geo_simd_fast_distance_batch(double lat1, double lng1, const double *lats, const double *lngs, double *out_dist, size_t count)
{
    if (count == 0) {
        return;
    }

    float64x2_t km_per_deg = vdupq_n_f64(GEO_INTERNAL_KM_PER_DEG);
    float64x2_t lat1_v = vdupq_n_f64(lat1);
    float64x2_t lng1_v = vdupq_n_f64(lng1);
    float64x2_t deg_to_rad = vdupq_n_f64(GEO_INTERNAL_DEG_TO_RAD);
    float64x2_t half = vdupq_n_f64(0.5);
    float64x2_t positive_180 = vdupq_n_f64(180.0);
    float64x2_t negative_180 = vdupq_n_f64(-180.0);
    float64x2_t full_circle = vdupq_n_f64(360.0);

    size_t i = 0;

    // Process 2 points at a time
    for (; i + 1 < count; i += 2) {
        float64x2_t lat2_v = vld1q_f64(&lats[i]);
        float64x2_t lng2_v = vld1q_f64(&lngs[i]);

        // Mid latitude for cos adjustment
        float64x2_t lat_mid = vmulq_f64(vaddq_f64(lat1_v, lat2_v), half);
        float64x2_t lat_mid_rad = vmulq_f64(lat_mid, deg_to_rad);
        float64x2_t cos_lat = neon_cos_reduced(lat_mid_rad);

        float64x2_t longitude_delta = vsubq_f64(lng2_v, lng1_v);

        longitude_delta = vbslq_f64(vcgtq_f64(longitude_delta, positive_180),
                                    vsubq_f64(longitude_delta, full_circle),
                                    longitude_delta);

        longitude_delta = vbslq_f64(vcltq_f64(longitude_delta, negative_180),
                                    vaddq_f64(longitude_delta, full_circle),
                                    longitude_delta);

        float64x2_t dlat = vmulq_f64(vsubq_f64(lat2_v, lat1_v), km_per_deg);
        float64x2_t dlng = vmulq_f64(vmulq_f64(longitude_delta, km_per_deg), cos_lat);

        // Distance = sqrt(dlat² + dlng²)
        float64x2_t dist_sq = vaddq_f64(vmulq_f64(dlat, dlat), vmulq_f64(dlng, dlng));
        float64x2_t dist = vsqrtq_f64(dist_sq);

        vst1q_f64(&out_dist[i], dist);
    }

    geo_scalar_fast_distance_batch(lat1, lng1, lats + i, lngs + i, out_dist + i, count - i);
}

// =========================================================
// SIMD Range Filtering - ARM64 NEON
// =========================================================

size_t geo_simd_filter_range(const uint64_t *z_codes, size_t count, uint64_t z_min, uint64_t z_max, uint8_t *out_mask)
{
    if (count == 0) {
        return 0;
    }

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

        out_mask[i] = mask0 ? 1 : 0;
        out_mask[i + 1] = mask1 ? 1 : 0;

        matches += (mask0 ? 1 : 0) + (mask1 ? 1 : 0);
    }

    matches += geo_scalar_filter_range(z_codes + i, count - i, z_min, z_max, out_mask + i);

    return matches;
}

// =========================================================
// SIMD Bounding Box Filtering - ARM64 NEON
// =========================================================

size_t geo_simd_filter_bbox(const double *lats, const double *lngs, size_t count,
                            double min_lat, double max_lat, double min_lng, double max_lng, uint8_t *out_mask)
{
    if (count == 0) {
        return 0;
    }

    bool wraps_antimeridian = min_lng > max_lng;

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
        uint64x2_t lng_ok = wraps_antimeridian ? vorrq_u64(lng_ge_min, lng_le_max) : vandq_u64(lng_ge_min, lng_le_max);

        // Combined
        uint64x2_t in_bbox = vandq_u64(lat_ok, lng_ok);

        uint64_t mask0 = vgetq_lane_u64(in_bbox, 0);
        uint64_t mask1 = vgetq_lane_u64(in_bbox, 1);

        out_mask[i] = mask0 ? 1 : 0;
        out_mask[i + 1] = mask1 ? 1 : 0;

        matches += (mask0 ? 1 : 0) + (mask1 ? 1 : 0);
    }

    matches += geo_scalar_filter_bbox(lats + i,
                                      lngs + i,
                                      count - i,
                                      min_lat,
                                      max_lat,
                                      min_lng,
                                      max_lng,
                                      out_mask + i);

    return matches;
}

// =========================================================
// SIMD Radius Filtering - ARM64 NEON
// =========================================================

void geo_simd_prepare_radius_query(double center_latitude,
                                   double center_longitude,
                                   double radius_km,
                                   GeoSimdRadiusQuery *query)
{
    geo_scalar_prepare_radius_query(center_latitude, center_longitude, radius_km, query);
}

static size_t filter_radius_neon(const double *lats,
                                 const double *lngs,
                                 size_t count,
                                 const GeoSimdRadiusQuery *query,
                                 uint8_t *out_mask,
                                 uint64_t *out_bits)
{
    if (count == 0) {
        return 0;
    }

    float64x2_t threshold = vdupq_n_f64(query->haversine_limit);
    size_t matches = 0;
    size_t i = 0;
    uint64_t packed_bits = 0;

    for (; i + 1 < count; i += 2) {
        float64x2_t haversine = neon_haversine_a(query,
                                                 vld1q_f64(&lats[i]),
                                                 vld1q_f64(&lngs[i]));
        uint64x2_t in_radius = vcleq_f64(haversine, threshold);
        uint8_t first_match = vgetq_lane_u64(in_radius, 0) != 0;
        uint8_t second_match = vgetq_lane_u64(in_radius, 1) != 0;

        if (out_mask) {
            out_mask[i] = first_match;
            out_mask[i + 1] = second_match;
        }

        if (out_bits) {
            packed_bits |= (uint64_t) first_match << (i & 63U);
            packed_bits |= (uint64_t) second_match << ((i + 1) & 63U);

            if ((i & 63U) == 62U) {
                out_bits[i >> 6] = packed_bits;
                packed_bits = 0;
            }
        }

        matches += first_match + second_match;
    }

    for (; i < count; ++i) {
        uint8_t matched = geo_scalar_radius_match(query, lats[i], lngs[i]);

        if (out_mask) {
            out_mask[i] = matched;
        }

        if (out_bits) {
            packed_bits |= (uint64_t) matched << (i & 63U);
        }

        matches += matched;
    }

    if (out_bits && (count & 63U)) {
        out_bits[count >> 6] = packed_bits;
    }

    return matches;
}

size_t geo_simd_filter_radius_prepared(const double *lats,
                                       const double *lngs,
                                       size_t count,
                                       const GeoSimdRadiusQuery *query,
                                       uint8_t *out_mask,
                                       uint64_t *out_bits)
{
    return filter_radius_neon(lats, lngs, count, query, out_mask, out_bits);
}

size_t geo_simd_filter_radius(const double *lats,
                              const double *lngs,
                              size_t count,
                              double center_lat,
                              double center_lng,
                              double radius_km,
                              uint8_t *out_mask)
{
    GeoSimdRadiusQuery query;

    geo_simd_prepare_radius_query(center_lat, center_lng, radius_km, &query);

    return geo_simd_filter_radius_prepared(lats, lngs, count, &query, out_mask, NULL);
}

size_t geo_simd_filter_radius_bits(const double *lats,
                                   const double *lngs,
                                   size_t count,
                                   double center_lat,
                                   double center_lng,
                                   double radius_km,
                                   uint64_t *out_bits)
{
    GeoSimdRadiusQuery query;

    geo_simd_prepare_radius_query(center_lat, center_lng, radius_km, &query);

    return geo_simd_filter_radius_prepared(lats, lngs, count, &query, NULL, out_bits);
}

static size_t filter_bbox_codes_neon(const uint64_t *z_codes,
                                     size_t count,
                                     double min_lat,
                                     double max_lat,
                                     double min_lng,
                                     double max_lng,
                                     uint8_t *out_mask,
                                     uint64_t *out_bits)
{
    uint32_t min_latitude = geo_internal_normalized_lat_lower(min_lat);
    uint32_t max_latitude = geo_internal_normalized_lat_upper(max_lat);
    uint32_t min_longitude = geo_internal_normalized_lng_lower(min_lng);
    uint32_t max_longitude = geo_internal_normalized_lng_upper(max_lng);
    bool wraps_antimeridian = min_lng > max_lng;
    uint64x2_t vector_minimum_latitude = vdupq_n_u64(min_latitude);
    uint64x2_t vector_maximum_latitude = vdupq_n_u64(max_latitude);
    uint64x2_t vector_minimum_longitude = vdupq_n_u64(min_longitude);
    uint64x2_t vector_maximum_longitude = vdupq_n_u64(max_longitude);
    size_t matches = 0;
    size_t i = 0;
    uint64_t packed_bits = 0;

    for (; i + 1 < count; i += 2) {
        uint64x2_t morton_codes = vld1q_u64(z_codes + i);
        uint64x2_t latitudes = neon_compact_bits(morton_codes);
        uint64x2_t longitudes = neon_compact_bits(vshrq_n_u64(morton_codes, 1));
        uint64x2_t latitude_matches = vandq_u64(vcgeq_u64(latitudes, vector_minimum_latitude),
                                                vcleq_u64(latitudes, vector_maximum_latitude));

        uint64x2_t above_minimum_longitude = vcgeq_u64(longitudes, vector_minimum_longitude);
        uint64x2_t below_maximum_longitude = vcleq_u64(longitudes, vector_maximum_longitude);
        uint64x2_t longitude_matches = wraps_antimeridian
                                           ? vorrq_u64(above_minimum_longitude, below_maximum_longitude)
                                           : vandq_u64(above_minimum_longitude, below_maximum_longitude);

        uint64x2_t in_bbox = vandq_u64(latitude_matches, longitude_matches);
        uint8_t first_match = vgetq_lane_u64(in_bbox, 0) != 0;
        uint8_t second_match = vgetq_lane_u64(in_bbox, 1) != 0;

        if (out_mask) {
            out_mask[i] = first_match;
            out_mask[i + 1] = second_match;
        }

        if (out_bits) {
            packed_bits |= (uint64_t) first_match << (i & 63U);
            packed_bits |= (uint64_t) second_match << ((i + 1) & 63U);

            if ((i & 63U) == 62U) {
                out_bits[i >> 6] = packed_bits;
                packed_bits = 0;
            }
        }

        matches += first_match + second_match;
    }

    for (; i < count; ++i) {
        uint32_t latitude = geo_internal_compact_bits(z_codes[i]);
        uint32_t longitude = geo_internal_compact_bits(z_codes[i] >> 1);
        bool longitude_matches = wraps_antimeridian ? longitude >= min_longitude || longitude <= max_longitude
                                                    : longitude >= min_longitude && longitude <= max_longitude;
        uint8_t matched = latitude >= min_latitude && latitude <= max_latitude && longitude_matches;

        if (out_mask) {
            out_mask[i] = matched;
        }

        if (out_bits) {
            packed_bits |= (uint64_t) matched << (i & 63U);
        }

        matches += matched;
    }

    if (out_bits && (count & 63U)) {
        out_bits[count >> 6] = packed_bits;
    }

    return matches;
}

size_t geo_simd_filter_bbox_codes(const uint64_t *z_codes,
                                  size_t count,
                                  double min_lat,
                                  double max_lat,
                                  double min_lng,
                                  double max_lng,
                                  uint8_t *out_mask)
{
    return filter_bbox_codes_neon(z_codes, count, min_lat, max_lat, min_lng, max_lng, out_mask, NULL);
}

size_t geo_simd_filter_bbox_codes_bits(const uint64_t *z_codes,
                                       size_t count,
                                       double min_lat,
                                       double max_lat,
                                       double min_lng,
                                       double max_lng,
                                       uint64_t *out_bits)
{
    return filter_bbox_codes_neon(z_codes, count, min_lat, max_lat, min_lng, max_lng, NULL, out_bits);
}

// =========================================================
// SIMD Bit Manipulation - ARM64 NEON
// =========================================================

void geo_simd_spread_bits_batch(const uint32_t *values, uint64_t *out, size_t count)
{
    size_t i = 0;

    for (; i + 1 < count; i += 2) {
        uint32x2_t packed_values = vld1_u32(&values[i]);

        vst1q_u64(&out[i], neon_spread_bits(vmovl_u32(packed_values)));
    }

    geo_scalar_spread_bits_batch(values + i, out + i, count - i);
}

void geo_simd_compact_bits_batch(const uint64_t *values, uint32_t *out, size_t count)
{
    size_t i = 0;

    for (; i + 1 < count; i += 2) {
        uint64x2_t packed_values = neon_compact_bits(vld1q_u64(&values[i]));

        vst1_u32(&out[i], vmovn_u64(packed_values));
    }

    geo_scalar_compact_bits_batch(values + i, out + i, count - i);
}

size_t geo_simd_exclude_id_bits(const GeoRecord *records,
                                size_t count,
                                uint64_t excluded_id,
                                uint64_t *candidate_bits)
{
    uint64x2_t excluded = vdupq_n_u64(excluded_id);
    size_t position = 0;

    for (; position + 2U <= count; position += 2U) {
        uint64x2x2_t interleaved = vld2q_u64((const uint64_t *) (records + position));
        uint64x2_t equal = vceqq_u64(interleaved.val[0], excluded);
        uint64_t equal_bits = (vgetq_lane_u64(equal, 0) >> 63U) |
                              ((vgetq_lane_u64(equal, 1) >> 63U) << 1U);

        candidate_bits[position >> 6U] &= ~(equal_bits << (position & 63U));
    }

    return geo_scalar_exclude_id_bits_from(records, position, count, excluded_id, candidate_bits);
}

// =========================================================
// Utility Functions - ARM64
// =========================================================

bool geo_simd_available(void)
{
    // Advanced SIMD is mandatory in the AArch64 architecture.
    return true;
}

const char *geo_simd_get_name(void)
{
    return GEO_SIMD_NAME;
}

size_t geo_simd_optimal_batch_size(void)
{
    return 256;
}

// =========================================================
// Benchmark Functions - ARM64
// =========================================================

GeoSimdBenchmark geo_simd_benchmark_encode(size_t count, int iterations)
{
    GeoSimdBenchmark result = { 0 };
    result.operation_name = "encode";
    result.operations = count * (size_t) (iterations > 0 ? iterations : 0);

    if (!count || iterations <= 0) {
        return result;
    }

    double *lats = (double *) malloc(count * sizeof(double));
    double *lngs = (double *) malloc(count * sizeof(double));
    uint64_t *z_codes = (uint64_t *) malloc(count * sizeof(uint64_t));

    if (!lats || !lngs || !z_codes) {
        free(lats);
        free(lngs);
        free(z_codes);
        return result;
    }

    // Initialize test data
    for (size_t i = 0; i < count; i++) {
        lats[i] = ((double) (i % 18000) / 100.0) - 90.0;
        lngs[i] = ((double) (i % 36000) / 100.0) - 180.0;
    }

    // Scalar benchmark
    double start = geo_get_time_ms();

    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            z_codes[i] = geo_encode(lats[i], lngs[i]);
        }

        g_benchmark_sink_u64 ^= z_codes[count / 2] ^ z_codes[count - 1];
    }
    result.scalar_time_ms = geo_get_time_ms() - start;

    // SIMD benchmark
    start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        geo_simd_encode_batch(lats, lngs, z_codes, count);
        g_benchmark_sink_u64 ^= z_codes[count / 2] ^ z_codes[count - 1];
    }
    result.simd_time_ms = geo_get_time_ms() - start;

    result.speedup = result.simd_time_ms ? result.scalar_time_ms / result.simd_time_ms : 0.0;

    free(lats);
    free(lngs);
    free(z_codes);
    return result;
}

GeoSimdBenchmark geo_simd_benchmark_decode(size_t count, int iterations)
{
    GeoSimdBenchmark result = { 0 };
    result.operation_name = "decode";
    result.operations = count * (size_t) (iterations > 0 ? iterations : 0);

    if (!count || iterations <= 0) {
        return result;
    }

    uint64_t *z_codes = (uint64_t *) malloc(count * sizeof(uint64_t));
    double *lats = (double *) malloc(count * sizeof(double));
    double *lngs = (double *) malloc(count * sizeof(double));

    if (!z_codes || !lats || !lngs) {
        free(z_codes);
        free(lats);
        free(lngs);
        return result;
    }

    // Initialize test data
    for (size_t i = 0; i < count; i++) {
        z_codes[i] = (uint64_t) i *12345ULL;
    }

    // Scalar benchmark
    double start = geo_get_time_ms();

    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            GeoPoint p = geo_decode(z_codes[i]);

            lats[i] = p.lat;
            lngs[i] = p.lng;
        }

        g_benchmark_sink_double += lats[count / 2] + lngs[count - 1];
    }
    result.scalar_time_ms = geo_get_time_ms() - start;

    // SIMD benchmark
    start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        geo_simd_decode_batch(z_codes, lats, lngs, count);
        g_benchmark_sink_double += lats[count / 2] + lngs[count - 1];
    }
    result.simd_time_ms = geo_get_time_ms() - start;

    result.speedup = result.simd_time_ms ? result.scalar_time_ms / result.simd_time_ms : 0.0;

    free(z_codes);
    free(lats);
    free(lngs);
    return result;
}

GeoSimdBenchmark geo_simd_benchmark_haversine(size_t count, int iterations)
{
    GeoSimdBenchmark result = { 0 };
    result.operation_name = "haversine";
    result.operations = count * (size_t) (iterations > 0 ? iterations : 0);

    if (!count || iterations <= 0) {
        return result;
    }

    double *lats = (double *) malloc(count * sizeof(double));
    double *lngs = (double *) malloc(count * sizeof(double));
    double *dists = (double *) malloc(count * sizeof(double));

    if (!lats || !lngs || !dists) {
        free(lats);
        free(lngs);
        free(dists);
        return result;
    }

    double center_lat = -23.5505;
    double center_lng = -46.6333;

    for (size_t i = 0; i < count; i++) {
        lats[i] = center_lat + ((double) (i % 1000) / 1000.0 - 0.5) * 10.0;
        lngs[i] = center_lng + ((double) (i % 1000) / 1000.0 - 0.5) * 10.0;
    }

    // Scalar benchmark
    double start = geo_get_time_ms();

    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            dists[i] = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);
        }

        g_benchmark_sink_double += dists[count / 2] + dists[count - 1];
    }
    result.scalar_time_ms = geo_get_time_ms() - start;

    // SIMD benchmark
    start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        geo_simd_haversine_batch(center_lat, center_lng, lats, lngs, dists, count);
        g_benchmark_sink_double += dists[count / 2] + dists[count - 1];
    }
    result.simd_time_ms = geo_get_time_ms() - start;

    result.speedup = result.simd_time_ms ? result.scalar_time_ms / result.simd_time_ms : 0.0;

    free(lats);
    free(lngs);
    free(dists);
    return result;
}

GeoSimdBenchmark geo_simd_benchmark_filter_radius(size_t count, int iterations)
{
    GeoSimdBenchmark result = { 0 };
    result.operation_name = "filter_radius";
    result.operations = count * (size_t) (iterations > 0 ? iterations : 0);

    if (!count || iterations <= 0) {
        return result;
    }

    double *lats = (double *) malloc(count * sizeof(double));
    double *lngs = (double *) malloc(count * sizeof(double));
    uint8_t *mask = (uint8_t *) malloc(count * sizeof(uint8_t));

    if (!lats || !lngs || !mask) {
        free(lats);
        free(lngs);
        free(mask);
        return result;
    }

    double center_lat = -23.5505;
    double center_lng = -46.6333;
    double radius_km = 50.0;

    for (size_t i = 0; i < count; i++) {
        lats[i] = ((double) (i % 18000) / 100.0) - 90.0;
        lngs[i] = ((double) (i % 36000) / 100.0) - 180.0;
    }

    // Scalar benchmark
    double start = geo_get_time_ms();

    for (int iter = 0; iter < iterations; iter++) {
        size_t matched = 0;

        for (size_t i = 0; i < count; i++) {
            double dist = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);

            mask[i] = (dist <= radius_km) ? 1 : 0;
            matched += mask[i];
        }

        g_benchmark_sink_u64 ^= matched;
    }
    result.scalar_time_ms = geo_get_time_ms() - start;

    // SIMD benchmark
    start = geo_get_time_ms();
    for (int iter = 0; iter < iterations; iter++) {
        g_benchmark_sink_u64 ^= geo_simd_filter_radius(lats, lngs, count, center_lat, center_lng, radius_km, mask);
    }
    result.simd_time_ms = geo_get_time_ms() - start;

    result.speedup = result.simd_time_ms ? result.scalar_time_ms / result.simd_time_ms : 0.0;

    free(lats);
    free(lngs);
    free(mask);
    return result;
}

#endif // GEO_SIMD_ARM64
