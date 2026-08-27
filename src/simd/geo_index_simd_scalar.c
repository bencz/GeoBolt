/*
 * GeoIndex SIMD Optimizations - Scalar Fallback Implementation
 * 
 * This file provides scalar implementations of the SIMD functions
 * for platforms where SIMD is not available or not detected.
 */

#include "geobolt/geo_index_simd.h"
#include "geobolt/geo_index.h"
#include "geo_index_internal.h"
#include "geo_index_simd_scalar_kernels.h"

#if !defined(GEO_SIMD_ARM64) && !defined(GEO_SIMD_X86_AVX2) && !defined(GEO_SIMD_X86_SSE4) && !defined(GEO_SIMD_X86_SSE2)

#include <stdlib.h>
#include <string.h>

static volatile uint64_t g_benchmark_sink_u64;
static volatile double g_benchmark_sink_double;

// =========================================================
// Scalar Batch Encoding
// =========================================================

void geo_simd_encode_batch(const double *lats, const double *lngs, uint64_t *out_z, size_t count)
{
    geo_scalar_encode_batch(lats, lngs, out_z, count);
}

void geo_simd_encode_records(const uint64_t *ids,
                             const double *lats,
                             const double *lngs,
                             GeoRecord *out_records,
                             size_t count)
{
    geo_scalar_encode_records(ids, lats, lngs, out_records, count);
}

bool geo_simd_validate_points(const double *lats, const double *lngs, size_t count)
{
    return geo_scalar_validate_points(lats, lngs, count);
}

// =========================================================
// Scalar Batch Decoding
// =========================================================

void geo_simd_decode_batch(const uint64_t *z_codes, double *out_lats, double *out_lngs, size_t count)
{
    geo_scalar_decode_batch(z_codes, out_lats, out_lngs, count);
}

void geo_simd_extract_interleaved_codes(const GeoRecord *records, uint64_t *out_codes, size_t count)
{
    geo_scalar_extract_interleaved_codes(records, out_codes, count);
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
    geo_scalar_decode_interleaved_records(records, out_lats, out_lngs, count);
}

// =========================================================
// Scalar Haversine Distance
// =========================================================

void geo_simd_haversine_batch(double lat1, double lng1, const double *lats, const double *lngs, double *out_dist, size_t count)
{
    geo_scalar_haversine_batch(lat1, lng1, lats, lngs, out_dist, count);
}

// =========================================================
// Scalar Fast Distance
// =========================================================

void geo_simd_fast_distance_batch(double lat1, double lng1, const double *lats, const double *lngs, double *out_dist, size_t count)
{
    geo_scalar_fast_distance_batch(lat1, lng1, lats, lngs, out_dist, count);
}

// =========================================================
// Scalar Range Filtering
// =========================================================

size_t geo_simd_filter_range(const uint64_t *z_codes, size_t count, uint64_t z_min, uint64_t z_max, uint8_t *out_mask)
{
    return geo_scalar_filter_range(z_codes, count, z_min, z_max, out_mask);
}

// =========================================================
// Scalar Bounding Box Filtering
// =========================================================

size_t geo_simd_filter_bbox(const double *lats, const double *lngs, size_t count,
                            double min_lat, double max_lat, double min_lng, double max_lng, uint8_t *out_mask)
{
    return geo_scalar_filter_bbox(lats, lngs, count, min_lat, max_lat, min_lng, max_lng, out_mask);
}

// =========================================================
// Scalar Radius Filtering
// =========================================================

void geo_simd_prepare_radius_query(double center_latitude,
                                   double center_longitude,
                                   double radius_km,
                                   GeoSimdRadiusQuery *query)
{
    geo_scalar_prepare_radius_query(center_latitude, center_longitude, radius_km, query);
}

size_t geo_simd_filter_radius_prepared(const double *lats,
                                       const double *lngs,
                                       size_t count,
                                       const GeoSimdRadiusQuery *query,
                                       uint8_t *out_mask,
                                       uint64_t *out_bits)
{
    return geo_scalar_filter_radius_prepared(lats, lngs, count, query, out_mask, out_bits);
}

size_t geo_simd_filter_radius(const double *lats, const double *lngs, size_t count,
                              double center_lat, double center_lng, double radius_km, uint8_t *out_mask)
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

size_t geo_simd_filter_bbox_codes(const uint64_t *z_codes,
                                  size_t count,
                                  double min_lat,
                                  double max_lat,
                                  double min_lng,
                                  double max_lng,
                                  uint8_t *out_mask)
{
    return geo_scalar_filter_bbox_codes(z_codes, count, min_lat, max_lat, min_lng, max_lng, out_mask, NULL);
}

size_t geo_simd_filter_bbox_codes_bits(const uint64_t *z_codes,
                                       size_t count,
                                       double min_lat,
                                       double max_lat,
                                       double min_lng,
                                       double max_lng,
                                       uint64_t *out_bits)
{
    return geo_scalar_filter_bbox_codes(z_codes, count, min_lat, max_lat, min_lng, max_lng, NULL, out_bits);
}

// =========================================================
// Scalar Bit Manipulation
// =========================================================

void geo_simd_spread_bits_batch(const uint32_t *values, uint64_t *out, size_t count)
{
    geo_scalar_spread_bits_batch(values, out, count);
}

void geo_simd_compact_bits_batch(const uint64_t *values, uint32_t *out, size_t count)
{
    geo_scalar_compact_bits_batch(values, out, count);
}

size_t geo_simd_exclude_id_bits(const GeoRecord *records,
                                size_t count,
                                uint64_t excluded_id,
                                uint64_t *candidate_bits)
{
    return geo_scalar_exclude_id_bits(records, count, excluded_id, candidate_bits);
}

// =========================================================
// Utility Functions - Scalar
// =========================================================

bool geo_simd_available(void)
{
    return false;
}

const char *geo_simd_get_name(void)
{
    return "Scalar (no SIMD)";
}

size_t geo_simd_optimal_batch_size(void)
{
    return 64;
}

// =========================================================
// Benchmark Functions - Scalar
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

    for (size_t i = 0; i < count; i++) {
        lats[i] = ((double) (i % 18000) / 100.0) - 90.0;
        lngs[i] = ((double) (i % 36000) / 100.0) - 180.0;
    }

    double start = geo_get_time_ms();

    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            z_codes[i] = geo_encode(lats[i], lngs[i]);
        }

        g_benchmark_sink_u64 ^= z_codes[count / 2] ^ z_codes[count - 1];
    }
    result.scalar_time_ms = geo_get_time_ms() - start;
    result.simd_time_ms = result.scalar_time_ms;        // Same as scalar
    result.speedup = 1.0;

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

    for (size_t i = 0; i < count; i++) {
        z_codes[i] = (uint64_t) i *12345ULL;
    }

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
    result.simd_time_ms = result.scalar_time_ms;
    result.speedup = 1.0;

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

    double start = geo_get_time_ms();

    for (int iter = 0; iter < iterations; iter++) {
        for (size_t i = 0; i < count; i++) {
            dists[i] = geo_haversine_km(center_lat, center_lng, lats[i], lngs[i]);
        }

        g_benchmark_sink_double += dists[count / 2] + dists[count - 1];
    }
    result.scalar_time_ms = geo_get_time_ms() - start;
    result.simd_time_ms = result.scalar_time_ms;
    result.speedup = 1.0;

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
    result.simd_time_ms = result.scalar_time_ms;
    result.speedup = 1.0;

    free(lats);
    free(lngs);
    free(mask);
    return result;
}

#endif // No SIMD
