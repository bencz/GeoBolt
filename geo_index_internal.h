/*
 * GeoIndex Internal Utilities
 * 
 * Shared inline functions and constants used by both scalar and SIMD
 * implementations. This file should only be included by .c implementation
 * files, not by external users.
 */

#ifndef GEO_INDEX_INTERNAL_H
#define GEO_INDEX_INTERNAL_H

#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

// =========================================================
// Shared Constants
// =========================================================

#define GEO_INTERNAL_DEG_TO_RAD  (M_PI / 180.0)
#define GEO_INTERNAL_RAD_TO_DEG  (180.0 / M_PI)
#define GEO_INTERNAL_KM_PER_DEG  111.32
#define GEO_INTERNAL_EARTH_RADIUS_KM 6371.0088

// Normalization constants for coordinate encoding
#define GEO_INTERNAL_COORD_MAX   4294967295.0
#define GEO_INTERNAL_LAT_RANGE   180.0
#define GEO_INTERNAL_LNG_RANGE   360.0

#define GEO_INTERNAL_LAT_SCALE   (GEO_INTERNAL_COORD_MAX / GEO_INTERNAL_LAT_RANGE)
#define GEO_INTERNAL_LNG_SCALE   (GEO_INTERNAL_COORD_MAX / GEO_INTERNAL_LNG_RANGE)
#define GEO_INTERNAL_LAT_OFFSET  90.0
#define GEO_INTERNAL_LNG_OFFSET  180.0

#define GEO_INTERNAL_LAT_DENORM_SCALE (GEO_INTERNAL_LAT_RANGE / GEO_INTERNAL_COORD_MAX)
#define GEO_INTERNAL_LNG_DENORM_SCALE (GEO_INTERNAL_LNG_RANGE / GEO_INTERNAL_COORD_MAX)

// =========================================================
// Morton Code Bit Manipulation (Inline)
// =========================================================

/*
 * Spread bits of a 32-bit value into a 64-bit value.
 * Each bit is spread to occupy every other bit position.
 * Used for Morton/Z-order encoding.
 */
static inline uint64_t geo_internal_spread_bits(uint32_t v) {
    uint64_t x = v;
    x = (x | (x << 32)) & 0x00000000FFFFFFFFULL;
    x = (x | (x << 16)) & 0x0000FFFF0000FFFFULL;
    x = (x | (x << 8))  & 0x00FF00FF00FF00FFULL;
    x = (x | (x << 4))  & 0x0F0F0F0F0F0F0F0FULL;
    x = (x | (x << 2))  & 0x3333333333333333ULL;
    x = (x | (x << 1))  & 0x5555555555555555ULL;
    return x;
}

/*
 * Compact bits of a 64-bit value into a 32-bit value.
 * Extracts bits at every other position.
 * Used for Morton/Z-order decoding.
 */
static inline uint32_t geo_internal_compact_bits(uint64_t v) {
    uint64_t x = v & 0x5555555555555555ULL;
    x = (x | (x >> 1))  & 0x3333333333333333ULL;
    x = (x | (x >> 2))  & 0x0F0F0F0F0F0F0F0FULL;
    x = (x | (x >> 4))  & 0x00FF00FF00FF00FFULL;
    x = (x | (x >> 8))  & 0x0000FFFF0000FFFFULL;
    x = (x | (x >> 16)) & 0x00000000FFFFFFFFULL;
    return (uint32_t)x;
}

// =========================================================
// Coordinate Normalization (Inline)
// =========================================================

/*
 * Normalize latitude to [0, UINT32_MAX] range.
 */
static inline uint32_t geo_internal_normalize_lat(double lat) {
    if (lat < -90.0) lat = -90.0;
    if (lat > 90.0) lat = 90.0;
    double normalized = (lat + GEO_INTERNAL_LAT_OFFSET) * GEO_INTERNAL_LAT_SCALE;
    return (uint32_t)normalized;
}

/*
 * Normalize longitude to [0, UINT32_MAX] range.
 */
static inline uint32_t geo_internal_normalize_lng(double lng) {
    if (lng < -180.0) lng = -180.0;
    if (lng > 180.0) lng = 180.0;
    double normalized = (lng + GEO_INTERNAL_LNG_OFFSET) * GEO_INTERNAL_LNG_SCALE;
    return (uint32_t)normalized;
}

/*
 * Denormalize uint32 to latitude.
 */
static inline double geo_internal_denormalize_lat(uint32_t v) {
    return ((double)v * GEO_INTERNAL_LAT_DENORM_SCALE) - GEO_INTERNAL_LAT_OFFSET;
}

/*
 * Denormalize uint32 to longitude.
 */
static inline double geo_internal_denormalize_lng(uint32_t v) {
    return ((double)v * GEO_INTERNAL_LNG_DENORM_SCALE) - GEO_INTERNAL_LNG_OFFSET;
}

// =========================================================
// Fast Encode/Decode (Inline)
// =========================================================

/*
 * Encode lat/lng to Morton code.
 */
static inline uint64_t geo_internal_encode(double lat, double lng) {
    uint32_t lat_norm = geo_internal_normalize_lat(lat);
    uint32_t lng_norm = geo_internal_normalize_lng(lng);
    return geo_internal_spread_bits(lat_norm) | (geo_internal_spread_bits(lng_norm) << 1);
}

/*
 * Decode Morton code to lat/lng.
 */
static inline void geo_internal_decode(uint64_t z, double *out_lat, double *out_lng) {
    *out_lat = geo_internal_denormalize_lat(geo_internal_compact_bits(z));
    *out_lng = geo_internal_denormalize_lng(geo_internal_compact_bits(z >> 1));
}

// =========================================================
// Min/Max/Clamp Utilities (Inline)
// =========================================================

static inline double geo_internal_min(double a, double b) {
    return (a < b) ? a : b;
}

static inline double geo_internal_max(double a, double b) {
    return (a > b) ? a : b;
}

static inline double geo_internal_clamp(double v, double lo, double hi) {
    return geo_internal_min(geo_internal_max(v, lo), hi);
}

#endif // GEO_INDEX_INTERNAL_H
