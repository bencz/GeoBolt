#ifndef GEO_INDEX_SIMD_SCALAR_KERNELS_H
#define GEO_INDEX_SIMD_SCALAR_KERNELS_H

#include "geo_index.h"
#include "geo_index_internal.h"
#include "geo_index_simd.h"

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

static inline void geo_scalar_encode_batch(const double *latitudes,
                                           const double *longitudes,
                                           uint64_t *codes,
                                           size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        codes[i] = geo_internal_encode(latitudes[i], longitudes[i]);
    }
}

static inline void geo_scalar_encode_records(const uint64_t *ids,
                                             const double *latitudes,
                                             const double *longitudes,
                                             GeoRecord *records,
                                             size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        records[i] = (GeoRecord) {
            .id = ids[i],
            .z = geo_internal_encode(latitudes[i], longitudes[i]),
        };
    }
}

static inline bool geo_scalar_validate_points(const double *latitudes, const double *longitudes, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        if (!geo_is_valid_point(latitudes[i], longitudes[i])) {
            return false;
        }
    }

    return true;
}

static inline void geo_scalar_decode_batch(const uint64_t *codes,
                                           double *latitudes,
                                           double *longitudes,
                                           size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        geo_internal_decode(codes[i], latitudes + i, longitudes + i);
    }
}

static inline void geo_scalar_extract_interleaved_codes(const GeoRecord *records, uint64_t *codes, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        codes[i] = records[i].z;
    }
}

static inline void geo_scalar_decode_interleaved_records(const GeoRecord *records,
                                                         double *latitudes,
                                                         double *longitudes,
                                                         size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        geo_internal_decode(records[i].z, latitudes + i, longitudes + i);
    }
}

static inline void geo_scalar_haversine_batch(double center_latitude,
                                              double center_longitude,
                                              const double *latitudes,
                                              const double *longitudes,
                                              double *distances,
                                              size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        distances[i] = geo_haversine_km(center_latitude, center_longitude, latitudes[i], longitudes[i]);
    }
}

static inline void geo_scalar_fast_distance_batch(double center_latitude,
                                                  double center_longitude,
                                                  const double *latitudes,
                                                  const double *longitudes,
                                                  double *distances,
                                                  size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        distances[i] = geo_fast_distance_km(center_latitude, center_longitude, latitudes[i], longitudes[i]);
    }
}

static inline size_t geo_scalar_filter_range(const uint64_t *codes,
                                             size_t count,
                                             uint64_t minimum,
                                             uint64_t maximum,
                                             uint8_t *mask)
{
    size_t matches = 0;

    for (size_t i = 0; i < count; ++i) {
        mask[i] = codes[i] >= minimum && codes[i] <= maximum;
        matches += mask[i];
    }

    return matches;
}

static inline size_t geo_scalar_filter_bbox(const double *latitudes,
                                            const double *longitudes,
                                            size_t count,
                                            double minimum_latitude,
                                            double maximum_latitude,
                                            double minimum_longitude,
                                            double maximum_longitude,
                                            uint8_t *mask)
{
    bool wraps = minimum_longitude > maximum_longitude;
    size_t matches = 0;

    for (size_t i = 0; i < count; ++i) {
        bool longitude_matches = wraps ? longitudes[i] >= minimum_longitude || longitudes[i] <= maximum_longitude
                                       : longitudes[i] >= minimum_longitude && longitudes[i] <= maximum_longitude;

        mask[i] = latitudes[i] >= minimum_latitude && latitudes[i] <= maximum_latitude && longitude_matches;
        matches += mask[i];
    }

    return matches;
}

static inline void geo_scalar_prepare_radius_query(double center_latitude,
                                                   double center_longitude,
                                                   double radius_km,
                                                   GeoSimdRadiusQuery *query)
{
    double angular = fmin(M_PI, radius_km / GEO_INTERNAL_EARTH_RADIUS_KM);
    double half_angle_sine = sin(angular * 0.5);

    query->center_latitude_radians = center_latitude * GEO_INTERNAL_DEG_TO_RAD;
    query->center_longitude_radians = center_longitude * GEO_INTERNAL_DEG_TO_RAD;
    query->cosine_center_latitude = cos(query->center_latitude_radians);
    query->sine_angular_radius = sin(angular);
    query->haversine_limit = half_angle_sine * half_angle_sine;
}

static inline bool geo_scalar_radius_match(const GeoSimdRadiusQuery *query, double latitude, double longitude)
{
    double latitude_radians = latitude * GEO_INTERNAL_DEG_TO_RAD;
    double longitude_radians = longitude * GEO_INTERNAL_DEG_TO_RAD;
    double sin_delta_latitude = sin((latitude_radians - query->center_latitude_radians) * 0.5);
    double sin_delta_longitude = sin((longitude_radians - query->center_longitude_radians) * 0.5);
    double haversine = sin_delta_latitude * sin_delta_latitude +
                       query->cosine_center_latitude * cos(latitude_radians) * sin_delta_longitude * sin_delta_longitude;

    return haversine <= query->haversine_limit;
}

static inline size_t geo_scalar_filter_radius_prepared(const double *latitudes,
                                                       const double *longitudes,
                                                       size_t count,
                                                       const GeoSimdRadiusQuery *query,
                                                       uint8_t *mask,
                                                       uint64_t *bit_mask)
{
    size_t matches = 0;
    uint64_t packed_bits = 0;

    for (size_t i = 0; i < count; ++i) {
        uint8_t matched = geo_scalar_radius_match(query, latitudes[i], longitudes[i]);

        if (mask) {
            mask[i] = matched;
        }

        if (bit_mask) {
            packed_bits |= (uint64_t) matched << (i & 63U);

            if ((i & 63U) == 63U) {
                bit_mask[i >> 6U] = packed_bits;
                packed_bits = 0;
            }
        }

        matches += matched;
    }

    if (bit_mask && (count & 63U)) {
        bit_mask[count >> 6U] = packed_bits;
    }

    return matches;
}

static inline size_t geo_scalar_filter_bbox_codes(const uint64_t *codes,
                                                  size_t count,
                                                  double minimum_latitude,
                                                  double maximum_latitude,
                                                  double minimum_longitude,
                                                  double maximum_longitude,
                                                  uint8_t *mask,
                                                  uint64_t *bit_mask)
{
    uint32_t minimum_latitude_bits = geo_internal_normalized_lat_lower(minimum_latitude);
    uint32_t maximum_latitude_bits = geo_internal_normalized_lat_upper(maximum_latitude);
    uint32_t minimum_longitude_bits = geo_internal_normalized_lng_lower(minimum_longitude);
    uint32_t maximum_longitude_bits = geo_internal_normalized_lng_upper(maximum_longitude);
    bool wraps = minimum_longitude > maximum_longitude;
    size_t matches = 0;
    uint64_t packed_bits = 0;

    for (size_t i = 0; i < count; ++i) {
        uint32_t latitude = geo_internal_compact_bits(codes[i]);
        uint32_t longitude = geo_internal_compact_bits(codes[i] >> 1U);
        bool longitude_matches = wraps ? longitude >= minimum_longitude_bits || longitude <= maximum_longitude_bits
                                       : longitude >= minimum_longitude_bits && longitude <= maximum_longitude_bits;
        uint8_t matched = latitude >= minimum_latitude_bits && latitude <= maximum_latitude_bits && longitude_matches;

        if (mask) {
            mask[i] = matched;
        }

        if (bit_mask) {
            packed_bits |= (uint64_t) matched << (i & 63U);

            if ((i & 63U) == 63U) {
                bit_mask[i >> 6U] = packed_bits;
                packed_bits = 0;
            }
        }

        matches += matched;
    }

    if (bit_mask && (count & 63U)) {
        bit_mask[count >> 6U] = packed_bits;
    }

    return matches;
}

static inline void geo_scalar_spread_bits_batch(const uint32_t *values, uint64_t *output, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        output[i] = geo_internal_spread_bits(values[i]);
    }
}

static inline void geo_scalar_compact_bits_batch(const uint64_t *values, uint32_t *output, size_t count)
{
    for (size_t i = 0; i < count; ++i) {
        output[i] = geo_internal_compact_bits(values[i]);
    }
}

static inline size_t geo_scalar_exclude_id_bits(const GeoRecord *records,
                                                size_t count,
                                                uint64_t excluded_id,
                                                uint64_t *candidate_bits)
{
    size_t word_count = geo_internal_bit_word_count(count);
    size_t matches = 0;

    for (size_t word = 0; word < word_count; ++word) {
        uint64_t candidates = candidate_bits[word];

        while (candidates) {
            unsigned lane = (unsigned) __builtin_ctzll(candidates);
            size_t record_index = word * 64U + lane;

            if (records[record_index].id == excluded_id) {
                candidate_bits[word] &= ~(UINT64_C(1) << lane);
            }

            candidates &= candidates - 1U;
        }

        matches += (size_t) __builtin_popcountll(candidate_bits[word]);
    }

    return matches;
}

static inline size_t geo_scalar_exclude_id_bits_from(const GeoRecord *records,
                                                     size_t first,
                                                     size_t count,
                                                     uint64_t excluded_id,
                                                     uint64_t *candidate_bits)
{
    for (size_t i = first; i < count; ++i) {
        size_t word = i >> 6U;
        uint64_t bit = UINT64_C(1) << (i & 63U);

        if ((candidate_bits[word] & bit) && records[i].id == excluded_id) {
            candidate_bits[word] &= ~bit;
        }
    }

    size_t matches = 0;
    size_t word_count = geo_internal_bit_word_count(count);

    for (size_t word = 0; word < word_count; ++word) {
        matches += (size_t) __builtin_popcountll(candidate_bits[word]);
    }

    return matches;
}

#endif
