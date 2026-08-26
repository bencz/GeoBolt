#include "geo_index_simd.h"
#include "geo_index.h"
#include "geo_index_internal.h"
#include "geo_index_simd_scalar_kernels.h"

#if defined(GEO_SIMD_X86_AVX2) || defined(GEO_SIMD_X86_SSE4) || defined(GEO_SIMD_X86_SSE2)

#include <immintrin.h>
#include <math.h>
#include <stdatomic.h>
#include <stdlib.h>

#if defined(__GNUC__) || defined(__clang__)
#define GEO_TARGET_AVX2 __attribute__((target("avx2,fma")))
#define GEO_TARGET_AVX512 __attribute__((target("avx512f,avx512dq,avx512vl,avx2,fma")))
#define GEO_CONSTRUCTOR __attribute__((constructor))
#else
#define GEO_TARGET_AVX2
#define GEO_TARGET_AVX512
#define GEO_CONSTRUCTOR
#endif

GEO_TARGET_AVX2 static inline __m256i spread4(__m256i x);

typedef void (*GeoEncodeKernel)(const double *, const double *, uint64_t *, size_t);
typedef void (*GeoEncodeRecordsKernel)(const uint64_t *, const double *, const double *, GeoRecord *, size_t);
typedef bool (*GeoValidateKernel)(const double *, const double *, size_t);
typedef void (*GeoDecodeKernel)(const uint64_t *, double *, double *, size_t);
typedef void (*GeoExtractKernel)(const GeoRecord *, uint64_t *, size_t);
typedef void (*GeoDecodeRecordsKernel)(const GeoRecord *, double *, double *, size_t);
typedef void (*GeoDistanceKernel)(double, double, const double *, const double *, double *, size_t);
typedef size_t (*GeoRangeKernel)(const uint64_t *, size_t, uint64_t, uint64_t, uint8_t *);
typedef size_t (*GeoBboxKernel)(const double *, const double *, size_t, double, double, double, double, uint8_t *);
typedef size_t (*GeoRadiusKernel)(const double *, const double *, size_t, const GeoSimdRadiusQuery *, uint8_t *, uint64_t *);
typedef size_t (*GeoBboxCodesKernel)(const uint64_t *, size_t, double, double, double, double, uint8_t *, uint64_t *);
typedef void (*GeoSpreadKernel)(const uint32_t *, uint64_t *, size_t);
typedef void (*GeoCompactKernel)(const uint64_t *, uint32_t *, size_t);
typedef size_t (*GeoExcludeIdKernel)(const GeoRecord *, size_t, uint64_t, uint64_t *);

typedef struct {
    GeoEncodeKernel encode;
    GeoEncodeRecordsKernel encode_records;
    GeoValidateKernel validate;
    GeoDecodeKernel decode_wide;
    GeoDecodeKernel decode_narrow;
    GeoExtractKernel extract_wide;
    GeoExtractKernel extract_narrow;
    GeoDecodeRecordsKernel decode_records_narrow;
    GeoDistanceKernel haversine;
    GeoDistanceKernel fast_distance;
    GeoRangeKernel filter_range;
    GeoBboxKernel filter_bbox;
    GeoRadiusKernel filter_radius;
    GeoBboxCodesKernel filter_bbox_codes;
    GeoSpreadKernel spread;
    GeoCompactKernel compact;
    GeoExcludeIdKernel exclude_id;
    const char *name;
    size_t optimal_batch_size;
    bool simd_available;
} GeoX86Backend;

static const GeoX86Backend g_scalar_backend;

static _Atomic(const GeoX86Backend *) g_backend = ATOMIC_VAR_INIT(&g_scalar_backend);
static volatile uint64_t g_benchmark_sink_u64;
static volatile double g_benchmark_sink_double;

static inline const GeoX86Backend *x86_backend(void)
{
    return atomic_load_explicit(&g_backend, memory_order_relaxed);
}

// =============================================================================
// Numerically stable vector trigonometric kernels
// =============================================================================

GEO_TARGET_AVX2 static inline __m256d sin_poly(__m256d x)
{
    const __m256d pi = _mm256_set1_pd(M_PI);
    const __m256d half_pi = _mm256_set1_pd(M_PI * 0.5);
    const __m256d neg_half_pi = _mm256_set1_pd(-M_PI * 0.5);

    __m256d high = _mm256_cmp_pd(x, half_pi, _CMP_GT_OQ);
    __m256d low = _mm256_cmp_pd(x, neg_half_pi, _CMP_LT_OQ);

    x = _mm256_blendv_pd(x, _mm256_sub_pd(pi, x), high);
    x = _mm256_blendv_pd(x, _mm256_sub_pd(_mm256_set1_pd(-M_PI), x), low);

    __m256d x2 = _mm256_mul_pd(x, x);
    __m256d p = _mm256_set1_pd(1.0 / 355687428096000.0);

    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(-1.0 / 1307674368000.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(1.0 / 6227020800.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(-1.0 / 39916800.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(1.0 / 362880.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(-1.0 / 5040.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(1.0 / 120.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(-1.0 / 6.0));

    return _mm256_fmadd_pd(_mm256_mul_pd(x, x2), p, x);
}

GEO_TARGET_AVX2 static inline __m256d cos_poly(__m256d x)
{
    __m256d x2 = _mm256_mul_pd(x, x);
    __m256d p = _mm256_set1_pd(1.0 / 20922789888000.0);

    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(-1.0 / 87178291200.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(1.0 / 479001600.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(-1.0 / 3628800.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(1.0 / 40320.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(-1.0 / 720.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(1.0 / 24.0));
    p = _mm256_fmadd_pd(p, x2, _mm256_set1_pd(-0.5));

    return _mm256_fmadd_pd(x2, p, _mm256_set1_pd(1.0));
}

GEO_TARGET_AVX2 static inline __m256d haversine_a4(const GeoSimdRadiusQuery *query, __m256d lats, __m256d lngs)
{
    __m256d radians_per_degree = _mm256_set1_pd(GEO_INTERNAL_DEG_TO_RAD);
    __m256d half = _mm256_set1_pd(0.5);
    __m256d center_latitude = _mm256_set1_pd(query->center_latitude_radians);
    __m256d center_longitude = _mm256_set1_pd(query->center_longitude_radians);
    __m256d cosine_center_latitude = _mm256_set1_pd(query->cosine_center_latitude);
    __m256d point_latitudes = _mm256_mul_pd(lats, radians_per_degree);
    __m256d point_longitudes = _mm256_mul_pd(lngs, radians_per_degree);

    __m256d sin_delta_latitude = sin_poly(_mm256_mul_pd(_mm256_sub_pd(point_latitudes, center_latitude), half));
    __m256d sin_delta_longitude = sin_poly(_mm256_mul_pd(_mm256_sub_pd(point_longitudes, center_longitude), half));

    __m256d haversine = _mm256_fmadd_pd(_mm256_mul_pd(cosine_center_latitude, cos_poly(point_latitudes)),
                                        _mm256_mul_pd(sin_delta_longitude, sin_delta_longitude),
                                        _mm256_mul_pd(sin_delta_latitude, sin_delta_latitude));

    return _mm256_max_pd(_mm256_setzero_pd(), _mm256_min_pd(_mm256_set1_pd(1.0), haversine));
}

// =============================================================================
// Morton batch encoding and decoding
// =============================================================================

GEO_TARGET_AVX2 static inline __m256i normalized_doubles_to_u64(__m256d values)
{
    __m256d unsigned_bias = _mm256_set1_pd(2147483648.0);
    __m256d uses_high_bit = _mm256_cmp_pd(values, unsigned_bias, _CMP_GE_OQ);
    __m256d adjusted = _mm256_blendv_pd(values, _mm256_sub_pd(values, unsigned_bias), uses_high_bit);
    __m128i packed = _mm256_cvttpd_epi32(adjusted);
    int high_bit_mask = _mm256_movemask_pd(uses_high_bit);

    __m128i high_bits = _mm_set_epi32((high_bit_mask & 8) ? INT32_MIN : 0,
                                      (high_bit_mask & 4) ? INT32_MIN : 0,
                                      (high_bit_mask & 2) ? INT32_MIN : 0,
                                      (high_bit_mask & 1) ? INT32_MIN : 0);

    return _mm256_cvtepu32_epi64(_mm_add_epi32(packed, high_bits));
}

GEO_TARGET_AVX2 static void encode_avx2(const double *lats, const double *lngs, uint64_t *out, size_t count)
{
    __m256d lat_scale = _mm256_set1_pd(GEO_INTERNAL_LAT_SCALE);
    __m256d lng_scale = _mm256_set1_pd(GEO_INTERNAL_LNG_SCALE);
    __m256d lat_offset = _mm256_set1_pd(GEO_INTERNAL_LAT_OFFSET);
    __m256d lng_offset = _mm256_set1_pd(GEO_INTERNAL_LNG_OFFSET);
    __m256d zero = _mm256_setzero_pd();
    __m256d maximum = _mm256_set1_pd(GEO_INTERNAL_COORD_MAX);
    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        __m256d normalized_latitudes = _mm256_mul_pd(_mm256_add_pd(_mm256_loadu_pd(lats + i), lat_offset), lat_scale);
        __m256d normalized_longitudes = _mm256_mul_pd(_mm256_add_pd(_mm256_loadu_pd(lngs + i), lng_offset), lng_scale);

        normalized_latitudes = _mm256_max_pd(zero, _mm256_min_pd(maximum, normalized_latitudes));
        normalized_longitudes = _mm256_max_pd(zero, _mm256_min_pd(maximum, normalized_longitudes));

        __m256i latitude_bits = spread4(normalized_doubles_to_u64(normalized_latitudes));
        __m256i longitude_bits = spread4(normalized_doubles_to_u64(normalized_longitudes));
        __m256i morton_codes = _mm256_or_si256(latitude_bits, _mm256_slli_epi64(longitude_bits, 1));

        _mm256_storeu_si256((__m256i *)(out + i), morton_codes);
    }

    for (; i < count; ++i) {
        out[i] = geo_internal_encode(lats[i], lngs[i]);
    }
}

GEO_TARGET_AVX2 static void encode_records_avx2(const uint64_t *ids,
                                               const double *lats,
                                               const double *lngs,
                                               GeoRecord *records,
                                               size_t count)
{
    __m256d lat_scale = _mm256_set1_pd(GEO_INTERNAL_LAT_SCALE);
    __m256d lng_scale = _mm256_set1_pd(GEO_INTERNAL_LNG_SCALE);
    __m256d lat_offset = _mm256_set1_pd(GEO_INTERNAL_LAT_OFFSET);
    __m256d lng_offset = _mm256_set1_pd(GEO_INTERNAL_LNG_OFFSET);
    __m256d zero = _mm256_setzero_pd();
    __m256d maximum = _mm256_set1_pd(GEO_INTERNAL_COORD_MAX);
    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        __m256d normalized_latitudes = _mm256_mul_pd(_mm256_add_pd(_mm256_loadu_pd(lats + i), lat_offset), lat_scale);
        __m256d normalized_longitudes = _mm256_mul_pd(_mm256_add_pd(_mm256_loadu_pd(lngs + i), lng_offset), lng_scale);

        normalized_latitudes = _mm256_max_pd(zero, _mm256_min_pd(maximum, normalized_latitudes));
        normalized_longitudes = _mm256_max_pd(zero, _mm256_min_pd(maximum, normalized_longitudes));

        __m256i latitude_bits = spread4(normalized_doubles_to_u64(normalized_latitudes));
        __m256i longitude_bits = spread4(normalized_doubles_to_u64(normalized_longitudes));
        __m256i morton_codes = _mm256_or_si256(latitude_bits, _mm256_slli_epi64(longitude_bits, 1));
        __m256i record_ids = _mm256_loadu_si256((const __m256i *) (ids + i));
        __m256i even_records = _mm256_unpacklo_epi64(record_ids, morton_codes);
        __m256i odd_records = _mm256_unpackhi_epi64(record_ids, morton_codes);
        __m256i first_records = _mm256_permute2x128_si256(even_records, odd_records, 0x20);
        __m256i second_records = _mm256_permute2x128_si256(even_records, odd_records, 0x31);

        _mm256_storeu_si256((__m256i *) (records + i), first_records);
        _mm256_storeu_si256((__m256i *) (records + i + 2), second_records);
    }

    for (; i < count; ++i) {
        records[i] = (GeoRecord) {
            .id = ids[i],
            .z = geo_internal_encode(lats[i], lngs[i]),
        };
    }
}

void geo_simd_encode_batch(const double *lats, const double *lngs, uint64_t *out, size_t count)
{
    x86_backend()->encode(lats, lngs, out, count);
}

void geo_simd_encode_records(const uint64_t *ids,
                             const double *lats,
                             const double *lngs,
                             GeoRecord *records,
                             size_t count)
{
    x86_backend()->encode_records(ids, lats, lngs, records, count);
}

GEO_TARGET_AVX2 static bool validate_points_avx2(const double *lats, const double *lngs, size_t count)
{
    __m256d minimum_latitude = _mm256_set1_pd(GEO_MIN_LAT);
    __m256d maximum_latitude = _mm256_set1_pd(GEO_MAX_LAT);
    __m256d minimum_longitude = _mm256_set1_pd(GEO_MIN_LNG);
    __m256d maximum_longitude = _mm256_set1_pd(GEO_MAX_LNG);
    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        __m256d latitudes = _mm256_loadu_pd(lats + i);
        __m256d longitudes = _mm256_loadu_pd(lngs + i);
        __m256d valid_latitudes = _mm256_and_pd(_mm256_cmp_pd(latitudes, minimum_latitude, _CMP_GE_OQ),
                                                _mm256_cmp_pd(latitudes, maximum_latitude, _CMP_LE_OQ));

        __m256d valid_longitudes = _mm256_and_pd(_mm256_cmp_pd(longitudes, minimum_longitude, _CMP_GE_OQ),
                                                 _mm256_cmp_pd(longitudes, maximum_longitude, _CMP_LE_OQ));

        if (_mm256_movemask_pd(_mm256_and_pd(valid_latitudes, valid_longitudes)) != 15) {
            return false;
        }
    }

    for (; i < count; ++i) {
        if (!geo_is_valid_point(lats[i], lngs[i])) {
            return false;
        }
    }

    return true;
}

bool geo_simd_validate_points(const double *lats, const double *lngs, size_t count)
{
    return x86_backend()->validate(lats, lngs, count);
}

GEO_TARGET_AVX2 static inline __m256i compact4(__m256i x)
{
    x = _mm256_and_si256(x, _mm256_set1_epi64x(0x5555555555555555LL));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 1)), _mm256_set1_epi64x(0x3333333333333333LL));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 2)), _mm256_set1_epi64x(0x0F0F0F0F0F0F0F0FLL));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 4)), _mm256_set1_epi64x(0x00FF00FF00FF00FFLL));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 8)), _mm256_set1_epi64x(0x0000FFFF0000FFFFLL));

    return _mm256_and_si256(_mm256_or_si256(x, _mm256_srli_epi64(x, 16)), _mm256_set1_epi64x(0xFFFFFFFFLL));
}

GEO_TARGET_AVX512 static inline __m512i compact8(__m512i x)
{
    x = _mm512_and_si512(x, _mm512_set1_epi64(0x5555555555555555LL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 1)), _mm512_set1_epi64(0x3333333333333333LL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 2)), _mm512_set1_epi64(0x0F0F0F0F0F0F0F0FLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 4)), _mm512_set1_epi64(0x00FF00FF00FF00FFLL));
    x = _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 8)), _mm512_set1_epi64(0x0000FFFF0000FFFFLL));

    return _mm512_and_si512(_mm512_or_si512(x, _mm512_srli_epi64(x, 16)), _mm512_set1_epi64(0xFFFFFFFFLL));
}

GEO_TARGET_AVX2 static inline void decode4_values(__m256i morton_codes, __m256d *latitudes, __m256d *longitudes)
{
    __m256i normalized_latitudes = compact4(morton_codes);
    __m256i normalized_longitudes = compact4(_mm256_srli_epi64(morton_codes, 1));
    __m128i latitude_low = _mm_shuffle_epi32(_mm256_castsi256_si128(normalized_latitudes), _MM_SHUFFLE(2, 0, 2, 0));
    __m128i latitude_high = _mm_shuffle_epi32(_mm256_extracti128_si256(normalized_latitudes, 1), _MM_SHUFFLE(2, 0, 2, 0));
    __m128i longitude_low = _mm_shuffle_epi32(_mm256_castsi256_si128(normalized_longitudes), _MM_SHUFFLE(2, 0, 2, 0));
    __m128i longitude_high = _mm_shuffle_epi32(_mm256_extracti128_si256(normalized_longitudes, 1), _MM_SHUFFLE(2, 0, 2, 0));
    __m128i packed_latitudes = _mm_unpacklo_epi64(latitude_low, latitude_high);
    __m128i packed_longitudes = _mm_unpacklo_epi64(longitude_low, longitude_high);
    __m256d latitude_values = _mm256_cvtepi32_pd(packed_latitudes);
    __m256d longitude_values = _mm256_cvtepi32_pd(packed_longitudes);
    __m256d unsigned_correction = _mm256_set1_pd(4294967296.0);
    __m256d zero = _mm256_setzero_pd();

    latitude_values = _mm256_add_pd(latitude_values,
                                    _mm256_and_pd(_mm256_cmp_pd(latitude_values, zero, _CMP_LT_OQ), unsigned_correction));

    longitude_values = _mm256_add_pd(longitude_values,
                                     _mm256_and_pd(_mm256_cmp_pd(longitude_values, zero, _CMP_LT_OQ), unsigned_correction));

    *latitudes = _mm256_fmadd_pd(latitude_values,
                                 _mm256_set1_pd(GEO_INTERNAL_LAT_DENORM_SCALE),
                                 _mm256_set1_pd(-90.0));

    *longitudes = _mm256_fmadd_pd(longitude_values,
                                  _mm256_set1_pd(GEO_INTERNAL_LNG_DENORM_SCALE),
                                  _mm256_set1_pd(-180.0));
}

GEO_TARGET_AVX2 static void decode_avx2(const uint64_t *codes, double *lats, double *lngs, size_t count)
{
    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        __m256i morton_codes = _mm256_loadu_si256((const __m256i *)(codes + i));
        __m256d latitude_values;
        __m256d longitude_values;

        decode4_values(morton_codes, &latitude_values, &longitude_values);

        _mm256_storeu_pd(lats + i, latitude_values);
        _mm256_storeu_pd(lngs + i, longitude_values);
    }

    for (; i < count; ++i) {
        geo_internal_decode(codes[i], lats + i, lngs + i);
    }
}

GEO_TARGET_AVX512 static void decode_avx512(const uint64_t *codes, double *lats, double *lngs, size_t count)
{
    __m512d latitude_scale = _mm512_set1_pd(GEO_INTERNAL_LAT_DENORM_SCALE);
    __m512d longitude_scale = _mm512_set1_pd(GEO_INTERNAL_LNG_DENORM_SCALE);
    __m512d latitude_offset = _mm512_set1_pd(-90.0);
    __m512d longitude_offset = _mm512_set1_pd(-180.0);
    size_t i = 0;

    for (; i + 7 < count; i += 8) {
        __m512i morton_codes = _mm512_loadu_si512((const void *) (codes + i));
        __m256i packed_latitudes = _mm512_cvtepi64_epi32(compact8(morton_codes));
        __m256i packed_longitudes = _mm512_cvtepi64_epi32(compact8(_mm512_srli_epi64(morton_codes, 1)));
        __m512d latitude_values = _mm512_cvtepu32_pd(packed_latitudes);
        __m512d longitude_values = _mm512_cvtepu32_pd(packed_longitudes);

        _mm512_storeu_pd(lats + i, _mm512_fmadd_pd(latitude_values, latitude_scale, latitude_offset));
        _mm512_storeu_pd(lngs + i, _mm512_fmadd_pd(longitude_values, longitude_scale, longitude_offset));
    }

    decode_avx2(codes + i, lats + i, lngs + i, count - i);
}

void geo_simd_decode_batch(const uint64_t *codes, double *lats, double *lngs, size_t count)
{
    x86_backend()->decode_wide(codes, lats, lngs, count);
}

void geo_simd_decode_batch_narrow(const uint64_t *codes, double *lats, double *lngs, size_t count)
{
    x86_backend()->decode_narrow(codes, lats, lngs, count);
}

GEO_TARGET_AVX2 static void extract_interleaved_codes_avx2(const GeoRecord *records, uint64_t *codes, size_t count)
{
    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        __m256i first_pairs = _mm256_loadu_si256((const __m256i *) (records + i));
        __m256i second_pairs = _mm256_loadu_si256((const __m256i *) (records + i + 2));
        __m256i interleaved_codes = _mm256_unpackhi_epi64(first_pairs, second_pairs);
        __m256i ordered_codes = _mm256_permute4x64_epi64(interleaved_codes, _MM_SHUFFLE(3, 1, 2, 0));

        _mm256_storeu_si256((__m256i *) (codes + i), ordered_codes);
    }

    for (; i < count; ++i) {
        codes[i] = records[i].z;
    }
}

GEO_TARGET_AVX512 static void extract_interleaved_codes_avx512(const GeoRecord *records, uint64_t *codes, size_t count)
{
    const __m512i code_indices = _mm512_set_epi64(15, 13, 11, 9, 7, 5, 3, 1);
    size_t i = 0;

    for (; i + 7 < count; i += 8) {
        __m512i first_records = _mm512_loadu_si512((const void *) (records + i));
        __m512i second_records = _mm512_loadu_si512((const void *) (records + i + 4));
        __m512i record_codes = _mm512_permutex2var_epi64(first_records, code_indices, second_records);

        _mm512_storeu_si512((void *) (codes + i), record_codes);
    }

    extract_interleaved_codes_avx2(records + i, codes + i, count - i);
}

void geo_simd_extract_interleaved_codes(const GeoRecord *records, uint64_t *codes, size_t count)
{
    x86_backend()->extract_wide(records, codes, count);
}

void geo_simd_extract_interleaved_codes_narrow(const GeoRecord *records, uint64_t *codes, size_t count)
{
    x86_backend()->extract_narrow(records, codes, count);
}

GEO_TARGET_AVX2 static void decode_interleaved_records_avx2(const GeoRecord *records,
                                                            double *lats,
                                                            double *lngs,
                                                            size_t count)
{
    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        __m256i first_pairs = _mm256_loadu_si256((const __m256i *) (records + i));
        __m256i second_pairs = _mm256_loadu_si256((const __m256i *) (records + i + 2));
        __m256i interleaved_codes = _mm256_unpackhi_epi64(first_pairs, second_pairs);
        __m256i ordered_codes = _mm256_permute4x64_epi64(interleaved_codes, _MM_SHUFFLE(3, 1, 2, 0));
        __m256d latitude_values;
        __m256d longitude_values;

        decode4_values(ordered_codes, &latitude_values, &longitude_values);
        _mm256_storeu_pd(lats + i, latitude_values);
        _mm256_storeu_pd(lngs + i, longitude_values);
    }

    for (; i < count; ++i) {
        geo_internal_decode(records[i].z, lats + i, lngs + i);
    }
}

void geo_simd_decode_interleaved_records_narrow(const GeoRecord *records,
                                                double *lats,
                                                double *lngs,
                                                size_t count)
{
    x86_backend()->decode_records_narrow(records, lats, lngs, count);
}

// =============================================================================
// Distance calculations
// =============================================================================

GEO_TARGET_AVX2 static void haversine_avx2(double lat1, double lng1, const double *lats, const double *lngs, double *out, size_t count)
{
    GeoSimdRadiusQuery query = {
        .center_latitude_radians = lat1 * GEO_INTERNAL_DEG_TO_RAD,
        .center_longitude_radians = lng1 * GEO_INTERNAL_DEG_TO_RAD,
        .cosine_center_latitude = cos(lat1 * GEO_INTERNAL_DEG_TO_RAD),
        .sine_angular_radius = 0.0,
        .haversine_limit = 0.0,
    };

    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        double haversine_lanes[4];

        _mm256_storeu_pd(haversine_lanes,
                         haversine_a4(&query, _mm256_loadu_pd(lats + i), _mm256_loadu_pd(lngs + i)));

        for (unsigned j = 0; j < 4; ++j) {
            double haversine = haversine_lanes[j];

            out[i + j] = 2.0 * GEO_INTERNAL_EARTH_RADIUS_KM * atan2(sqrt(haversine), sqrt(1.0 - haversine));
        }
    }

    for (; i < count; ++i) {
        out[i] = geo_haversine_km(lat1, lng1, lats[i], lngs[i]);
    }
}

void geo_simd_haversine_batch(double lat1, double lng1, const double *lats, const double *lngs, double *out, size_t count)
{
    x86_backend()->haversine(lat1, lng1, lats, lngs, out, count);
}

GEO_TARGET_AVX2 static void fast_distance_avx2(double lat1, double lng1, const double *lats, const double *lngs, double *out, size_t count)
{
    __m256d center_latitude = _mm256_set1_pd(lat1);
    __m256d center_longitude = _mm256_set1_pd(lng1);
    __m256d kilometers_per_degree = _mm256_set1_pd(GEO_INTERNAL_KM_PER_DEG);
    __m256d half = _mm256_set1_pd(0.5);
    __m256d radians_per_degree = _mm256_set1_pd(GEO_INTERNAL_DEG_TO_RAD);
    __m256d positive_180 = _mm256_set1_pd(180.0);
    __m256d negative_180 = _mm256_set1_pd(-180.0);
    __m256d full_circle = _mm256_set1_pd(360.0);
    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        __m256d point_latitudes = _mm256_loadu_pd(lats + i);
        __m256d point_longitudes = _mm256_loadu_pd(lngs + i);
        __m256d longitude_delta = _mm256_sub_pd(point_longitudes, center_longitude);

        longitude_delta = _mm256_blendv_pd(longitude_delta,
                                           _mm256_sub_pd(longitude_delta, full_circle),
                                           _mm256_cmp_pd(longitude_delta, positive_180, _CMP_GT_OQ));

        longitude_delta = _mm256_blendv_pd(longitude_delta,
                                           _mm256_add_pd(longitude_delta, full_circle),
                                           _mm256_cmp_pd(longitude_delta, negative_180, _CMP_LT_OQ));

        __m256d latitude_distance = _mm256_mul_pd(_mm256_sub_pd(point_latitudes, center_latitude), kilometers_per_degree);
        __m256d middle_latitude = _mm256_mul_pd(_mm256_mul_pd(_mm256_add_pd(center_latitude, point_latitudes), half),
                                                radians_per_degree);
        __m256d longitude_distance = _mm256_mul_pd(_mm256_mul_pd(longitude_delta, kilometers_per_degree),
                                                   cos_poly(middle_latitude));

        __m256d squared_distance = _mm256_fmadd_pd(latitude_distance,
                                                   latitude_distance,
                                                   _mm256_mul_pd(longitude_distance, longitude_distance));

        _mm256_storeu_pd(out + i, _mm256_sqrt_pd(squared_distance));
    }

    for (; i < count; ++i) {
        out[i] = geo_fast_distance_km(lat1, lng1, lats[i], lngs[i]);
    }
}

void geo_simd_fast_distance_batch(double lat1, double lng1, const double *lats, const double *lngs, double *out, size_t count)
{
    x86_backend()->fast_distance(lat1, lng1, lats, lngs, out, count);
}

// =============================================================================
// Vector filters
// =============================================================================

GEO_TARGET_AVX2 static size_t filter_range_avx2(const uint64_t *codes,
                                                size_t count,
                                                uint64_t minimum,
                                                uint64_t maximum,
                                                uint8_t *mask)
{
    __m256i sign = _mm256_set1_epi64x((long long) UINT64_C(0x8000000000000000));
    __m256i vector_minimum = _mm256_xor_si256(_mm256_set1_epi64x((long long)minimum), sign);
    __m256i vector_maximum = _mm256_xor_si256(_mm256_set1_epi64x((long long)maximum), sign);
    size_t matches = 0;
    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        __m256i values = _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)(codes + i)), sign);
        __m256i below_minimum = _mm256_cmpgt_epi64(vector_minimum, values);
        __m256i above_maximum = _mm256_cmpgt_epi64(values, vector_maximum);
        int matched_bits = (~_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_or_si256(below_minimum, above_maximum)))) & 15;

        for (unsigned j = 0; j < 4; ++j) {
            mask[i + j] = (uint8_t)((matched_bits >> j) & 1);
        }

        matches += (size_t)__builtin_popcount((unsigned)matched_bits);
    }

    for (; i < count; ++i) {
        mask[i] = codes[i] >= minimum && codes[i] <= maximum;
        matches += mask[i];
    }

    return matches;
}

size_t geo_simd_filter_range(const uint64_t *codes, size_t count, uint64_t minimum, uint64_t maximum, uint8_t *mask)
{
    return x86_backend()->filter_range(codes, count, minimum, maximum, mask);
}

GEO_TARGET_AVX2 static size_t filter_bbox_avx2(const double *lats,
                                               const double *lngs,
                                               size_t count,
                                               double min_lat,
                                               double max_lat,
                                               double min_lng,
                                               double max_lng,
                                               uint8_t *mask)
{
    __m256d minimum_latitude = _mm256_set1_pd(min_lat);
    __m256d maximum_latitude = _mm256_set1_pd(max_lat);
    __m256d minimum_longitude = _mm256_set1_pd(min_lng);
    __m256d maximum_longitude = _mm256_set1_pd(max_lng);
    bool wraps = min_lng > max_lng;
    size_t matches = 0;
    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        __m256d latitudes = _mm256_loadu_pd(lats + i);
        __m256d longitudes = _mm256_loadu_pd(lngs + i);

        __m256d latitude_matches = _mm256_and_pd(_mm256_cmp_pd(latitudes, minimum_latitude, _CMP_GE_OQ),
                                                  _mm256_cmp_pd(latitudes, maximum_latitude, _CMP_LE_OQ));

        __m256d above_minimum_longitude = _mm256_cmp_pd(longitudes, minimum_longitude, _CMP_GE_OQ);
        __m256d below_maximum_longitude = _mm256_cmp_pd(longitudes, maximum_longitude, _CMP_LE_OQ);
        __m256d longitude_matches = wraps ? _mm256_or_pd(above_minimum_longitude, below_maximum_longitude)
                                         : _mm256_and_pd(above_minimum_longitude, below_maximum_longitude);

        int matched_bits = _mm256_movemask_pd(_mm256_and_pd(latitude_matches, longitude_matches));

        for (unsigned j = 0; j < 4; ++j) {
            mask[i + j] = (uint8_t)((matched_bits >> j) & 1);
        }

        matches += (size_t)__builtin_popcount((unsigned)matched_bits);
    }

    for (; i < count; ++i) {
        bool longitude_matches = wraps ? (lngs[i] >= min_lng || lngs[i] <= max_lng)
                                       : (lngs[i] >= min_lng && lngs[i] <= max_lng);

        mask[i] = lats[i] >= min_lat && lats[i] <= max_lat && longitude_matches;
        matches += mask[i];
    }

    return matches;
}

size_t geo_simd_filter_bbox(const double *lats, const double *lngs, size_t count,
                            double min_lat, double max_lat, double min_lng, double max_lng, uint8_t *mask)
{
    return x86_backend()->filter_bbox(lats, lngs, count, min_lat, max_lat, min_lng, max_lng, mask);
}

void geo_simd_prepare_radius_query(double center_latitude,
                                   double center_longitude,
                                   double radius_km,
                                   GeoSimdRadiusQuery *query)
{
    geo_scalar_prepare_radius_query(center_latitude, center_longitude, radius_km, query);
}

GEO_TARGET_AVX2 static size_t filter_radius_avx2(const double *lats,
                                                 const double *lngs,
                                                 size_t count,
                                                 const GeoSimdRadiusQuery *query,
                                                 uint8_t *mask,
                                                 uint64_t *bit_mask)
{
    __m256d limit = _mm256_set1_pd(query->haversine_limit);
    size_t matches = 0;
    size_t i = 0;
    uint64_t packed_bits = 0;

    for (; i + 3 < count; i += 4) {
        __m256d haversine = haversine_a4(query, _mm256_loadu_pd(lats + i), _mm256_loadu_pd(lngs + i));
        int matched_bits = _mm256_movemask_pd(_mm256_cmp_pd(haversine, limit, _CMP_LE_OQ));

        if (mask) {
            for (unsigned j = 0; j < 4; ++j) {
                mask[i + j] = (uint8_t)((matched_bits >> j) & 1);
            }
        }

        if (bit_mask) {
            packed_bits |= (uint64_t) (unsigned) matched_bits << (i & 63U);

            if ((i & 63U) == 60U) {
                bit_mask[i >> 6] = packed_bits;
                packed_bits = 0;
            }
        }

        matches += (size_t)__builtin_popcount((unsigned)matched_bits);
    }

    for (; i < count; ++i) {
        uint8_t matched = geo_scalar_radius_match(query, lats[i], lngs[i]);

        if (mask) {
            mask[i] = matched;
        }

        if (bit_mask) {
            packed_bits |= (uint64_t) matched << (i & 63U);
        }

        matches += matched;
    }

    if (bit_mask && (count & 63U)) {
        bit_mask[count >> 6] = packed_bits;
    }

    return matches;
}

size_t geo_simd_filter_radius_prepared(const double *lats,
                                       const double *lngs,
                                       size_t count,
                                       const GeoSimdRadiusQuery *query,
                                       uint8_t *mask,
                                       uint64_t *bit_mask)
{
    return x86_backend()->filter_radius(lats, lngs, count, query, mask, bit_mask);
}

size_t geo_simd_filter_radius(const double *lats, const double *lngs, size_t count, double lat, double lng, double radius, uint8_t *mask)
{
    GeoSimdRadiusQuery query;

    geo_simd_prepare_radius_query(lat, lng, radius, &query);

    return geo_simd_filter_radius_prepared(lats, lngs, count, &query, mask, NULL);
}

size_t geo_simd_filter_radius_bits(const double *lats,
                                   const double *lngs,
                                   size_t count,
                                   double lat,
                                   double lng,
                                   double radius,
                                   uint64_t *bit_mask)
{
    GeoSimdRadiusQuery query;

    geo_simd_prepare_radius_query(lat, lng, radius, &query);

    return geo_simd_filter_radius_prepared(lats, lngs, count, &query, NULL, bit_mask);
}

GEO_TARGET_AVX2 static size_t filter_bbox_codes_avx2(const uint64_t *codes,
                                                     size_t count,
                                                     double min_lat,
                                                     double max_lat,
                                                     double min_lng,
                                                     double max_lng,
                                                     uint8_t *mask,
                                                     uint64_t *bit_mask)
{
    uint32_t min_latitude = geo_internal_normalized_lat_lower(min_lat);
    uint32_t max_latitude = geo_internal_normalized_lat_upper(max_lat);
    uint32_t min_longitude = geo_internal_normalized_lng_lower(min_lng);
    uint32_t max_longitude = geo_internal_normalized_lng_upper(max_lng);
    __m256i vector_minimum_latitude = _mm256_set1_epi64x(min_latitude);
    __m256i vector_maximum_latitude = _mm256_set1_epi64x(max_latitude);
    __m256i vector_minimum_longitude = _mm256_set1_epi64x(min_longitude);
    __m256i vector_maximum_longitude = _mm256_set1_epi64x(max_longitude);
    __m256i all_bits = _mm256_cmpeq_epi64(_mm256_setzero_si256(), _mm256_setzero_si256());
    bool wraps = min_lng > max_lng;
    size_t matches = 0;
    size_t i = 0;
    uint64_t packed_bits = 0;

    for (; i + 3 < count; i += 4) {
        __m256i morton_codes = _mm256_loadu_si256((const __m256i *) (codes + i));
        __m256i latitudes = compact4(morton_codes);
        __m256i longitudes = compact4(_mm256_srli_epi64(morton_codes, 1));

        __m256i latitude_above_minimum = _mm256_xor_si256(_mm256_cmpgt_epi64(vector_minimum_latitude, latitudes), all_bits);
        __m256i latitude_below_maximum = _mm256_xor_si256(_mm256_cmpgt_epi64(latitudes, vector_maximum_latitude), all_bits);
        __m256i latitude_matches = _mm256_and_si256(latitude_above_minimum, latitude_below_maximum);

        __m256i longitude_above_minimum = _mm256_xor_si256(_mm256_cmpgt_epi64(vector_minimum_longitude, longitudes), all_bits);
        __m256i longitude_below_maximum = _mm256_xor_si256(_mm256_cmpgt_epi64(longitudes, vector_maximum_longitude), all_bits);
        __m256i longitude_matches = wraps ? _mm256_or_si256(longitude_above_minimum, longitude_below_maximum)
                                         : _mm256_and_si256(longitude_above_minimum, longitude_below_maximum);

        int matched_bits = _mm256_movemask_pd(_mm256_castsi256_pd(_mm256_and_si256(latitude_matches, longitude_matches)));

        if (mask) {
            for (unsigned lane = 0; lane < 4; ++lane) {
                mask[i + lane] = (uint8_t) ((matched_bits >> lane) & 1);
            }
        }

        if (bit_mask) {
            packed_bits |= (uint64_t) (unsigned) matched_bits << (i & 63U);

            if ((i & 63U) == 60U) {
                bit_mask[i >> 6] = packed_bits;
                packed_bits = 0;
            }
        }

        matches += (size_t) __builtin_popcount((unsigned) matched_bits);
    }

    for (; i < count; ++i) {
        uint32_t latitude = geo_internal_compact_bits(codes[i]);
        uint32_t longitude = geo_internal_compact_bits(codes[i] >> 1);
        bool longitude_matches = wraps ? longitude >= min_longitude || longitude <= max_longitude
                                       : longitude >= min_longitude && longitude <= max_longitude;
        uint8_t matched = latitude >= min_latitude && latitude <= max_latitude && longitude_matches;

        if (mask) {
            mask[i] = matched;
        }

        if (bit_mask) {
            packed_bits |= (uint64_t) matched << (i & 63U);
        }

        matches += matched;
    }

    if (bit_mask && (count & 63U)) {
        bit_mask[count >> 6] = packed_bits;
    }

    return matches;
}

GEO_TARGET_AVX512 static size_t filter_bbox_codes_avx512(const uint64_t *codes,
                                                         size_t count,
                                                         double min_lat,
                                                         double max_lat,
                                                         double min_lng,
                                                         double max_lng,
                                                         uint8_t *mask,
                                                         uint64_t *bit_mask)
{
    uint32_t min_latitude = geo_internal_normalized_lat_lower(min_lat);
    uint32_t max_latitude = geo_internal_normalized_lat_upper(max_lat);
    uint32_t min_longitude = geo_internal_normalized_lng_lower(min_lng);
    uint32_t max_longitude = geo_internal_normalized_lng_upper(max_lng);
    __m512i vector_minimum_latitude = _mm512_set1_epi64(min_latitude);
    __m512i vector_maximum_latitude = _mm512_set1_epi64(max_latitude);
    __m512i vector_minimum_longitude = _mm512_set1_epi64(min_longitude);
    __m512i vector_maximum_longitude = _mm512_set1_epi64(max_longitude);
    bool wraps = min_lng > max_lng;
    size_t matches = 0;
    size_t i = 0;
    uint64_t packed_bits = 0;

    for (; i + 7 < count; i += 8) {
        __m512i morton_codes = _mm512_loadu_si512((const void *) (codes + i));
        __m512i latitudes = compact8(morton_codes);
        __m512i longitudes = compact8(_mm512_srli_epi64(morton_codes, 1));
        __mmask8 latitude_matches = _mm512_cmp_epu64_mask(latitudes, vector_minimum_latitude, _MM_CMPINT_GE) &
                                     _mm512_cmp_epu64_mask(latitudes, vector_maximum_latitude, _MM_CMPINT_LE);
        __mmask8 above_minimum_longitude = _mm512_cmp_epu64_mask(longitudes,
                                                                 vector_minimum_longitude,
                                                                 _MM_CMPINT_GE);
        __mmask8 below_maximum_longitude = _mm512_cmp_epu64_mask(longitudes,
                                                                 vector_maximum_longitude,
                                                                 _MM_CMPINT_LE);
        __mmask8 longitude_matches = wraps ? above_minimum_longitude | below_maximum_longitude
                                           : above_minimum_longitude & below_maximum_longitude;
        unsigned matched_bits = (unsigned) (latitude_matches & longitude_matches);

        if (mask) {
            for (unsigned lane = 0; lane < 8; ++lane) {
                mask[i + lane] = (uint8_t) ((matched_bits >> lane) & 1U);
            }
        }

        if (bit_mask) {
            packed_bits |= (uint64_t) matched_bits << (i & 63U);

            if ((i & 63U) == 56U) {
                bit_mask[i >> 6] = packed_bits;
                packed_bits = 0;
            }
        }

        matches += (size_t) __builtin_popcount(matched_bits);
    }

    for (; i < count; ++i) {
        uint32_t latitude = geo_internal_compact_bits(codes[i]);
        uint32_t longitude = geo_internal_compact_bits(codes[i] >> 1);
        bool longitude_matches = wraps ? longitude >= min_longitude || longitude <= max_longitude
                                       : longitude >= min_longitude && longitude <= max_longitude;
        uint8_t matched = latitude >= min_latitude && latitude <= max_latitude && longitude_matches;

        if (mask) {
            mask[i] = matched;
        }

        if (bit_mask) {
            packed_bits |= (uint64_t) matched << (i & 63U);
        }

        matches += matched;
    }

    if (bit_mask && (count & 63U)) {
        bit_mask[count >> 6] = packed_bits;
    }

    return matches;
}

size_t geo_simd_filter_bbox_codes(const uint64_t *codes,
                                  size_t count,
                                  double min_lat,
                                  double max_lat,
                                  double min_lng,
                                  double max_lng,
                                  uint8_t *mask)
{
    return x86_backend()->filter_bbox_codes(codes, count, min_lat, max_lat, min_lng, max_lng, mask, NULL);
}

size_t geo_simd_filter_bbox_codes_bits(const uint64_t *codes,
                                       size_t count,
                                       double min_lat,
                                       double max_lat,
                                       double min_lng,
                                       double max_lng,
                                       uint64_t *bit_mask)
{
    return x86_backend()->filter_bbox_codes(codes, count, min_lat, max_lat, min_lng, max_lng, NULL, bit_mask);
}

// =============================================================================
// Vector Morton bit manipulation
// =============================================================================

GEO_TARGET_AVX2 static inline __m256i spread4(__m256i x)
{
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 16)), _mm256_set1_epi64x(0x0000FFFF0000FFFFLL));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 8)), _mm256_set1_epi64x(0x00FF00FF00FF00FFLL));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 4)), _mm256_set1_epi64x(0x0F0F0F0F0F0F0F0FLL));
    x = _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 2)), _mm256_set1_epi64x(0x3333333333333333LL));

    return _mm256_and_si256(_mm256_or_si256(x, _mm256_slli_epi64(x, 1)), _mm256_set1_epi64x(0x5555555555555555LL));
}

GEO_TARGET_AVX2 static void spread_batch_avx2(const uint32_t *values, uint64_t *out, size_t count)
{
    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        __m128i packed_values = _mm_loadu_si128((const __m128i *)(values + i));

        _mm256_storeu_si256((__m256i *)(out + i), spread4(_mm256_cvtepu32_epi64(packed_values)));
    }

    for (; i < count; ++i) {
        out[i] = geo_internal_spread_bits(values[i]);
    }
}

void geo_simd_spread_bits_batch(const uint32_t *values, uint64_t *out, size_t count)
{
    x86_backend()->spread(values, out, count);
}

GEO_TARGET_AVX2 static void compact_batch_avx2(const uint64_t *values, uint32_t *out, size_t count)
{
    size_t i = 0;

    for (; i + 3 < count; i += 4) {
        uint64_t temp[4];

        _mm256_storeu_si256((__m256i *)temp, compact4(_mm256_loadu_si256((const __m256i *)(values + i))));

        for (unsigned j = 0; j < 4; ++j) {
            out[i + j] = (uint32_t) temp[j];
        }
    }

    for (; i < count; ++i) {
        out[i] = geo_internal_compact_bits(values[i]);
    }
}

void geo_simd_compact_bits_batch(const uint64_t *values, uint32_t *out, size_t count)
{
    x86_backend()->compact(values, out, count);
}

GEO_TARGET_AVX2 static size_t exclude_id_bits_avx2(const GeoRecord *records,
                                                   size_t count,
                                                   uint64_t excluded_id,
                                                   uint64_t *candidate_bits)
{
    __m256i excluded = _mm256_set1_epi64x((long long) excluded_id);
    size_t position = 0;

    for (; position + 4U <= count; position += 4U) {
        __m256i first_records = _mm256_loadu_si256((const __m256i *) (records + position));
        __m256i second_records = _mm256_loadu_si256((const __m256i *) (records + position + 2U));
        __m256i first_ids = _mm256_permute4x64_epi64(first_records, 0x88);
        __m256i second_ids = _mm256_permute4x64_epi64(second_records, 0x88);
        __m256i ids = _mm256_permute2x128_si256(first_ids, second_ids, 0x20);
        __m256i equal = _mm256_cmpeq_epi64(ids, excluded);
        uint64_t equal_bits = (uint64_t) _mm256_movemask_pd(_mm256_castsi256_pd(equal));

        candidate_bits[position >> 6U] &= ~(equal_bits << (position & 63U));
    }

    return geo_scalar_exclude_id_bits_from(records, position, count, excluded_id, candidate_bits);
}

GEO_TARGET_AVX512 static size_t exclude_id_bits_avx512(const GeoRecord *records,
                                                       size_t count,
                                                       uint64_t excluded_id,
                                                       uint64_t *candidate_bits)
{
    const __m512i indices = _mm512_set_epi64(14, 12, 10, 8, 6, 4, 2, 0);
    __m512i excluded = _mm512_set1_epi64((long long) excluded_id);
    size_t position = 0;

    for (; position + 8U <= count; position += 8U) {
        __m512i first_records = _mm512_loadu_si512((const void *) (records + position));
        __m512i second_records = _mm512_loadu_si512((const void *) (records + position + 4U));
        __m512i ids = _mm512_permutex2var_epi64(first_records, indices, second_records);
        uint64_t equal_bits = (uint64_t) _mm512_cmpeq_epi64_mask(ids, excluded);

        candidate_bits[position >> 6U] &= ~(equal_bits << (position & 63U));
    }

    return geo_scalar_exclude_id_bits_from(records, position, count, excluded_id, candidate_bits);
}

size_t geo_simd_exclude_id_bits(const GeoRecord *records,
                                size_t count,
                                uint64_t excluded_id,
                                uint64_t *candidate_bits)
{
    return x86_backend()->exclude_id(records, count, excluded_id, candidate_bits);
}

// =============================================================================
// Immutable runtime dispatch
// =============================================================================

static const GeoX86Backend g_scalar_backend = {
    .encode = geo_scalar_encode_batch,
    .encode_records = geo_scalar_encode_records,
    .validate = geo_scalar_validate_points,
    .decode_wide = geo_scalar_decode_batch,
    .decode_narrow = geo_scalar_decode_batch,
    .extract_wide = geo_scalar_extract_interleaved_codes,
    .extract_narrow = geo_scalar_extract_interleaved_codes,
    .decode_records_narrow = geo_scalar_decode_interleaved_records,
    .haversine = geo_scalar_haversine_batch,
    .fast_distance = geo_scalar_fast_distance_batch,
    .filter_range = geo_scalar_filter_range,
    .filter_bbox = geo_scalar_filter_bbox,
    .filter_radius = geo_scalar_filter_radius_prepared,
    .filter_bbox_codes = geo_scalar_filter_bbox_codes,
    .spread = geo_scalar_spread_bits_batch,
    .compact = geo_scalar_compact_bits_batch,
    .exclude_id = geo_scalar_exclude_id_bits,
    .name = "x86-64 scalar fallback",
    .optimal_batch_size = 128,
    .simd_available = false,
};

static const GeoX86Backend g_avx2_backend = {
    .encode = encode_avx2,
    .encode_records = encode_records_avx2,
    .validate = validate_points_avx2,
    .decode_wide = decode_avx2,
    .decode_narrow = decode_avx2,
    .extract_wide = extract_interleaved_codes_avx2,
    .extract_narrow = extract_interleaved_codes_avx2,
    .decode_records_narrow = decode_interleaved_records_avx2,
    .haversine = haversine_avx2,
    .fast_distance = fast_distance_avx2,
    .filter_range = filter_range_avx2,
    .filter_bbox = filter_bbox_avx2,
    .filter_radius = filter_radius_avx2,
    .filter_bbox_codes = filter_bbox_codes_avx2,
    .spread = spread_batch_avx2,
    .compact = compact_batch_avx2,
    .exclude_id = exclude_id_bits_avx2,
    .name = "x86-64 AVX2+FMA (runtime)",
    .optimal_batch_size = 512,
    .simd_available = true,
};

static const GeoX86Backend g_avx512_backend = {
    .encode = encode_avx2,
    .encode_records = encode_records_avx2,
    .validate = validate_points_avx2,
    .decode_wide = decode_avx512,
    .decode_narrow = decode_avx2,
    .extract_wide = extract_interleaved_codes_avx512,
    .extract_narrow = extract_interleaved_codes_avx2,
    .decode_records_narrow = decode_interleaved_records_avx2,
    .haversine = haversine_avx2,
    .fast_distance = fast_distance_avx2,
    .filter_range = filter_range_avx2,
    .filter_bbox = filter_bbox_avx2,
    .filter_radius = filter_radius_avx2,
    .filter_bbox_codes = filter_bbox_codes_avx512,
    .spread = spread_batch_avx2,
    .compact = compact_batch_avx2,
    .exclude_id = exclude_id_bits_avx512,
    .name = "x86-64 AVX-512F+DQ+VL (runtime)",
    .optimal_batch_size = 512,
    .simd_available = true,
};

static const GeoX86Backend *detect_x86_backend(void)
{
#if defined(__GNUC__) || defined(__clang__)
    __builtin_cpu_init();

    bool has_avx2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    bool has_avx512 = has_avx2 &&
                      __builtin_cpu_supports("avx512f") &&
                      __builtin_cpu_supports("avx512dq") &&
                      __builtin_cpu_supports("avx512vl");

    if (has_avx512) {
        return &g_avx512_backend;
    }

    if (has_avx2) {
        return &g_avx2_backend;
    }
#endif

    return &g_scalar_backend;
}

GEO_CONSTRUCTOR static void initialize_x86_backend(void)
{
    const GeoX86Backend *selected = detect_x86_backend();

    atomic_store_explicit(&g_backend, selected, memory_order_relaxed);
}

// =============================================================================
// Public backend information
// =============================================================================

bool geo_simd_available(void)
{
    return x86_backend()->simd_available;
}

const char *geo_simd_get_name(void)
{
    return x86_backend()->name;
}

size_t geo_simd_optimal_batch_size(void)
{
    return x86_backend()->optimal_batch_size;
}

// =============================================================================
// Optimizer-resistant microbenchmarks
// =============================================================================

static void consume_u64(const uint64_t *values, size_t count)
{
    if (count) {
        g_benchmark_sink_u64 ^= values[count / 2] ^ values[count - 1];
    }
}

static void consume_double(const double *values, size_t count)
{
    if (count) {
        g_benchmark_sink_double += values[count / 2] + values[count - 1];
    }
}

GeoSimdBenchmark geo_simd_benchmark_encode(size_t count, int iterations)
{
    GeoSimdBenchmark result = {
        .scalar_time_ms = 0.0,
        .simd_time_ms = 0.0,
        .speedup = 0.0,
        .operations = count * (size_t)(iterations > 0 ? iterations : 0),
        .operation_name = "encode",
    };

    if (!count || iterations <= 0) {
        return result;
    }

    double *latitudes = malloc(count * sizeof(*latitudes));
    double *longitudes = malloc(count * sizeof(*longitudes));
    uint64_t *codes = malloc(count * sizeof(*codes));

    if (!latitudes || !longitudes || !codes) {
        free(latitudes);
        free(longitudes);
        free(codes);

        return result;
    }

    for (size_t i = 0; i < count; ++i) {
        latitudes[i] = (double)(i % 18000) / 100.0 - 90.0;
        longitudes[i] = (double)(i % 36000) / 100.0 - 180.0;
    }

    double start = geo_get_time_ms();

    for (int iteration = 0; iteration < iterations; ++iteration) {
        for (size_t i = 0; i < count; ++i) {
            codes[i] = geo_encode(latitudes[i], longitudes[i]);
        }

        consume_u64(codes, count);
    }

    result.scalar_time_ms = geo_get_time_ms() - start;
    start = geo_get_time_ms();

    for (int iteration = 0; iteration < iterations; ++iteration) {
        geo_simd_encode_batch(latitudes, longitudes, codes, count);
        consume_u64(codes, count);
    }

    result.simd_time_ms = geo_get_time_ms() - start;
    result.speedup = result.simd_time_ms ? result.scalar_time_ms / result.simd_time_ms : 0.0;

    free(latitudes);
    free(longitudes);
    free(codes);

    return result;
}

GeoSimdBenchmark geo_simd_benchmark_decode(size_t count, int iterations)
{
    GeoSimdBenchmark result = {
        .scalar_time_ms = 0.0,
        .simd_time_ms = 0.0,
        .speedup = 0.0,
        .operations = count * (size_t)(iterations > 0 ? iterations : 0),
        .operation_name = "decode",
    };

    if (!count || iterations <= 0) {
        return result;
    }

    uint64_t *codes = malloc(count * sizeof(*codes));
    double *latitudes = malloc(count * sizeof(*latitudes));
    double *longitudes = malloc(count * sizeof(*longitudes));

    if (!codes || !latitudes || !longitudes) {
        free(codes);
        free(latitudes);
        free(longitudes);

        return result;
    }

    for (size_t i = 0; i < count; ++i) {
        codes[i] = i * UINT64_C(12345);
    }

    double start = geo_get_time_ms();

    for (int iteration = 0; iteration < iterations; ++iteration) {
        for (size_t i = 0; i < count; ++i) {
            GeoPoint point = geo_decode(codes[i]);

            latitudes[i] = point.lat;
            longitudes[i] = point.lng;
        }

        consume_double(latitudes, count);
    }

    result.scalar_time_ms = geo_get_time_ms() - start;
    start = geo_get_time_ms();

    for (int iteration = 0; iteration < iterations; ++iteration) {
        geo_simd_decode_batch(codes, latitudes, longitudes, count);
        consume_double(latitudes, count);
    }

    result.simd_time_ms = geo_get_time_ms() - start;
    result.speedup = result.simd_time_ms ? result.scalar_time_ms / result.simd_time_ms : 0.0;

    free(codes);
    free(latitudes);
    free(longitudes);

    return result;
}

GeoSimdBenchmark geo_simd_benchmark_haversine(size_t count, int iterations)
{
    GeoSimdBenchmark result = {
        .scalar_time_ms = 0.0,
        .simd_time_ms = 0.0,
        .speedup = 0.0,
        .operations = count * (size_t)(iterations > 0 ? iterations : 0),
        .operation_name = "haversine",
    };

    if (!count || iterations <= 0) {
        return result;
    }

    double *latitudes = malloc(count * sizeof(*latitudes));
    double *longitudes = malloc(count * sizeof(*longitudes));
    double *distances = malloc(count * sizeof(*distances));

    if (!latitudes || !longitudes || !distances) {
        free(latitudes);
        free(longitudes);
        free(distances);

        return result;
    }

    for (size_t i = 0; i < count; ++i) {
        latitudes[i] = (double)(i % 18000) / 100.0 - 90.0;
        longitudes[i] = (double)(i % 36000) / 100.0 - 180.0;
    }

    double start = geo_get_time_ms();

    for (int iteration = 0; iteration < iterations; ++iteration) {
        for (size_t i = 0; i < count; ++i) {
            distances[i] = geo_haversine_km(-23.55, -46.63, latitudes[i], longitudes[i]);
        }

        consume_double(distances, count);
    }

    result.scalar_time_ms = geo_get_time_ms() - start;
    start = geo_get_time_ms();

    for (int iteration = 0; iteration < iterations; ++iteration) {
        geo_simd_haversine_batch(-23.55, -46.63, latitudes, longitudes, distances, count);
        consume_double(distances, count);
    }

    result.simd_time_ms = geo_get_time_ms() - start;
    result.speedup = result.simd_time_ms ? result.scalar_time_ms / result.simd_time_ms : 0.0;

    free(latitudes);
    free(longitudes);
    free(distances);

    return result;
}

GeoSimdBenchmark geo_simd_benchmark_filter_radius(size_t count, int iterations)
{
    GeoSimdBenchmark result = {
        .scalar_time_ms = 0.0,
        .simd_time_ms = 0.0,
        .speedup = 0.0,
        .operations = count * (size_t)(iterations > 0 ? iterations : 0),
        .operation_name = "filter_radius",
    };

    if (!count || iterations <= 0) {
        return result;
    }

    double *latitudes = malloc(count * sizeof(*latitudes));
    double *longitudes = malloc(count * sizeof(*longitudes));
    uint8_t *mask = malloc(count);

    if (!latitudes || !longitudes || !mask) {
        free(latitudes);
        free(longitudes);
        free(mask);

        return result;
    }

    for (size_t i = 0; i < count; ++i) {
        latitudes[i] = -23.55 + ((double)(i % 1001) / 1000.0 - 0.5) * 2.0;
        longitudes[i] = -46.63 + ((double)((i * 37) % 1001) / 1000.0 - 0.5) * 2.0;
    }

    double start = geo_get_time_ms();

    for (int iteration = 0; iteration < iterations; ++iteration) {
        size_t matched = 0;

        for (size_t i = 0; i < count; ++i) {
            mask[i] = geo_haversine_km(-23.55, -46.63, latitudes[i], longitudes[i]) <= 50.0;
            matched += mask[i];
        }

        g_benchmark_sink_u64 ^= matched;
    }

    result.scalar_time_ms = geo_get_time_ms() - start;
    start = geo_get_time_ms();

    for (int iteration = 0; iteration < iterations; ++iteration) {
        g_benchmark_sink_u64 ^= geo_simd_filter_radius(latitudes, longitudes, count, -23.55, -46.63, 50.0, mask);
    }

    result.simd_time_ms = geo_get_time_ms() - start;
    result.speedup = result.simd_time_ms ? result.scalar_time_ms / result.simd_time_ms : 0.0;

    free(latitudes);
    free(longitudes);
    free(mask);

    return result;
}

#endif
