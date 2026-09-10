#include "geobolt/geo_index.h"
#include "geo_index_internal.h"
#include "geo_index_io.h"
#include "geo_index_persistence.h"
#include "geo_index_private.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#define GEO_HAS_MMAP 1
#else
#define GEO_HAS_MMAP 0
#endif

#define GEO_QUERY_BATCH 512
#define GEO_RADIX_THRESHOLD 2048
#define GEO_RADIX_BITS 11
#define GEO_RADIX_BUCKETS (1U << GEO_RADIX_BITS)
#define GEO_KNN_LEAF_RECORDS 64
#define GEO_COVER_LOCAL_CELLS 128
#define GEO_COVER_TARGET_RECORDS_PER_CELL UINT64_C(64)
#ifndef GEO_COVER_MAX_DENSITY_REFINEMENT
#define GEO_COVER_MAX_DENSITY_REFINEMENT 6
#endif
#define GEO_COVER_CONTAINED UINT8_C(1)
#define GEO_COVER_EXHAUSTED UINT8_C(2)

typedef struct {
    uint32_t min_lat;
    uint32_t max_lat;
    uint32_t min_lng;
    uint32_t max_lng;
} GeoNormalizedBox;

typedef struct {
    ZRange *out;
    uint8_t *needs_filter;
    int count;
    int capacity;
    bool overflow;
} GeoRangeBuilder;

typedef struct {
    uint64_t morton_prefix;
    uint32_t min_latitude;
    uint32_t min_longitude;
    uint8_t depth;
    uint8_t flags;
} GeoCoverCell;

typedef struct {
    GeoRecord record;
    double distance;
} GeoRankedRecord;

typedef struct {
    uint64_t morton_prefix;
    size_t begin;
    size_t end;
    size_t source_index;
    double lower_bound;
    uint64_t size;
    uint32_t min_latitude;
    uint32_t min_longitude;
} GeoKnnCell;

typedef struct {
    GeoKnnCell *cells;
    size_t count;
    size_t capacity;
} GeoKnnQueue;

struct GeoKnnWorkspace {
    GeoRankedRecord *heap;
    size_t heap_capacity;
    GeoKnnQueue queue;
};

// =============================================================================
// Morton prefix directory
// =============================================================================

static void prefix_directory_destroy(GeoIndex *index)
{
    geo_density_index_destroy(index);

    free(index->prefix_offsets);

    index->prefix_offsets = NULL;
    index->prefix_bits = 0;
    index->maximum_density_refinement = 0;
}

static bool prefix_directory_build(GeoIndex *index)
{
    prefix_directory_destroy(index);

    uint8_t bits = geo_persisted_prefix_bits(index->count);

    if (!bits) {
        return geo_density_index_build(index);
    }

    size_t bucket_count = (size_t) 1 << bits;
    size_t *offsets = malloc((bucket_count + 1) * sizeof(*offsets));

    if (!offsets) {
        return false;
    }

    size_t next_bucket = 0;
    unsigned shift = 64U - bits;

    for (size_t position = 0; position < index->count; ++position) {
        size_t record_bucket = (size_t) (index->records[position].z >> shift);

        while (next_bucket <= record_bucket) {
            offsets[next_bucket++] = position;
        }
    }

    while (next_bucket <= bucket_count) {
        offsets[next_bucket++] = index->count;
    }

    index->prefix_offsets = offsets;
    index->prefix_bits = bits;

    return geo_density_index_build(index);
}

// =============================================================================
// Checked size arithmetic
// =============================================================================

static bool size_mul(size_t a, size_t b, size_t *out)
{
    if (a && b > SIZE_MAX / a) {
        return false;
    }

    *out = a * b;

    return true;
}

static bool size_add(size_t a, size_t b, size_t *out)
{
    if (b > SIZE_MAX - a) {
        return false;
    }

    *out = a + b;

    return true;
}

// =============================================================================
// Morton encoding and coordinate conversion
// =============================================================================

uint64_t geo_spread_bits(uint32_t v)
{
    return geo_internal_spread_bits(v);
}

uint32_t geo_compact_bits(uint64_t v)
{
    return geo_internal_compact_bits(v);
}

uint32_t geo_normalize_lat(double lat)
{
    return geo_internal_normalize_lat(lat);
}

uint32_t geo_normalize_lng(double lng)
{
    return geo_internal_normalize_lng(lng);
}

double geo_denormalize_lat(uint32_t v)
{
    return geo_internal_denormalize_lat(v);
}

double geo_denormalize_lng(uint32_t v)
{
    return geo_internal_denormalize_lng(v);
}

uint64_t geo_encode(double lat, double lng)
{
    return geo_internal_encode(lat, lng);
}

GeoPoint geo_decode(uint64_t z)
{
    GeoPoint point;

    geo_internal_decode(z, &point.lat, &point.lng);

    return point;
}

double geo_to_radians(double degrees)
{
    return degrees * GEO_INTERNAL_DEG_TO_RAD;
}

double geo_to_degrees(double radians)
{
    return radians * GEO_INTERNAL_RAD_TO_DEG;
}

// =============================================================================
// Distance calculations and spherical bounding boxes
// =============================================================================

static double haversine_a(double lat1_rad, double cos_lat1, double lng1_rad, double lat2, double lng2)
{
    double lat2_rad = geo_to_radians(lat2);
    double sin_delta_lat = sin((lat2_rad - lat1_rad) * 0.5);
    double sin_delta_lng = sin((geo_to_radians(lng2) - lng1_rad) * 0.5);
    double haversine = sin_delta_lat * sin_delta_lat + cos_lat1 * cos(lat2_rad) * sin_delta_lng * sin_delta_lng;

    return fmax(0.0, fmin(1.0, haversine));
}

double geo_haversine_km(double lat1, double lng1, double lat2, double lng2)
{
    double lat1_rad = geo_to_radians(lat1);
    double haversine = haversine_a(lat1_rad, cos(lat1_rad), geo_to_radians(lng1), lat2, lng2);

    return 2.0 * GEO_EARTH_RADIUS_KM * atan2(sqrt(haversine), sqrt(1.0 - haversine));
}

double geo_haversine_m(double lat1, double lng1, double lat2, double lng2)
{
    return geo_haversine_km(lat1, lng1, lat2, lng2) * 1000.0;
}

double geo_fast_distance_km(double lat1, double lng1, double lat2, double lng2)
{
    double dlat = (lat2 - lat1) * GEO_INTERNAL_KM_PER_DEG;
    double dlng = geo_wrap_lng(lng2 - lng1) * GEO_INTERNAL_KM_PER_DEG * cos(geo_to_radians((lat1 + lat2) * 0.5));

    return hypot(dlat, dlng);
}

static void bounding_box_prepared(double longitude,
                                  double angular_radius,
                                  double latitude_radians,
                                  double cosine_latitude,
                                  double sine_angular_radius,
                                  double *minimum_latitude,
                                  double *maximum_latitude,
                                  double *minimum_longitude,
                                  double *maximum_longitude)
{
    if (angular_radius >= M_PI) {
        *minimum_latitude = GEO_MIN_LAT;
        *maximum_latitude = GEO_MAX_LAT;
        *minimum_longitude = GEO_MIN_LNG;
        *maximum_longitude = GEO_MAX_LNG;

        return;
    }

    double low = latitude_radians - angular_radius;
    double high = latitude_radians + angular_radius;
    const double half_pi = M_PI * 0.5;

    *minimum_latitude = geo_to_degrees(fmax(low, -half_pi));
    *maximum_latitude = geo_to_degrees(fmin(high, half_pi));

    if (low <= -half_pi || high >= half_pi) {
        *minimum_longitude = GEO_MIN_LNG;
        *maximum_longitude = GEO_MAX_LNG;

        return;
    }

    double ratio = fmax(-1.0, fmin(1.0, sine_angular_radius / cosine_latitude));
    double delta = geo_to_degrees(asin(ratio));

    *minimum_longitude = geo_wrap_lng(longitude - delta);
    *maximum_longitude = geo_wrap_lng(longitude + delta);
}

void geo_bounding_box(double lat, double lng, double radius_km, double *min_lat, double *max_lat, double *min_lng, double *max_lng)
{
    if (!min_lat || !max_lat || !min_lng || !max_lng) {
        return;
    }

    if (!geo_is_valid_point(lat, lng) || !isfinite(radius_km) || radius_km < 0.0) {
        *min_lat = *max_lat = *min_lng = *max_lng = NAN;

        return;
    }

    double angular_radius = radius_km / GEO_EARTH_RADIUS_KM;
    double latitude_radians = geo_to_radians(lat);

    bounding_box_prepared(lng,
                          angular_radius,
                          latitude_radians,
                          cos(latitude_radians),
                          sin(angular_radius),
                          min_lat,
                          max_lat,
                          min_lng,
                          max_lng);
}

// =============================================================================
// Mutable index lifecycle and ingestion
// =============================================================================

GeoIndex *geo_index_create(size_t capacity)
{
    if (!capacity) {
        capacity = 1024;
    }

    size_t bytes;

    if (!size_mul(capacity, sizeof(GeoRecord), &bytes)) {
        return NULL;
    }

    GeoIndex *index = calloc(1, sizeof(*index));

    if (!index) {
        return NULL;
    }

    index->records = malloc(bytes);

    if (!index->records) {
        free(index);

        return NULL;
    }

    index->capacity = capacity;

    return index;
}

void geo_index_destroy(GeoIndex *index)
{
    if (!index) {
        return;
    }

    prefix_directory_destroy(index);

#if GEO_HAS_MMAP
    if (index->mapping) {
        munmap(index->mapping, index->mapping_size);
    } else
#endif
    {
        free(index->records);
    }

    free(index);
}

bool geo_index_reserve(GeoIndex *index, size_t capacity)
{
    if (!index || index->read_only) {
        return false;
    }

    if (capacity <= index->capacity) {
        return true;
    }

    size_t bytes;

    if (!size_mul(capacity, sizeof(GeoRecord), &bytes)) {
        return false;
    }

    GeoRecord *records = realloc(index->records, bytes);

    if (!records) {
        return false;
    }

    index->records = records;
    index->capacity = capacity;

    return true;
}

static bool index_grow_for(GeoIndex *index, size_t extra)
{
    size_t required;

    if (!size_add(index->count, extra, &required)) {
        return false;
    }

    if (required <= index->capacity) {
        return true;
    }

    size_t capacity = index->capacity ? index->capacity : 1024;

    while (capacity < required) {
        if (capacity > SIZE_MAX / 2) {
            capacity = required;

            break;
        }

        capacity *= 2;
    }

    return geo_index_reserve(index, capacity);
}

bool geo_index_add(GeoIndex *index, uint64_t id, double lat, double lng)
{
    if (!index || index->read_only || !geo_is_valid_point(lat, lng) || !index_grow_for(index, 1)) {
        return false;
    }

    index->records[index->count++] = (GeoRecord) {
        .id = id,
        .z = geo_encode(lat, lng),
    };

    prefix_directory_destroy(index);
    index->sorted = false;

    return true;
}

bool geo_index_add_batch(GeoIndex *index, const uint64_t *ids, const double *lats, const double *lngs, size_t count)
{
    if (!index || index->read_only || (count && (!ids || !lats || !lngs))) {
        return false;
    }

    if (!count) {
        return true;
    }

    if (!geo_simd_validate_points(lats, lngs, count)) {
        return false;
    }

    return geo_index_add_batch_unchecked(index, ids, lats, lngs, count);
}

bool geo_index_add_batch_unchecked(GeoIndex *index,
                                   const uint64_t *ids,
                                   const double *latitudes,
                                   const double *longitudes,
                                   size_t count)
{
    if (!index || index->read_only || (count && (!ids || !latitudes || !longitudes))) {
        return false;
    }

    if (!count) {
        return true;
    }

    if (!index_grow_for(index, count)) {
        return false;
    }

    geo_simd_encode_records(ids,
                            latitudes,
                            longitudes,
                            index->records + index->count,
                            count);

    index->count += count;
    prefix_directory_destroy(index);
    index->sorted = false;

    return true;
}

bool geo_index_add_records(GeoIndex *index, const GeoRecord *records, size_t count)
{
    if (!index || index->read_only || (count && !records)) {
        return false;
    }

    if (!count) {
        return true;
    }

    if (!index_grow_for(index, count)) {
        return false;
    }

    memcpy(index->records + index->count, records, count * sizeof(*records));
    index->count += count;
    prefix_directory_destroy(index);
    index->sorted = false;

    return true;
}

// =============================================================================
// Index construction: radix sort by Morton key
// =============================================================================

int geo_compare_records_by_z(const void *a, const void *b)
{
    uint64_t first_z = ((const GeoRecord *)a)->z;
    uint64_t second_z = ((const GeoRecord *)b)->z;

    return (first_z > second_z) - (first_z < second_z);
}

int geo_compare_records_by_id(const void *a, const void *b)
{
    uint64_t first_id = ((const GeoRecord *)a)->id;
    uint64_t second_id = ((const GeoRecord *)b)->id;

    return (first_id > second_id) - (first_id < second_id);
}

static bool radix_sort_records(GeoRecord *records, size_t count)
{
    size_t bytes;

    if (!size_mul(count, sizeof(*records), &bytes)) {
        return false;
    }

    GeoRecord *scratch = malloc(bytes);

    if (!scratch) {
        return false;
    }

    uint64_t varying_bits = 0;
    uint64_t first_z = records[0].z;

    for (size_t i = 1; i < count; ++i) {
        varying_bits |= first_z ^ records[i].z;
    }

    GeoRecord *source = records;
    GeoRecord *destination = scratch;
    unsigned passes = 0;

    for (unsigned shift = 0; shift < 64; shift += GEO_RADIX_BITS) {
        unsigned pass_bits = 64U - shift < GEO_RADIX_BITS ? 64U - shift : GEO_RADIX_BITS;
        size_t bucket_count = (size_t) 1 << pass_bits;
        uint64_t bucket_mask = bucket_count - 1;

        if (!((varying_bits >> shift) & bucket_mask)) {
            continue;
        }

        size_t offsets[GEO_RADIX_BUCKETS] = { 0 };

        for (size_t i = 0; i < count; ++i) {
            offsets[(source[i].z >> shift) & bucket_mask]++;
        }

        size_t sum = 0;

        for (size_t i = 0; i < bucket_count; ++i) {
            size_t records_in_bucket = offsets[i];

            offsets[i] = sum;
            sum += records_in_bucket;
        }

        for (size_t i = 0; i < count; ++i) {
            size_t bucket = (size_t) ((source[i].z >> shift) & bucket_mask);

            destination[offsets[bucket]++] = source[i];
        }

        GeoRecord *temporary = source;

        source = destination;
        destination = temporary;
        passes++;
    }

    if (passes & 1U) {
        memcpy(records, source, bytes);
    }

    free(scratch);

    return true;
}

void geo_index_sort_transient(GeoIndex *index)
{
    if (!index || index->sorted || index->read_only) {
        return;
    }

    if (index->count >= GEO_RADIX_THRESHOLD) {
        if (!radix_sort_records(index->records, index->count)) {
            qsort(index->records, index->count, sizeof(GeoRecord), geo_compare_records_by_z);
        }
    } else if (index->count > 1) {
        qsort(index->records, index->count, sizeof(GeoRecord), geo_compare_records_by_z);
    }

    index->sorted = true;
}

bool geo_index_build(GeoIndex *index)
{
    if (!index || index->read_only) {
        return false;
    }

    if (!index->sorted) {
        geo_index_sort_transient(index);
    }

    return geo_index_finalize_sorted(index);
}

bool geo_index_finalize_sorted(GeoIndex *index)
{
    if (!index || !index->sorted || index->read_only) {
        return false;
    }

    return prefix_directory_build(index);
}

void geo_index_clear(GeoIndex *index)
{
    if (!index || index->read_only) {
        return;
    }

    prefix_directory_destroy(index);
    index->count = 0;
    index->sorted = false;
}

bool geo_index_is_read_only(const GeoIndex *index)
{
    return index && index->read_only;
}

// =============================================================================
// Persistent read-only indexes
// =============================================================================

bool geo_index_save(const GeoIndex *index, const char *path)
{
    size_t density_bytes = geo_density_index_serialized_size(index);
    size_t record_bytes;

    if (!index || !path || !index->sorted || !density_bytes ||
        !size_mul(index->count, sizeof(*index->records), &record_bytes)) {
        return false;
    }

    char *temporary_path = NULL;
    FILE *file = geo_io_create_atomic_file(path, &temporary_path);

    if (!file) {
        return false;
    }

    GeoFileHeader header = {
        .magic = { 0 },
        .version = GEO_FILE_VERSION,
        .record_size = sizeof(GeoRecord),
        .endian_marker = GEO_FILE_ENDIAN_MARKER,
        .prefix_bits = index->prefix_bits,
        .count = index->count,
        .records_checksum = geo_persisted_checksum_update(geo_persisted_checksum_initial(),
                                                          index->records,
                                                          record_bytes),
        .prefix_checksum = geo_persisted_checksum_initial(),
        .density_checksum = 0,
        .density_bytes = density_bytes,
    };

    memcpy(header.magic, GEO_FILE_MAGIC, sizeof(header.magic));

    bool succeeded = fwrite(&header, sizeof(header), 1, file) == 1;

    if (succeeded && index->count) {
        succeeded = fwrite(index->records, sizeof(GeoRecord), index->count, file) == index->count;
    }

    if (succeeded && index->prefix_bits) {
        size_t prefix_count = ((size_t) 1 << index->prefix_bits) + 1;

        succeeded = geo_io_write_u64_offsets(file,
                                             index->prefix_offsets,
                                             prefix_count,
                                             &header.prefix_checksum);
    }

    void *density_buffer = succeeded ? malloc(density_bytes) : NULL;

    if (succeeded) {
        succeeded = density_buffer &&
                    geo_density_index_serialize(index, density_buffer, density_bytes) &&
                    fwrite(density_buffer, 1, density_bytes, file) == density_bytes;
    }

    if (succeeded) {
        header.density_checksum = geo_persisted_checksum_update(geo_persisted_checksum_initial(),
                                                                density_buffer,
                                                                density_bytes);
        succeeded = fseeko(file, 0, SEEK_SET) == 0 && fwrite(&header, sizeof(header), 1, file) == 1;
    }

    free(density_buffer);

    if (succeeded) {
        succeeded = geo_io_publish_atomic_file(file, temporary_path, path);
    } else {
        geo_io_discard_atomic_file(file, temporary_path);
    }

    free(temporary_path);

    return succeeded;
}

GeoIndex *geo_index_open_mmap(const char *path)
{
    const GeoIndexOpenOptions options = {
        .verify_checksums = true,
        .validate_sorted_order = true,
        .memory_advice = GEO_INDEX_ADVICE_NORMAL,
    };

    return geo_index_open_mmap_with_options(path, &options);
}

GeoIndex *geo_index_open_mmap_with_options(const char *path, const GeoIndexOpenOptions *options)
{
#if GEO_HAS_MMAP
    if (!path || (options && options->memory_advice > GEO_INDEX_ADVICE_HUGE_PAGE)) {
        return NULL;
    }

    int file_descriptor = open(path, O_RDONLY);

    if (file_descriptor < 0) {
        return NULL;
    }

    struct stat st;

    if (fstat(file_descriptor, &st) || st.st_size < (off_t)sizeof(GeoFileHeader)) {
        close(file_descriptor);

        return NULL;
    }

    size_t map_size = (size_t) st.st_size;
    void *mapping = mmap(NULL, map_size, PROT_READ, MAP_PRIVATE, file_descriptor, 0);

    close(file_descriptor);

    if (mapping == MAP_FAILED) {
        return NULL;
    }

    const GeoFileHeader *header = mapping;
    size_t record_bytes = 0;
    size_t expected_size = 0;
    size_t prefix_bytes = 0;
    uint8_t prefix_bits = (uint8_t) header->prefix_bits;
    bool prefix_layout_valid = header->prefix_bits <= GEO_PERSISTED_PREFIX_MAX_BITS;

    if (prefix_layout_valid && prefix_bits) {
        size_t prefix_count = ((size_t) 1 << prefix_bits) + 1;

        if (!size_mul(prefix_count, sizeof(uint64_t), &prefix_bytes)) {
            prefix_layout_valid = false;
        }
    }

    bool valid = !memcmp(header->magic, GEO_FILE_MAGIC, sizeof(header->magic)) &&
                 header->version == GEO_FILE_VERSION &&
                 header->record_size == sizeof(GeoRecord) &&
                 header->endian_marker == GEO_FILE_ENDIAN_MARKER &&
                 prefix_layout_valid &&
                 header->count <= SIZE_MAX &&
                 header->density_bytes <= SIZE_MAX &&
                 size_mul((size_t)header->count, sizeof(GeoRecord), &record_bytes) &&
                 size_add(sizeof(*header), record_bytes, &expected_size) &&
                 size_add(expected_size, prefix_bytes, &expected_size) &&
                 size_add(expected_size, (size_t) header->density_bytes, &expected_size) &&
                 expected_size == map_size;

    const unsigned char *mapped_bytes = mapping;
    const GeoRecord *mapped_records = (const GeoRecord *) (mapped_bytes + sizeof(*header));
    const uint64_t *mapped_prefix_offsets = (const uint64_t *) (mapped_bytes + sizeof(*header) + record_bytes);
    const void *density_data = mapped_bytes + sizeof(*header) + record_bytes + prefix_bytes;

    if (valid && (!options || options->verify_checksums)) {
        uint64_t records_checksum = geo_persisted_checksum_update(geo_persisted_checksum_initial(),
                                                                  mapped_records,
                                                                  record_bytes);
        uint64_t prefix_checksum = geo_persisted_checksum_update(geo_persisted_checksum_initial(),
                                                                 mapped_prefix_offsets,
                                                                 prefix_bytes);
        uint64_t density_checksum = geo_persisted_checksum_update(geo_persisted_checksum_initial(),
                                                                  density_data,
                                                                  (size_t) header->density_bytes);

        valid = records_checksum == header->records_checksum &&
                prefix_checksum == header->prefix_checksum &&
                density_checksum == header->density_checksum;
    }

    for (size_t i = 1; valid && (!options || options->validate_sorted_order) && i < (size_t) header->count; ++i) {
        const GeoRecord *previous = mapped_records + i - 1U;
        const GeoRecord *current = mapped_records + i;

        valid = previous->z <= current->z;
    }

    if (!valid) {
        munmap(mapping, map_size);

        return NULL;
    }

    GeoIndex *index = calloc(1, sizeof(*index));

    if (!index) {
        munmap(mapping, map_size);

        return NULL;
    }

    index->records = (GeoRecord *) mapped_records;
    index->count = (size_t)header->count;
    index->sorted = true;
    index->read_only = true;
    index->mapping = mapping;
    index->mapping_size = map_size;
    index->content_checksum = geo_persisted_checksum_initial();
    index->content_checksum = geo_persisted_checksum_update(index->content_checksum,
                                                            &header->records_checksum,
                                                            sizeof(header->records_checksum));
    index->content_checksum = geo_persisted_checksum_update(index->content_checksum,
                                                            &header->prefix_checksum,
                                                            sizeof(header->prefix_checksum));
    index->content_checksum = geo_persisted_checksum_update(index->content_checksum,
                                                            &header->density_checksum,
                                                            sizeof(header->density_checksum));

    if (prefix_bits) {
        size_t prefix_count = ((size_t) 1 << prefix_bits) + 1;

        if (prefix_count < 2U || prefix_count > SIZE_MAX / sizeof(*index->prefix_offsets)) {
            geo_index_destroy(index);

            return NULL;
        }

        index->prefix_offsets = malloc(prefix_count * sizeof(*index->prefix_offsets));

        if (!index->prefix_offsets) {
            geo_index_destroy(index);

            return NULL;
        }

        for (size_t i = 0; i < prefix_count; ++i) {
            if (mapped_prefix_offsets[i] > SIZE_MAX) {
                geo_index_destroy(index);

                return NULL;
            }

            index->prefix_offsets[i] = (size_t) mapped_prefix_offsets[i];
        }

        index->prefix_bits = prefix_bits;

        bool offsets_valid = index->prefix_offsets[0] == 0 &&
                             index->prefix_offsets[prefix_count - 1] == index->count;

        for (size_t i = 1; offsets_valid && i < prefix_count; ++i) {
            offsets_valid = index->prefix_offsets[i - 1] <= index->prefix_offsets[i] &&
                            index->prefix_offsets[i] <= index->count;
        }

        if (!offsets_valid) {
            geo_index_destroy(index);

            return NULL;
        }
    }

    size_t consumed = 0;

    if (!geo_density_index_attach(index,
                                  density_data,
                                  (size_t) header->density_bytes,
                                  &consumed) ||
        consumed != (size_t) header->density_bytes) {
        geo_index_destroy(index);

        return NULL;
    }

    if (options && options->memory_advice != GEO_INDEX_ADVICE_NORMAL &&
        !geo_index_advise(index, options->memory_advice)) {
        geo_index_destroy(index);

        return NULL;
    }

    return index;
#else
    (void) path;

    return NULL;
#endif
}

bool geo_index_advise(GeoIndex *index, GeoIndexMemoryAdvice advice)
{
#if GEO_HAS_MMAP
    if (!index || !index->mapping || advice > GEO_INDEX_ADVICE_HUGE_PAGE) {
        return false;
    }

    int native_advice;

    switch (advice) {
        case GEO_INDEX_ADVICE_NORMAL:
            native_advice = POSIX_MADV_NORMAL;
            break;
        case GEO_INDEX_ADVICE_RANDOM:
            native_advice = POSIX_MADV_RANDOM;
            break;
        case GEO_INDEX_ADVICE_SEQUENTIAL:
            native_advice = POSIX_MADV_SEQUENTIAL;
            break;
        case GEO_INDEX_ADVICE_WILL_NEED:
            native_advice = POSIX_MADV_WILLNEED;
            break;
        case GEO_INDEX_ADVICE_HUGE_PAGE:
#if defined(MADV_HUGEPAGE)
            native_advice = MADV_HUGEPAGE;
            break;
#else
            return false;
#endif
        default:
            return false;
    }

    return posix_madvise(index->mapping, index->mapping_size, native_advice) == 0;
#else
    (void) index;
    (void) advice;

    return false;
#endif
}

// =============================================================================
// Binary search over the sorted Morton key array
// =============================================================================

size_t geo_lower_bound(const GeoRecord *records, size_t count, uint64_t key)
{
    size_t first = 0;

    while (count) {
        size_t step = count >> 1;
        size_t middle = first + step;

        if (records[middle].z < key) {
            first = middle + 1;
            count -= step + 1;
        } else {
            count = step;
        }
    }

    return first;
}

size_t geo_upper_bound(const GeoRecord *records, size_t count, uint64_t key)
{
    size_t first = 0;

    while (count) {
        size_t step = count >> 1;
        size_t middle = first + step;

        if (records[middle].z <= key) {
            first = middle + 1;
            count -= step + 1;
        } else {
            count = step;
        }
    }

    return first;
}

static size_t index_lower_bound(const GeoIndex *index, uint64_t key)
{
    if (!index->prefix_bits) {
        return geo_lower_bound(index->records, index->count, key);
    }

    size_t bucket = (size_t) (key >> (64U - index->prefix_bits));
    size_t begin = index->prefix_offsets[bucket];
    size_t end = index->prefix_offsets[bucket + 1];

    return begin + geo_lower_bound(index->records + begin, end - begin, key);
}

static size_t index_upper_bound(const GeoIndex *index, uint64_t key)
{
    if (!index->prefix_bits) {
        return geo_upper_bound(index->records, index->count, key);
    }

    size_t bucket = (size_t) (key >> (64U - index->prefix_bits));
    size_t begin = index->prefix_offsets[bucket];
    size_t end = index->prefix_offsets[bucket + 1];

    return begin + geo_upper_bound(index->records + begin, end - begin, key);
}

// =============================================================================
// Adaptive hierarchical Morton cover
// =============================================================================

static void range_append(GeoRangeBuilder *builder, uint64_t range_min, uint64_t range_max, bool needs_filter)
{
    if (builder->overflow) {
        return;
    }

    if (builder->count > 0) {
        ZRange *previous = &builder->out[builder->count - 1];
        bool previous_needs_filter = !builder->needs_filter || builder->needs_filter[builder->count - 1];
        bool overlaps = range_min <= previous->max;
        bool adjacent = previous->max != UINT64_MAX && range_min == previous->max + 1;
        bool compatible_filter = !builder->needs_filter || previous_needs_filter == needs_filter;

        if (overlaps || (adjacent && compatible_filter)) {
            if (range_max > previous->max) {
                previous->max = range_max;
            }

            if (overlaps && builder->needs_filter) {
                builder->needs_filter[builder->count - 1] = (uint8_t) (previous_needs_filter || needs_filter);
            }

            return;
        }
    }

    if (builder->count == builder->capacity) {
        builder->overflow = true;
        return;
    }

    builder->out[builder->count] = (ZRange) {
        .min = range_min,
        .max = range_max,
    };

    if (builder->needs_filter) {
        builder->needs_filter[builder->count] = (uint8_t) needs_filter;
    }

    builder->count++;
}

static bool cover_classify_cell(GeoCoverCell *cell, const GeoNormalizedBox *boxes, int box_count)
{
    uint64_t cell_size = UINT64_C(1) << (32U - cell->depth);
    uint64_t cell_max_latitude = (uint64_t) cell->min_latitude + cell_size - 1;
    uint64_t cell_max_longitude = (uint64_t) cell->min_longitude + cell_size - 1;
    bool intersects = false;
    bool fully_contained = false;

    for (int box_index = 0; box_index < box_count; ++box_index) {
        const GeoNormalizedBox *box = &boxes[box_index];

        bool cell_intersects_box = cell_max_latitude >= box->min_lat &&
                                   cell->min_latitude <= box->max_lat &&
                                   cell_max_longitude >= box->min_lng &&
                                   cell->min_longitude <= box->max_lng;

        bool box_contains_cell = cell->min_latitude >= box->min_lat &&
                                 cell_max_latitude <= box->max_lat &&
                                 cell->min_longitude >= box->min_lng &&
                                 cell_max_longitude <= box->max_lng;

        intersects = intersects || cell_intersects_box;
        fully_contained = fully_contained || box_contains_cell;
    }

    if (fully_contained) {
        cell->flags |= GEO_COVER_CONTAINED;
    }

    return intersects;
}

static size_t cover_make_children(const GeoCoverCell *parent,
                                  const GeoNormalizedBox *boxes,
                                  int box_count,
                                  GeoCoverCell children[4])
{
    uint64_t parent_size = UINT64_C(1) << (32U - parent->depth);
    uint64_t child_size = parent_size >> 1;
    size_t child_count = 0;

    for (unsigned quadrant = 0; quadrant < 4; ++quadrant) {
        GeoCoverCell child = {
            .morton_prefix = (parent->morton_prefix << 2) | quadrant,
            .min_latitude = (uint32_t) ((uint64_t) parent->min_latitude + ((quadrant & 1U) ? child_size : 0)),
            .min_longitude = (uint32_t) ((uint64_t) parent->min_longitude + ((quadrant & 2U) ? child_size : 0)),
            .depth = (uint8_t) (parent->depth + 1),
            .flags = 0,
        };

        if (cover_classify_cell(&child, boxes, box_count)) {
            children[child_count++] = child;
        }
    }

    return child_count;
}

static void cover_cell_range(const GeoCoverCell *cell, uint64_t *range_min, uint64_t *range_max)
{
    unsigned remaining = 64U - (unsigned) cell->depth * 2U;

    if (remaining == 64U) {
        *range_min = 0;
        *range_max = UINT64_MAX;
    } else if (!remaining) {
        *range_min = cell->morton_prefix;
        *range_max = cell->morton_prefix;
    } else {
        *range_min = cell->morton_prefix << remaining;
        *range_max = *range_min | (UINT64_MAX >> (64U - remaining));
    }
}

static int cover_boxes_internal(const GeoNormalizedBox *boxes,
                                int box_count,
                                ZRange *out,
                                uint8_t *needs_filter,
                                int capacity,
                                int requested_depth)
{
    if (!boxes || box_count < 1 || !out || capacity < 1) {
        return 0;
    }

    if (requested_depth < 0) {
        requested_depth = 0;
    } else if (requested_depth > 32) {
        requested_depth = 32;
    }

    GeoCoverCell local_cells[GEO_COVER_LOCAL_CELLS];
    GeoCoverCell *cells = local_cells;

    if (capacity > GEO_COVER_LOCAL_CELLS) {
        if ((size_t) capacity > SIZE_MAX / sizeof(*cells)) {
            return 0;
        }

        cells = malloc((size_t) capacity * sizeof(*cells));

        if (!cells) {
            return 0;
        }
    }

    cells[0] = (GeoCoverCell) {
        .morton_prefix = 0,
        .min_latitude = 0,
        .min_longitude = 0,
        .depth = 0,
        .flags = 0,
    };

    if (!cover_classify_cell(cells, boxes, box_count)) {
        if (cells != local_cells) {
            free(cells);
        }

        return 0;
    }

    size_t cell_count = 1;

    for (;;) {
        size_t selected = SIZE_MAX;
        uint64_t selected_size = 0;

        for (size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
            const GeoCoverCell *cell = cells + cell_index;

            if ((cell->flags & (GEO_COVER_CONTAINED | GEO_COVER_EXHAUSTED)) ||
                cell->depth >= requested_depth ||
                cell->depth == 32) {
                continue;
            }

            uint64_t cell_size = UINT64_C(1) << (32U - cell->depth);

            if (cell_size <= selected_size) {
                continue;
            }

            selected = cell_index;
            selected_size = cell_size;
        }

        if (selected == SIZE_MAX) {
            break;
        }

        GeoCoverCell selected_children[4];
        size_t selected_child_count = cover_make_children(cells + selected,
                                                          boxes,
                                                          box_count,
                                                          selected_children);

        if (!selected_child_count || cell_count + selected_child_count - 1 > (size_t) capacity) {
            cells[selected].flags |= GEO_COVER_EXHAUSTED;

            continue;
        }

        size_t tail_count = cell_count - selected - 1;

        if (selected_child_count > 1 && tail_count) {
            memmove(cells + selected + selected_child_count,
                    cells + selected + 1,
                    tail_count * sizeof(*cells));
        }

        memcpy(cells + selected, selected_children, selected_child_count * sizeof(*cells));
        cell_count += selected_child_count - 1;
    }

    GeoRangeBuilder builder = {
        .out = out,
        .needs_filter = needs_filter,
        .count = 0,
        .capacity = capacity,
        .overflow = false,
    };

    for (size_t cell_index = 0; cell_index < cell_count; ++cell_index) {
        uint64_t range_min;
        uint64_t range_max;

        cover_cell_range(cells + cell_index, &range_min, &range_max);
        range_append(&builder,
                     range_min,
                     range_max,
                     !(cells[cell_index].flags & GEO_COVER_CONTAINED));
    }

    if (cells != local_cells) {
        free(cells);
    }

    return builder.overflow ? 0 : builder.count;
}

static int cover_boxes(const GeoNormalizedBox *boxes, int box_count, ZRange *out, int capacity, int requested_depth)
{
    return cover_boxes_internal(boxes, box_count, out, NULL, capacity, requested_depth);
}

static int make_boxes(double min_lat, double max_lat, double min_lng, double max_lng, GeoNormalizedBox boxes[2])
{
    uint32_t normalized_min_lat = geo_normalize_lat(min_lat);
    uint32_t normalized_max_lat = geo_normalize_lat(max_lat);

    if (min_lng <= max_lng) {
        boxes[0] = (GeoNormalizedBox) {
            .min_lat = normalized_min_lat,
            .max_lat = normalized_max_lat,
            .min_lng = geo_normalize_lng(min_lng),
            .max_lng = geo_normalize_lng(max_lng),
        };

        return 1;
    }

    boxes[0] = (GeoNormalizedBox) {
        .min_lat = normalized_min_lat,
        .max_lat = normalized_max_lat,
        .min_lng = geo_normalize_lng(min_lng),
        .max_lng = UINT32_MAX,
    };

    boxes[1] = (GeoNormalizedBox) {
        .min_lat = normalized_min_lat,
        .max_lat = normalized_max_lat,
        .min_lng = 0,
        .max_lng = geo_normalize_lng(max_lng),
    };

    return 2;
}

static int radius_ranges(double lat, double lng, double radius, ZRange *out, int capacity, int depth)
{
    double min_lat;
    double max_lat;
    double min_lng;
    double max_lng;

    geo_bounding_box(lat, lng, radius, &min_lat, &max_lat, &min_lng, &max_lng);

    GeoNormalizedBox boxes[2];
    int count = make_boxes(min_lat, max_lat, min_lng, max_lng, boxes);

    return cover_boxes(boxes, count, out, capacity, depth);
}

static int radius_boxes_prepared(const GeoSimdRadiusQuery *query,
                                 double longitude,
                                 double radius,
                                 GeoNormalizedBox boxes[2])
{
    double minimum_latitude;
    double maximum_latitude;
    double minimum_longitude;
    double maximum_longitude;
    double angular_radius = fmin(M_PI, radius / GEO_EARTH_RADIUS_KM);
    bounding_box_prepared(longitude,
                          angular_radius,
                          query->center_latitude_radians,
                          query->cosine_center_latitude,
                          query->sine_angular_radius,
                          &minimum_latitude,
                          &maximum_latitude,
                          &minimum_longitude,
                          &maximum_longitude);

    return make_boxes(minimum_latitude,
                      maximum_latitude,
                      minimum_longitude,
                      maximum_longitude,
                      boxes);
}

static int radius_ranges_prepared(const GeoSimdRadiusQuery *query,
                                  double longitude,
                                  double radius,
                                  ZRange *out,
                                  int capacity,
                                  int depth)
{
    GeoNormalizedBox boxes[2];
    int box_count = radius_boxes_prepared(query, longitude, radius, boxes);

    return cover_boxes(boxes, box_count, out, capacity, depth);
}

static int cover_depth_for_record_count(uint64_t record_count)
{
    uint64_t covered_records = GEO_COVER_TARGET_RECORDS_PER_CELL;
    int depth = 0;

    while (record_count > covered_records && depth < 32) {
        if (covered_records > UINT64_MAX / 4U) {
            return 32;
        }

        covered_records *= 4U;
        depth++;
    }

    return depth;
}

unsigned geo_index_density_refinement(const GeoIndex *index, double latitude, double longitude)
{
#if defined(GEO_DISABLE_DENSITY_COVER)
    (void) index;
    (void) latitude;
    (void) longitude;

    return 0;
#else
    if (!index ||
        !index->maximum_density_refinement ||
        !index->prefix_bits ||
        !index->prefix_offsets ||
        !index->count ||
        !geo_is_valid_point(latitude, longitude)) {
        return 0;
    }

    return geo_density_index_query(index, geo_encode(latitude, longitude), NULL);
#endif
}

static int cover_depth_for_query(uint64_t record_count, unsigned density_refinement)
{
    int depth = cover_depth_for_record_count(record_count);

    for (unsigned refinement = 0; refinement < density_refinement && depth < 32; ++refinement) {
        depth++;
    }

    return depth;
}

int geo_build_ranges(double lat, double lng, double radius, ZRange *out, int capacity)
{
    if (!geo_is_valid_point(lat, lng) || !isfinite(radius) || radius < 0.0) {
        return 0;
    }

    return radius_ranges(lat, lng, radius, out, capacity, 24);
}

int geo_build_ranges_adaptive(double lat, double lng, double radius, ZRange *out, int capacity, int precision)
{
    if (!geo_is_valid_point(lat, lng) || !isfinite(radius) || radius < 0.0) {
        return 0;
    }

    return radius_ranges(lat, lng, radius, out, capacity, precision);
}

// =============================================================================
// Search-result storage and distance ordering
// =============================================================================

GeoSearchResult *geo_result_create(size_t capacity)
{
    if (!capacity) {
        capacity = 64;
    }

    GeoSearchResult *result = calloc(1, sizeof(*result));

    if (!result) {
        return NULL;
    }

    if (!geo_result_reserve(result, capacity)) {
        free(result);

        return NULL;
    }

    return result;
}

void geo_result_destroy(GeoSearchResult *result)
{
    if (result) {
        free(result->results);
        free(result);
    }
}

void geo_result_clear(GeoSearchResult *result)
{
    if (result) {
        result->count = 0;
    }
}

bool geo_result_reserve(GeoSearchResult *result, size_t capacity)
{
    if (!result) {
        return false;
    }

    if (capacity <= result->capacity) {
        return true;
    }

    size_t bytes;

    if (!size_mul(capacity, sizeof(GeoRecord), &bytes)) {
        return false;
    }

    GeoRecord *records = realloc(result->results, bytes);

    if (!records) {
        return false;
    }

    result->results = records;
    result->capacity = capacity;

    return true;
}

bool geo_result_add(GeoSearchResult *result, const GeoRecord *record)
{
    if (!result || !record) {
        return false;
    }

    if (result->count == result->capacity) {
        size_t capacity = result->capacity ? result->capacity : 64;

        if (capacity > SIZE_MAX / 2 || !geo_result_reserve(result, capacity * 2)) {
            return false;
        }
    }

    result->results[result->count++] = *record;

    return true;
}

GeoIdResult *geo_id_result_create(size_t capacity)
{
    if (!capacity) {
        capacity = 64;
    }

    GeoIdResult *result = calloc(1, sizeof(*result));

    if (!result) {
        return NULL;
    }

    if (!geo_id_result_reserve(result, capacity)) {
        free(result);

        return NULL;
    }

    return result;
}

void geo_id_result_destroy(GeoIdResult *result)
{
    if (result) {
        free(result->ids);
        free(result);
    }
}

void geo_id_result_clear(GeoIdResult *result)
{
    if (result) {
        result->count = 0;
    }
}

bool geo_id_result_reserve(GeoIdResult *result, size_t capacity)
{
    if (!result) {
        return false;
    }

    if (capacity <= result->capacity) {
        return true;
    }

    size_t bytes;

    if (!size_mul(capacity, sizeof(*result->ids), &bytes)) {
        return false;
    }

    uint64_t *ids = realloc(result->ids, bytes);

    if (!ids) {
        return false;
    }

    result->ids = ids;
    result->capacity = capacity;

    return true;
}

static int compare_ranked(const void *a, const void *b)
{
    const GeoRankedRecord *first = a;
    const GeoRankedRecord *second = b;

    if (first->distance < second->distance) {
        return -1;
    }

    if (first->distance > second->distance) {
        return 1;
    }

    return (first->record.id > second->record.id) - (first->record.id < second->record.id);
}

void geo_result_sort_by_distance(GeoSearchResult *result, double lat, double lng)
{
    if (!result || result->count < 2 || !geo_is_valid_point(lat, lng)) {
        return;
    }

    GeoRankedRecord *ranked = malloc(result->count * sizeof(*ranked));

    if (!ranked) {
        return;
    }

    for (size_t i = 0; i < result->count; ++i) {
        GeoPoint point = geo_decode(result->results[i].z);

        ranked[i] = (GeoRankedRecord) {
            .record = result->results[i],
            .distance = geo_haversine_km(lat, lng, point.lat, point.lng),
        };
    }

    qsort(ranked, result->count, sizeof(*ranked), compare_ranked);

    for (size_t i = 0; i < result->count; ++i) {
        result->results[i] = ranked[i].record;
    }

    free(ranked);
}

// =============================================================================
// Radius and bounding-box search engines
// =============================================================================

static bool valid_index(const GeoIndex *index)
{
    return index && index->sorted && (index->records || !index->count);
}

static void reset_stats(GeoSearchStats *stats)
{
    if (stats) {
        memset(stats, 0, sizeof(*stats));
    }
}

static void append_records_from_bits(const GeoRecord *records,
                                     size_t record_count,
                                     const uint64_t *match_bits,
                                     GeoSearchResult *result)
{
    size_t word_count = geo_internal_bit_word_count(record_count);

    for (size_t word_index = 0; word_index < word_count; ++word_index) {
        size_t word_offset = word_index * 64U;
        size_t records_in_word = record_count - word_offset;

        if (records_in_word > 64U) {
            records_in_word = 64U;
        }

        uint64_t valid_bits = records_in_word == 64U
                                  ? UINT64_MAX
                                  : (UINT64_C(1) << records_in_word) - 1;
        uint64_t bits = match_bits[word_index] & valid_bits;

        if (bits == valid_bits) {
            memcpy(result->results + result->count,
                   records + word_offset,
                   records_in_word * sizeof(*records));
            result->count += records_in_word;

            continue;
        }

        while (bits) {
            unsigned matched_lane = (unsigned) __builtin_ctzll(bits);

            result->results[result->count++] = records[word_offset + matched_lane];
            bits &= bits - 1;
        }
    }
}

static bool append_ids_from_bits(const GeoRecord *records,
                                 size_t record_count,
                                 const uint64_t *match_bits,
                                 GeoIdResult *result,
                                 bool allow_growth)
{
    size_t word_count = geo_internal_bit_word_count(record_count);

    for (size_t word_index = 0; word_index < word_count; ++word_index) {
        size_t word_offset = word_index * 64U;
        size_t records_in_word = record_count - word_offset;

        if (records_in_word > 64U) {
            records_in_word = 64U;
        }

        uint64_t valid_bits = records_in_word == 64U
                                  ? UINT64_MAX
                                  : (UINT64_C(1) << records_in_word) - 1;
        uint64_t bits = match_bits[word_index] & valid_bits;
        size_t matches = (size_t) __builtin_popcountll(bits);
        size_t required;

        if (!size_add(result->count, matches, &required)) {
            return false;
        }

        if (required > result->capacity) {
            if (!allow_growth) {
                return false;
            }

            size_t capacity = result->capacity ? result->capacity : 64;

            while (capacity < required) {
                if (capacity > SIZE_MAX / 2U) {
                    capacity = required;
                    break;
                }

                capacity *= 2U;
            }

            if (!geo_id_result_reserve(result, capacity)) {
                return false;
            }
        }

        while (bits) {
            unsigned matched_lane = (unsigned) __builtin_ctzll(bits);

            result->ids[result->count++] = records[word_offset + matched_lane].id;
            bits &= bits - 1;
        }
    }

    return true;
}

static size_t filter_record_match_bits(const GeoRecord *records,
                                       size_t count,
                                       uint64_t *match_bits,
                                       GeoRecordBitFilter bit_filter,
                                       GeoRecordFilter filter,
                                       void *filter_context)
{
    if (bit_filter) {
        return bit_filter(records, count, match_bits, filter_context);
    }

    size_t matched = 0;
    size_t word_count = geo_internal_bit_word_count(count);

    for (size_t word_index = 0; word_index < word_count; ++word_index) {
        uint64_t candidates = match_bits[word_index];

        while (candidates) {
            unsigned lane = (unsigned) __builtin_ctzll(candidates);
            size_t record_index = word_index * 64U + lane;

            if (!filter(records + record_index, filter_context)) {
                match_bits[word_index] &= ~(UINT64_C(1) << lane);
            }

            candidates &= candidates - 1U;
        }

        matched += (size_t) __builtin_popcountll(match_bits[word_index]);
    }

    return matched;
}

static void set_all_record_match_bits(uint64_t *match_bits, size_t count)
{
    size_t full_words = count / 64U;
    size_t remaining_bits = count % 64U;

    for (size_t word = 0; word < full_words; ++word) {
        match_bits[word] = UINT64_MAX;
    }

    if (remaining_bits) {
        match_bits[full_words] = (UINT64_C(1) << remaining_bits) - 1U;
    }
}

static bool scan_radius(const GeoIndex *index,
                        const ZRange *ranges,
                        const size_t *cached_begins,
                        const size_t *cached_ends,
                        int range_count,
                        const GeoSimdRadiusQuery *radius_query,
                        GeoQueryStrategy strategy,
                        GeoSearchResult *result,
                        GeoIdResult *id_result,
                        bool allow_id_growth,
                        size_t *matched_count,
                        GeoSearchStats *stats,
                        GeoRecordBitFilter bit_filter,
                        GeoRecordFilter filter,
                        void *filter_context)
{
    size_t begins[GEO_QUERY_MAX_RANGES];
    size_t ends[GEO_QUERY_MAX_RANGES];
    size_t total = 0;

    for (int i = 0; i < range_count; ++i) {
        begins[i] = cached_begins ? cached_begins[i] : index_lower_bound(index, ranges[i].min);
        ends[i] = cached_ends ? cached_ends[i] : index_upper_bound(index, ranges[i].max);

        if (!size_add(total, ends[i] - begins[i], &total)) {
            return false;
        }
    }

    if (result) {
        size_t required;

        if (!size_add(result->count, total, &required) ||
            (required > result->capacity && !geo_result_reserve(result, required))) {
            return false;
        }
    }

    if (strategy == GEO_QUERY_STRATEGY_COPY && !filter) {
        if (id_result) {
            size_t required;

            if (!size_add(id_result->count, total, &required) ||
                (required > id_result->capacity &&
                 (!allow_id_growth || !geo_id_result_reserve(id_result, required)))) {
                return false;
            }
        }

        for (int range = 0; range < range_count; ++range) {
            size_t range_records = ends[range] - begins[range];

            if (result && range_records) {
                memcpy(result->results + result->count,
                       index->records + begins[range],
                       range_records * sizeof(*index->records));
                result->count += range_records;
            }

            if (id_result) {
                for (size_t position = begins[range]; position < ends[range]; ++position) {
                    id_result->ids[id_result->count++] = index->records[position].id;
                }
            }
        }

        if (matched_count) {
            *matched_count = total;
        }

        if (stats) {
            stats->records_scanned = 0;
            stats->records_matched = total;
            stats->ranges_checked = (uint64_t) range_count;
        }

        return true;
    }

    if (strategy == GEO_QUERY_STRATEGY_COPY) {
        uint64_t match_bits[(GEO_QUERY_BATCH + 63U) / 64U];
        size_t matched = 0;

        for (int range = 0; range < range_count; ++range) {
            for (size_t position = begins[range]; position < ends[range];) {
                size_t batch_count = ends[range] - position;

                if (batch_count > GEO_QUERY_BATCH) {
                    batch_count = GEO_QUERY_BATCH;
                }

                set_all_record_match_bits(match_bits, batch_count);
                size_t batch_matches = filter_record_match_bits(index->records + position,
                                                                 batch_count,
                                                                 match_bits,
                                                                 bit_filter,
                                                                 filter,
                                                                 filter_context);

                matched += batch_matches;

                if (result) {
                    append_records_from_bits(index->records + position, batch_count, match_bits, result);
                }

                if (id_result && !append_ids_from_bits(index->records + position,
                                                       batch_count,
                                                       match_bits,
                                                       id_result,
                                                       allow_id_growth)) {
                    return false;
                }

                position += batch_count;
            }
        }

        if (matched_count) {
            *matched_count = matched;
        }

        if (stats) {
            stats->records_scanned = total;
            stats->records_matched = matched;
            stats->ranges_checked = (uint64_t) range_count;
        }

        return true;
    }

    double latitudes[GEO_QUERY_BATCH];
    double longitudes[GEO_QUERY_BATCH];
    uint64_t match_bits[(GEO_QUERY_BATCH + 63U) / 64U];
    size_t matched = 0;

    for (int range = 0; range < range_count; ++range) {
        for (size_t position = begins[range]; position < ends[range];) {
            size_t batch_count = ends[range] - position;

            if (batch_count > GEO_QUERY_BATCH) {
                batch_count = GEO_QUERY_BATCH;
            }

            geo_simd_decode_interleaved_records_narrow(index->records + position,
                                                        latitudes,
                                                        longitudes,
                                                        batch_count);
            size_t batch_matches = geo_simd_filter_radius_prepared(latitudes,
                                                                   longitudes,
                                                                   batch_count,
                                                                   radius_query,
                                                                   NULL,
                                                                   (result || id_result || filter) ? match_bits : NULL);

            if (filter) {
                batch_matches = filter_record_match_bits(index->records + position,
                                                         batch_count,
                                                         match_bits,
                                                         bit_filter,
                                                         filter,
                                                         filter_context);
            }

            matched += batch_matches;

            if (result) {
                append_records_from_bits(index->records + position, batch_count, match_bits, result);
            }

            if (id_result && !append_ids_from_bits(index->records + position,
                                                   batch_count,
                                                   match_bits,
                                                   id_result,
                                                   allow_id_growth)) {
                return false;
            }

            position += batch_count;
        }
    }

    if (matched_count) {
        *matched_count = matched;
    }

    if (stats) {
        stats->records_scanned = total;
        stats->records_matched = matched;
        stats->ranges_checked = (uint64_t) range_count;
    }

    return true;
}

static bool scan_bbox(const GeoIndex *index,
                      const ZRange *ranges,
                      const uint8_t *needs_filter,
                      int range_count,
                      double min_lat,
                      double max_lat,
                      double min_lng,
                      double max_lng,
                      GeoSearchResult *result,
                      size_t *matched_count,
                      GeoSearchStats *stats,
                      GeoRecordBitFilter bit_filter,
                      GeoRecordFilter filter,
                      void *filter_context)
{
    size_t begins[GEO_QUERY_MAX_RANGES];
    size_t ends[GEO_QUERY_MAX_RANGES];
    size_t total = 0;

    for (int range = 0; range < range_count; ++range) {
        begins[range] = index_lower_bound(index, ranges[range].min);
        ends[range] = index_upper_bound(index, ranges[range].max);

        if (!size_add(total, ends[range] - begins[range], &total)) {
            return false;
        }
    }

    if (result) {
        size_t required;

        if (!size_add(result->count, total, &required) ||
            (required > result->capacity && !geo_result_reserve(result, required))) {
            return false;
        }
    }

    uint64_t codes[GEO_QUERY_BATCH];
    uint64_t match_bits[(GEO_QUERY_BATCH + 63U) / 64U];
    size_t matched = 0;

    for (int range = 0; range < range_count; ++range) {
        size_t range_count_records = ends[range] - begins[range];

        if (!needs_filter[range]) {
            if (!filter) {
                if (result && range_count_records) {
                    memcpy(result->results + result->count,
                           index->records + begins[range],
                           range_count_records * sizeof(*index->records));
                    result->count += range_count_records;
                }

                matched += range_count_records;
                continue;
            }

            for (size_t position = begins[range]; position < ends[range];) {
                size_t batch_count = ends[range] - position;

                if (batch_count > GEO_QUERY_BATCH) {
                    batch_count = GEO_QUERY_BATCH;
                }

                set_all_record_match_bits(match_bits, batch_count);
                size_t batch_matches = filter_record_match_bits(index->records + position,
                                                                 batch_count,
                                                                 match_bits,
                                                                 bit_filter,
                                                                 filter,
                                                                 filter_context);

                matched += batch_matches;

                if (result) {
                    append_records_from_bits(index->records + position, batch_count, match_bits, result);
                }

                position += batch_count;
            }

            continue;
        }

        for (size_t position = begins[range]; position < ends[range];) {
            size_t batch_count = ends[range] - position;

            if (batch_count > GEO_QUERY_BATCH) {
                batch_count = GEO_QUERY_BATCH;
            }

            geo_simd_extract_interleaved_codes(index->records + position, codes, batch_count);

            size_t batch_matches = (result || filter)
                                       ? geo_simd_filter_bbox_codes_bits(codes,
                                                                         batch_count,
                                                                         min_lat,
                                                                         max_lat,
                                                                         min_lng,
                                                                         max_lng,
                                                                         match_bits)
                                       : geo_simd_filter_bbox_codes(codes,
                                                                    batch_count,
                                                                    min_lat,
                                                                    max_lat,
                                                                    min_lng,
                                                                    max_lng,
                                                                    NULL);

            if (filter) {
                batch_matches = filter_record_match_bits(index->records + position,
                                                         batch_count,
                                                         match_bits,
                                                         bit_filter,
                                                         filter,
                                                         filter_context);
            }

            matched += batch_matches;

            if (result) {
                append_records_from_bits(index->records + position, batch_count, match_bits, result);
            }

            position += batch_count;
        }
    }

    if (matched_count) {
        *matched_count = matched;
    }

    if (stats) {
        stats->records_scanned = total;
        stats->records_matched = matched;
        stats->ranges_checked = (uint64_t) range_count;
    }

    return true;
}

bool geo_index_prepare_radius_query(double latitude,
                                    double longitude,
                                    double radius_km,
                                    uint64_t record_count,
                                    unsigned density_refinement,
                                    GeoRadiusQueryPlan *plan)
{
    if (!plan || !geo_is_valid_point(latitude, longitude) || !isfinite(radius_km) || radius_km < 0.0) {
        return false;
    }

    plan->latitude = latitude;
    plan->longitude = longitude;
    plan->radius_km = radius_km;
    plan->bounds_index = NULL;
    plan->bounds_cached = false;
    geo_simd_prepare_radius_query(latitude, longitude, radius_km, &plan->radius_query);

    if (radius_km / GEO_EARTH_RADIUS_KM >= M_PI) {
        plan->ranges[0] = (ZRange) {
            .min = 0,
            .max = UINT64_MAX,
        };
        plan->candidate_records = record_count;
        plan->estimated_output_bytes = record_count <= UINT64_MAX / sizeof(GeoRecord)
                                           ? record_count * sizeof(GeoRecord)
                                           : UINT64_MAX;
        plan->estimated_cost = plan->estimated_output_bytes;
        plan->range_count = 1;
        plan->cover_depth = 0;
        plan->strategy = GEO_QUERY_STRATEGY_COPY;

        return true;
    }

    int cover_depth = cover_depth_for_query(record_count, density_refinement);

    plan->range_count = radius_ranges_prepared(&plan->radius_query,
                                               longitude,
                                               radius_km,
                                               plan->ranges,
                                               GEO_QUERY_MAX_RANGES,
                                               cover_depth);
    plan->candidate_records = 0;
    plan->estimated_output_bytes = 0;
    plan->estimated_cost = 0;
    plan->cover_depth = (uint8_t) cover_depth;
    plan->strategy = GEO_QUERY_STRATEGY_FILTER;

    if (plan->range_count <= 0) {
        return false;
    }

    return true;
}

uint64_t geo_index_radius_plan_candidate_count(const GeoIndex *index, const GeoRadiusQueryPlan *plan)
{
    if (!valid_index(index) || !plan || plan->range_count <= 0) {
        return UINT64_MAX;
    }

    if (plan->bounds_cached && plan->bounds_index == index) {
        return plan->candidate_records;
    }

    uint64_t candidate_count = 0;

    for (int range = 0; range < plan->range_count; ++range) {
        size_t begin = index_lower_bound(index, plan->ranges[range].min);
        size_t end = index_upper_bound(index, plan->ranges[range].max);
        uint64_t range_records = end - begin;

        if (range_records > UINT64_MAX - candidate_count) {
            return UINT64_MAX;
        }

        candidate_count += range_records;
    }

    return candidate_count;
}

static uint64_t radius_plan_saturating_multiply(uint64_t value, uint64_t multiplier)
{
    return value && multiplier > UINT64_MAX / value ? UINT64_MAX : value * multiplier;
}

static uint64_t radius_plan_saturating_add(uint64_t first, uint64_t second)
{
    return second > UINT64_MAX - first ? UINT64_MAX : first + second;
}

static void radius_plan_set_cost(GeoRadiusQueryPlan *plan, size_t output_record_size)
{
    uint64_t output_bytes = radius_plan_saturating_multiply(plan->candidate_records, output_record_size);

    // A circle occupies pi/4 of its bounding box. Using 3/4 is deliberately conservative
    // and keeps the model integer-only in the query hot path.
    plan->estimated_output_bytes = output_bytes - output_bytes / 4U;

    uint64_t filter_cost = radius_plan_saturating_multiply(plan->candidate_records, 10U);
    uint64_t range_cost = radius_plan_saturating_multiply((uint64_t) plan->range_count, 96U);
    uint64_t output_cost = plan->estimated_output_bytes / 16U;

    plan->estimated_cost = radius_plan_saturating_add(radius_plan_saturating_add(filter_cost, range_cost), output_cost);
}

bool geo_index_prepare_radius_query_for_index(const GeoIndex *index,
                                              double latitude,
                                              double longitude,
                                              double radius_km,
                                              size_t output_record_size,
                                              GeoRadiusQueryPlan *plan)
{
    if (!valid_index(index) || !plan ||
        !geo_is_valid_point(latitude, longitude) ||
        !isfinite(radius_km) ||
        radius_km < 0.0) {
        return false;
    }

    if (radius_km / GEO_EARTH_RADIUS_KM >= M_PI) {
        return geo_index_prepare_radius_query(latitude,
                                              longitude,
                                              radius_km,
                                              index->count,
                                              0,
                                              plan);
    }

    GeoRadiusQueryPlan candidate = {
        .latitude = latitude,
        .longitude = longitude,
        .radius_km = radius_km,
        .strategy = GEO_QUERY_STRATEGY_FILTER,
    };

    geo_simd_prepare_radius_query(latitude, longitude, radius_km, &candidate.radius_query);

    int minimum_depth = cover_depth_for_record_count(index->count);
    unsigned refinement = geo_density_index_query(index, geo_encode(latitude, longitude), NULL);

    if (!refinement) {
        return geo_index_prepare_radius_query(latitude,
                                              longitude,
                                              radius_km,
                                              index->count,
                                              0,
                                              plan);
    }

    int maximum_depth = minimum_depth + (int) refinement;

    if (maximum_depth > 32) {
        maximum_depth = 32;
    }

    bool selected = false;
    GeoNormalizedBox boxes[2];
    int box_count = radius_boxes_prepared(&candidate.radius_query, longitude, radius_km, boxes);

    for (int depth = maximum_depth; depth >= minimum_depth; --depth) {
        candidate.range_count = cover_boxes(boxes,
                                            box_count,
                                            candidate.ranges,
                                            GEO_QUERY_MAX_RANGES,
                                            depth);

        if (candidate.range_count <= 0) {
            continue;
        }

        candidate.cover_depth = (uint8_t) depth;
        candidate.bounds_index = index;
        candidate.bounds_cached = true;
        candidate.candidate_records = 0;

        for (int range = 0; range < candidate.range_count; ++range) {
            candidate.range_begins[range] = index_lower_bound(index, candidate.ranges[range].min);
            candidate.range_ends[range] = index_upper_bound(index, candidate.ranges[range].max);
            candidate.candidate_records = radius_plan_saturating_add(candidate.candidate_records,
                                                                     candidate.range_ends[range] -
                                                                         candidate.range_begins[range]);
        }

        radius_plan_set_cost(&candidate, output_record_size);

        if (!selected || candidate.estimated_cost < plan->estimated_cost) {
            *plan = candidate;
            selected = true;
        }

        GeoRadiusQueryPlan lower_bound = candidate;

        lower_bound.range_count = 1;
        radius_plan_set_cost(&lower_bound, output_record_size);

        /*
         * Coarsening a conservative Morton cover cannot reduce its candidate set. Once even
         * a hypothetical one-range coarser cover cannot beat the best cost, every remaining
         * depth is dominated and its recursive cover construction can be skipped.
         */
        if (selected && lower_bound.estimated_cost >= plan->estimated_cost) {
            break;
        }
    }

    return selected;
}

bool geo_index_search_radius_plan_append(const GeoIndex *index,
                                         const GeoRadiusQueryPlan *plan,
                                         GeoSearchResult *result,
                                         size_t *matched_count,
                                         GeoSearchStats *stats)
{
    return geo_index_search_radius_plan_append_filtered(index, plan, result, matched_count, stats, NULL, NULL, NULL);
}

bool geo_index_search_radius_plan_append_filtered(const GeoIndex *index,
                                                  const GeoRadiusQueryPlan *plan,
                                                  GeoSearchResult *result,
                                                  size_t *matched_count,
                                                  GeoSearchStats *stats,
                                                  GeoRecordBitFilter bit_filter,
                                                  GeoRecordFilter filter,
                                                  void *filter_context)
{
    reset_stats(stats);

    if (!valid_index(index) || !plan || plan->range_count <= 0 || (!result && !matched_count)) {
        return false;
    }

    double start_time = stats ? geo_get_time_ms() : 0.0;
    bool succeeded = scan_radius(index,
                                 plan->ranges,
                                 plan->bounds_cached && plan->bounds_index == index ? plan->range_begins : NULL,
                                 plan->bounds_cached && plan->bounds_index == index ? plan->range_ends : NULL,
                                 plan->range_count,
                                 &plan->radius_query,
                                 plan->strategy,
                                 result,
                                 NULL,
                                 false,
                                 matched_count,
                                 stats,
                                 bit_filter,
                                 filter,
                                 filter_context);

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start_time;
    }

    return succeeded;
}

bool geo_index_search_radius_plan_ids_append_filtered(const GeoIndex *index,
                                                      const GeoRadiusQueryPlan *plan,
                                                      GeoIdResult *result,
                                                      bool allow_growth,
                                                      size_t *matched_count,
                                                      GeoSearchStats *stats,
                                                      GeoRecordBitFilter bit_filter,
                                                      GeoRecordFilter filter,
                                                      void *filter_context)
{
    reset_stats(stats);

    if (!valid_index(index) || !plan || plan->range_count <= 0 || !result || (result->capacity && !result->ids)) {
        return false;
    }

    double start_time = stats ? geo_get_time_ms() : 0.0;
    bool succeeded = scan_radius(index,
                                 plan->ranges,
                                 plan->bounds_cached && plan->bounds_index == index ? plan->range_begins : NULL,
                                 plan->bounds_cached && plan->bounds_index == index ? plan->range_ends : NULL,
                                 plan->range_count,
                                 &plan->radius_query,
                                 plan->strategy,
                                 NULL,
                                 result,
                                 allow_growth,
                                 matched_count,
                                 stats,
                                 bit_filter,
                                 filter,
                                 filter_context);

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start_time;
    }

    return succeeded;
}

static bool search_radius_execute(const GeoIndex *index,
                                  double latitude,
                                  double longitude,
                                  double radius_km,
                                  GeoSearchResult *result,
                                  GeoSearchStats *stats,
                                  bool clear_result)
{
    reset_stats(stats);

    if (!valid_index(index) || !result) {
        return false;
    }

    double start_time = stats ? geo_get_time_ms() : 0.0;
    GeoRadiusQueryPlan plan;

    if (!geo_index_prepare_radius_query_for_index(index,
                                                  latitude,
                                                  longitude,
                                                  radius_km,
                                                  sizeof(GeoRecord),
                                                  &plan)) {
        return false;
    }

    if (clear_result) {
        geo_result_clear(result);
    }

    bool succeeded = geo_index_search_radius_plan_append(index, &plan, result, NULL, stats);

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start_time;
    }

    return succeeded;
}

bool geo_search_radius_reuse(const GeoIndex *index, double lat, double lng, double radius, GeoSearchResult *result, GeoSearchStats *stats)
{
    return search_radius_execute(index, lat, lng, radius, result, stats, true);
}

bool geo_index_search_radius_append(const GeoIndex *index,
                                    double lat,
                                    double lng,
                                    double radius_km,
                                    GeoSearchResult *result,
                                    GeoSearchStats *stats)
{
    return search_radius_execute(index, lat, lng, radius_km, result, stats, false);
}

bool geo_search_radius_count(const GeoIndex *index,
                             double lat,
                             double lng,
                             double radius,
                             size_t *count,
                             GeoSearchStats *stats)
{
    if (!valid_index(index) || !count) {
        return false;
    }

    reset_stats(stats);
    double start_time = stats ? geo_get_time_ms() : 0.0;
    *count = 0;
    GeoRadiusQueryPlan plan;

    bool succeeded = geo_index_prepare_radius_query_for_index(index, lat, lng, radius, 0, &plan) &&
                     geo_index_search_radius_plan_append(index, &plan, NULL, count, stats);

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start_time;
    }

    return succeeded;
}

static bool search_radius_ids_execute(const GeoIndex *index,
                                      double latitude,
                                      double longitude,
                                      double radius_km,
                                      GeoIdResult *result,
                                      bool allow_growth,
                                      GeoSearchStats *stats)
{
    reset_stats(stats);

    if (!valid_index(index) || !result || (result->capacity && !result->ids)) {
        return false;
    }

    double start_time = stats ? geo_get_time_ms() : 0.0;
    GeoRadiusQueryPlan plan;
    if (!geo_index_prepare_radius_query_for_index(index,
                                                  latitude,
                                                  longitude,
                                                  radius_km,
                                                  sizeof(uint64_t),
                                                  &plan)) {
        return false;
    }

    result->count = 0;
    size_t matched = 0;
    bool succeeded = scan_radius(index,
                                 plan.ranges,
                                 plan.bounds_cached && plan.bounds_index == index ? plan.range_begins : NULL,
                                 plan.bounds_cached && plan.bounds_index == index ? plan.range_ends : NULL,
                                 plan.range_count,
                                 &plan.radius_query,
                                 plan.strategy,
                                 NULL,
                                 result,
                                 allow_growth,
                                 &matched,
                                 stats,
                                 NULL,
                                 NULL,
                                 NULL);

    if (!succeeded) {
        result->count = 0;
    }

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start_time;
    }

    return succeeded && result->count == matched;
}

bool geo_search_radius_ids_reuse(const GeoIndex *index,
                                 double lat,
                                 double lng,
                                 double radius_km,
                                 GeoIdResult *result,
                                 GeoSearchStats *stats)
{
    return search_radius_ids_execute(index, lat, lng, radius_km, result, true, stats);
}

bool geo_search_radius_ids_into(const GeoIndex *index,
                                double lat,
                                double lng,
                                double radius_km,
                                uint64_t *ids,
                                size_t capacity,
                                size_t *count,
                                GeoSearchStats *stats)
{
    if (!count || (capacity && !ids)) {
        return false;
    }

    GeoIdResult result = {
        .ids = ids,
        .count = 0,
        .capacity = capacity,
    };

    if (!search_radius_ids_execute(index, lat, lng, radius_km, &result, false, stats)) {
        return false;
    }

    *count = result.count;

    return true;
}

GeoSearchResult *geo_search_radius(const GeoIndex *index, double lat, double lng, double radius, GeoSearchStats *stats)
{
    GeoSearchResult *result = geo_result_create(64);

    if (!result) {
        return NULL;
    }

    if (!geo_search_radius_reuse(index, lat, lng, radius, result, stats)) {
        geo_result_destroy(result);

        return NULL;
    }

    return result;
}

bool geo_index_prepare_bbox_query(double min_latitude,
                                  double max_latitude,
                                  double min_longitude,
                                  double max_longitude,
                                  GeoBboxQueryPlan *plan)
{
    if (!plan || !geo_is_valid_lat(min_latitude) || !geo_is_valid_lat(max_latitude) ||
        !geo_is_valid_lng(min_longitude) || !geo_is_valid_lng(max_longitude) ||
        min_latitude > max_latitude) {
        return false;
    }

    GeoNormalizedBox boxes[2];
    int box_count = make_boxes(min_latitude, max_latitude, min_longitude, max_longitude, boxes);

    plan->min_latitude = min_latitude;
    plan->max_latitude = max_latitude;
    plan->min_longitude = min_longitude;
    plan->max_longitude = max_longitude;
    plan->range_count = cover_boxes_internal(boxes,
                                             box_count,
                                             plan->ranges,
                                             plan->needs_filter,
                                             GEO_QUERY_MAX_RANGES,
                                             24);

    return plan->range_count > 0;
}

bool geo_index_search_bbox_plan_append(const GeoIndex *index,
                                       const GeoBboxQueryPlan *plan,
                                       GeoSearchResult *result,
                                       size_t *matched_count,
                                       GeoSearchStats *stats)
{
    return geo_index_search_bbox_plan_append_filtered(index, plan, result, matched_count, stats, NULL, NULL, NULL);
}

bool geo_index_search_bbox_plan_append_filtered(const GeoIndex *index,
                                                const GeoBboxQueryPlan *plan,
                                                GeoSearchResult *result,
                                                size_t *matched_count,
                                                GeoSearchStats *stats,
                                                GeoRecordBitFilter bit_filter,
                                                GeoRecordFilter filter,
                                                void *filter_context)
{
    reset_stats(stats);

    if (!valid_index(index) || !plan || plan->range_count <= 0 || (!result && !matched_count)) {
        return false;
    }

    double start_time = stats ? geo_get_time_ms() : 0.0;
    bool succeeded = scan_bbox(index,
                               plan->ranges,
                               plan->needs_filter,
                               plan->range_count,
                               plan->min_latitude,
                               plan->max_latitude,
                               plan->min_longitude,
                               plan->max_longitude,
                               result,
                               matched_count,
                               stats,
                               bit_filter,
                               filter,
                               filter_context);

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start_time;
    }

    return succeeded;
}

static bool search_bbox_execute(const GeoIndex *index,
                                double min_latitude,
                                double max_latitude,
                                double min_longitude,
                                double max_longitude,
                                GeoSearchResult *result,
                                GeoSearchStats *stats,
                                bool clear_result)
{
    reset_stats(stats);

    if (!valid_index(index) || !result) {
        return false;
    }

    double start_time = stats ? geo_get_time_ms() : 0.0;
    GeoBboxQueryPlan plan;

    if (!geo_index_prepare_bbox_query(min_latitude,
                                      max_latitude,
                                      min_longitude,
                                      max_longitude,
                                      &plan)) {
        return false;
    }

    if (clear_result) {
        geo_result_clear(result);
    }

    bool succeeded = geo_index_search_bbox_plan_append(index, &plan, result, NULL, stats);

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start_time;
    }

    return succeeded;
}

bool geo_search_bbox_reuse(const GeoIndex *index,
                           double min_lat,
                           double max_lat,
                           double min_lng,
                           double max_lng,
                           GeoSearchResult *result,
                           GeoSearchStats *stats)
{
    return search_bbox_execute(index, min_lat, max_lat, min_lng, max_lng, result, stats, true);
}

bool geo_index_search_bbox_append(const GeoIndex *index,
                                  double min_lat,
                                  double max_lat,
                                  double min_lng,
                                  double max_lng,
                                  GeoSearchResult *result,
                                  GeoSearchStats *stats)
{
    return search_bbox_execute(index, min_lat, max_lat, min_lng, max_lng, result, stats, false);
}

bool geo_search_bbox_count(const GeoIndex *index,
                           double min_lat,
                           double max_lat,
                           double min_lng,
                           double max_lng,
                           size_t *count,
                           GeoSearchStats *stats)
{
    if (!valid_index(index) || !count) {
        return false;
    }

    reset_stats(stats);
    double start_time = stats ? geo_get_time_ms() : 0.0;
    GeoBboxQueryPlan plan;
    *count = 0;

    bool succeeded = geo_index_prepare_bbox_query(min_lat, max_lat, min_lng, max_lng, &plan) &&
                     geo_index_search_bbox_plan_append(index, &plan, NULL, count, stats);

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start_time;
    }

    return succeeded;
}

GeoSearchResult *geo_search_bbox(const GeoIndex *index,
                                 double min_lat,
                                 double max_lat,
                                 double min_lng,
                                 double max_lng,
                                 GeoSearchStats *stats)
{
    GeoSearchResult *result = geo_result_create(64);

    if (!result) {
        return NULL;
    }

    if (!geo_search_bbox_reuse(index, min_lat, max_lat, min_lng, max_lng, result, stats)) {
        geo_result_destroy(result);

        return NULL;
    }

    return result;
}

// =============================================================================
// K-nearest-neighbor best-first traversal
// =============================================================================

static double knn_longitude_delta(double longitude, double minimum, double maximum)
{
    if (longitude >= minimum && longitude <= maximum) {
        return 0.0;
    }

    double to_minimum = fabs(geo_wrap_lng(longitude - minimum));
    double to_maximum = fabs(geo_wrap_lng(longitude - maximum));

    return geo_to_radians(fmin(to_minimum, to_maximum));
}

static double knn_cell_lower_bound(double latitude,
                                   double longitude,
                                   uint32_t minimum_latitude,
                                   uint32_t minimum_longitude,
                                   uint64_t cell_size)
{
    uint32_t maximum_latitude = (uint32_t) ((uint64_t) minimum_latitude + cell_size - 1);
    uint32_t maximum_longitude = (uint32_t) ((uint64_t) minimum_longitude + cell_size - 1);
    double minimum_latitude_radians = geo_to_radians(geo_denormalize_lat(minimum_latitude));
    double maximum_latitude_radians = geo_to_radians(geo_denormalize_lat(maximum_latitude));
    double minimum_longitude_degrees = geo_denormalize_lng(minimum_longitude);
    double maximum_longitude_degrees = geo_denormalize_lng(maximum_longitude);
    double query_latitude = geo_to_radians(latitude);
    double longitude_delta = knn_longitude_delta(longitude, minimum_longitude_degrees, maximum_longitude_degrees);
    double latitude_factor = sin(query_latitude);
    double longitude_factor = cos(query_latitude) * cos(longitude_delta);
    double maximum_dot = fmax(latitude_factor * sin(minimum_latitude_radians) +
                                  longitude_factor * cos(minimum_latitude_radians),
                              latitude_factor * sin(maximum_latitude_radians) +
                                  longitude_factor * cos(maximum_latitude_radians));

    double stationary_latitude = atan2(latitude_factor, longitude_factor);

    if (stationary_latitude >= minimum_latitude_radians && stationary_latitude <= maximum_latitude_radians) {
        double stationary_dot = latitude_factor * sin(stationary_latitude) +
                                longitude_factor * cos(stationary_latitude);

        maximum_dot = fmax(maximum_dot, stationary_dot);
    }

    maximum_dot = fmax(-1.0, fmin(1.0, maximum_dot));

    return GEO_EARTH_RADIUS_KM * acos(maximum_dot);
}

static bool knn_queue_reserve(GeoKnnQueue *queue, size_t capacity)
{
    if (capacity <= queue->capacity) {
        return true;
    }

    size_t bytes;

    if (!size_mul(capacity, sizeof(*queue->cells), &bytes)) {
        return false;
    }

    GeoKnnCell *cells = realloc(queue->cells, bytes);

    if (!cells) {
        return false;
    }

    queue->cells = cells;
    queue->capacity = capacity;

    return true;
}

bool geo_knn_workspace_reserve(GeoKnnWorkspace *workspace, size_t neighbor_capacity, size_t cell_capacity)
{
    if (!workspace) {
        return false;
    }

    if (neighbor_capacity > workspace->heap_capacity) {
        size_t bytes;

        if (!size_mul(neighbor_capacity, sizeof(*workspace->heap), &bytes)) {
            return false;
        }

        GeoRankedRecord *heap = realloc(workspace->heap, bytes);

        if (!heap) {
            return false;
        }

        workspace->heap = heap;
        workspace->heap_capacity = neighbor_capacity;
    }

    return knn_queue_reserve(&workspace->queue, cell_capacity);
}

GeoKnnWorkspace *geo_knn_workspace_create(size_t neighbor_capacity, size_t cell_capacity)
{
    GeoKnnWorkspace *workspace = calloc(1, sizeof(*workspace));

    if (!workspace) {
        return NULL;
    }

    if (!geo_knn_workspace_reserve(workspace, neighbor_capacity, cell_capacity)) {
        geo_knn_workspace_destroy(workspace);

        return NULL;
    }

    return workspace;
}

void geo_knn_workspace_destroy(GeoKnnWorkspace *workspace)
{
    if (workspace) {
        free(workspace->queue.cells);
        free(workspace->heap);
        free(workspace);
    }
}

static bool knn_queue_push(GeoKnnQueue *queue, GeoKnnCell cell)
{
    if (queue->count == queue->capacity) {
        size_t capacity = queue->capacity ? queue->capacity * 2 : 64;

        if (capacity < queue->capacity || !knn_queue_reserve(queue, capacity)) {
            return false;
        }
    }

    size_t position = queue->count++;

    while (position) {
        size_t parent = (position - 1) >> 1;

        if (queue->cells[parent].lower_bound <= cell.lower_bound) {
            break;
        }

        queue->cells[position] = queue->cells[parent];
        position = parent;
    }

    queue->cells[position] = cell;

    return true;
}

static GeoKnnCell knn_queue_pop(GeoKnnQueue *queue)
{
    GeoKnnCell nearest = queue->cells[0];
    GeoKnnCell replacement = queue->cells[--queue->count];
    size_t position = 0;

    while (position * 2 + 1 < queue->count) {
        size_t left = position * 2 + 1;
        size_t right = left + 1;
        size_t smallest = right < queue->count &&
                                  queue->cells[right].lower_bound < queue->cells[left].lower_bound
                              ? right
                              : left;

        if (replacement.lower_bound <= queue->cells[smallest].lower_bound) {
            break;
        }

        queue->cells[position] = queue->cells[smallest];
        position = smallest;
    }

    if (queue->count) {
        queue->cells[position] = replacement;
    }

    return nearest;
}

static void heap_up(GeoRankedRecord *heap, size_t pos)
{
    while (pos) {
        size_t parent = (pos - 1) >> 1;

        if (heap[parent].distance >= heap[pos].distance) {
            break;
        }

        GeoRankedRecord temporary = heap[parent];

        heap[parent] = heap[pos];
        heap[pos] = temporary;
        pos = parent;
    }
}

static void heap_down(GeoRankedRecord *heap, size_t count)
{
    size_t pos = 0;

    for (;;) {
        size_t left = pos * 2 + 1;

        if (left >= count) {
            break;
        }

        size_t right = left + 1;
        size_t largest = right < count && heap[right].distance > heap[left].distance ? right : left;

        if (heap[pos].distance >= heap[largest].distance) {
            break;
        }

        GeoRankedRecord temporary = heap[pos];

        heap[pos] = heap[largest];
        heap[largest] = temporary;
        pos = largest;
    }
}

static void knn_consider_record(GeoRankedRecord *heap,
                                size_t capacity,
                                size_t *count,
                                const GeoRecord *record,
                                double latitude,
                                double longitude,
                                double max_radius)
{
    GeoPoint point = geo_decode(record->z);
    double distance = geo_haversine_km(latitude, longitude, point.lat, point.lng);

    if (distance > max_radius) {
        return;
    }

    GeoRankedRecord item = {
        .record = *record,
        .distance = distance,
    };

    if (*count < capacity) {
        heap[*count] = item;
        heap_up(heap, (*count)++);
    } else if (capacity && distance < heap[0].distance) {
        heap[0] = item;
        heap_down(heap, *count);
    }
}

GeoSearchResult *geo_search_knn(const GeoIndex *index, double lat, double lng, size_t k, double max_radius, GeoSearchStats *stats)
{
    return geo_index_search_knn_filtered(index, lat, lng, k, max_radius, stats, NULL, NULL);
}

bool geo_index_search_knn_sources(const GeoKnnSource *sources,
                                  size_t source_count,
                                  double lat,
                                  double lng,
                                  size_t k,
                                  double max_radius,
                                  GeoSearchResult *result,
                                  GeoKnnWorkspace *workspace,
                                  GeoSearchStats *stats)
{
    reset_stats(stats);

    if ((!sources && source_count) || !geo_is_valid_point(lat, lng) || !k || !isfinite(max_radius) || max_radius < 0.0 ||
        !result || !workspace) {
        return false;
    }

    size_t total_records = 0;

    for (size_t source_index = 0; source_index < source_count; ++source_index) {
        if (!valid_index(sources[source_index].index) ||
            sources[source_index].index->count > SIZE_MAX - total_records) {
            return false;
        }

        total_records += sources[source_index].index->count;
    }

    size_t capacity = k < total_records ? k : total_records;

    if (!geo_knn_workspace_reserve(workspace, capacity, 64) || !geo_result_reserve(result, capacity ? capacity : 1)) {
        return false;
    }

    geo_result_clear(result);
    workspace->queue.count = 0;

    double start_time = stats ? geo_get_time_ms() : 0.0;
    GeoKnnQueue *queue = &workspace->queue;
    GeoRankedRecord *heap = workspace->heap;
    size_t heap_count = 0;
    uint64_t records_scanned = 0;
    uint64_t cells_visited = 0;
    bool succeeded = true;

    for (size_t source_index = 0; succeeded && source_index < source_count; ++source_index) {
        const GeoIndex *index = sources[source_index].index;

        if (!index->count) {
            continue;
        }

        succeeded = knn_queue_push(queue,
                                   (GeoKnnCell) {
                                       .min_latitude = 0,
                                       .min_longitude = 0,
                                       .size = UINT64_C(1) << 32,
                                       .morton_prefix = 0,
                                       .begin = 0,
                                       .end = index->count,
                                       .source_index = source_index,
                                       .lower_bound = 0.0,
                                   });
    }

    while (succeeded && queue->count) {
        GeoKnnCell cell = knn_queue_pop(queue);
        const GeoKnnSource *source = sources + cell.source_index;
        const GeoIndex *index = source->index;
        double distance_limit = heap_count == capacity && capacity ? fmin(max_radius, heap[0].distance) : max_radius;

        if (cell.lower_bound > distance_limit + 1e-9) {
            break;
        }

        cells_visited++;

        size_t record_count = cell.end - cell.begin;

        if (record_count <= GEO_KNN_LEAF_RECORDS || cell.size == 1) {
            records_scanned += record_count;

            for (size_t position = cell.begin; position < cell.end; ++position) {
                if (source->filter && !source->filter(index->records + position, source->filter_context)) {
                    continue;
                }

                knn_consider_record(heap,
                                    capacity,
                                    &heap_count,
                                    index->records + position,
                                    lat,
                                    lng,
                                    max_radius);
            }

            continue;
        }

        uint64_t child_size = cell.size >> 1;
        unsigned remaining_bits = (unsigned) __builtin_ctzll(child_size) * 2U;
        size_t boundaries[5];

        boundaries[0] = cell.begin;
        boundaries[4] = cell.end;

        for (unsigned quadrant = 1; quadrant < 4; ++quadrant) {
            uint64_t child_prefix = (cell.morton_prefix << 2) | quadrant;
            uint64_t child_minimum = remaining_bits ? child_prefix << remaining_bits : child_prefix;

            boundaries[quadrant] = index_lower_bound(index, child_minimum);
        }

        for (unsigned quadrant = 0; quadrant < 4; ++quadrant) {
            if (boundaries[quadrant] == boundaries[quadrant + 1]) {
                continue;
            }

            uint32_t child_minimum_latitude =
                (uint32_t) ((uint64_t) cell.min_latitude + ((quadrant & 1U) ? child_size : 0));

            uint32_t child_minimum_longitude =
                (uint32_t) ((uint64_t) cell.min_longitude + ((quadrant & 2U) ? child_size : 0));
            double lower_bound = knn_cell_lower_bound(lat,
                                                      lng,
                                                      child_minimum_latitude,
                                                      child_minimum_longitude,
                                                      child_size);

            distance_limit = heap_count == capacity && capacity ? fmin(max_radius, heap[0].distance) : max_radius;

            if (lower_bound > distance_limit + 1e-9) {
                continue;
            }

            succeeded = knn_queue_push(queue,
                                       (GeoKnnCell) {
                                           .min_latitude = child_minimum_latitude,
                                           .min_longitude = child_minimum_longitude,
                                           .size = child_size,
                                           .morton_prefix = (cell.morton_prefix << 2) | quadrant,
                                           .begin = boundaries[quadrant],
                                           .end = boundaries[quadrant + 1],
                                           .source_index = cell.source_index,
                                           .lower_bound = lower_bound,
                                       });

            if (!succeeded) {
                break;
            }
        }
    }

    if (!succeeded) {
        queue->count = 0;
        geo_result_clear(result);

        return false;
    }

    qsort(heap, heap_count, sizeof(*heap), compare_ranked);

    for (size_t i = 0; i < heap_count; ++i) {
        result->results[result->count++] = heap[i].record;
    }

    if (stats) {
        stats->records_scanned = records_scanned;
        stats->records_matched = heap_count;
        stats->ranges_checked = cells_visited;
        stats->search_time_ms = geo_get_time_ms() - start_time;
    }

    queue->count = 0;

    return true;
}

bool geo_index_search_knn_reuse_filtered(const GeoIndex *index,
                                         double lat,
                                         double lng,
                                         size_t k,
                                         double max_radius,
                                         GeoSearchResult *result,
                                         GeoKnnWorkspace *workspace,
                                         GeoSearchStats *stats,
                                         GeoRecordFilter filter,
                                         void *filter_context)
{
    const GeoKnnSource source = {
        .index = index,
        .filter = filter,
        .filter_context = filter_context,
    };

    return geo_index_search_knn_sources(&source,
                                        1,
                                        lat,
                                        lng,
                                        k,
                                        max_radius,
                                        result,
                                        workspace,
                                        stats);
}

bool geo_search_knn_reuse(const GeoIndex *index,
                          double lat,
                          double lng,
                          size_t k,
                          double max_radius,
                          GeoSearchResult *result,
                          GeoKnnWorkspace *workspace,
                          GeoSearchStats *stats)
{
    return geo_index_search_knn_reuse_filtered(index,
                                               lat,
                                               lng,
                                               k,
                                               max_radius,
                                               result,
                                               workspace,
                                               stats,
                                               NULL,
                                               NULL);
}

GeoSearchResult *geo_index_search_knn_filtered(const GeoIndex *index,
                                               double lat,
                                               double lng,
                                               size_t k,
                                               double max_radius,
                                               GeoSearchStats *stats,
                                               GeoRecordFilter filter,
                                               void *filter_context)
{
    reset_stats(stats);

    if (!valid_index(index) || !geo_is_valid_point(lat, lng) || !k || !isfinite(max_radius) || max_radius < 0.0) {
        return NULL;
    }

    size_t capacity = k < index->count ? k : index->count;
    GeoSearchResult *result = geo_result_create(capacity ? capacity : 1);
    GeoKnnWorkspace *workspace = geo_knn_workspace_create(capacity, 64);

    if (!result || !workspace ||
        !geo_index_search_knn_reuse_filtered(index,
                                             lat,
                                             lng,
                                             k,
                                             max_radius,
                                             result,
                                             workspace,
                                             stats,
                                             filter,
                                             filter_context)) {
        geo_knn_workspace_destroy(workspace);
        geo_result_destroy(result);

        return NULL;
    }

    geo_knn_workspace_destroy(workspace);

    return result;
}

// =============================================================================
// Validation and utility functions
// =============================================================================

double geo_get_time_ms(void)
{
    struct timespec ts;

    return clock_gettime(CLOCK_MONOTONIC, &ts) ? 0.0 : (double) ts.tv_sec * 1000.0 + (double) ts.tv_nsec / 1000000.0;
}

bool geo_is_valid_lat(double lat)
{
    return isfinite(lat) && lat >= GEO_MIN_LAT && lat <= GEO_MAX_LAT;
}

bool geo_is_valid_lng(double lng)
{
    return isfinite(lng) && lng >= GEO_MIN_LNG && lng <= GEO_MAX_LNG;
}

bool geo_is_valid_point(double lat, double lng)
{
    return geo_is_valid_lat(lat) && geo_is_valid_lng(lng);
}

double geo_clamp_lat(double lat)
{
    return fmax(GEO_MIN_LAT, fmin(GEO_MAX_LAT, lat));
}

double geo_clamp_lng(double lng)
{
    return fmax(GEO_MIN_LNG, fmin(GEO_MAX_LNG, lng));
}

double geo_wrap_lng(double lng)
{
    if (!isfinite(lng)) {
        return lng;
    }

    double wrapped = fmod(lng + 180.0, 360.0);

    if (wrapped < 0.0) {
        wrapped += 360.0;
    }

    wrapped -= 180.0;

    return wrapped == -180.0 && lng > 0.0 ? 180.0 : wrapped;
}
