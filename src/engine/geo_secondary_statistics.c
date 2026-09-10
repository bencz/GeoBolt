#include "geo_secondary_statistics.h"

#include "geo_db_format.h"

#include <stdlib.h>
#include <string.h>

#define GEO_SECONDARY_HISTOGRAM_MAGIC UINT32_C(0x31484247)
#define GEO_SECONDARY_HISTOGRAM_METADATA_MAGIC UINT32_C(0x314d4247)
#define GEO_SECONDARY_HISTOGRAM_VERSION 1U
#define GEO_SECONDARY_HISTOGRAM_VALUE_SIZE 32U
#define GEO_SECONDARY_HISTOGRAM_METADATA_HEADER_SIZE 40U
#define GEO_SECONDARY_HISTOGRAM_DIRECTORY_ENTRY_SIZE 8U

static const unsigned char GEO_SECONDARY_HISTOGRAM_PREFIX[] = {
    'i', 'n', 'd', 'e', 'x', '-', 'h', 'i', 's', 't', '/',
};
static const unsigned char GEO_SECONDARY_HISTOGRAM_METADATA_PREFIX[] = {
    'i', 'n', 'd', 'e', 'x', '-', 'h', 'i', 's', 't', '-', 'm', 'e', 't', 'a', '/',
};

static uint16_t statistics_load_u16_le(const unsigned char *input)
{
    return (uint16_t) input[0] | (uint16_t) ((uint16_t) input[1] << 8U);
}

static uint32_t statistics_load_u32_le(const unsigned char *input)
{
    return (uint32_t) input[0] |
           (uint32_t) input[1] << 8U |
           (uint32_t) input[2] << 16U |
           (uint32_t) input[3] << 24U;
}

static uint64_t statistics_load_u64_le(const unsigned char *input)
{
    return (uint64_t) statistics_load_u32_le(input) | (uint64_t) statistics_load_u32_le(input + 4U) << 32U;
}

static void statistics_store_u16_le(unsigned char *output, uint16_t value)
{
    output[0] = (unsigned char) value;
    output[1] = (unsigned char) (value >> 8U);
}

static void statistics_store_u32_le(unsigned char *output, uint32_t value)
{
    output[0] = (unsigned char) value;
    output[1] = (unsigned char) (value >> 8U);
    output[2] = (unsigned char) (value >> 16U);
    output[3] = (unsigned char) (value >> 24U);
}

static void statistics_store_u64_le(unsigned char *output, uint64_t value)
{
    statistics_store_u32_le(output, (uint32_t) value);
    statistics_store_u32_le(output + 4U, (uint32_t) (value >> 32U));
}

static uint64_t statistics_checksum(const unsigned char *data, size_t size, size_t checksum_offset)
{
    uint64_t checksum = UINT64_C(1469598103934665603);

    for (size_t index = 0U; index < size; ++index) {
        unsigned char byte = index >= checksum_offset && index < checksum_offset + sizeof(uint64_t) ? 0U : data[index];

        checksum ^= byte;
        checksum *= UINT64_C(1099511628211);
    }

    return checksum;
}

static void statistics_histogram_metadata_key(
    unsigned char output[sizeof(GEO_SECONDARY_HISTOGRAM_METADATA_PREFIX) + 8U],
    uint64_t index_id)
{
    memcpy(output, GEO_SECONDARY_HISTOGRAM_METADATA_PREFIX, sizeof(GEO_SECONDARY_HISTOGRAM_METADATA_PREFIX));
    geo_db_store_u64_be(output + sizeof(GEO_SECONDARY_HISTOGRAM_METADATA_PREFIX), index_id);
}

static void statistics_histogram_key(unsigned char output[sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX) + 10U],
                                     uint64_t index_id,
                                     uint16_t bin)
{
    memcpy(output, GEO_SECONDARY_HISTOGRAM_PREFIX, sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX));
    geo_db_store_u64_be(output + sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX), index_id);
    output[sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX) + 8U] = (unsigned char) (bin >> 8U);
    output[sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX) + 9U] = (unsigned char) bin;
}

static void statistics_histogram_range_key(unsigned char output[sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX) + 8U],
                                           uint64_t index_id)
{
    memcpy(output, GEO_SECONDARY_HISTOGRAM_PREFIX, sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX));
    geo_db_store_u64_be(output + sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX), index_id);
}

void geo_secondary_histogram_destroy(GeoSecondaryHistogram *histogram)
{
    if (!histogram) {
        return;
    }

    free(histogram->boundaries);
    *histogram = (GeoSecondaryHistogram) { 0 };
}

int geo_secondary_compare_encoded(const unsigned char *first,
                                  size_t first_size,
                                  const unsigned char *second,
                                  size_t second_size)
{
    size_t common_size = first_size < second_size ? first_size : second_size;
    int comparison = common_size ? memcmp(first, second, common_size) : 0;

    if (comparison != 0) {
        return comparison;
    }

    return (first_size > second_size) - (first_size < second_size);
}

bool geo_secondary_encoded_key_valid(GeoDatabaseIndexType type,
                                     const unsigned char *value,
                                     size_t value_size)
{
    if (!value) {
        return false;
    }
    if (type == GEO_DATABASE_INDEX_BOOL) {
        return value_size == 1U && value[0] <= 1U;
    }
    if (type >= GEO_DATABASE_INDEX_INT64 && type <= GEO_DATABASE_INDEX_DATETIME) {
        return value_size == 8U;
    }
    if ((type != GEO_DATABASE_INDEX_STRING && type != GEO_DATABASE_INDEX_BYTES) ||
        value_size < 2U || value_size > GEO_DATABASE_MAX_INDEXED_VALUE_SIZE * 2U + 2U ||
        value[value_size - 2U] != 0U || value[value_size - 1U] != 0U) {
        return false;
    }

    for (size_t index = 0U; index + 2U < value_size; ++index) {
        if (value[index] == 0U) {
            if (value[index + 1U] != UINT8_MAX) {
                return false;
            }
            index++;
        }
    }

    return true;
}

static bool statistics_boundary_reserve(unsigned char **boundaries, size_t *capacity, size_t required)
{
    if (required <= *capacity) {
        return true;
    }

    size_t updated_capacity = *capacity ? *capacity : 256U;

    while (updated_capacity < required) {
        if (updated_capacity > SIZE_MAX / 2U) {
            return false;
        }

        updated_capacity *= 2U;
    }

    unsigned char *updated = realloc(*boundaries, updated_capacity);

    if (!updated) {
        return false;
    }

    *boundaries = updated;
    *capacity = updated_capacity;
    return true;
}

static GeoDatabaseStatus statistics_count_distinct_limit(GeoRocksDatabase *rocks,
                                                         uint64_t index_id,
                                                         GeoDatabaseIndexType type,
                                                         uint32_t *distinct_count)
{
    unsigned char seek_key[8];
    unsigned char previous[GEO_DATABASE_MAX_INDEXED_VALUE_SIZE * 2U + 2U];
    size_t previous_size = 0U;
    GeoRocksStatus rocks_status;
    GeoRocksIterator *iterator = geo_rocks_iterator_create(rocks,
                                                           GEO_ROCKS_CF_SECONDARY_INDEX,
                                                           NULL,
                                                           &rocks_status);

    if (!iterator) {
        return geo_db_status_from_rocks(&rocks_status);
    }

    geo_db_store_u64_be(seek_key, index_id);
    geo_rocks_iterator_seek(iterator, seek_key, sizeof(seek_key));
    GeoDatabaseStatus status = GEO_DATABASE_OK;
    uint32_t count = 0U;

    while (geo_rocks_iterator_valid(iterator) && count <= GEO_SECONDARY_HISTOGRAM_MAX_BINS) {
        size_t key_size = 0U;
        const unsigned char *key = geo_rocks_iterator_key(iterator, &key_size);

        if (!key || key_size < 16U || geo_db_load_u64_be(key) != index_id) {
            break;
        }

        const unsigned char *value = key + 8U;
        size_t value_size = key_size - 16U;

        if (!geo_secondary_encoded_key_valid(type, value, value_size)) {
            status = GEO_DATABASE_CORRUPTION;
            break;
        }

        if (count == 0U || geo_secondary_compare_encoded(previous, previous_size, value, value_size) != 0) {
            memcpy(previous, value, value_size);
            previous_size = value_size;
            count++;
        }

        geo_rocks_iterator_next(iterator);
    }

    if (status == GEO_DATABASE_OK && !geo_rocks_iterator_status(iterator, &rocks_status)) {
        status = geo_db_status_from_rocks(&rocks_status);
    }

    geo_rocks_iterator_destroy(iterator);
    *distinct_count = count;
    return status;
}

static bool statistics_finalize_bin(GeoSecondaryHistogram *histogram,
                                    size_t *boundary_capacity,
                                    const unsigned char *upper_boundary,
                                    size_t upper_boundary_size,
                                    uint64_t count,
                                    uint32_t distinct_values)
{
    if (histogram->bin_count >= GEO_SECONDARY_HISTOGRAM_MAX_BINS || upper_boundary_size > UINT32_MAX ||
        histogram->boundary_size > UINT32_MAX - upper_boundary_size || !count || !distinct_values) {
        return false;
    }

    size_t required = (size_t) histogram->boundary_size + upper_boundary_size;

    if (!statistics_boundary_reserve(&histogram->boundaries, boundary_capacity, required)) {
        return false;
    }

    uint32_t bin = histogram->bin_count;

    histogram->boundary_offsets[bin] = histogram->boundary_size;
    memcpy(histogram->boundaries + histogram->boundary_size, upper_boundary, upper_boundary_size);
    histogram->boundary_size += (uint32_t) upper_boundary_size;
    histogram->boundary_offsets[bin + 1U] = histogram->boundary_size;
    histogram->counts[bin] = count;
    histogram->distinct_values[bin] = distinct_values;
    histogram->bin_count++;
    return true;
}

static GeoDatabaseStatus statistics_build_bins(GeoRocksDatabase *rocks,
                                               uint64_t index_id,
                                               GeoDatabaseIndexType type,
                                               uint64_t entry_count,
                                               bool exact_distinct,
                                               GeoSecondaryHistogram *histogram)
{
    unsigned char seek_key[8];
    unsigned char group_value[GEO_DATABASE_MAX_INDEXED_VALUE_SIZE * 2U + 2U];
    size_t group_value_size = 0U;
    uint64_t group_count = 0U;
    uint64_t bin_count = 0U;
    uint32_t bin_distinct = 0U;
    uint64_t target = exact_distinct ? 1U : entry_count / GEO_SECONDARY_HISTOGRAM_MAX_BINS;

    if (!exact_distinct && entry_count % GEO_SECONDARY_HISTOGRAM_MAX_BINS) {
        target++;
    }
    size_t boundary_capacity = 0U;
    GeoRocksStatus rocks_status;
    GeoRocksIterator *iterator = geo_rocks_iterator_create(rocks,
                                                           GEO_ROCKS_CF_SECONDARY_INDEX,
                                                           NULL,
                                                           &rocks_status);

    if (!iterator) {
        return geo_db_status_from_rocks(&rocks_status);
    }

    geo_db_store_u64_be(seek_key, index_id);
    geo_rocks_iterator_seek(iterator, seek_key, sizeof(seek_key));
    GeoDatabaseStatus status = GEO_DATABASE_OK;

    while (geo_rocks_iterator_valid(iterator)) {
        size_t key_size = 0U;
        const unsigned char *key = geo_rocks_iterator_key(iterator, &key_size);

        if (!key || key_size < 16U || geo_db_load_u64_be(key) != index_id) {
            break;
        }

        const unsigned char *value = key + 8U;
        size_t value_size = key_size - 16U;

        if (!geo_secondary_encoded_key_valid(type, value, value_size)) {
            status = GEO_DATABASE_CORRUPTION;
            break;
        }

        bool same_group = group_count &&
                          geo_secondary_compare_encoded(group_value, group_value_size, value, value_size) == 0;

        if (!same_group && group_count) {
            if (group_count > UINT64_MAX - bin_count) {
                status = GEO_DATABASE_CORRUPTION;
                break;
            }

            bin_count += group_count;

            if (bin_distinct != UINT32_MAX) {
                bin_distinct++;
            }
            bool should_finalize = histogram->bin_count + 1U < GEO_SECONDARY_HISTOGRAM_MAX_BINS &&
                                   (exact_distinct || bin_count >= target);

            if (should_finalize &&
                !statistics_finalize_bin(histogram,
                                         &boundary_capacity,
                                         group_value,
                                         group_value_size,
                                         bin_count,
                                         bin_distinct)) {
                status = GEO_DATABASE_OUT_OF_MEMORY;
                break;
            }

            if (should_finalize) {
                bin_count = 0U;
                bin_distinct = 0U;
            }

            group_count = 0U;
        }

        if (!same_group) {
            memcpy(group_value, value, value_size);
            group_value_size = value_size;
        }

        if (group_count == UINT64_MAX) {
            status = GEO_DATABASE_CORRUPTION;
            break;
        }

        group_count++;
        geo_rocks_iterator_next(iterator);
    }

    if (status == GEO_DATABASE_OK && !geo_rocks_iterator_status(iterator, &rocks_status)) {
        status = geo_db_status_from_rocks(&rocks_status);
    }

    if (status == GEO_DATABASE_OK && group_count) {
        if (group_count > UINT64_MAX - bin_count) {
            status = GEO_DATABASE_CORRUPTION;
        }

        if (status == GEO_DATABASE_OK) {
            bin_count += group_count;
        }

        if (bin_distinct != UINT32_MAX) {
            bin_distinct++;
        }

        if (status == GEO_DATABASE_OK &&
            !statistics_finalize_bin(histogram,
                                     &boundary_capacity,
                                     group_value,
                                     group_value_size,
                                     bin_count,
                                     bin_distinct)) {
            status = GEO_DATABASE_OUT_OF_MEMORY;
        }
    }

    geo_rocks_iterator_destroy(iterator);
    return status;
}

GeoDatabaseStatus geo_secondary_histogram_build(GeoRocksDatabase *rocks,
                                                uint64_t index_id,
                                                GeoDatabaseIndexType type,
                                                uint64_t entry_count,
                                                GeoSecondaryHistogram *histogram)
{
    if (!rocks || !index_id || !histogram) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    geo_secondary_histogram_destroy(histogram);

    if (!entry_count) {
        return GEO_DATABASE_OK;
    }

    uint32_t distinct_limit = 0U;
    GeoDatabaseStatus status = statistics_count_distinct_limit(rocks, index_id, type, &distinct_limit);

    if (status == GEO_DATABASE_OK) {
        status = statistics_build_bins(rocks,
                                       index_id,
                                       type,
                                       entry_count,
                                       distinct_limit <= GEO_SECONDARY_HISTOGRAM_MAX_BINS,
                                       histogram);
    }

    uint64_t counted = 0U;

    for (uint32_t bin = 0U; status == GEO_DATABASE_OK && bin < histogram->bin_count; ++bin) {
        if (histogram->counts[bin] > UINT64_MAX - counted) {
            status = GEO_DATABASE_CORRUPTION;
            break;
        }

        counted += histogram->counts[bin];
    }

    if (status == GEO_DATABASE_OK && (counted != entry_count || histogram->bin_count == 0U)) {
        status = GEO_DATABASE_CORRUPTION;
    }

    if (status != GEO_DATABASE_OK) {
        geo_secondary_histogram_destroy(histogram);
    }

    return status;
}

static size_t statistics_histogram_metadata_size(const GeoSecondaryHistogram *histogram)
{
    if (!histogram || histogram->bin_count > GEO_SECONDARY_HISTOGRAM_MAX_BINS) {
        return 0U;
    }

    return (size_t) histogram->bin_count * GEO_SECONDARY_HISTOGRAM_DIRECTORY_ENTRY_SIZE + histogram->boundary_size;
}

static bool statistics_histogram_encode_metadata(const GeoSecondaryHistogram *histogram,
                                                 unsigned char *output,
                                                 size_t output_size)
{
    size_t required = statistics_histogram_metadata_size(histogram);

    if ((!output && output_size) || output_size != required) {
        return false;
    }

    size_t boundary_output = (size_t) histogram->bin_count * GEO_SECONDARY_HISTOGRAM_DIRECTORY_ENTRY_SIZE;

    for (uint32_t bin = 0U; bin < histogram->bin_count; ++bin) {
        uint32_t boundary_length = histogram->boundary_offsets[bin + 1U] - histogram->boundary_offsets[bin];

        statistics_store_u32_le(output + (size_t) bin * GEO_SECONDARY_HISTOGRAM_DIRECTORY_ENTRY_SIZE,
                                boundary_length);
        statistics_store_u32_le(output + (size_t) bin * GEO_SECONDARY_HISTOGRAM_DIRECTORY_ENTRY_SIZE + 4U,
                                histogram->distinct_values[bin]);
        memcpy(output + boundary_output,
               histogram->boundaries + histogram->boundary_offsets[bin],
               boundary_length);
        boundary_output += boundary_length;
    }

    return boundary_output == required;
}

static bool statistics_histogram_decode_metadata(GeoSecondaryHistogram *histogram,
                                                 GeoDatabaseIndexType type,
                                                 uint32_t bin_count,
                                                 uint32_t boundary_size,
                                                 const unsigned char *metadata,
                                                 size_t metadata_size)
{
    if (!histogram || bin_count > GEO_SECONDARY_HISTOGRAM_MAX_BINS ||
        metadata_size != (size_t) bin_count * GEO_SECONDARY_HISTOGRAM_DIRECTORY_ENTRY_SIZE + boundary_size ||
        (metadata_size && !metadata)) {
        return false;
    }

    GeoSecondaryHistogram decoded = { .bin_count = bin_count, .boundary_size = boundary_size };

    if (boundary_size) {
        decoded.boundaries = malloc(boundary_size);

        if (!decoded.boundaries) {
            return false;
        }
    }

    size_t source_offset = (size_t) bin_count * GEO_SECONDARY_HISTOGRAM_DIRECTORY_ENTRY_SIZE;
    uint32_t destination_offset = 0U;

    for (uint32_t bin = 0U; bin < bin_count; ++bin) {
        uint32_t boundary_length = statistics_load_u32_le(
            metadata + (size_t) bin * GEO_SECONDARY_HISTOGRAM_DIRECTORY_ENTRY_SIZE);
        uint32_t distinct_values = statistics_load_u32_le(
            metadata + (size_t) bin * GEO_SECONDARY_HISTOGRAM_DIRECTORY_ENTRY_SIZE + 4U);

        if (!boundary_length || !distinct_values || boundary_length > boundary_size - destination_offset) {
            geo_secondary_histogram_destroy(&decoded);
            return false;
        }

        decoded.boundary_offsets[bin] = destination_offset;
        memcpy(decoded.boundaries + destination_offset, metadata + source_offset, boundary_length);
        decoded.distinct_values[bin] = distinct_values;
        destination_offset += boundary_length;
        source_offset += boundary_length;
        decoded.boundary_offsets[bin + 1U] = destination_offset;

        if (!geo_secondary_encoded_key_valid(type,
                                             decoded.boundaries + decoded.boundary_offsets[bin],
                                             boundary_length) ||
            (bin && geo_secondary_compare_encoded(decoded.boundaries + decoded.boundary_offsets[bin - 1U],
                                                  decoded.boundary_offsets[bin] - decoded.boundary_offsets[bin - 1U],
                                                  decoded.boundaries + decoded.boundary_offsets[bin],
                                                  boundary_length) >= 0)) {
            geo_secondary_histogram_destroy(&decoded);
            return false;
        }
    }

    if (source_offset != metadata_size || destination_offset != boundary_size) {
        geo_secondary_histogram_destroy(&decoded);
        return false;
    }

    geo_secondary_histogram_destroy(histogram);
    *histogram = decoded;
    return true;
}

GeoDatabaseStatus geo_secondary_histogram_put_metadata(GeoRocksBatch *batch,
                                                       uint64_t index_id,
                                                       GeoDatabaseIndexType type,
                                                       const GeoSecondaryHistogram *histogram)
{
    if (!batch || !index_id || !histogram || histogram->bin_count > GEO_SECONDARY_HISTOGRAM_MAX_BINS) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    size_t payload_size = statistics_histogram_metadata_size(histogram);

    if (payload_size > SIZE_MAX - GEO_SECONDARY_HISTOGRAM_METADATA_HEADER_SIZE) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    size_t value_size = GEO_SECONDARY_HISTOGRAM_METADATA_HEADER_SIZE + payload_size;
    unsigned char *value = calloc(1U, value_size);

    if (!value) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    statistics_store_u32_le(value, GEO_SECONDARY_HISTOGRAM_METADATA_MAGIC);
    statistics_store_u16_le(value + 4U, GEO_SECONDARY_HISTOGRAM_VERSION);
    statistics_store_u16_le(value + 6U, GEO_SECONDARY_HISTOGRAM_METADATA_HEADER_SIZE);
    statistics_store_u32_le(value + 8U, (uint32_t) type);
    statistics_store_u32_le(value + 12U, histogram->bin_count);
    statistics_store_u64_le(value + 16U, index_id);
    statistics_store_u32_le(value + 24U, histogram->boundary_size);

    if (!statistics_histogram_encode_metadata(histogram,
                                              value + GEO_SECONDARY_HISTOGRAM_METADATA_HEADER_SIZE,
                                              payload_size)) {
        free(value);
        return GEO_DATABASE_CORRUPTION;
    }

    statistics_store_u64_le(value + 32U, statistics_checksum(value, value_size, 32U));

    unsigned char key[sizeof(GEO_SECONDARY_HISTOGRAM_METADATA_PREFIX) + 8U];
    GeoRocksStatus rocks_status;

    statistics_histogram_metadata_key(key, index_id);
    bool succeeded = geo_rocks_batch_put(batch,
                                         GEO_ROCKS_CF_CATALOG,
                                         key,
                                         sizeof(key),
                                         value,
                                         value_size,
                                         &rocks_status);
    free(value);
    return succeeded ? GEO_DATABASE_OK : geo_db_status_from_rocks(&rocks_status);
}

GeoDatabaseStatus geo_secondary_histogram_load_metadata(GeoRocksDatabase *rocks,
                                                        uint64_t index_id,
                                                        GeoDatabaseIndexType type,
                                                        uint32_t bin_count,
                                                        uint32_t boundary_size,
                                                        GeoSecondaryHistogram *histogram)
{
    if (!rocks || !index_id || !histogram || !bin_count || !boundary_size ||
        bin_count > GEO_SECONDARY_HISTOGRAM_MAX_BINS) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    unsigned char key[sizeof(GEO_SECONDARY_HISTOGRAM_METADATA_PREFIX) + 8U];
    GeoRocksBuffer value = { 0 };
    GeoRocksStatus rocks_status;

    statistics_histogram_metadata_key(key, index_id);

    if (!geo_rocks_get(rocks, GEO_ROCKS_CF_CATALOG, NULL, key, sizeof(key), &value, &rocks_status)) {
        geo_rocks_buffer_release(&value);
        return rocks_status.code == GEO_ROCKS_NOT_FOUND ? GEO_DATABASE_CORRUPTION
                                                        : geo_db_status_from_rocks(&rocks_status);
    }

    const unsigned char *bytes = value.data;
    uint64_t payload_size = (uint64_t) bin_count * GEO_SECONDARY_HISTOGRAM_DIRECTORY_ENTRY_SIZE + boundary_size;
    uint64_t expected_size = GEO_SECONDARY_HISTOGRAM_METADATA_HEADER_SIZE + payload_size;
    bool valid = value.size == expected_size &&
                 statistics_load_u32_le(bytes) == GEO_SECONDARY_HISTOGRAM_METADATA_MAGIC &&
                 statistics_load_u16_le(bytes + 4U) == GEO_SECONDARY_HISTOGRAM_VERSION &&
                 statistics_load_u16_le(bytes + 6U) == GEO_SECONDARY_HISTOGRAM_METADATA_HEADER_SIZE &&
                 statistics_load_u32_le(bytes + 8U) == (uint32_t) type &&
                 statistics_load_u32_le(bytes + 12U) == bin_count &&
                 statistics_load_u64_le(bytes + 16U) == index_id &&
                 statistics_load_u32_le(bytes + 24U) == boundary_size &&
                 statistics_load_u32_le(bytes + 28U) == 0U &&
                 statistics_load_u64_le(bytes + 32U) == statistics_checksum(bytes, value.size, 32U);

    if (valid) {
        valid = statistics_histogram_decode_metadata(histogram,
                                                     type,
                                                     bin_count,
                                                     boundary_size,
                                                     bytes + GEO_SECONDARY_HISTOGRAM_METADATA_HEADER_SIZE,
                                                     (size_t) payload_size);
    }

    geo_rocks_buffer_release(&value);
    return valid ? GEO_DATABASE_OK : GEO_DATABASE_CORRUPTION;
}

GeoDatabaseStatus geo_secondary_histogram_put_count(GeoRocksBatch *batch,
                                                    uint64_t index_id,
                                                    uint32_t bin,
                                                    uint64_t count)
{
    if (!batch || !index_id || bin >= GEO_SECONDARY_HISTOGRAM_MAX_BINS) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    unsigned char key[sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX) + 10U];
    unsigned char value[GEO_SECONDARY_HISTOGRAM_VALUE_SIZE] = { 0 };
    GeoRocksStatus rocks_status;

    statistics_histogram_key(key, index_id, (uint16_t) bin);
    statistics_store_u32_le(value, GEO_SECONDARY_HISTOGRAM_MAGIC);
    statistics_store_u16_le(value + 4U, GEO_SECONDARY_HISTOGRAM_VERSION);
    statistics_store_u16_le(value + 6U, GEO_SECONDARY_HISTOGRAM_VALUE_SIZE);
    statistics_store_u64_le(value + 8U, index_id);
    statistics_store_u64_le(value + 16U, count);
    statistics_store_u64_le(value + 24U, statistics_checksum(value, sizeof(value), 24U));

    return geo_rocks_batch_put(batch,
                               GEO_ROCKS_CF_CATALOG,
                               key,
                               sizeof(key),
                               value,
                               sizeof(value),
                               &rocks_status)
               ? GEO_DATABASE_OK
               : geo_db_status_from_rocks(&rocks_status);
}

GeoDatabaseStatus geo_secondary_histogram_put_all_counts(GeoRocksBatch *batch,
                                                         uint64_t index_id,
                                                         const GeoSecondaryHistogram *histogram)
{
    if (!histogram) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    for (uint32_t bin = 0U; bin < histogram->bin_count; ++bin) {
        GeoDatabaseStatus status = geo_secondary_histogram_put_count(batch,
                                                                     index_id,
                                                                     bin,
                                                                     histogram->counts[bin]);

        if (status != GEO_DATABASE_OK) {
            return status;
        }
    }

    return GEO_DATABASE_OK;
}

GeoDatabaseStatus geo_secondary_histogram_delete(GeoRocksBatch *batch, uint64_t index_id)
{
    if (!batch || !index_id || index_id == UINT64_MAX) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    unsigned char begin[sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX) + 8U];
    unsigned char end[sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX) + 8U];
    unsigned char metadata_key[sizeof(GEO_SECONDARY_HISTOGRAM_METADATA_PREFIX) + 8U];
    GeoRocksStatus rocks_status;

    statistics_histogram_range_key(begin, index_id);
    statistics_histogram_range_key(end, index_id + 1U);
    statistics_histogram_metadata_key(metadata_key, index_id);

    if (!geo_rocks_batch_delete_range(batch,
                                      GEO_ROCKS_CF_CATALOG,
                                      begin,
                                      sizeof(begin),
                                      end,
                                      sizeof(end),
                                      &rocks_status) ||
        !geo_rocks_batch_delete(batch,
                                GEO_ROCKS_CF_CATALOG,
                                metadata_key,
                                sizeof(metadata_key),
                                &rocks_status)) {
        return geo_db_status_from_rocks(&rocks_status);
    }

    return GEO_DATABASE_OK;
}

GeoDatabaseStatus geo_secondary_histogram_load_counts(GeoRocksDatabase *rocks,
                                                      uint64_t index_id,
                                                      uint64_t entry_count,
                                                      GeoSecondaryHistogram *histogram)
{
    if (!rocks || !index_id || !histogram || !histogram->bin_count) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    uint64_t total = 0U;

    for (uint32_t bin = 0U; bin < histogram->bin_count; ++bin) {
        unsigned char key[sizeof(GEO_SECONDARY_HISTOGRAM_PREFIX) + 10U];
        GeoRocksBuffer value = { 0 };
        GeoRocksStatus rocks_status;

        statistics_histogram_key(key, index_id, (uint16_t) bin);

        if (!geo_rocks_get(rocks, GEO_ROCKS_CF_CATALOG, NULL, key, sizeof(key), &value, &rocks_status)) {
            geo_rocks_buffer_release(&value);
            return rocks_status.code == GEO_ROCKS_NOT_FOUND ? GEO_DATABASE_CORRUPTION
                                                            : geo_db_status_from_rocks(&rocks_status);
        }

        const unsigned char *bytes = value.data;
        bool valid = value.size == GEO_SECONDARY_HISTOGRAM_VALUE_SIZE &&
                     statistics_load_u32_le(bytes) == GEO_SECONDARY_HISTOGRAM_MAGIC &&
                     statistics_load_u16_le(bytes + 4U) == GEO_SECONDARY_HISTOGRAM_VERSION &&
                     statistics_load_u16_le(bytes + 6U) == GEO_SECONDARY_HISTOGRAM_VALUE_SIZE &&
                     statistics_load_u64_le(bytes + 8U) == index_id &&
                     statistics_load_u64_le(bytes + 24U) == statistics_checksum(bytes, value.size, 24U);
        uint64_t count = valid ? statistics_load_u64_le(bytes + 16U) : 0U;

        geo_rocks_buffer_release(&value);

        if (!valid || count > UINT64_MAX - total) {
            return GEO_DATABASE_CORRUPTION;
        }

        histogram->counts[bin] = count;
        total += count;
    }

    return total == entry_count ? GEO_DATABASE_OK : GEO_DATABASE_CORRUPTION;
}

uint32_t geo_secondary_histogram_find_bin(const GeoSecondaryHistogram *histogram,
                                         const unsigned char *encoded_value,
                                         size_t encoded_size)
{
    uint32_t low = 0U;
    uint32_t high = histogram->bin_count;

    while (low < high) {
        uint32_t middle = low + (high - low) / 2U;
        uint32_t begin = histogram->boundary_offsets[middle];
        uint32_t end = histogram->boundary_offsets[middle + 1U];
        int comparison = geo_secondary_compare_encoded(encoded_value,
                                                       encoded_size,
                                                       histogram->boundaries + begin,
                                                       end - begin);

        if (comparison <= 0) {
            high = middle;
        } else {
            low = middle + 1U;
        }
    }

    return low < histogram->bin_count ? low : histogram->bin_count - 1U;
}

static uint64_t statistics_saturating_add(uint64_t first, uint64_t second)
{
    return second > UINT64_MAX - first ? UINT64_MAX : first + second;
}

static uint64_t statistics_equality_estimate(const GeoSecondaryHistogram *histogram,
                                             const unsigned char *value,
                                             size_t value_size)
{
    if (!histogram->bin_count) {
        return 0U;
    }

    uint32_t bin = geo_secondary_histogram_find_bin(histogram, value, value_size);
    uint32_t begin = histogram->boundary_offsets[bin];
    uint32_t end = histogram->boundary_offsets[bin + 1U];
    int upper_comparison = geo_secondary_compare_encoded(value,
                                                         value_size,
                                                         histogram->boundaries + begin,
                                                         end - begin);

    if (histogram->distinct_values[bin] == 1U && upper_comparison != 0) {
        return histogram->counts[bin];
    }

    uint64_t estimate = histogram->counts[bin] / histogram->distinct_values[bin];
    return estimate ? estimate : 1U;
}

static uint64_t statistics_less_estimate(const GeoSecondaryHistogram *histogram,
                                         const unsigned char *value,
                                         size_t value_size,
                                         bool inclusive)
{
    if (!histogram->bin_count) {
        return 0U;
    }

    uint32_t bin = geo_secondary_histogram_find_bin(histogram, value, value_size);
    uint64_t estimate = 0U;

    for (uint32_t index = 0U; index < bin; ++index) {
        estimate = statistics_saturating_add(estimate, histogram->counts[index]);
    }

    uint32_t begin = histogram->boundary_offsets[bin];
    uint32_t end = histogram->boundary_offsets[bin + 1U];
    int upper_comparison = geo_secondary_compare_encoded(value,
                                                         value_size,
                                                         histogram->boundaries + begin,
                                                         end - begin);

    if (upper_comparison == 0) {
        uint64_t equality = statistics_equality_estimate(histogram, value, value_size);
        uint64_t partial = inclusive ? histogram->counts[bin] : histogram->counts[bin] - equality;
        return statistics_saturating_add(estimate, partial);
    }

    if (upper_comparison > 0 && bin + 1U == histogram->bin_count) {
        return statistics_saturating_add(estimate, histogram->counts[bin]);
    }

    if (histogram->distinct_values[bin] == 1U) {
        return estimate;
    }

    return statistics_saturating_add(estimate, histogram->counts[bin] / 2U);
}

static uint64_t statistics_total(const GeoSecondaryHistogram *histogram)
{
    uint64_t total = 0U;

    for (uint32_t bin = 0U; bin < histogram->bin_count; ++bin) {
        total = statistics_saturating_add(total, histogram->counts[bin]);
    }

    return total;
}

uint64_t geo_secondary_histogram_estimate(const GeoSecondaryHistogram *histogram,
                                         GeoDatabaseIndexOperator operation,
                                         const unsigned char *lower,
                                         size_t lower_size,
                                         const unsigned char *upper,
                                         size_t upper_size)
{
    if (!histogram || !lower || !histogram->bin_count) {
        return 0U;
    }

    uint64_t total = statistics_total(histogram);
    uint64_t less = statistics_less_estimate(histogram, lower, lower_size, false);
    uint64_t less_equal = statistics_less_estimate(histogram, lower, lower_size, true);
    uint64_t estimate;

    switch (operation) {
        case GEO_DATABASE_INDEX_EQUAL:
            estimate = statistics_equality_estimate(histogram, lower, lower_size);
            break;
        case GEO_DATABASE_INDEX_LESS:
            estimate = less;
            break;
        case GEO_DATABASE_INDEX_LESS_EQUAL:
            estimate = less_equal;
            break;
        case GEO_DATABASE_INDEX_GREATER:
            estimate = total - (less_equal < total ? less_equal : total);
            break;
        case GEO_DATABASE_INDEX_GREATER_EQUAL:
            estimate = total - (less < total ? less : total);
            break;
        case GEO_DATABASE_INDEX_BETWEEN: {
            uint64_t through_upper = statistics_less_estimate(histogram, upper, upper_size, true);
            estimate = through_upper > less ? through_upper - less : 0U;
            break;
        }
        default:
            return total;
    }

    return estimate < total ? estimate : total;
}
