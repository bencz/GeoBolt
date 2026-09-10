#include "geo_secondary_index.h"

#include "geobolt/geodoc.h"
#include "geo_secondary_statistics.h"

#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#define GEO_SECONDARY_DEFINITION_MAGIC UINT32_C(0x31494247)
#define GEO_SECONDARY_DEFINITION_VERSION 2U
#define GEO_SECONDARY_DEFINITION_HEADER_SIZE 64U
#define GEO_SECONDARY_BUILD_BATCH_ENTRIES 4096U
#define GEO_SECONDARY_STACK_KEY_SIZE 64U

static const unsigned char GEO_SECONDARY_CATALOG_PREFIX[] = { 'i', 'n', 'd', 'e', 'x', '/' };

typedef enum {
    GEO_SECONDARY_BUILDING = 1,
    GEO_SECONDARY_READY = 2,
} GeoSecondaryState;

typedef struct {
    GeoDatabaseIndexInfo info;
    GeoSecondaryHistogram histogram;
    size_t histogram_offset;
    GeoSecondaryState state;
} GeoSecondaryEntry;

typedef struct {
    unsigned char fixed[64];
    unsigned char *data;
    unsigned char *heap;
    size_t size;
    bool present;
} GeoSecondaryEncodedValue;

struct GeoSecondaryCatalog {
    pthread_rwlock_t lock;
    GeoSecondaryEntry *entries;
    size_t count;
    size_t capacity;
    size_t histogram_bin_count;
    bool lock_initialized;
};

static uint16_t secondary_load_u16_le(const unsigned char *input)
{
    return (uint16_t) input[0] | (uint16_t) ((uint16_t) input[1] << 8U);
}

static uint32_t secondary_load_u32_le(const unsigned char *input)
{
    return (uint32_t) input[0] |
           (uint32_t) input[1] << 8U |
           (uint32_t) input[2] << 16U |
           (uint32_t) input[3] << 24U;
}

static uint64_t secondary_load_u64_le(const unsigned char *input)
{
    return (uint64_t) secondary_load_u32_le(input) | (uint64_t) secondary_load_u32_le(input + 4U) << 32U;
}

static void secondary_store_u16_le(unsigned char *output, uint16_t value)
{
    output[0] = (unsigned char) value;
    output[1] = (unsigned char) (value >> 8U);
}

static void secondary_store_u32_le(unsigned char *output, uint32_t value)
{
    output[0] = (unsigned char) value;
    output[1] = (unsigned char) (value >> 8U);
    output[2] = (unsigned char) (value >> 16U);
    output[3] = (unsigned char) (value >> 24U);
}

static void secondary_store_u64_le(unsigned char *output, uint64_t value)
{
    secondary_store_u32_le(output, (uint32_t) value);
    secondary_store_u32_le(output + 4U, (uint32_t) (value >> 32U));
}

static uint64_t secondary_checksum(const unsigned char *data, size_t size)
{
    uint64_t checksum = UINT64_C(1469598103934665603);

    for (size_t index = 0U; index < size; ++index) {
        unsigned char byte = index >= 40U && index < 48U ? 0U : data[index];

        checksum ^= byte;
        checksum *= UINT64_C(1099511628211);
    }

    return checksum;
}

static bool secondary_type_valid(GeoDatabaseIndexType type)
{
    return type >= GEO_DATABASE_INDEX_BOOL && type <= GEO_DATABASE_INDEX_BYTES;
}

static bool secondary_pointer_valid(const char *json_pointer)
{
    if (!json_pointer || (json_pointer[0] != '\0' && json_pointer[0] != '/')) {
        return false;
    }

    for (const char *cursor = json_pointer; *cursor; ++cursor) {
        if (*cursor == '~' && (cursor[1] != '0' && cursor[1] != '1')) {
            return false;
        }
    }

    return true;
}

static bool secondary_query_accepts(GeoDatabaseIndexOperator operation,
                                    int lower_comparison,
                                    int upper_comparison);

static void secondary_catalog_key(unsigned char output[sizeof(GEO_SECONDARY_CATALOG_PREFIX) + 8U], uint64_t index_id)
{
    memcpy(output, GEO_SECONDARY_CATALOG_PREFIX, sizeof(GEO_SECONDARY_CATALOG_PREFIX));
    geo_db_store_u64_be(output + sizeof(GEO_SECONDARY_CATALOG_PREFIX), index_id);
}

static GeoDatabaseStatus secondary_definition_put(GeoRocksBatch *batch, const GeoSecondaryEntry *entry)
{
    size_t name_size = strlen(entry->info.name);
    size_t pointer_size = strlen(entry->info.json_pointer);

    if (name_size > SIZE_MAX - GEO_SECONDARY_DEFINITION_HEADER_SIZE ||
        pointer_size > SIZE_MAX - GEO_SECONDARY_DEFINITION_HEADER_SIZE - name_size) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    size_t value_size = GEO_SECONDARY_DEFINITION_HEADER_SIZE + name_size + pointer_size;
    unsigned char *value = calloc(1U, value_size);

    if (!value) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    secondary_store_u32_le(value, GEO_SECONDARY_DEFINITION_MAGIC);
    secondary_store_u16_le(value + 4U, GEO_SECONDARY_DEFINITION_VERSION);
    secondary_store_u16_le(value + 6U, GEO_SECONDARY_DEFINITION_HEADER_SIZE);
    secondary_store_u32_le(value + 8U, (uint32_t) entry->state);
    secondary_store_u32_le(value + 12U, (uint32_t) entry->info.type);
    secondary_store_u64_le(value + 16U, entry->info.index_id);
    secondary_store_u64_le(value + 24U, entry->info.entry_count);
    secondary_store_u32_le(value + 32U, (uint32_t) name_size);
    secondary_store_u32_le(value + 36U, (uint32_t) pointer_size);
    secondary_store_u32_le(value + 48U, entry->histogram.bin_count);
    secondary_store_u32_le(value + 52U, entry->histogram.boundary_size);
    memcpy(value + GEO_SECONDARY_DEFINITION_HEADER_SIZE, entry->info.name, name_size);
    memcpy(value + GEO_SECONDARY_DEFINITION_HEADER_SIZE + name_size, entry->info.json_pointer, pointer_size);

    secondary_store_u64_le(value + 40U, secondary_checksum(value, value_size));

    unsigned char key[sizeof(GEO_SECONDARY_CATALOG_PREFIX) + 8U];
    GeoRocksStatus rocks_status;

    secondary_catalog_key(key, entry->info.index_id);
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

static bool secondary_definition_decode(const void *key,
                                        size_t key_size,
                                        const void *value,
                                        size_t value_size,
                                        GeoSecondaryEntry *entry)
{
    if (!key || key_size != sizeof(GEO_SECONDARY_CATALOG_PREFIX) + 8U || !value ||
        value_size < GEO_SECONDARY_DEFINITION_HEADER_SIZE || !entry) {
        return false;
    }

    const unsigned char *key_bytes = key;
    const unsigned char *bytes = value;
    uint32_t name_size = secondary_load_u32_le(bytes + 32U);
    uint32_t pointer_size = secondary_load_u32_le(bytes + 36U);
    uint32_t histogram_bin_count = secondary_load_u32_le(bytes + 48U);
    uint32_t histogram_boundary_size = secondary_load_u32_le(bytes + 52U);
    uint64_t index_id = secondary_load_u64_le(bytes + 16U);
    uint64_t expected_size = GEO_SECONDARY_DEFINITION_HEADER_SIZE + (uint64_t) name_size + pointer_size;
    GeoSecondaryState state = (GeoSecondaryState) secondary_load_u32_le(bytes + 8U);
    GeoDatabaseIndexType type = (GeoDatabaseIndexType) secondary_load_u32_le(bytes + 12U);

    if (memcmp(key_bytes, GEO_SECONDARY_CATALOG_PREFIX, sizeof(GEO_SECONDARY_CATALOG_PREFIX)) != 0 ||
        geo_db_load_u64_be(key_bytes + sizeof(GEO_SECONDARY_CATALOG_PREFIX)) != index_id ||
        secondary_load_u32_le(bytes) != GEO_SECONDARY_DEFINITION_MAGIC ||
        secondary_load_u16_le(bytes + 4U) != GEO_SECONDARY_DEFINITION_VERSION ||
        secondary_load_u16_le(bytes + 6U) != GEO_SECONDARY_DEFINITION_HEADER_SIZE ||
        (state != GEO_SECONDARY_BUILDING && state != GEO_SECONDARY_READY) || !secondary_type_valid(type) || !index_id ||
        name_size == 0U || name_size > GEO_DATABASE_MAX_INDEX_NAME_SIZE || pointer_size > GEO_DATABASE_MAX_INDEX_POINTER_SIZE ||
        histogram_bin_count > GEO_SECONDARY_HISTOGRAM_MAX_BINS || expected_size != value_size ||
        ((histogram_bin_count == 0U) != (histogram_boundary_size == 0U)) ||
        secondary_load_u32_le(bytes + 56U) != 0U || secondary_load_u32_le(bytes + 60U) != 0U ||
        secondary_load_u64_le(bytes + 40U) != secondary_checksum(bytes, value_size)) {
        return false;
    }

    *entry = (GeoSecondaryEntry) {
        .info = {
            .index_id = index_id,
            .entry_count = secondary_load_u64_le(bytes + 24U),
            .type = type,
        },
        .histogram = {
            .boundary_size = histogram_boundary_size,
            .bin_count = histogram_bin_count,
        },
        .state = state,
    };
    memcpy(entry->info.name, bytes + GEO_SECONDARY_DEFINITION_HEADER_SIZE, name_size);
    memcpy(entry->info.json_pointer,
           bytes + GEO_SECONDARY_DEFINITION_HEADER_SIZE + name_size,
           pointer_size);

    bool valid = secondary_pointer_valid(entry->info.json_pointer) &&
                 (state == GEO_SECONDARY_READY || (!entry->info.entry_count && !entry->histogram.bin_count));

    if (!valid) {
        geo_secondary_histogram_destroy(&entry->histogram);
    }

    return valid;
}

static bool secondary_catalog_reserve(GeoSecondaryCatalog *catalog, size_t required)
{
    if (required <= catalog->capacity) {
        return true;
    }

    size_t capacity = catalog->capacity ? catalog->capacity : 8U;

    while (capacity < required) {
        if (capacity > SIZE_MAX / 2U) {
            return false;
        }

        capacity *= 2U;
    }

    if (capacity > SIZE_MAX / sizeof(*catalog->entries)) {
        return false;
    }

    GeoSecondaryEntry *entries = realloc(catalog->entries, capacity * sizeof(*entries));

    if (!entries) {
        return false;
    }

    catalog->entries = entries;
    catalog->capacity = capacity;
    return true;
}

static size_t secondary_catalog_find_name_unlocked(const GeoSecondaryCatalog *catalog, const char *name)
{
    for (size_t index = 0U; index < catalog->count; ++index) {
        if (strcmp(catalog->entries[index].info.name, name) == 0) {
            return index;
        }
    }

    return SIZE_MAX;
}

static bool secondary_catalog_recompute_histogram_offsets(GeoSecondaryCatalog *catalog)
{
    size_t total = 0U;

    for (size_t index = 0U; index < catalog->count; ++index) {
        if (catalog->entries[index].histogram.bin_count > SIZE_MAX - total) {
            return false;
        }

        catalog->entries[index].histogram_offset = total;
        total += catalog->entries[index].histogram.bin_count;
    }

    catalog->histogram_bin_count = total;
    return true;
}

GeoSecondaryCatalog *geo_secondary_catalog_open(GeoRocksDatabase *rocks,
                                                uint64_t next_secondary_index_id,
                                                GeoDatabaseStatus *status)
{
    if (!rocks || next_secondary_index_id == 0U) {
        if (status) {
            *status = GEO_DATABASE_INVALID_ARGUMENT;
        }
        return NULL;
    }

    GeoSecondaryCatalog *catalog = calloc(1U, sizeof(*catalog));

    if (!catalog) {
        if (status) {
            *status = GEO_DATABASE_OUT_OF_MEMORY;
        }
        return NULL;
    }

    if (pthread_rwlock_init(&catalog->lock, NULL) != 0) {
        free(catalog);
        if (status) {
            *status = GEO_DATABASE_FAILED_STATE;
        }
        return NULL;
    }
    catalog->lock_initialized = true;

    GeoRocksStatus rocks_status;
    GeoRocksIterator *iterator = geo_rocks_iterator_create(rocks, GEO_ROCKS_CF_CATALOG, NULL, &rocks_status);

    if (!iterator) {
        geo_secondary_catalog_destroy(catalog);
        if (status) {
            *status = geo_db_status_from_rocks(&rocks_status);
        }
        return NULL;
    }

    geo_rocks_iterator_seek(iterator, GEO_SECONDARY_CATALOG_PREFIX, sizeof(GEO_SECONDARY_CATALOG_PREFIX));
    GeoDatabaseStatus load_status = GEO_DATABASE_OK;

    while (geo_rocks_iterator_valid(iterator)) {
        size_t key_size = 0U;
        size_t value_size = 0U;
        const void *key = geo_rocks_iterator_key(iterator, &key_size);
        const void *value = geo_rocks_iterator_value(iterator, &value_size);

        if (!key || key_size < sizeof(GEO_SECONDARY_CATALOG_PREFIX) ||
            memcmp(key, GEO_SECONDARY_CATALOG_PREFIX, sizeof(GEO_SECONDARY_CATALOG_PREFIX)) != 0) {
            break;
        }

        GeoSecondaryEntry entry = { 0 };

        if (!secondary_definition_decode(key, key_size, value, value_size, &entry)) {
            load_status = GEO_DATABASE_CORRUPTION;
            break;
        }

        bool valid_entry = entry.info.index_id < next_secondary_index_id &&
                           secondary_catalog_find_name_unlocked(catalog, entry.info.name) == SIZE_MAX &&
                           entry.histogram.bin_count <= SIZE_MAX - catalog->histogram_bin_count &&
                           secondary_catalog_reserve(catalog, catalog->count + 1U);

        if (valid_entry && entry.state == GEO_SECONDARY_READY && entry.histogram.bin_count) {
            uint32_t bin_count = entry.histogram.bin_count;
            uint32_t boundary_size = entry.histogram.boundary_size;

            load_status = geo_secondary_histogram_load_metadata(rocks,
                                                                entry.info.index_id,
                                                                entry.info.type,
                                                                bin_count,
                                                                boundary_size,
                                                                &entry.histogram);

            if (load_status == GEO_DATABASE_OK) {
                load_status = geo_secondary_histogram_load_counts(rocks,
                                                                  entry.info.index_id,
                                                                  entry.info.entry_count,
                                                                  &entry.histogram);
            }
            valid_entry = load_status == GEO_DATABASE_OK;
        }

        if (!valid_entry) {
            geo_secondary_histogram_destroy(&entry.histogram);

            if (load_status == GEO_DATABASE_OK) {
                load_status = GEO_DATABASE_CORRUPTION;
            }
            break;
        }

        entry.histogram_offset = catalog->histogram_bin_count;
        catalog->histogram_bin_count += entry.histogram.bin_count;
        catalog->entries[catalog->count++] = entry;
        geo_rocks_iterator_next(iterator);
    }

    if (load_status == GEO_DATABASE_OK && !geo_rocks_iterator_status(iterator, &rocks_status)) {
        load_status = geo_db_status_from_rocks(&rocks_status);
    }

    geo_rocks_iterator_destroy(iterator);

    if (load_status != GEO_DATABASE_OK) {
        geo_secondary_catalog_destroy(catalog);
        if (status) {
            *status = load_status;
        }
        return NULL;
    }

    if (status) {
        *status = GEO_DATABASE_OK;
    }
    return catalog;
}

void geo_secondary_catalog_destroy(GeoSecondaryCatalog *catalog)
{
    if (!catalog) {
        return;
    }

    if (catalog->lock_initialized) {
        (void) pthread_rwlock_destroy(&catalog->lock);
    }

    for (size_t index = 0U; index < catalog->count; ++index) {
        geo_secondary_histogram_destroy(&catalog->entries[index].histogram);
    }

    free(catalog->entries);
    free(catalog);
}

bool geo_secondary_catalog_has_indexes(const GeoSecondaryCatalog *catalog)
{
    if (!catalog) {
        return false;
    }

    for (size_t index = 0U; index < catalog->count; ++index) {
        if (catalog->entries[index].state == GEO_SECONDARY_READY) {
            return true;
        }
    }

    return false;
}

static void secondary_encoded_release(GeoSecondaryEncodedValue *encoded)
{
    free(encoded->heap);

    *encoded = (GeoSecondaryEncodedValue) { 0 };
}

static void secondary_encode_u64(GeoSecondaryEncodedValue *encoded, uint64_t value)
{
    geo_db_store_u64_be(encoded->fixed, value);
    encoded->data = encoded->fixed;
    encoded->size = sizeof(uint64_t);
    encoded->present = true;
}

static GeoDatabaseStatus secondary_encode_bytes(GeoSecondaryEncodedValue *encoded, const void *data, size_t size)
{
    if (size > GEO_DATABASE_MAX_INDEXED_VALUE_SIZE || size > (SIZE_MAX - 2U) / 2U) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    size_t capacity = size * 2U + 2U;
    unsigned char *output = encoded->fixed;

    if (capacity > sizeof(encoded->fixed)) {
        encoded->heap = malloc(capacity);
        output = encoded->heap;
    }

    if (!output) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    const unsigned char *input = data;
    size_t written = 0U;

    for (size_t index = 0U; index < size; ++index) {
        output[written++] = input[index];

        if (input[index] == 0U) {
            output[written++] = UINT8_MAX;
        }
    }

    output[written++] = 0U;
    output[written++] = 0U;
    encoded->data = output;
    encoded->size = written;
    encoded->present = true;
    return GEO_DATABASE_OK;
}

static GeoDatabaseStatus secondary_encode_doc_value(const GeoSecondaryEntry *entry,
                                                    GeoDocView document,
                                                    GeoSecondaryEncodedValue *encoded)
{
    *encoded = (GeoSecondaryEncodedValue) { 0 };

    if (!document.size) {
        return GEO_DATABASE_OK;
    }

    GeoDocValue value;
    GeoDocStatus doc_status = geo_doc_find_pointer(document, entry->info.json_pointer, &value);

    if (doc_status == GEO_DOC_NOT_FOUND) {
        return GEO_DATABASE_OK;
    }
    if (doc_status != GEO_DOC_OK) {
        if (doc_status == GEO_DOC_OUT_OF_MEMORY) {
            return GEO_DATABASE_OUT_OF_MEMORY;
        }

        return doc_status == GEO_DOC_INVALID_FORMAT ? GEO_DATABASE_CORRUPTION : GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoDocType actual_type = geo_doc_value_type(value);

    if (actual_type == GEO_DOC_NULL) {
        return GEO_DATABASE_OK;
    }

    switch (entry->info.type) {
        case GEO_DATABASE_INDEX_BOOL: {
            bool result;

            if (actual_type != GEO_DOC_BOOL || geo_doc_value_bool(value, &result) != GEO_DOC_OK) {
                return GEO_DATABASE_INVALID_ARGUMENT;
            }
            encoded->fixed[0] = result ? 1U : 0U;
            encoded->data = encoded->fixed;
            encoded->size = 1U;
            encoded->present = true;
            return GEO_DATABASE_OK;
        }
        case GEO_DATABASE_INDEX_INT64:
        case GEO_DATABASE_INDEX_DATETIME: {
            int64_t result;

            if (actual_type != GEO_DOC_INT64 || geo_doc_value_int64(value, &result) != GEO_DOC_OK) {
                return GEO_DATABASE_INVALID_ARGUMENT;
            }
            secondary_encode_u64(encoded, (uint64_t) result ^ (UINT64_C(1) << 63U));
            return GEO_DATABASE_OK;
        }
        case GEO_DATABASE_INDEX_UINT64: {
            uint64_t result;

            if (actual_type != GEO_DOC_UINT64 || geo_doc_value_uint64(value, &result) != GEO_DOC_OK) {
                return GEO_DATABASE_INVALID_ARGUMENT;
            }
            secondary_encode_u64(encoded, result);
            return GEO_DATABASE_OK;
        }
        case GEO_DATABASE_INDEX_DOUBLE: {
            double result;
            uint64_t bits;

            if (actual_type != GEO_DOC_DOUBLE || geo_doc_value_double(value, &result) != GEO_DOC_OK || !isfinite(result)) {
                return GEO_DATABASE_INVALID_ARGUMENT;
            }

            if (result == 0.0) {
                result = 0.0;
            }
            memcpy(&bits, &result, sizeof(bits));
            bits = bits & (UINT64_C(1) << 63U) ? ~bits : bits ^ (UINT64_C(1) << 63U);
            secondary_encode_u64(encoded, bits);
            return GEO_DATABASE_OK;
        }
        case GEO_DATABASE_INDEX_STRING:
        case GEO_DATABASE_INDEX_BYTES: {
            const void *data;
            size_t size;
            GeoDocType expected = entry->info.type == GEO_DATABASE_INDEX_STRING ? GEO_DOC_STRING : GEO_DOC_BYTES;

            if (actual_type != expected || geo_doc_value_data(value, &data, &size) != GEO_DOC_OK) {
                return GEO_DATABASE_INVALID_ARGUMENT;
            }
            return secondary_encode_bytes(encoded, data, size);
        }
    }

    return GEO_DATABASE_INVALID_ARGUMENT;
}

static GeoDatabaseStatus secondary_encode_query_value(const GeoDatabaseIndexValue *value,
                                                      GeoSecondaryEncodedValue *encoded)
{
    *encoded = (GeoSecondaryEncodedValue) { 0 };

    if (!value || !secondary_type_valid(value->type)) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    switch (value->type) {
        case GEO_DATABASE_INDEX_BOOL:
            encoded->fixed[0] = value->as.boolean ? 1U : 0U;
            encoded->data = encoded->fixed;
            encoded->size = 1U;
            encoded->present = true;
            return GEO_DATABASE_OK;
        case GEO_DATABASE_INDEX_INT64:
            secondary_encode_u64(encoded, (uint64_t) value->as.signed_integer ^ (UINT64_C(1) << 63U));
            return GEO_DATABASE_OK;
        case GEO_DATABASE_INDEX_UINT64:
            secondary_encode_u64(encoded, value->as.unsigned_integer);
            return GEO_DATABASE_OK;
        case GEO_DATABASE_INDEX_DATETIME:
            secondary_encode_u64(encoded, (uint64_t) value->as.datetime ^ (UINT64_C(1) << 63U));
            return GEO_DATABASE_OK;
        case GEO_DATABASE_INDEX_DOUBLE: {
            double number = value->as.floating_point;
            uint64_t bits;

            if (!isfinite(number)) {
                return GEO_DATABASE_INVALID_ARGUMENT;
            }
            if (number == 0.0) {
                number = 0.0;
            }
            memcpy(&bits, &number, sizeof(bits));
            bits = bits & (UINT64_C(1) << 63U) ? ~bits : bits ^ (UINT64_C(1) << 63U);
            secondary_encode_u64(encoded, bits);
            return GEO_DATABASE_OK;
        }
        case GEO_DATABASE_INDEX_STRING:
        case GEO_DATABASE_INDEX_BYTES:
            if (value->as.bytes.size && !value->as.bytes.data) {
                return GEO_DATABASE_INVALID_ARGUMENT;
            }
            return secondary_encode_bytes(encoded, value->as.bytes.data, value->as.bytes.size);
    }

    return GEO_DATABASE_INVALID_ARGUMENT;
}

static unsigned char *secondary_make_key(uint64_t index_id,
                                         const GeoSecondaryEncodedValue *value,
                                         uint64_t object_id,
                                         unsigned char stack_storage[GEO_SECONDARY_STACK_KEY_SIZE],
                                         size_t *key_size,
                                         bool *allocated)
{
    if (value->size > SIZE_MAX - 16U) {
        return NULL;
    }

    size_t size = 8U + value->size + 8U;
    unsigned char *key = size <= GEO_SECONDARY_STACK_KEY_SIZE ? stack_storage : malloc(size);

    if (!key) {
        return NULL;
    }

    geo_db_store_u64_be(key, index_id);
    memcpy(key + 8U, value->data, value->size);
    geo_db_store_u64_be(key + 8U + value->size, object_id);
    *key_size = size;
    *allocated = key != stack_storage;
    return key;
}

static void secondary_index_bounds(uint64_t index_id, unsigned char begin[8], unsigned char end[8])
{
    geo_db_store_u64_be(begin, index_id);
    geo_db_store_u64_be(end, index_id + 1U);
}

static GeoDatabaseStatus secondary_delete_index_range(GeoRocksDatabase *rocks, uint64_t index_id, bool synchronize)
{
    unsigned char begin[8];
    unsigned char end[8];
    GeoRocksStatus rocks_status;
    GeoRocksBatch *batch = geo_rocks_batch_create(rocks, 32U, &rocks_status);

    secondary_index_bounds(index_id, begin, end);
    bool succeeded = batch &&
                     geo_rocks_batch_delete_range(batch,
                                                  GEO_ROCKS_CF_SECONDARY_INDEX,
                                                  begin,
                                                  sizeof(begin),
                                                  end,
                                                  sizeof(end),
                                                  &rocks_status) &&
                     geo_rocks_write(rocks, batch, synchronize, &rocks_status);
    geo_rocks_batch_destroy(batch);
    return succeeded ? GEO_DATABASE_OK : geo_db_status_from_rocks(&rocks_status);
}

static GeoDatabaseStatus secondary_build_index(GeoRocksDatabase *rocks, GeoSecondaryEntry *entry)
{
    GeoDatabaseStatus status = secondary_delete_index_range(rocks, entry->info.index_id, true);

    if (status != GEO_DATABASE_OK) {
        return status;
    }

    GeoRocksStatus rocks_status;
    GeoRocksIterator *iterator = geo_rocks_iterator_create(rocks, GEO_ROCKS_CF_OBJECTS, NULL, &rocks_status);

    if (!iterator) {
        return geo_db_status_from_rocks(&rocks_status);
    }

    GeoRocksBatch *batch = geo_rocks_batch_create(rocks, 1024U * 1024U, &rocks_status);

    if (!batch) {
        geo_rocks_iterator_destroy(iterator);
        return geo_db_status_from_rocks(&rocks_status);
    }

    geo_rocks_iterator_seek_first(iterator);
    size_t batch_entries = 0U;
    uint64_t entry_count = 0U;

    while (status == GEO_DATABASE_OK && geo_rocks_iterator_valid(iterator)) {
        size_t key_size = 0U;
        size_t value_size = 0U;
        const void *object_key = geo_rocks_iterator_key(iterator, &key_size);
        const void *object_value = geo_rocks_iterator_value(iterator, &value_size);
        uint64_t sequence;
        uint64_t morton_code;
        GeoDocView document;

        if (!object_key || key_size != 8U ||
            !geo_db_object_decode_view(object_value, value_size, &sequence, &morton_code, &document)) {
            status = GEO_DATABASE_CORRUPTION;
            break;
        }

        GeoSecondaryEncodedValue encoded;
        status = secondary_encode_doc_value(entry, document, &encoded);

        if (status == GEO_DATABASE_OK && encoded.present) {
            unsigned char stack_key[GEO_SECONDARY_STACK_KEY_SIZE];
            size_t secondary_key_size;
            bool key_allocated = false;
            unsigned char *secondary_key = secondary_make_key(entry->info.index_id,
                                                               &encoded,
                                                               geo_db_load_u64_be(object_key),
                                                               stack_key,
                                                               &secondary_key_size,
                                                               &key_allocated);

            if (!secondary_key) {
                status = GEO_DATABASE_OUT_OF_MEMORY;
            } else if (!geo_rocks_batch_put(batch,
                                            GEO_ROCKS_CF_SECONDARY_INDEX,
                                            secondary_key,
                                            secondary_key_size,
                                            NULL,
                                            0U,
                                            &rocks_status)) {
                status = geo_db_status_from_rocks(&rocks_status);
            }

            if (key_allocated) {
                free(secondary_key);
            }

            if (status == GEO_DATABASE_OK && entry_count == UINT64_MAX) {
                status = GEO_DATABASE_CORRUPTION;
            }

            if (status == GEO_DATABASE_OK) {
                batch_entries++;
                entry_count++;
            }
        }

        secondary_encoded_release(&encoded);
        geo_rocks_iterator_next(iterator);

        if (status == GEO_DATABASE_OK && batch_entries == GEO_SECONDARY_BUILD_BATCH_ENTRIES) {
            if (!geo_rocks_write(rocks, batch, false, &rocks_status)) {
                status = geo_db_status_from_rocks(&rocks_status);
            }

            geo_rocks_batch_destroy(batch);
            batch = status == GEO_DATABASE_OK ? geo_rocks_batch_create(rocks, 1024U * 1024U, &rocks_status) : NULL;
            batch_entries = 0U;

            if (status == GEO_DATABASE_OK && !batch) {
                status = geo_db_status_from_rocks(&rocks_status);
            }
        }
    }

    if (status == GEO_DATABASE_OK && !geo_rocks_iterator_status(iterator, &rocks_status)) {
        status = geo_db_status_from_rocks(&rocks_status);
    }
    if (status == GEO_DATABASE_OK && batch_entries && !geo_rocks_write(rocks, batch, false, &rocks_status)) {
        status = geo_db_status_from_rocks(&rocks_status);
    }

    geo_rocks_batch_destroy(batch);
    geo_rocks_iterator_destroy(iterator);

    if (status == GEO_DATABASE_OK) {
        entry->info.entry_count = entry_count;
        status = geo_secondary_histogram_build(rocks,
                                               entry->info.index_id,
                                               entry->info.type,
                                               entry_count,
                                               &entry->histogram);
    }
    return status;
}

static GeoDatabaseStatus secondary_store_definition(GeoRocksDatabase *rocks,
                                                    const GeoSecondaryEntry *entry,
                                                    bool synchronize,
                                                    bool store_histogram)
{
    GeoRocksStatus rocks_status;
    GeoRocksBatch *batch = geo_rocks_batch_create(rocks, 4096U, &rocks_status);
    GeoDatabaseStatus status = batch ? secondary_definition_put(batch, entry) : geo_db_status_from_rocks(&rocks_status);

    if (status == GEO_DATABASE_OK && store_histogram && entry->histogram.bin_count) {
        status = geo_secondary_histogram_put_metadata(batch,
                                                      entry->info.index_id,
                                                      entry->info.type,
                                                      &entry->histogram);
    }

    if (status == GEO_DATABASE_OK && store_histogram && entry->histogram.bin_count) {
        status = geo_secondary_histogram_put_all_counts(batch, entry->info.index_id, &entry->histogram);
    }

    if (status == GEO_DATABASE_OK && !geo_rocks_write(rocks, batch, synchronize, &rocks_status)) {
        status = geo_db_status_from_rocks(&rocks_status);
    }

    geo_rocks_batch_destroy(batch);
    return status;
}

GeoDatabaseStatus geo_secondary_catalog_recover(GeoRocksDatabase *rocks, GeoSecondaryCatalog *catalog)
{
    if (!rocks || !catalog) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    for (size_t index = 0U; index < catalog->count;) {
        GeoSecondaryEntry *entry = catalog->entries + index;

        if (entry->state != GEO_SECONDARY_BUILDING) {
            index++;
            continue;
        }

        GeoDatabaseStatus status = secondary_build_index(rocks, entry);

        /* A crash can interrupt rollback of a type-incompatible build. Such a definition
         * is not a committed index and must not prevent access to valid canonical objects. */
        if (status == GEO_DATABASE_INVALID_ARGUMENT) {
            status = geo_secondary_catalog_drop(rocks, catalog, entry->info.name);

            if (status != GEO_DATABASE_OK) {
                return status;
            }

            continue;
        }

        if (status != GEO_DATABASE_OK) {
            return status;
        }

        entry->state = GEO_SECONDARY_READY;
        status = secondary_store_definition(rocks, entry, true, true);

        if (status != GEO_DATABASE_OK) {
            entry->state = GEO_SECONDARY_BUILDING;
            return status;
        }

        index++;
    }

    return secondary_catalog_recompute_histogram_offsets(catalog) ? GEO_DATABASE_OK : GEO_DATABASE_CORRUPTION;
}

GeoDatabaseStatus geo_secondary_catalog_create(GeoRocksDatabase *rocks,
                                               GeoSecondaryCatalog *catalog,
                                               GeoDbCatalog *database_catalog,
                                               const char *name,
                                               const char *json_pointer,
                                               GeoDatabaseIndexType type)
{
    if (!rocks || !catalog || !database_catalog || !name || !name[0] || !secondary_type_valid(type) ||
        !secondary_pointer_valid(json_pointer)) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    size_t name_size = strlen(name);
    size_t pointer_size = strlen(json_pointer);

    if (name_size > GEO_DATABASE_MAX_INDEX_NAME_SIZE || pointer_size > GEO_DATABASE_MAX_INDEX_POINTER_SIZE ||
        database_catalog->next_secondary_index_id == 0U || database_catalog->next_secondary_index_id == UINT64_MAX) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    if (pthread_rwlock_wrlock(&catalog->lock) != 0) {
        return GEO_DATABASE_FAILED_STATE;
    }

    bool already_exists = secondary_catalog_find_name_unlocked(catalog, name) != SIZE_MAX;
    bool reserved = already_exists || secondary_catalog_reserve(catalog, catalog->count + 1U);
    (void) pthread_rwlock_unlock(&catalog->lock);

    if (already_exists) {
        return GEO_DATABASE_ALREADY_EXISTS;
    }
    if (!reserved) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    GeoSecondaryEntry entry = {
        .info = {
            .index_id = database_catalog->next_secondary_index_id,
            .type = type,
        },
        .state = GEO_SECONDARY_BUILDING,
    };
    memcpy(entry.info.name, name, name_size + 1U);
    memcpy(entry.info.json_pointer, json_pointer, pointer_size + 1U);

    GeoDbCatalog updated_catalog = *database_catalog;
    updated_catalog.next_secondary_index_id++;
    GeoRocksStatus rocks_status;
    GeoRocksBatch *batch = geo_rocks_batch_create(rocks, 2048U, &rocks_status);
    GeoDatabaseStatus status = batch ? secondary_definition_put(batch, &entry) : geo_db_status_from_rocks(&rocks_status);

    if (status == GEO_DATABASE_OK && !geo_db_catalog_put(batch, &updated_catalog, &status)) {
        status = GEO_DATABASE_IO_ERROR;
    }
    if (status == GEO_DATABASE_OK && !geo_rocks_write(rocks, batch, true, &rocks_status)) {
        status = geo_db_status_from_rocks(&rocks_status);
    }
    geo_rocks_batch_destroy(batch);

    if (status != GEO_DATABASE_OK) {
        geo_secondary_histogram_destroy(&entry.histogram);
        return status;
    }

    *database_catalog = updated_catalog;

    if (pthread_rwlock_wrlock(&catalog->lock) != 0) {
        return GEO_DATABASE_FAILED_STATE;
    }

    size_t catalog_index = catalog->count;
    catalog->entries[catalog->count++] = entry;
    (void) pthread_rwlock_unlock(&catalog->lock);
    status = secondary_build_index(rocks, &entry);

    if (status == GEO_DATABASE_OK) {
        if (pthread_rwlock_wrlock(&catalog->lock) != 0) {
            status = GEO_DATABASE_FAILED_STATE;
        } else {
            /* Reserve publication state before READY is durable. No fallible allocation or
             * lock acquisition may leave a durable READY index unmaintained in memory. */
            if (entry.histogram.bin_count > SIZE_MAX - catalog->histogram_bin_count) {
                status = GEO_DATABASE_OUT_OF_MEMORY;
            } else {
                entry.state = GEO_SECONDARY_READY;
                status = secondary_store_definition(rocks, &entry, true, true);
            }

            if (status == GEO_DATABASE_OK) {
                entry.histogram_offset = catalog->histogram_bin_count;
                catalog->histogram_bin_count += entry.histogram.bin_count;
                catalog->entries[catalog_index] = entry;
            }

            (void) pthread_rwlock_unlock(&catalog->lock);
        }
    }

    if (status == GEO_DATABASE_OK) {
        return status;
    }

    geo_secondary_histogram_destroy(&entry.histogram);
    GeoDatabaseStatus rollback_status = geo_secondary_catalog_drop(rocks, catalog, name);

    /* If rollback itself cannot be made durable, callers must quarantine writes: READY
     * publication may have reached storage even when its acknowledgement failed. */
    return rollback_status == GEO_DATABASE_OK ? status : GEO_DATABASE_FAILED_STATE;
}

GeoDatabaseStatus geo_secondary_catalog_drop(GeoRocksDatabase *rocks,
                                             GeoSecondaryCatalog *catalog,
                                             const char *name)
{
    if (!rocks || !catalog || !name || !name[0]) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    size_t index = secondary_catalog_find_name_unlocked(catalog, name);

    if (index == SIZE_MAX) {
        return GEO_DATABASE_NOT_FOUND;
    }

    unsigned char catalog_key[sizeof(GEO_SECONDARY_CATALOG_PREFIX) + 8U];
    unsigned char begin[8];
    unsigned char end[8];
    uint64_t index_id = catalog->entries[index].info.index_id;
    GeoRocksStatus rocks_status;
    GeoRocksBatch *batch = geo_rocks_batch_create(rocks, 128U, &rocks_status);
    GeoDatabaseStatus status = batch ? GEO_DATABASE_OK : geo_db_status_from_rocks(&rocks_status);

    secondary_catalog_key(catalog_key, index_id);
    secondary_index_bounds(index_id, begin, end);

    if (status == GEO_DATABASE_OK &&
        !geo_rocks_batch_delete(batch,
                                GEO_ROCKS_CF_CATALOG,
                                catalog_key,
                                sizeof(catalog_key),
                                &rocks_status)) {
        status = geo_db_status_from_rocks(&rocks_status);
    }
    if (status == GEO_DATABASE_OK &&
        !geo_rocks_batch_delete_range(batch,
                                      GEO_ROCKS_CF_SECONDARY_INDEX,
                                      begin,
                                      sizeof(begin),
                                      end,
                                      sizeof(end),
                                      &rocks_status)) {
        status = geo_db_status_from_rocks(&rocks_status);
    }
    if (status == GEO_DATABASE_OK) {
        status = geo_secondary_histogram_delete(batch, index_id);
    }
    if (status == GEO_DATABASE_OK && !geo_rocks_write(rocks, batch, true, &rocks_status)) {
        status = geo_db_status_from_rocks(&rocks_status);
    }

    geo_rocks_batch_destroy(batch);

    if (status != GEO_DATABASE_OK) {
        return status;
    }

    if (pthread_rwlock_wrlock(&catalog->lock) != 0) {
        return GEO_DATABASE_FAILED_STATE;
    }

    geo_secondary_histogram_destroy(&catalog->entries[index].histogram);
    memmove(catalog->entries + index,
            catalog->entries + index + 1U,
            (catalog->count - index - 1U) * sizeof(*catalog->entries));
    catalog->count--;
    bool offsets_valid = secondary_catalog_recompute_histogram_offsets(catalog);
    (void) pthread_rwlock_unlock(&catalog->lock);
    return offsets_valid ? GEO_DATABASE_OK : GEO_DATABASE_FAILED_STATE;
}

GeoDatabaseStatus geo_secondary_catalog_list(const GeoSecondaryCatalog *catalog,
                                             GeoDatabaseIndexInfo *indexes,
                                             size_t capacity,
                                             size_t *count)
{
    if (!catalog || !count || (capacity && !indexes)) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoSecondaryCatalog *mutable_catalog = (GeoSecondaryCatalog *) catalog;

    if (pthread_rwlock_rdlock(&mutable_catalog->lock) != 0) {
        return GEO_DATABASE_FAILED_STATE;
    }

    size_t ready_count = 0U;

    for (size_t index = 0U; index < catalog->count; ++index) {
        if (catalog->entries[index].state != GEO_SECONDARY_READY) {
            continue;
        }

        if (ready_count < capacity) {
            indexes[ready_count] = catalog->entries[index].info;
        }
        ready_count++;
    }

    *count = ready_count;
    (void) pthread_rwlock_unlock(&mutable_catalog->lock);
    return ready_count <= capacity ? GEO_DATABASE_OK : GEO_DATABASE_OUT_OF_MEMORY;
}

GeoDatabaseStatus geo_secondary_write_begin(const GeoSecondaryCatalog *catalog, GeoSecondaryWriteState *state)
{
    if (!catalog || !state) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    *state = (GeoSecondaryWriteState) { 0 };
    state->index_count = catalog->count;
    state->histogram_bin_count = catalog->histogram_bin_count;

    if (!state->index_count && !state->histogram_bin_count) {
        return GEO_DATABASE_OK;
    }

    if (state->histogram_bin_count > SIZE_MAX - state->index_count) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    size_t delta_count = state->index_count + state->histogram_bin_count;

    if (delta_count > SIZE_MAX / sizeof(*state->entry_deltas)) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    state->entry_deltas = calloc(delta_count, sizeof(*state->entry_deltas));
    state->histogram_deltas = state->entry_deltas ? state->entry_deltas + state->index_count : NULL;
    return state->entry_deltas ? GEO_DATABASE_OK : GEO_DATABASE_OUT_OF_MEMORY;
}

static bool secondary_delta_add(int64_t *delta, int increment)
{
    if ((increment > 0 && *delta == INT64_MAX) || (increment < 0 && *delta == INT64_MIN)) {
        return false;
    }

    *delta += increment;
    return true;
}

GeoDatabaseStatus geo_secondary_write_object(const GeoSecondaryCatalog *catalog,
                                             GeoRocksBatch *batch,
                                             GeoSecondaryWriteState *state,
                                             uint64_t object_id,
                                             GeoDocView old_document,
                                             GeoDocView new_document)
{
    if (!catalog || !batch || !state || state->index_count != catalog->count ||
        state->histogram_bin_count != catalog->histogram_bin_count || !object_id ||
        (old_document.size && !old_document.data) || (new_document.size && !new_document.data)) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    for (size_t index = 0U; index < catalog->count; ++index) {
        const GeoSecondaryEntry *entry = catalog->entries + index;

        if (entry->state != GEO_SECONDARY_READY) {
            continue;
        }

        GeoSecondaryEncodedValue old_value;
        GeoSecondaryEncodedValue new_value;
        GeoDatabaseStatus status = secondary_encode_doc_value(entry, old_document, &old_value);

        if (status == GEO_DATABASE_OK) {
            status = secondary_encode_doc_value(entry, new_document, &new_value);
        } else {
            new_value = (GeoSecondaryEncodedValue) { 0 };
        }

        bool unchanged = status == GEO_DATABASE_OK && old_value.present == new_value.present &&
                         (!old_value.present ||
                          (old_value.size == new_value.size && memcmp(old_value.data, new_value.data, old_value.size) == 0));

        if (status == GEO_DATABASE_OK && !unchanged && old_value.present) {
            unsigned char stack_key[GEO_SECONDARY_STACK_KEY_SIZE];
            size_t key_size;
            bool key_allocated = false;
            unsigned char *key = secondary_make_key(entry->info.index_id,
                                                     &old_value,
                                                     object_id,
                                                     stack_key,
                                                     &key_size,
                                                     &key_allocated);
            GeoRocksStatus rocks_status;

            if (!key) {
                status = GEO_DATABASE_OUT_OF_MEMORY;
            } else if (!geo_rocks_batch_delete(batch, GEO_ROCKS_CF_SECONDARY_INDEX, key, key_size, &rocks_status)) {
                status = geo_db_status_from_rocks(&rocks_status);
            }
            if (key_allocated) {
                free(key);
            }
        }

        if (status == GEO_DATABASE_OK && !unchanged && new_value.present) {
            unsigned char stack_key[GEO_SECONDARY_STACK_KEY_SIZE];
            size_t key_size;
            bool key_allocated = false;
            unsigned char *key = secondary_make_key(entry->info.index_id,
                                                     &new_value,
                                                     object_id,
                                                     stack_key,
                                                     &key_size,
                                                     &key_allocated);
            GeoRocksStatus rocks_status;

            if (!key) {
                status = GEO_DATABASE_OUT_OF_MEMORY;
            } else if (!geo_rocks_batch_put(batch, GEO_ROCKS_CF_SECONDARY_INDEX, key, key_size, NULL, 0U, &rocks_status)) {
                status = geo_db_status_from_rocks(&rocks_status);
            }
            if (key_allocated) {
                free(key);
            }
        }

        if (status == GEO_DATABASE_OK && !unchanged) {
            int entry_delta = (int) new_value.present - (int) old_value.present;

            if (entry_delta && !secondary_delta_add(state->entry_deltas + index, entry_delta)) {
                status = GEO_DATABASE_CORRUPTION;
            }
        }

        if (status == GEO_DATABASE_OK && !unchanged && entry->histogram.bin_count) {
            if (old_value.present) {
                uint32_t old_bin = geo_secondary_histogram_find_bin(&entry->histogram,
                                                                    old_value.data,
                                                                    old_value.size);

                if (!secondary_delta_add(state->histogram_deltas + entry->histogram_offset + old_bin, -1)) {
                    status = GEO_DATABASE_CORRUPTION;
                }
            }

            if (status == GEO_DATABASE_OK && new_value.present) {
                uint32_t new_bin = geo_secondary_histogram_find_bin(&entry->histogram,
                                                                    new_value.data,
                                                                    new_value.size);

                if (!secondary_delta_add(state->histogram_deltas + entry->histogram_offset + new_bin, 1)) {
                    status = GEO_DATABASE_CORRUPTION;
                }
            }
        }

        secondary_encoded_release(&new_value);
        secondary_encoded_release(&old_value);

        if (status != GEO_DATABASE_OK) {
            return status;
        }
    }

    return GEO_DATABASE_OK;
}

static bool secondary_apply_entry_delta(uint64_t count, int64_t delta, uint64_t *result)
{
    if (delta < 0) {
        uint64_t magnitude = (uint64_t) (-(delta + 1)) + 1U;

        if (magnitude > count) {
            return false;
        }
        *result = count - magnitude;
        return true;
    }

    if ((uint64_t) delta > UINT64_MAX - count) {
        return false;
    }

    *result = count + (uint64_t) delta;
    return true;
}

GeoDatabaseStatus geo_secondary_write_statistics(const GeoSecondaryCatalog *catalog,
                                                 GeoRocksBatch *batch,
                                                 const GeoSecondaryWriteState *state)
{
    if (!catalog || !batch || !state || state->index_count != catalog->count ||
        state->histogram_bin_count != catalog->histogram_bin_count) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    for (size_t index = 0U; index < catalog->count; ++index) {
        if (!state->entry_deltas[index]) {
            continue;
        }

        GeoSecondaryEntry updated = catalog->entries[index];

        if (!secondary_apply_entry_delta(updated.info.entry_count,
                                         state->entry_deltas[index],
                                         &updated.info.entry_count)) {
            return GEO_DATABASE_CORRUPTION;
        }

        GeoDatabaseStatus status = secondary_definition_put(batch, &updated);

        if (status != GEO_DATABASE_OK) {
            return status;
        }
    }

    for (size_t index = 0U; index < catalog->count; ++index) {
        const GeoSecondaryEntry *entry = catalog->entries + index;

        for (uint32_t bin = 0U; bin < entry->histogram.bin_count; ++bin) {
            int64_t delta = state->histogram_deltas[entry->histogram_offset + bin];

            if (!delta) {
                continue;
            }

            uint64_t updated_count;

            if (!secondary_apply_entry_delta(entry->histogram.counts[bin], delta, &updated_count)) {
                return GEO_DATABASE_CORRUPTION;
            }

            GeoDatabaseStatus status = geo_secondary_histogram_put_count(batch,
                                                                         entry->info.index_id,
                                                                         bin,
                                                                         updated_count);

            if (status != GEO_DATABASE_OK) {
                return status;
            }
        }
    }

    return GEO_DATABASE_OK;
}

bool geo_secondary_write_commit(GeoSecondaryCatalog *catalog, const GeoSecondaryWriteState *state)
{
    if (!catalog || !state || state->index_count != catalog->count ||
        state->histogram_bin_count != catalog->histogram_bin_count) {
        return false;
    }

    if (pthread_rwlock_wrlock(&catalog->lock) != 0) {
        return false;
    }

    for (size_t index = 0U; index < catalog->count; ++index) {
        uint64_t updated_count;

        if (!secondary_apply_entry_delta(catalog->entries[index].info.entry_count,
                                         state->entry_deltas[index],
                                         &updated_count)) {
            (void) pthread_rwlock_unlock(&catalog->lock);
            return false;
        }

        const GeoSecondaryEntry *entry = catalog->entries + index;

        for (uint32_t bin = 0U; bin < entry->histogram.bin_count; ++bin) {
            if (!secondary_apply_entry_delta(entry->histogram.counts[bin],
                                             state->histogram_deltas[entry->histogram_offset + bin],
                                             &updated_count)) {
                (void) pthread_rwlock_unlock(&catalog->lock);
                return false;
            }
        }
    }

    for (size_t index = 0U; index < catalog->count; ++index) {
        bool applied = secondary_apply_entry_delta(catalog->entries[index].info.entry_count,
                                                   state->entry_deltas[index],
                                                   &catalog->entries[index].info.entry_count);

        if (!applied) {
            (void) pthread_rwlock_unlock(&catalog->lock);
            return false;
        }

        GeoSecondaryEntry *entry = catalog->entries + index;

        for (uint32_t bin = 0U; bin < entry->histogram.bin_count; ++bin) {
            applied = secondary_apply_entry_delta(entry->histogram.counts[bin],
                                                  state->histogram_deltas[entry->histogram_offset + bin],
                                                  &entry->histogram.counts[bin]);

            if (!applied) {
                (void) pthread_rwlock_unlock(&catalog->lock);
                return false;
            }
        }
    }

    return pthread_rwlock_unlock(&catalog->lock) == 0;
}

void geo_secondary_write_destroy(GeoSecondaryWriteState *state)
{
    if (!state) {
        return;
    }

    free(state->entry_deltas);
    *state = (GeoSecondaryWriteState) { 0 };
}

static bool secondary_query_accepts(GeoDatabaseIndexOperator operation, int lower_comparison, int upper_comparison)
{
    switch (operation) {
        case GEO_DATABASE_INDEX_EQUAL:
            return lower_comparison == 0;
        case GEO_DATABASE_INDEX_LESS:
            return lower_comparison < 0;
        case GEO_DATABASE_INDEX_LESS_EQUAL:
            return lower_comparison <= 0;
        case GEO_DATABASE_INDEX_GREATER:
            return lower_comparison > 0;
        case GEO_DATABASE_INDEX_GREATER_EQUAL:
            return lower_comparison >= 0;
        case GEO_DATABASE_INDEX_BETWEEN:
            return lower_comparison >= 0 && upper_comparison <= 0;
    }

    return false;
}

static bool secondary_query_past_end(GeoDatabaseIndexOperator operation, int lower_comparison, int upper_comparison)
{
    if (operation == GEO_DATABASE_INDEX_EQUAL) {
        return lower_comparison > 0;
    }
    if (operation == GEO_DATABASE_INDEX_LESS) {
        return lower_comparison >= 0;
    }
    if (operation == GEO_DATABASE_INDEX_LESS_EQUAL) {
        return lower_comparison > 0;
    }
    return operation == GEO_DATABASE_INDEX_BETWEEN && upper_comparison > 0;
}

static bool secondary_result_append(GeoIdResult *result, uint64_t object_id)
{
    if (result->count == result->capacity) {
        size_t capacity = result->capacity ? result->capacity : 64U;

        if (capacity > SIZE_MAX / 2U) {
            return false;
        }

        capacity *= 2U;

        if (!geo_id_result_reserve(result, capacity)) {
            return false;
        }
    }

    result->ids[result->count++] = object_id;
    return true;
}

static GeoDatabaseStatus secondary_select_entry_locked(const GeoSecondaryCatalog *catalog,
                                                       const GeoDatabaseIndexPredicate *predicate,
                                                       const GeoSecondaryEntry **selected_entry)
{
    if (!catalog || !predicate || !predicate->index_name || !predicate->index_name[0] || !selected_entry ||
        predicate->operation < GEO_DATABASE_INDEX_EQUAL || predicate->operation > GEO_DATABASE_INDEX_BETWEEN) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoSecondaryCatalog *mutable_catalog = (GeoSecondaryCatalog *) catalog;

    if (pthread_rwlock_rdlock(&mutable_catalog->lock) != 0) {
        return GEO_DATABASE_FAILED_STATE;
    }

    size_t index = secondary_catalog_find_name_unlocked(catalog, predicate->index_name);

    if (index == SIZE_MAX || catalog->entries[index].state != GEO_SECONDARY_READY) {
        (void) pthread_rwlock_unlock(&mutable_catalog->lock);
        return GEO_DATABASE_NOT_FOUND;
    }

    *selected_entry = catalog->entries + index;

    if (predicate->lower.type != (*selected_entry)->info.type ||
        (predicate->operation == GEO_DATABASE_INDEX_BETWEEN && predicate->upper.type != (*selected_entry)->info.type)) {
        (void) pthread_rwlock_unlock(&mutable_catalog->lock);
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    return GEO_DATABASE_OK;
}

static bool secondary_select_entry_unlock(const GeoSecondaryCatalog *catalog)
{
    return pthread_rwlock_unlock(&((GeoSecondaryCatalog *) catalog)->lock) == 0;
}

GeoDatabaseStatus geo_secondary_prepare_predicates(const GeoSecondaryCatalog *catalog,
                                                   const GeoDatabaseIndexPredicate *predicates,
                                                   size_t predicate_count,
                                                   GeoSecondaryPreparedPredicate *prepared)
{
    if (!catalog || !predicates || !predicate_count || !prepared) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoSecondaryCatalog *mutable_catalog = (GeoSecondaryCatalog *) catalog;

    if (pthread_rwlock_rdlock(&mutable_catalog->lock) != 0) {
        return GEO_DATABASE_FAILED_STATE;
    }

    GeoDatabaseStatus status = GEO_DATABASE_OK;

    for (size_t predicate_index = 0U; predicate_index < predicate_count; ++predicate_index) {
        const GeoDatabaseIndexPredicate *predicate = predicates + predicate_index;

        if (!predicate->index_name || !predicate->index_name[0] ||
            predicate->operation < GEO_DATABASE_INDEX_EQUAL || predicate->operation > GEO_DATABASE_INDEX_BETWEEN) {
            status = GEO_DATABASE_INVALID_ARGUMENT;
            break;
        }

        size_t entry_index = secondary_catalog_find_name_unlocked(catalog, predicate->index_name);

        if (entry_index == SIZE_MAX || catalog->entries[entry_index].state != GEO_SECONDARY_READY) {
            status = GEO_DATABASE_NOT_FOUND;
            break;
        }

        const GeoSecondaryEntry *entry = catalog->entries + entry_index;

        if (predicate->lower.type != entry->info.type ||
            (predicate->operation == GEO_DATABASE_INDEX_BETWEEN && predicate->upper.type != entry->info.type)) {
            status = GEO_DATABASE_INVALID_ARGUMENT;
            break;
        }

        GeoSecondaryEncodedValue lower;
        GeoSecondaryEncodedValue upper = { 0 };

        status = secondary_encode_query_value(&predicate->lower, &lower);

        if (status == GEO_DATABASE_OK && predicate->operation == GEO_DATABASE_INDEX_BETWEEN) {
            status = secondary_encode_query_value(&predicate->upper, &upper);

            if (status == GEO_DATABASE_OK &&
                geo_secondary_compare_encoded(lower.data, lower.size, upper.data, upper.size) > 0) {
                status = GEO_DATABASE_INVALID_ARGUMENT;
            }
        }

        secondary_encoded_release(&upper);
        secondary_encoded_release(&lower);

        if (status != GEO_DATABASE_OK) {
            break;
        }

        prepared[predicate_index] = (GeoSecondaryPreparedPredicate) {
            .predicate = predicate,
            .type = entry->info.type,
        };
        memcpy(prepared[predicate_index].json_pointer,
               entry->info.json_pointer,
               strlen(entry->info.json_pointer) + 1U);
    }

    if (pthread_rwlock_unlock(&mutable_catalog->lock) != 0) {
        return GEO_DATABASE_FAILED_STATE;
    }

    return status;
}

static int secondary_compare_bytes(const void *first, size_t first_size, const void *second, size_t second_size)
{
    size_t common_size = first_size < second_size ? first_size : second_size;
    int comparison = common_size ? memcmp(first, second, common_size) : 0;

    return comparison ? comparison : (first_size > second_size) - (first_size < second_size);
}

static GeoDatabaseStatus secondary_compare_doc_value(GeoDocValue document_value,
                                                     GeoDatabaseIndexType type,
                                                     const GeoDatabaseIndexValue *query_value,
                                                     int *comparison)
{
    GeoDocType actual_type = geo_doc_value_type(document_value);

    switch (type) {
        case GEO_DATABASE_INDEX_BOOL: {
            bool value;

            if (actual_type != GEO_DOC_BOOL || geo_doc_value_bool(document_value, &value) != GEO_DOC_OK) {
                return GEO_DATABASE_CORRUPTION;
            }

            *comparison = (value > query_value->as.boolean) - (value < query_value->as.boolean);
            return GEO_DATABASE_OK;
        }
        case GEO_DATABASE_INDEX_INT64:
        case GEO_DATABASE_INDEX_DATETIME: {
            int64_t value;
            int64_t query = type == GEO_DATABASE_INDEX_DATETIME ? query_value->as.datetime
                                                                : query_value->as.signed_integer;

            if (actual_type != GEO_DOC_INT64 || geo_doc_value_int64(document_value, &value) != GEO_DOC_OK) {
                return GEO_DATABASE_CORRUPTION;
            }

            *comparison = (value > query) - (value < query);
            return GEO_DATABASE_OK;
        }
        case GEO_DATABASE_INDEX_UINT64: {
            uint64_t value;

            if (actual_type != GEO_DOC_UINT64 || geo_doc_value_uint64(document_value, &value) != GEO_DOC_OK) {
                return GEO_DATABASE_CORRUPTION;
            }

            *comparison = (value > query_value->as.unsigned_integer) - (value < query_value->as.unsigned_integer);
            return GEO_DATABASE_OK;
        }
        case GEO_DATABASE_INDEX_DOUBLE: {
            double value;

            if (actual_type != GEO_DOC_DOUBLE || geo_doc_value_double(document_value, &value) != GEO_DOC_OK || !isfinite(value)) {
                return GEO_DATABASE_CORRUPTION;
            }

            double query = query_value->as.floating_point;

            *comparison = (value > query) - (value < query);
            return GEO_DATABASE_OK;
        }
        case GEO_DATABASE_INDEX_STRING:
        case GEO_DATABASE_INDEX_BYTES: {
            const void *value;
            size_t value_size;
            GeoDocType expected_type = type == GEO_DATABASE_INDEX_STRING ? GEO_DOC_STRING : GEO_DOC_BYTES;

            if (actual_type != expected_type || geo_doc_value_data(document_value, &value, &value_size) != GEO_DOC_OK) {
                return GEO_DATABASE_CORRUPTION;
            }

            *comparison = secondary_compare_bytes(value,
                                                  value_size,
                                                  query_value->as.bytes.data,
                                                  query_value->as.bytes.size);
            return GEO_DATABASE_OK;
        }
    }

    return GEO_DATABASE_INVALID_ARGUMENT;
}

GeoDatabaseStatus geo_secondary_document_matches(const GeoSecondaryPreparedPredicate *predicates,
                                                 size_t predicate_count,
                                                 GeoDocView document,
                                                 bool *matches)
{
    if (!predicates || !predicate_count || !matches || (document.size && !document.data)) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    *matches = false;

    if (!document.size) {
        return GEO_DATABASE_OK;
    }

    for (size_t index = 0U; index < predicate_count; ++index) {
        const GeoSecondaryPreparedPredicate *prepared = predicates + index;
        const GeoDatabaseIndexPredicate *predicate = prepared->predicate;
        GeoDocValue document_value;
        GeoDocStatus doc_status = geo_doc_find_pointer(document, prepared->json_pointer, &document_value);

        if (doc_status == GEO_DOC_NOT_FOUND ||
            (doc_status == GEO_DOC_OK && geo_doc_value_type(document_value) == GEO_DOC_NULL)) {
            return GEO_DATABASE_OK;
        }
        if (doc_status != GEO_DOC_OK) {
            return GEO_DATABASE_CORRUPTION;
        }

        int lower_comparison;
        GeoDatabaseStatus status = secondary_compare_doc_value(document_value,
                                                               prepared->type,
                                                               &predicate->lower,
                                                               &lower_comparison);
        int upper_comparison = 0;

        if (status == GEO_DATABASE_OK && predicate->operation == GEO_DATABASE_INDEX_BETWEEN) {
            status = secondary_compare_doc_value(document_value,
                                                 prepared->type,
                                                 &predicate->upper,
                                                 &upper_comparison);
        }

        if (status != GEO_DATABASE_OK) {
            return status;
        }
        if (!secondary_query_accepts(predicate->operation, lower_comparison, upper_comparison)) {
            return GEO_DATABASE_OK;
        }
    }

    *matches = true;
    return GEO_DATABASE_OK;
}

GeoDatabaseStatus geo_secondary_query_filtered(const GeoSecondaryCatalog *catalog,
                                               GeoRocksDatabase *rocks,
                                               const GeoDatabaseIndexPredicate *predicate,
                                               GeoIdResult *result,
                                               GeoDatabaseIndexQueryStats *stats,
                                               bool (*filter)(uint64_t object_id, void *context),
                                               void *filter_context)
{
    if (!rocks || !result) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    const GeoSecondaryEntry *entry;
    GeoDatabaseStatus status = secondary_select_entry_locked(catalog, predicate, &entry);

    if (status != GEO_DATABASE_OK) {
        return status;
    }

    uint64_t index_id = entry->info.index_id;
    GeoDatabaseIndexType index_type = entry->info.type;

    if (!secondary_select_entry_unlock(catalog)) {
        return GEO_DATABASE_FAILED_STATE;
    }

    GeoSecondaryEncodedValue lower;
    GeoSecondaryEncodedValue upper = { 0 };
    status = secondary_encode_query_value(&predicate->lower, &lower);

    if (status == GEO_DATABASE_OK && predicate->operation == GEO_DATABASE_INDEX_BETWEEN) {
        status = secondary_encode_query_value(&predicate->upper, &upper);

        if (status == GEO_DATABASE_OK &&
            geo_secondary_compare_encoded(lower.data, lower.size, upper.data, upper.size) > 0) {
            status = GEO_DATABASE_INVALID_ARGUMENT;
        }
    }

    GeoRocksStatus rocks_status;
    GeoRocksIterator *iterator = status == GEO_DATABASE_OK
                                     ? geo_rocks_iterator_create(rocks, GEO_ROCKS_CF_SECONDARY_INDEX, NULL, &rocks_status)
                                     : NULL;

    if (status == GEO_DATABASE_OK && !iterator) {
        status = geo_db_status_from_rocks(&rocks_status);
    }

    unsigned char seek_storage[GEO_SECONDARY_STACK_KEY_SIZE];
    unsigned char *seek_key = seek_storage;
    unsigned char *seek_heap = NULL;

    if (status == GEO_DATABASE_OK) {
        bool starts_at_index = predicate->operation == GEO_DATABASE_INDEX_LESS ||
                               predicate->operation == GEO_DATABASE_INDEX_LESS_EQUAL;
        size_t seek_size = starts_at_index ? 8U : 8U + lower.size;

        if (seek_size > sizeof(seek_storage)) {
            seek_heap = malloc(seek_size);
            seek_key = seek_heap;
        }

        if (!seek_key) {
            status = GEO_DATABASE_OUT_OF_MEMORY;
        } else {
            geo_db_store_u64_be(seek_key, index_id);

            if (!starts_at_index) {
                memcpy(seek_key + 8U, lower.data, lower.size);
            }
            geo_rocks_iterator_seek(iterator, seek_key, seek_size);
        }
    }

    geo_id_result_clear(result);
    GeoDatabaseIndexQueryStats query_stats = { 0 };

    while (status == GEO_DATABASE_OK && geo_rocks_iterator_valid(iterator)) {
        size_t key_size = 0U;
        const unsigned char *key = geo_rocks_iterator_key(iterator, &key_size);

        if (!key || key_size < 16U || geo_db_load_u64_be(key) != index_id) {
            break;
        }

        const unsigned char *encoded_value = key + 8U;
        size_t encoded_size = key_size - 16U;

        if (!geo_secondary_encoded_key_valid(index_type, encoded_value, encoded_size)) {
            status = GEO_DATABASE_CORRUPTION;
            break;
        }

        int lower_comparison = geo_secondary_compare_encoded(encoded_value, encoded_size, lower.data, lower.size);
        int upper_comparison = predicate->operation == GEO_DATABASE_INDEX_BETWEEN
                                   ? geo_secondary_compare_encoded(encoded_value, encoded_size, upper.data, upper.size)
                                   : 0;

        if (secondary_query_past_end(predicate->operation, lower_comparison, upper_comparison)) {
            break;
        }

        query_stats.scanned_entries++;

        if (secondary_query_accepts(predicate->operation, lower_comparison, upper_comparison)) {
            uint64_t object_id = geo_db_load_u64_be(key + key_size - 8U);

            if (!object_id) {
                status = GEO_DATABASE_CORRUPTION;
                break;
            }

            bool filter_matches = !filter || filter(object_id, filter_context);

            if (filter_matches && !secondary_result_append(result, object_id)) {
                status = GEO_DATABASE_OUT_OF_MEMORY;
                break;
            }

            query_stats.matched_entries += filter_matches;
        }

        geo_rocks_iterator_next(iterator);
    }

    if (status == GEO_DATABASE_OK && iterator && !geo_rocks_iterator_status(iterator, &rocks_status)) {
        status = geo_db_status_from_rocks(&rocks_status);
    }

    if (stats) {
        *stats = query_stats;
    }

    free(seek_heap);
    geo_rocks_iterator_destroy(iterator);
    secondary_encoded_release(&upper);
    secondary_encoded_release(&lower);
    return status;
}

GeoDatabaseStatus geo_secondary_query(const GeoSecondaryCatalog *catalog,
                                      GeoRocksDatabase *rocks,
                                      const GeoDatabaseIndexPredicate *predicate,
                                      GeoIdResult *result,
                                      GeoDatabaseIndexQueryStats *stats)
{
    return geo_secondary_query_filtered(catalog, rocks, predicate, result, stats, NULL, NULL);
}

GeoDatabaseStatus geo_secondary_estimate(const GeoSecondaryCatalog *catalog,
                                         const GeoDatabaseIndexPredicate *predicate,
                                         uint64_t *estimated_matches)
{
    if (!estimated_matches) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    const GeoSecondaryEntry *entry;
    GeoDatabaseStatus status = secondary_select_entry_locked(catalog, predicate, &entry);

    if (status != GEO_DATABASE_OK) {
        return status;
    }

    GeoSecondaryEncodedValue lower;
    GeoSecondaryEncodedValue upper = { 0 };

    status = secondary_encode_query_value(&predicate->lower, &lower);

    if (status == GEO_DATABASE_OK && predicate->operation == GEO_DATABASE_INDEX_BETWEEN) {
        status = secondary_encode_query_value(&predicate->upper, &upper);

        if (status == GEO_DATABASE_OK &&
            geo_secondary_compare_encoded(lower.data, lower.size, upper.data, upper.size) > 0) {
            status = GEO_DATABASE_INVALID_ARGUMENT;
        }
    }

    if (status == GEO_DATABASE_OK) {
        *estimated_matches = entry->histogram.bin_count
                                 ? geo_secondary_histogram_estimate(&entry->histogram,
                                                                    predicate->operation,
                                                                    lower.data,
                                                                    lower.size,
                                                                    upper.data,
                                                                    upper.size)
                                 : entry->info.entry_count;
    }

    secondary_encoded_release(&upper);
    secondary_encoded_release(&lower);

    if (!secondary_select_entry_unlock(catalog)) {
        return GEO_DATABASE_FAILED_STATE;
    }

    return status;
}
