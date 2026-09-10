#include "geo_spatial_memtable.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define GEO_SPATIAL_BUCKET_BITS 12U
#define GEO_SPATIAL_BUCKET_COUNT (1U << GEO_SPATIAL_BUCKET_BITS)
#define GEO_SPATIAL_CELL_BITS (GEO_SPATIAL_BUCKET_BITS / 2U)
#define GEO_SPATIAL_CELL_COUNT (1U << GEO_SPATIAL_CELL_BITS)
#define GEO_SPATIAL_CELL_SHIFT (32U - GEO_SPATIAL_CELL_BITS)
#define GEO_SPATIAL_NO_ENTRY UINT32_MAX

typedef struct {
    uint64_t object_id;
    uint64_t morton_code;
    uint32_t next_in_bucket;
    uint32_t operation;
} GeoSpatialEntry;

typedef struct {
    uint64_t object_id;
    uint32_t entry_index;
    uint32_t reserved;
} GeoSpatialHashSlot;

typedef struct {
    GeoSpatialEntry *entries;
    GeoSpatialHashSlot *hash_slots;
    uint32_t *bucket_heads;
    size_t count;
    size_t unique_count;
    size_t capacity;
    size_t hash_capacity;
    uint64_t first_sequence;
    uint64_t next_sequence;
} GeoSpatialTable;

struct GeoSpatialMemtable {
    GeoSpatialTable *active;
    GeoSpatialTable *frozen;
    size_t maximum_operations;
};

_Static_assert(sizeof(GeoSpatialEntry) == 24U, "Spatial entry must not carry redundant sequence or pointer fields");
_Static_assert(sizeof(GeoSpatialHashSlot) == 16U, "Spatial hash slot must remain 16 bytes");

static uint64_t spatial_hash_id(uint64_t object_id)
{
    object_id ^= object_id >> 30U;
    object_id *= UINT64_C(0xbf58476d1ce4e5b9);
    object_id ^= object_id >> 27U;
    object_id *= UINT64_C(0x94d049bb133111eb);
    return object_id ^ (object_id >> 31U);
}

static bool spatial_hash_capacity(size_t maximum_operations, size_t *capacity)
{
    if (!maximum_operations || maximum_operations > UINT32_MAX || maximum_operations > SIZE_MAX / 2U) {
        return false;
    }

    size_t required = maximum_operations * 2U;
    size_t selected = 2U;

    while (selected < required) {
        if (selected > SIZE_MAX / 2U) {
            return false;
        }

        selected *= 2U;
    }

    *capacity = selected;
    return true;
}

static GeoSpatialTable *spatial_table_create(size_t maximum_operations)
{
    size_t hash_capacity;

    if (!spatial_hash_capacity(maximum_operations, &hash_capacity) ||
        maximum_operations > SIZE_MAX / sizeof(GeoSpatialEntry) ||
        hash_capacity > SIZE_MAX / sizeof(GeoSpatialHashSlot)) {
        return NULL;
    }

    GeoSpatialTable *table = calloc(1U, sizeof(*table));

    if (!table) {
        return NULL;
    }

    table->entries = malloc(maximum_operations * sizeof(*table->entries));
    table->hash_slots = calloc(hash_capacity, sizeof(*table->hash_slots));
    table->bucket_heads = malloc(GEO_SPATIAL_BUCKET_COUNT * sizeof(*table->bucket_heads));

    if (!table->entries || !table->hash_slots || !table->bucket_heads) {
        free(table->bucket_heads);
        free(table->hash_slots);
        free(table->entries);
        free(table);
        return NULL;
    }

    for (size_t bucket = 0U; bucket < GEO_SPATIAL_BUCKET_COUNT; ++bucket) {
        table->bucket_heads[bucket] = GEO_SPATIAL_NO_ENTRY;
    }

    table->capacity = maximum_operations;
    table->hash_capacity = hash_capacity;
    return table;
}

static void spatial_table_destroy(GeoSpatialTable *table)
{
    if (!table) {
        return;
    }

    free(table->bucket_heads);
    free(table->hash_slots);
    free(table->entries);
    free(table);
}

static GeoSpatialHashSlot *spatial_table_find_slot(GeoSpatialTable *table, uint64_t object_id)
{
    size_t mask = table->hash_capacity - 1U;
    size_t slot = (size_t) spatial_hash_id(object_id) & mask;

    while (table->hash_slots[slot].object_id && table->hash_slots[slot].object_id != object_id) {
        slot = (slot + 1U) & mask;
    }

    return table->hash_slots + slot;
}

static const GeoSpatialHashSlot *spatial_table_lookup(const GeoSpatialTable *table, uint64_t object_id)
{
    if (!table || !table->count) {
        return NULL;
    }

    size_t mask = table->hash_capacity - 1U;
    size_t slot = (size_t) spatial_hash_id(object_id) & mask;

    while (table->hash_slots[slot].object_id) {
        if (table->hash_slots[slot].object_id == object_id) {
            return table->hash_slots + slot;
        }

        slot = (slot + 1U) & mask;
    }

    return NULL;
}

GeoSpatialMemtable *geo_spatial_memtable_create(size_t maximum_operations)
{
    GeoSpatialMemtable *memtable = calloc(1U, sizeof(*memtable));

    if (!memtable) {
        return NULL;
    }

    memtable->active = spatial_table_create(maximum_operations);

    if (!memtable->active) {
        free(memtable);
        return NULL;
    }

    memtable->maximum_operations = maximum_operations;
    return memtable;
}

void geo_spatial_memtable_destroy(GeoSpatialMemtable *memtable)
{
    if (!memtable) {
        return;
    }

    spatial_table_destroy(memtable->frozen);
    spatial_table_destroy(memtable->active);
    free(memtable);
}

bool geo_spatial_memtable_can_append(const GeoSpatialMemtable *memtable, size_t operation_count)
{
    return memtable && operation_count <= memtable->active->capacity - memtable->active->count;
}

bool geo_spatial_memtable_append(GeoSpatialMemtable *memtable,
                                 const GeoDatabaseMutation *operations,
                                 size_t operation_count,
                                 uint64_t first_sequence)
{
    if (!memtable || (!operations && operation_count) || !operation_count ||
        !geo_spatial_memtable_can_append(memtable, operation_count) ||
        operation_count > UINT64_MAX - first_sequence) {
        return false;
    }

    GeoSpatialTable *table = memtable->active;

    if (!table->count) {
        table->first_sequence = first_sequence;
    } else if (table->next_sequence != first_sequence) {
        return false;
    }

    for (size_t operation_index = 0U; operation_index < operation_count; ++operation_index) {
        if (!operations[operation_index].object_id ||
            (operations[operation_index].operation != GEO_DATABASE_UPSERT &&
             operations[operation_index].operation != GEO_DATABASE_DELETE)) {
            return false;
        }
    }

    for (size_t operation_index = 0U; operation_index < operation_count; ++operation_index) {
        const GeoDatabaseMutation *operation = operations + operation_index;
        size_t entry_index = table->count++;
        GeoSpatialEntry *entry = table->entries + entry_index;
        uint32_t bucket = (uint32_t) (operation->morton_code >> (64U - GEO_SPATIAL_BUCKET_BITS));

        *entry = (GeoSpatialEntry) {
            .object_id = operation->object_id,
            .morton_code = operation->morton_code,
            .next_in_bucket = operation->operation == GEO_DATABASE_UPSERT
                                  ? table->bucket_heads[bucket]
                                  : GEO_SPATIAL_NO_ENTRY,
            .operation = operation->operation,
        };

        if (operation->operation == GEO_DATABASE_UPSERT) {
            table->bucket_heads[bucket] = (uint32_t) entry_index;
        }

        GeoSpatialHashSlot *slot = spatial_table_find_slot(table, operation->object_id);

        if (!slot->object_id) {
            slot->object_id = operation->object_id;
            table->unique_count++;
        }

        slot->entry_index = (uint32_t) entry_index;
    }

    table->next_sequence = first_sequence + operation_count;
    return true;
}

bool geo_spatial_memtable_should_rotate(const GeoSpatialMemtable *memtable)
{
    return memtable && memtable->active->count >= memtable->maximum_operations;
}

bool geo_spatial_memtable_rotate(GeoSpatialMemtable *memtable)
{
    if (!memtable || memtable->frozen || !memtable->active->count) {
        return false;
    }

    GeoSpatialTable *replacement = spatial_table_create(memtable->maximum_operations);

    if (!replacement) {
        return false;
    }

    memtable->frozen = memtable->active;
    memtable->active = replacement;
    return true;
}

bool geo_spatial_memtable_has_active(const GeoSpatialMemtable *memtable)
{
    return memtable && memtable->active->count != 0U;
}

bool geo_spatial_memtable_has_frozen(const GeoSpatialMemtable *memtable)
{
    return memtable && memtable->frozen != NULL;
}

bool geo_spatial_memtable_export_frozen(const GeoSpatialMemtable *memtable,
                                        GeoDatabaseMutation **operations,
                                        size_t *operation_count,
                                        uint64_t *first_sequence,
                                        uint64_t *next_sequence)
{
    if (!memtable || !memtable->frozen || !operations || !operation_count || !first_sequence || !next_sequence) {
        return false;
    }

    const GeoSpatialTable *table = memtable->frozen;
    GeoDatabaseMutation *output = malloc(table->unique_count * sizeof(*output));

    if (!output) {
        return false;
    }

    size_t count = 0U;

    for (size_t slot = 0U; slot < table->hash_capacity; ++slot) {
        if (!table->hash_slots[slot].object_id) {
            continue;
        }

        const GeoSpatialEntry *entry = table->entries + table->hash_slots[slot].entry_index;
        output[count++] = (GeoDatabaseMutation) {
            .object_id = entry->object_id,
            .morton_code = entry->morton_code,
            .operation = entry->operation,
        };
    }

    *operations = output;
    *operation_count = count;
    *first_sequence = table->first_sequence;
    *next_sequence = table->next_sequence;
    return count == table->unique_count;
}

void geo_spatial_memtable_release_frozen(GeoSpatialMemtable *memtable)
{
    if (!memtable) {
        return;
    }

    spatial_table_destroy(memtable->frozen);
    memtable->frozen = NULL;
}

bool geo_spatial_memtable_allows_persisted_record(const GeoRecord *record, void *context)
{
    const GeoSpatialMemtable *memtable = context;

    return record && memtable &&
           !spatial_table_lookup(memtable->active, record->id) &&
           !spatial_table_lookup(memtable->frozen, record->id);
}

static uint32_t spatial_bucket_from_cells(uint32_t latitude_cell, uint32_t longitude_cell)
{
    uint32_t latitude = latitude_cell << GEO_SPATIAL_CELL_SHIFT;
    uint32_t longitude = longitude_cell << GEO_SPATIAL_CELL_SHIFT;
    uint64_t morton = geo_spread_bits(latitude) | geo_spread_bits(longitude) << 1U;

    return (uint32_t) (morton >> (64U - GEO_SPATIAL_BUCKET_BITS));
}

static void spatial_select_longitude_cells(uint64_t selected[GEO_SPATIAL_BUCKET_COUNT / 64U],
                                           uint32_t minimum_latitude_cell,
                                           uint32_t maximum_latitude_cell,
                                           double minimum_longitude,
                                           double maximum_longitude)
{
    uint32_t minimum_longitude_cell = geo_normalize_lng(minimum_longitude) >> GEO_SPATIAL_CELL_SHIFT;
    uint32_t maximum_longitude_cell = geo_normalize_lng(maximum_longitude) >> GEO_SPATIAL_CELL_SHIFT;

    for (uint32_t latitude_cell = minimum_latitude_cell;
         latitude_cell <= maximum_latitude_cell;
         ++latitude_cell) {
        for (uint32_t longitude_cell = minimum_longitude_cell;
             longitude_cell <= maximum_longitude_cell;
             ++longitude_cell) {
            uint32_t bucket = spatial_bucket_from_cells(latitude_cell, longitude_cell);
            selected[bucket / 64U] |= UINT64_C(1) << (bucket % 64U);
        }
    }
}

static bool spatial_entry_matches(const GeoSpatialEntry *entry,
                                  double latitude,
                                  double longitude,
                                  double radius_km,
                                  double minimum_latitude,
                                  double maximum_latitude,
                                  double minimum_longitude,
                                  double maximum_longitude)
{
    GeoPoint point = geo_decode(entry->morton_code);
    bool longitude_matches = minimum_longitude <= maximum_longitude
                               ? point.lng >= minimum_longitude && point.lng <= maximum_longitude
                               : point.lng >= minimum_longitude || point.lng <= maximum_longitude;

    return point.lat >= minimum_latitude && point.lat <= maximum_latitude && longitude_matches &&
           geo_haversine_km(latitude, longitude, point.lat, point.lng) <= radius_km;
}

static bool spatial_search_table(const GeoSpatialTable *table,
                                 const GeoSpatialTable *newer,
                                 const uint64_t selected[GEO_SPATIAL_BUCKET_COUNT / 64U],
                                 double latitude,
                                 double longitude,
                                 double radius_km,
                                 double minimum_latitude,
                                 double maximum_latitude,
                                 double minimum_longitude,
                                 double maximum_longitude,
                                 GeoSearchResult *result,
                                 size_t *matched_count,
                                 GeoRecordFilter filter,
                                 void *filter_context,
                                 GeoSearchStats *stats)
{
    if (!table) {
        return true;
    }

    for (size_t word = 0U; word < GEO_SPATIAL_BUCKET_COUNT / 64U; ++word) {
        uint64_t buckets = selected[word];

        while (buckets) {
            unsigned lane = (unsigned) __builtin_ctzll(buckets);
            uint32_t bucket = (uint32_t) (word * 64U + lane);

            for (uint32_t entry_index = table->bucket_heads[bucket];
                 entry_index != GEO_SPATIAL_NO_ENTRY;
                 entry_index = table->entries[entry_index].next_in_bucket) {
                const GeoSpatialEntry *entry = table->entries + entry_index;
                const GeoSpatialHashSlot *latest = spatial_table_lookup(table, entry->object_id);

                if (!latest || latest->entry_index != entry_index || spatial_table_lookup(newer, entry->object_id)) {
                    continue;
                }

                GeoRecord record = { .id = entry->object_id, .z = entry->morton_code };

                if (stats) {
                    stats->records_scanned++;
                }

                if ((filter && !filter(&record, filter_context)) ||
                    !spatial_entry_matches(entry,
                                           latitude,
                                           longitude,
                                           radius_km,
                                           minimum_latitude,
                                           maximum_latitude,
                                           minimum_longitude,
                                           maximum_longitude)) {
                    continue;
                }

                if ((result && !geo_result_add(result, &record)) || *matched_count == SIZE_MAX) {
                    return false;
                }

                (*matched_count)++;

                if (stats) {
                    stats->records_matched++;
                }
            }

            buckets &= buckets - 1U;
        }
    }

    return true;
}

bool geo_spatial_memtable_search_radius(const GeoSpatialMemtable *memtable,
                                        double latitude,
                                        double longitude,
                                        double radius_km,
                                        GeoSearchResult *result,
                                        size_t *count,
                                        GeoRecordFilter filter,
                                        void *filter_context,
                                        GeoSearchStats *stats)
{
    if (!memtable || (!result && !count) || (result && count)) {
        return false;
    }

    size_t initial_result_count = result ? result->count : 0U;
    size_t initial_count = count ? *count : 0U;

    double minimum_latitude;
    double maximum_latitude;
    double minimum_longitude;
    double maximum_longitude;

    geo_bounding_box(latitude,
                     longitude,
                     radius_km,
                     &minimum_latitude,
                     &maximum_latitude,
                     &minimum_longitude,
                     &maximum_longitude);

    if (!isfinite(minimum_latitude)) {
        return false;
    }

    uint64_t selected[GEO_SPATIAL_BUCKET_COUNT / 64U] = { 0 };
    uint32_t minimum_latitude_cell = geo_normalize_lat(minimum_latitude) >> GEO_SPATIAL_CELL_SHIFT;
    uint32_t maximum_latitude_cell = geo_normalize_lat(maximum_latitude) >> GEO_SPATIAL_CELL_SHIFT;

    if (minimum_longitude <= maximum_longitude) {
        spatial_select_longitude_cells(selected,
                                       minimum_latitude_cell,
                                       maximum_latitude_cell,
                                       minimum_longitude,
                                       maximum_longitude);
    } else {
        spatial_select_longitude_cells(selected,
                                       minimum_latitude_cell,
                                       maximum_latitude_cell,
                                       minimum_longitude,
                                       GEO_MAX_LNG);
        spatial_select_longitude_cells(selected,
                                       minimum_latitude_cell,
                                       maximum_latitude_cell,
                                       GEO_MIN_LNG,
                                       maximum_longitude);
    }

    size_t matched_count = 0U;
    bool succeeded = spatial_search_table(memtable->active,
                                          NULL,
                                          selected,
                                          latitude,
                                          longitude,
                                          radius_km,
                                          minimum_latitude,
                                          maximum_latitude,
                                          minimum_longitude,
                                          maximum_longitude,
                                          result,
                                          &matched_count,
                                          filter,
                                          filter_context,
                                          stats) &&
                     spatial_search_table(memtable->frozen,
                                          memtable->active,
                                          selected,
                                          latitude,
                                          longitude,
                                          radius_km,
                                          minimum_latitude,
                                          maximum_latitude,
                                          minimum_longitude,
                                          maximum_longitude,
                                          result,
                                          &matched_count,
                                          filter,
                                          filter_context,
                                          stats);

    if (succeeded && count && matched_count <= SIZE_MAX - *count) {
        *count += matched_count;
    } else if (succeeded && count) {
        succeeded = false;
    }

    if (!succeeded) {
        if (result) {
            result->count = initial_result_count;
        }

        if (count) {
            *count = initial_count;
        }
    }

    return succeeded;
}

size_t geo_spatial_memtable_active_operations(const GeoSpatialMemtable *memtable)
{
    return memtable ? memtable->active->count : 0U;
}

size_t geo_spatial_memtable_frozen_operations(const GeoSpatialMemtable *memtable)
{
    return memtable && memtable->frozen ? memtable->frozen->count : 0U;
}
