#include "geo_index_private.h"

#include <stdlib.h>
#include <string.h>

#define GEO_DENSITY_THRESHOLD_MULTIPLIER 4U
#define GEO_DENSITY_DEFERRED_ROOT_MAX_RECORDS 131072U

typedef struct {
    uint64_t range_min;
    size_t begin;
    size_t end;
} GeoDensityParent;

typedef struct {
    GeoPersistedDensityCell *cells;
    size_t count;
    size_t capacity;
} GeoDensityWriterLevel;

struct GeoDensityWriter {
    GeoDensityWriterLevel levels[GEO_PERSISTED_DENSITY_MAX_LEVELS];
    uint64_t *node_counts;
    uint64_t *node_begins;
    uint64_t *deferred_morton;
    uint32_t *touched_nodes;
    size_t level_offsets[GEO_PERSISTED_DENSITY_MAX_LEVELS];
    size_t level_node_counts[GEO_PERSISTED_DENSITY_MAX_LEVELS];
    size_t touched_counts[GEO_PERSISTED_DENSITY_MAX_LEVELS];
    uint64_t expected_records;
    uint64_t first_position;
    uint64_t next_position;
    uint64_t final_position;
    uint64_t baseline;
    uint64_t current_root;
    uint64_t current_root_count;
    uint64_t current_root_begin;
    size_t deferred_count;
    size_t deferred_capacity;
    size_t total_nodes;
    uint8_t prefix_bits;
    uint8_t level_count;
    const size_t *root_offsets;
    bool root_active;
    bool root_hot;
    bool defer_cold_roots;
    bool partition_fragment;
    bool failed;
};

struct GeoDensityIndex {
    GeoPersistedDensityHeader header;
    GeoPersistedDensityCell *cells;
    size_t cell_capacity;
    bool cells_mapped;
};

static bool density_size_add(size_t first, size_t second, size_t *result)
{
    if (second > SIZE_MAX - first) {
        return false;
    }

    *result = first + second;

    return true;
}

static bool density_size_multiply(size_t first, size_t second, size_t *result)
{
    if (first && second > SIZE_MAX / first) {
        return false;
    }

    *result = first * second;

    return true;
}

void geo_density_index_destroy(GeoIndex *index)
{
    if (!index || !index->density) {
        return;
    }

    if (!index->density->cells_mapped) {
        free(index->density->cells);
    }

    free(index->density);
    index->density = NULL;
    index->maximum_density_refinement = 0;
}

static bool density_cells_reserve(GeoDensityIndex *density, size_t required)
{
    if (required <= density->cell_capacity) {
        return true;
    }

    size_t capacity = density->cell_capacity ? density->cell_capacity : 64;

    while (capacity < required) {
        if (capacity > SIZE_MAX / 2U) {
            capacity = required;
            break;
        }

        capacity *= 2U;
    }

    size_t bytes;

    if (!density_size_multiply(capacity, sizeof(*density->cells), &bytes)) {
        return false;
    }

    GeoPersistedDensityCell *cells = realloc(density->cells, bytes);

    if (!cells) {
        return false;
    }

    density->cells = cells;
    density->cell_capacity = capacity;

    return true;
}

static bool density_parent_reserve(GeoDensityParent **parents, size_t *capacity, size_t required)
{
    if (required <= *capacity) {
        return true;
    }

    size_t new_capacity = *capacity ? *capacity : 16;

    while (new_capacity < required) {
        if (new_capacity > SIZE_MAX / 2U) {
            new_capacity = required;
            break;
        }

        new_capacity *= 2U;
    }

    size_t bytes;

    if (!density_size_multiply(new_capacity, sizeof(**parents), &bytes)) {
        return false;
    }

    GeoDensityParent *resized = realloc(*parents, bytes);

    if (!resized) {
        return false;
    }

    *parents = resized;
    *capacity = new_capacity;

    return true;
}

static bool density_is_hot(size_t records, size_t baseline)
{
    return baseline <= SIZE_MAX / GEO_DENSITY_THRESHOLD_MULTIPLIER &&
           records > baseline * GEO_DENSITY_THRESHOLD_MULTIPLIER;
}

static bool density_writer_is_hot(uint64_t records, uint64_t baseline)
{
    return baseline <= UINT64_MAX / GEO_DENSITY_THRESHOLD_MULTIPLIER &&
           records > baseline * GEO_DENSITY_THRESHOLD_MULTIPLIER;
}

static bool density_writer_level_reserve(GeoDensityWriterLevel *level, size_t required)
{
    if (required <= level->capacity) {
        return true;
    }

    size_t capacity = level->capacity ? level->capacity * 2U : 64U;

    if (capacity < level->capacity || capacity < required) {
        capacity = required;
    }

    size_t bytes;

    if (!density_size_multiply(capacity, sizeof(*level->cells), &bytes)) {
        return false;
    }

    GeoPersistedDensityCell *cells = realloc(level->cells, bytes);

    if (!cells) {
        return false;
    }

    level->cells = cells;
    level->capacity = capacity;

    return true;
}

static void density_writer_add_hot_record(GeoDensityWriter *writer, uint64_t morton, uint64_t absolute_position)
{
    for (size_t level_index = 0; level_index < writer->level_count; ++level_index) {
        unsigned local_bits = (unsigned) (level_index + 1U) * 2U;
        unsigned shift = 64U - (unsigned) writer->prefix_bits - local_bits;
        size_t local_index = (size_t) ((morton >> shift) & ((UINT64_C(1) << local_bits) - 1U));
        size_t offset = writer->level_offsets[level_index];
        size_t node_index = offset + local_index;

        if (!writer->node_counts[node_index]) {
            size_t touched_position = offset + writer->touched_counts[level_index]++;

            writer->touched_nodes[touched_position] = (uint32_t) local_index;
            writer->node_begins[node_index] = absolute_position;
        }

        writer->node_counts[node_index]++;
    }
}

static bool density_writer_flush_root(GeoDensityWriter *writer)
{
    if (!writer->root_active) {
        return true;
    }

    bool root_hot = writer->partition_fragment || writer->defer_cold_roots
                        ? writer->root_hot
                        : density_writer_is_hot(writer->current_root_count, writer->baseline);

    for (size_t level_index = 0; root_hot && level_index < writer->level_count; ++level_index) {
        size_t offset = writer->level_offsets[level_index];
        GeoDensityWriterLevel *output = writer->levels + level_index;

        for (size_t touched = 0; touched < writer->touched_counts[level_index]; ++touched) {
            size_t local_index = writer->touched_nodes[offset + touched];
            size_t node_index = offset + local_index;
            uint64_t node_count = writer->node_counts[node_index];
            bool parent_hot = level_index == 0 || writer->partition_fragment;

            if (level_index && !writer->partition_fragment) {
                size_t parent_index = writer->level_offsets[level_index - 1U] + (local_index >> 2U);

                parent_hot = density_writer_is_hot(writer->node_counts[parent_index], writer->baseline);
            }

            if (!parent_hot) {
                continue;
            }

            if (!density_writer_level_reserve(output, output->count + 1U)) {
                return false;
            }

            unsigned cell_bits = (unsigned) writer->prefix_bits + (unsigned) (level_index + 1U) * 2U;
            unsigned shift = 64U - cell_bits;
            uint64_t prefix = (writer->current_root << ((level_index + 1U) * 2U)) | local_index;

            output->cells[output->count++] = (GeoPersistedDensityCell) {
                .range_min = shift ? prefix << shift : prefix,
                .begin = writer->node_begins[node_index],
                .end = writer->node_begins[node_index] + node_count,
            };
        }
    }

    for (size_t level_index = 0; level_index < writer->level_count; ++level_index) {
        size_t offset = writer->level_offsets[level_index];

        for (size_t touched = 0; touched < writer->touched_counts[level_index]; ++touched) {
            size_t node_index = offset + writer->touched_nodes[offset + touched];

            writer->node_counts[node_index] = 0;
        }

        writer->touched_counts[level_index] = 0;
    }

    writer->current_root_count = 0;
    writer->deferred_count = 0;
    writer->root_hot = false;

    return true;
}

GeoDensityWriter *geo_density_writer_create(uint64_t record_count, uint8_t prefix_bits)
{
    if (prefix_bits > GEO_PERSISTED_PREFIX_MAX_BITS) {
        return NULL;
    }

    GeoDensityWriter *writer = calloc(1, sizeof(*writer));

    if (!writer) {
        return NULL;
    }

    writer->expected_records = record_count;
    writer->final_position = record_count;
    writer->prefix_bits = prefix_bits;

    if (prefix_bits) {
        uint64_t bucket_count = UINT64_C(1) << prefix_bits;

        writer->baseline = record_count / bucket_count + (record_count % bucket_count != 0);

        if (writer->baseline <= (UINT64_MAX - 1U) / GEO_DENSITY_THRESHOLD_MULTIPLIER) {
            uint64_t hot_record_count = writer->baseline * GEO_DENSITY_THRESHOLD_MULTIPLIER + 1U;

            if (hot_record_count <= GEO_DENSITY_DEFERRED_ROOT_MAX_RECORDS &&
                hot_record_count <= SIZE_MAX / sizeof(*writer->deferred_morton)) {
                writer->deferred_capacity = (size_t) hot_record_count;
                writer->deferred_morton = malloc(writer->deferred_capacity * sizeof(*writer->deferred_morton));
                writer->defer_cold_roots = writer->deferred_morton != NULL;
            }
        }
    }

    size_t level_nodes = 4U;

    while (prefix_bits &&
           writer->level_count < GEO_PERSISTED_DENSITY_MAX_LEVELS &&
           (unsigned) prefix_bits + (unsigned) (writer->level_count + 1U) * 2U <= 64U) {
        size_t level = writer->level_count++;

        writer->level_offsets[level] = writer->total_nodes;
        writer->level_node_counts[level] = level_nodes;

        if (level_nodes > SIZE_MAX - writer->total_nodes) {
            geo_density_writer_destroy(writer);

            return NULL;
        }

        writer->total_nodes += level_nodes;

        if (level_nodes > SIZE_MAX / 4U) {
            break;
        }

        level_nodes *= 4U;
    }

    if (writer->total_nodes) {
        writer->node_counts = calloc(writer->total_nodes, sizeof(*writer->node_counts));
        writer->node_begins = malloc(writer->total_nodes * sizeof(*writer->node_begins));
        writer->touched_nodes = malloc(writer->total_nodes * sizeof(*writer->touched_nodes));
    }

    if (writer->total_nodes && (!writer->node_counts || !writer->node_begins || !writer->touched_nodes)) {
        geo_density_writer_destroy(writer);

        return NULL;
    }

    return writer;
}

GeoDensityWriter *geo_density_writer_create_partition(uint64_t record_count,
                                                       uint8_t prefix_bits,
                                                       uint64_t first_position,
                                                       uint64_t partition_records,
                                                       const size_t *root_offsets)
{
    if ((prefix_bits && !root_offsets) || first_position > record_count || partition_records > record_count - first_position) {
        return NULL;
    }

    GeoDensityWriter *writer = geo_density_writer_create(record_count, prefix_bits);

    if (!writer) {
        return NULL;
    }

    free(writer->deferred_morton);
    writer->deferred_morton = NULL;
    writer->deferred_capacity = 0;
    writer->defer_cold_roots = false;
    writer->next_position = first_position;
    writer->first_position = first_position;
    writer->final_position = first_position + partition_records;
    writer->root_offsets = root_offsets;
    writer->partition_fragment = true;

    return writer;
}

void geo_density_writer_destroy(GeoDensityWriter *writer)
{
    if (!writer) {
        return;
    }

    for (size_t level = 0; level < GEO_PERSISTED_DENSITY_MAX_LEVELS; ++level) {
        free(writer->levels[level].cells);
    }

    free(writer->touched_nodes);
    free(writer->deferred_morton);
    free(writer->node_begins);
    free(writer->node_counts);
    free(writer);
}

bool geo_density_writer_add(GeoDensityWriter *writer,
                            const GeoRecord *records,
                            size_t count,
                            uint64_t first_position)
{
    if (!writer || (!records && count) || writer->failed || first_position != writer->next_position ||
        writer->next_position > writer->final_position || count > writer->final_position - writer->next_position) {
        return false;
    }

    for (size_t position = 0; position < count; ++position) {
        uint64_t absolute_position = first_position + position;
        uint64_t root = writer->prefix_bits ? records[position].z >> (64U - writer->prefix_bits) : 0;

        if (!writer->root_active || root != writer->current_root) {
            if (!density_writer_flush_root(writer)) {
                writer->failed = true;

                return false;
            }

            writer->root_active = true;
            writer->current_root = root;
            writer->current_root_begin = absolute_position;

            if (writer->partition_fragment) {
                uint64_t root_records = writer->prefix_bits
                                            ? writer->root_offsets[root + 1U] - writer->root_offsets[root]
                                            : writer->expected_records;

                writer->root_hot = density_writer_is_hot(root_records, writer->baseline);
            }
        }

        writer->current_root_count++;

        if (writer->partition_fragment) {
            if (writer->root_hot) {
                density_writer_add_hot_record(writer, records[position].z, absolute_position);
            }

            continue;
        }

        if (writer->defer_cold_roots && !writer->root_hot) {
            writer->deferred_morton[writer->deferred_count++] = records[position].z;

            if (writer->deferred_count == writer->deferred_capacity) {
                writer->root_hot = true;

                for (size_t deferred = 0; deferred < writer->deferred_count; ++deferred) {
                    density_writer_add_hot_record(writer,
                                                  writer->deferred_morton[deferred],
                                                  writer->current_root_begin + deferred);
                }

                writer->deferred_count = 0;
            }

            continue;
        }

        density_writer_add_hot_record(writer, records[position].z, absolute_position);
    }

    writer->next_position += count;

    return true;
}

static bool density_writer_append_partition_cells(GeoDensityWriter *writer,
                                                  const GeoDensityWriter *partition,
                                                  size_t level_index)
{
    GeoDensityWriterLevel *output = writer->levels + level_index;
    const GeoDensityWriterLevel *input = partition->levels + level_index;

    for (size_t cell_index = 0; cell_index < input->count; ++cell_index) {
        const GeoPersistedDensityCell *cell = input->cells + cell_index;

        if (output->count && output->cells[output->count - 1U].range_min == cell->range_min) {
            GeoPersistedDensityCell *previous = output->cells + output->count - 1U;

            if (previous->end != cell->begin) {
                return false;
            }

            previous->end = cell->end;
            continue;
        }

        if (!density_writer_level_reserve(output, output->count + 1U)) {
            return false;
        }

        output->cells[output->count++] = *cell;
    }

    return true;
}

static void density_writer_prune_level(GeoDensityWriter *writer, size_t level_index)
{
    if (!level_index) {
        return;
    }

    GeoDensityWriterLevel *level = writer->levels + level_index;
    const GeoDensityWriterLevel *parents = writer->levels + level_index - 1U;
    unsigned parent_bits = (unsigned) writer->prefix_bits + (unsigned) level_index * 2U;
    unsigned parent_shift = 64U - parent_bits;
    size_t parent_index = 0;
    size_t output_count = 0;

    for (size_t cell_index = 0; cell_index < level->count; ++cell_index) {
        GeoPersistedDensityCell cell = level->cells[cell_index];
        uint64_t parent_range = parent_shift ? (cell.range_min >> parent_shift) << parent_shift : cell.range_min;

        while (parent_index < parents->count && parents->cells[parent_index].range_min < parent_range) {
            parent_index++;
        }

        if (parent_index == parents->count || parents->cells[parent_index].range_min != parent_range) {
            continue;
        }

        uint64_t parent_records = parents->cells[parent_index].end - parents->cells[parent_index].begin;

        if (density_writer_is_hot(parent_records, writer->baseline)) {
            level->cells[output_count++] = cell;
        }
    }

    level->count = output_count;
}

bool geo_density_writer_merge_partitions(GeoDensityWriter *writer,
                                         GeoDensityWriter *const *partitions,
                                         size_t partition_count)
{
    if (!writer || (!partitions && partition_count) || writer->partition_fragment || writer->next_position || writer->root_active) {
        return false;
    }

    uint64_t next_position = 0;

    for (size_t partition_index = 0; partition_index < partition_count; ++partition_index) {
        GeoDensityWriter *partition = partitions[partition_index];

        if (!partition || !partition->partition_fragment || partition->expected_records != writer->expected_records ||
            partition->prefix_bits != writer->prefix_bits || partition->next_position != partition->final_position ||
            partition->first_position != next_position || partition->final_position < partition->first_position) {
            return false;
        }

        if (!density_writer_flush_root(partition)) {
            return false;
        }

        for (size_t level = 0; level < writer->level_count; ++level) {
            if (!density_writer_append_partition_cells(writer, partition, level)) {
                return false;
            }
        }

        next_position = partition->final_position;
    }

    if (next_position != writer->expected_records) {
        return false;
    }

    for (size_t level = 1; level < writer->level_count; ++level) {
        density_writer_prune_level(writer, level);
    }

    writer->next_position = writer->expected_records;

    return true;
}

bool geo_density_writer_append_file(GeoDensityWriter *writer,
                                    FILE *file,
                                    uint64_t *serialized_bytes,
                                    uint64_t *checksum)
{
    if (!writer || !file || !serialized_bytes || !checksum || writer->failed ||
        writer->next_position != writer->final_position || !density_writer_flush_root(writer)) {
        return false;
    }

    GeoPersistedDensityHeader header = {
        .version = GEO_PERSISTED_DENSITY_VERSION,
        .cell_size = sizeof(GeoPersistedDensityCell),
    };
    uint64_t cell_count = 0;

    for (size_t level = 0; level < writer->level_count && writer->levels[level].count; ++level) {
        header.level_first[level] = cell_count;
        header.level_prefix_bits[level] = writer->prefix_bits + (uint8_t) ((level + 1U) * 2U);
        cell_count += writer->levels[level].count;
        header.level_count++;
    }

    header.level_first[header.level_count] = cell_count;
    header.cell_count = cell_count;
    uint64_t value_checksum = geo_persisted_checksum_update(geo_persisted_checksum_initial(),
                                                             &header,
                                                             sizeof(header));
    bool succeeded = fwrite(&header, sizeof(header), 1, file) == 1;

    for (size_t level = 0; succeeded && level < header.level_count; ++level) {
        GeoDensityWriterLevel *source = writer->levels + level;
        size_t bytes = source->count * sizeof(*source->cells);

        succeeded = fwrite(source->cells, 1, bytes, file) == bytes;
        value_checksum = geo_persisted_checksum_update(value_checksum, source->cells, bytes);
    }

    if (succeeded) {
        *serialized_bytes = sizeof(header) + cell_count * sizeof(GeoPersistedDensityCell);
        *checksum = value_checksum;
    }

    return succeeded;
}

static bool density_collect_root_parents(const GeoIndex *index,
                                         size_t baseline,
                                         GeoDensityParent **parents,
                                         size_t *parent_count,
                                         size_t *parent_capacity)
{
    size_t bucket_count = (size_t) 1 << index->prefix_bits;
    unsigned shift = 64U - index->prefix_bits;

    for (size_t bucket = 0; bucket < bucket_count; ++bucket) {
        size_t begin = index->prefix_offsets[bucket];
        size_t end = index->prefix_offsets[bucket + 1];

        if (!density_is_hot(end - begin, baseline)) {
            continue;
        }

        if (!density_parent_reserve(parents, parent_capacity, *parent_count + 1)) {
            return false;
        }

        (*parents)[(*parent_count)++] = (GeoDensityParent) {
            .range_min = (uint64_t) bucket << shift,
            .begin = begin,
            .end = end,
        };
    }

    return true;
}

static bool density_append_child(GeoDensityIndex *density,
                                 GeoDensityParent **next_parents,
                                 size_t *next_count,
                                 size_t *next_capacity,
                                 uint64_t range_min,
                                 size_t begin,
                                 size_t end,
                                 size_t baseline)
{
    size_t cell_index = (size_t) density->header.cell_count;

    if (!density_cells_reserve(density, cell_index + 1)) {
        return false;
    }

    density->cells[cell_index] = (GeoPersistedDensityCell) {
        .range_min = range_min,
        .begin = begin,
        .end = end,
    };
    density->header.cell_count++;

    if (!density_is_hot(end - begin, baseline)) {
        return true;
    }

    if (!density_parent_reserve(next_parents, next_capacity, *next_count + 1)) {
        return false;
    }

    (*next_parents)[(*next_count)++] = (GeoDensityParent) {
        .range_min = range_min,
        .begin = begin,
        .end = end,
    };

    return true;
}

static bool density_subdivide_level(const GeoIndex *index,
                                    GeoDensityIndex *density,
                                    const GeoDensityParent *parents,
                                    size_t parent_count,
                                    uint8_t prefix_bits,
                                    size_t baseline,
                                    GeoDensityParent **next_parents,
                                    size_t *next_count,
                                    size_t *next_capacity)
{
    unsigned shift = 64U - prefix_bits;

    for (size_t parent_index = 0; parent_index < parent_count; ++parent_index) {
        const GeoDensityParent *parent = parents + parent_index;
        size_t position = parent->begin;

        while (position < parent->end) {
            uint64_t range_min = shift ? (index->records[position].z >> shift) << shift : index->records[position].z;
            size_t child_end = position + 1;

            while (child_end < parent->end) {
                uint64_t next_range_min = shift
                                              ? (index->records[child_end].z >> shift) << shift
                                              : index->records[child_end].z;

                if (next_range_min != range_min) {
                    break;
                }

                child_end++;
            }

            if (!density_append_child(density,
                                      next_parents,
                                      next_count,
                                      next_capacity,
                                      range_min,
                                      position,
                                      child_end,
                                      baseline)) {
                return false;
            }

            position = child_end;
        }
    }

    return true;
}

bool geo_density_index_build(GeoIndex *index)
{
    if (!index) {
        return false;
    }

    geo_density_index_destroy(index);

    GeoDensityIndex *density = calloc(1, sizeof(*density));

    if (!density) {
        return false;
    }

    density->header.version = GEO_PERSISTED_DENSITY_VERSION;
    density->header.cell_size = sizeof(GeoPersistedDensityCell);

    if (!index->prefix_bits || !index->prefix_offsets || !index->count) {
        index->density = density;

        return true;
    }

    size_t bucket_count = (size_t) 1 << index->prefix_bits;
    size_t baseline = index->count / bucket_count + (index->count % bucket_count != 0);
    GeoDensityParent *parents = NULL;
    GeoDensityParent *next_parents = NULL;
    size_t parent_count = 0;
    size_t parent_capacity = 0;
    size_t next_capacity = 0;
    bool succeeded = density_collect_root_parents(index,
                                                  baseline,
                                                  &parents,
                                                  &parent_count,
                                                  &parent_capacity);

    for (unsigned level = 0;
         succeeded && parent_count && level < GEO_PERSISTED_DENSITY_MAX_LEVELS;
         ++level) {
        unsigned prefix_bits = (unsigned) index->prefix_bits + (level + 1U) * 2U;

        if (prefix_bits > 64U) {
            prefix_bits = 64U;
        }

        density->header.level_first[level] = density->header.cell_count;
        density->header.level_prefix_bits[level] = (uint8_t) prefix_bits;
        size_t next_count = 0;

        succeeded = density_subdivide_level(index,
                                            density,
                                            parents,
                                            parent_count,
                                            (uint8_t) prefix_bits,
                                            baseline,
                                            &next_parents,
                                            &next_count,
                                            &next_capacity);

        if (!succeeded) {
            break;
        }

        density->header.level_count++;
        GeoDensityParent *temporary = parents;

        parents = next_parents;
        next_parents = temporary;
        parent_count = next_count;
        parent_capacity = next_capacity;
        next_capacity = 0;
        free(next_parents);
        next_parents = NULL;

        if (prefix_bits == 64U) {
            break;
        }
    }

    density->header.level_first[density->header.level_count] = density->header.cell_count;
    free(next_parents);
    free(parents);

    if (!succeeded) {
        if (!density->cells_mapped) {
            free(density->cells);
        }

        free(density);

        return false;
    }

    index->density = density;
    index->maximum_density_refinement = (uint8_t) density->header.level_count;

    return true;
}

size_t geo_density_index_serialized_size(const GeoIndex *index)
{
    if (!index || !index->density || index->density->header.cell_count > SIZE_MAX) {
        return 0;
    }

    size_t cell_bytes;
    size_t total;

    if (!density_size_multiply((size_t) index->density->header.cell_count,
                               sizeof(*index->density->cells),
                               &cell_bytes) ||
        !density_size_add(sizeof(index->density->header), cell_bytes, &total)) {
        return 0;
    }

    return total;
}

bool geo_density_index_serialize(const GeoIndex *index, void *destination, size_t size)
{
    size_t required = geo_density_index_serialized_size(index);

    if (!required || !destination || size < required) {
        return false;
    }

    memcpy(destination, &index->density->header, sizeof(index->density->header));

    if (index->density->header.cell_count) {
        memcpy((unsigned char *) destination + sizeof(index->density->header),
               index->density->cells,
               (size_t) index->density->header.cell_count * sizeof(*index->density->cells));
    }

    return true;
}

static bool density_header_validate(const GeoIndex *index,
                                    const GeoPersistedDensityHeader *header,
                                    size_t available,
                                    size_t *serialized_size)
{
    if (header->version != GEO_PERSISTED_DENSITY_VERSION ||
        header->cell_size != sizeof(GeoPersistedDensityCell) ||
        header->level_count > GEO_PERSISTED_DENSITY_MAX_LEVELS ||
        header->reserved != 0 ||
        header->cell_count > SIZE_MAX ||
        header->level_first[0] != 0 ||
        header->level_first[header->level_count] != header->cell_count) {
        return false;
    }

    for (size_t level = 0; level < header->level_count; ++level) {
        unsigned expected_bits = (unsigned) index->prefix_bits + (unsigned) (level + 1U) * 2U;

        if (expected_bits > 64U) {
            expected_bits = 64U;
        }

        if (header->level_first[level] > header->level_first[level + 1] ||
            header->level_prefix_bits[level] != expected_bits) {
            return false;
        }
    }

    size_t cell_bytes;

    return density_size_multiply((size_t) header->cell_count, sizeof(GeoPersistedDensityCell), &cell_bytes) &&
           density_size_add(sizeof(*header), cell_bytes, serialized_size) &&
           *serialized_size <= available;
}

bool geo_density_index_attach(GeoIndex *index, const void *data, size_t available, size_t *consumed)
{
    if (!index || !data || available < sizeof(GeoPersistedDensityHeader)) {
        return false;
    }

    const GeoPersistedDensityHeader *header = data;
    size_t serialized_size = 0;

    if (!density_header_validate(index, header, available, &serialized_size)) {
        return false;
    }

    const GeoPersistedDensityCell *cells = (const GeoPersistedDensityCell *) ((const unsigned char *) data + sizeof(*header));

    for (size_t level = 0; level < header->level_count; ++level) {
        size_t first = (size_t) header->level_first[level];
        size_t end = (size_t) header->level_first[level + 1];
        uint64_t previous_range_min = 0;

        for (size_t cell = first; cell < end; ++cell) {
            bool ordered = cell == first || previous_range_min < cells[cell].range_min;

            if (!ordered || cells[cell].begin >= cells[cell].end || cells[cell].end > index->count) {
                return false;
            }

            previous_range_min = cells[cell].range_min;
        }
    }

    GeoDensityIndex *density = calloc(1, sizeof(*density));

    if (!density) {
        return false;
    }

    density->header = *header;
    density->cells = (GeoPersistedDensityCell *) cells;
    density->cell_capacity = (size_t) header->cell_count;
    density->cells_mapped = true;
    geo_density_index_destroy(index);
    index->density = density;
    index->maximum_density_refinement = (uint8_t) header->level_count;

    if (consumed) {
        *consumed = serialized_size;
    }

    return true;
}

static const GeoPersistedDensityCell *density_find_cell(const GeoDensityIndex *density,
                                                        size_t level,
                                                        uint64_t range_min)
{
    size_t first = (size_t) density->header.level_first[level];
    size_t count = (size_t) (density->header.level_first[level + 1] - density->header.level_first[level]);

    while (count) {
        size_t step = count >> 1;
        size_t middle = first + step;

        if (density->cells[middle].range_min < range_min) {
            first = middle + 1;
            count -= step + 1;
        } else {
            count = step;
        }
    }

    size_t level_end = (size_t) density->header.level_first[level + 1];

    if (first < level_end && density->cells[first].range_min == range_min) {
        return density->cells + first;
    }

    return NULL;
}

unsigned geo_density_index_query(const GeoIndex *index, uint64_t morton, size_t *local_records)
{
    if (!index || !index->count) {
        return 0;
    }

    if (local_records) {
        if (index->prefix_bits && index->prefix_offsets) {
            size_t bucket = (size_t) (morton >> (64U - index->prefix_bits));

            *local_records = index->prefix_offsets[bucket + 1] - index->prefix_offsets[bucket];
        } else {
            *local_records = index->count;
        }
    }

    if (!index->density) {
        return 0;
    }

    unsigned refinement = 0;

    for (size_t level = 0; level < index->density->header.level_count; ++level) {
        unsigned prefix_bits = index->density->header.level_prefix_bits[level];
        unsigned shift = 64U - prefix_bits;
        uint64_t range_min = shift ? (morton >> shift) << shift : morton;
        const GeoPersistedDensityCell *cell = density_find_cell(index->density, level, range_min);

        if (!cell) {
            break;
        }

        refinement = (unsigned) level + 1U;

        if (local_records) {
            *local_records = (size_t) (cell->end - cell->begin);
        }
    }

    return refinement;
}
