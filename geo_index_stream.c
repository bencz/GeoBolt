#include "geo_index.h"
#include "geo_index_io.h"
#include "geo_index_persistence.h"
#include "geo_index_private.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <unistd.h>
#define GEO_STREAM_SUPPORTED 1
#else
#define GEO_STREAM_SUPPORTED 0
#endif

#define GEO_STREAM_DEFAULT_CHUNK_RECORDS 1048576
#define GEO_STREAM_INPUT_BUFFER_RECORDS 4096
#define GEO_STREAM_OUTPUT_BUFFER_RECORDS 4096
#define GEO_STREAM_MERGE_FAN_IN 32

typedef struct {
    FILE *file;
    uint64_t count;
    uint64_t remaining;
    uint32_t level;
} GeoStreamRun;

typedef struct {
    GeoRecord record;
    size_t run_index;
} GeoMergeNode;

typedef struct {
    GeoStreamRun *run;
    GeoRecord *records;
    size_t position;
    size_t buffered;
} GeoMergeCursor;

typedef struct {
    FILE *file;
    GeoRecord *buffer;
    size_t buffered;
    size_t *prefix_offsets;
    size_t prefix_count;
    size_t next_prefix;
    uint8_t prefix_bits;
    GeoDensityWriter *density_writer;
    uint64_t records_written;
    uint64_t records_checksum;
    bool checksum_enabled;
} GeoMergeOutput;

struct GeoStreamBuilder {
    char *output_path;
    char *temporary_directory;
    GeoIndex *chunk;
    GeoParallelSorter *sorter;
    GeoStreamRun *runs;
    size_t run_count;
    size_t run_capacity;
    size_t runs_created;
    size_t intermediate_merges;
    size_t peak_open_runs;
    uint64_t total_records;
    double chunk_build_time_ms;
    double intermediate_merge_time_ms;
    bool failed;
    bool finished;
};

static bool stream_compact_tail(GeoStreamBuilder *builder);

static void stream_advise_sequential(FILE *file)
{
#if GEO_STREAM_SUPPORTED && defined(POSIX_FADV_SEQUENTIAL)
    (void) posix_fadvise(fileno(file), 0, 0, POSIX_FADV_SEQUENTIAL);
#else
    (void) file;
#endif
}

static void stream_preallocate(FILE *file, uint64_t record_count, size_t extra_bytes)
{
#if defined(__linux__)
    if (extra_bytes > (size_t) INT64_MAX) {
        return;
    }

    uint64_t maximum_records = ((uint64_t) INT64_MAX - extra_bytes) / sizeof(GeoRecord);

    if (record_count <= maximum_records) {
        off_t file_size = (off_t) (record_count * sizeof(GeoRecord) + extra_bytes);

        (void) posix_fallocate(fileno(file), 0, file_size);
    }
#else
    (void) file;
    (void) record_count;
    (void) extra_bytes;
#endif
}

static char *stream_copy_string(const char *source)
{
    size_t length = strlen(source);
    char *copy = malloc(length + 1);

    if (copy) {
        memcpy(copy, source, length + 1);
    }

    return copy;
}

static FILE *stream_create_temporary_run(const char *directory)
{
#if GEO_STREAM_SUPPORTED
    static const char suffix[] = "/geobolt-run-XXXXXX";
    size_t directory_length = strlen(directory);
    size_t template_size = directory_length + sizeof(suffix);
    char *path_template = malloc(template_size);

    if (!path_template) {
        return NULL;
    }

    memcpy(path_template, directory, directory_length);
    memcpy(path_template + directory_length, suffix, sizeof(suffix));

    int descriptor = mkstemp(path_template);

    if (descriptor < 0) {
        free(path_template);

        return NULL;
    }

    unlink(path_template);
    free(path_template);

    FILE *file = fdopen(descriptor, "w+b");

    if (!file) {
        close(descriptor);
    }

    return file;
#else
    (void) directory;

    return NULL;
#endif
}

static bool stream_reserve_runs(GeoStreamBuilder *builder, size_t capacity)
{
    if (capacity <= builder->run_capacity) {
        return true;
    }

    if (capacity > SIZE_MAX / sizeof(*builder->runs)) {
        return false;
    }

    GeoStreamRun *runs = realloc(builder->runs, capacity * sizeof(*builder->runs));

    if (!runs) {
        return false;
    }

    builder->runs = runs;
    builder->run_capacity = capacity;

    return true;
}

static bool stream_append_run(GeoStreamBuilder *builder, FILE *file, uint64_t count, uint32_t level)
{
    if (builder->run_count == builder->run_capacity) {
        size_t capacity = builder->run_capacity ? builder->run_capacity * 2 : 8;

        if (capacity < builder->run_capacity || !stream_reserve_runs(builder, capacity)) {
            return false;
        }
    }

    builder->runs[builder->run_count++] = (GeoStreamRun) {
        .file = file,
        .count = count,
        .remaining = 0,
        .level = level,
    };

    if (builder->run_count > builder->peak_open_runs) {
        builder->peak_open_runs = builder->run_count;
    }

    return true;
}

static bool stream_sort_chunk(GeoStreamBuilder *builder)
{
    return geo_index_sort_transient_with_sorter(builder->chunk, builder->sorter);
}

static bool stream_flush_chunk(GeoStreamBuilder *builder)
{
    if (!builder->chunk->count) {
        return true;
    }

    double start = geo_get_time_ms();

    if (!stream_sort_chunk(builder)) {
        return false;
    }

    FILE *run_file = stream_create_temporary_run(builder->temporary_directory);

    if (!run_file) {
        return false;
    }

    size_t count = builder->chunk->count;

    stream_advise_sequential(run_file);
    stream_preallocate(run_file, count, 0);

    bool succeeded = fwrite(builder->chunk->records, sizeof(GeoRecord), count, run_file) == count &&
                     fflush(run_file) == 0 &&
                     stream_append_run(builder, run_file, count, 0);

    if (!succeeded) {
        fclose(run_file);

        return false;
    }

    builder->chunk_build_time_ms += geo_get_time_ms() - start;
    builder->runs_created++;
    geo_index_clear(builder->chunk);

    return stream_compact_tail(builder);
}

GeoStreamBuilder *geo_stream_builder_create(const char *output_path,
                                            const char *temporary_directory,
                                            size_t chunk_capacity)
{
    return geo_stream_builder_create_parallel(output_path, temporary_directory, chunk_capacity, 1);
}

GeoStreamBuilder *geo_stream_builder_create_parallel(const char *output_path,
                                                     const char *temporary_directory,
                                                     size_t chunk_capacity,
                                                     size_t sort_threads)
{
#if GEO_STREAM_SUPPORTED
    if (!output_path || !output_path[0]) {
        return NULL;
    }

    if (!temporary_directory || !temporary_directory[0]) {
        temporary_directory = "/tmp";
    }

    if (!chunk_capacity) {
        chunk_capacity = GEO_STREAM_DEFAULT_CHUNK_RECORDS;
    }

    GeoStreamBuilder *builder = calloc(1, sizeof(*builder));

    if (!builder) {
        return NULL;
    }

    builder->output_path = stream_copy_string(output_path);
    builder->temporary_directory = stream_copy_string(temporary_directory);
    builder->chunk = geo_index_create(chunk_capacity);
    builder->sorter = geo_parallel_sorter_create(chunk_capacity, sort_threads);

    if (!builder->output_path || !builder->temporary_directory || !builder->chunk || !builder->sorter) {
        geo_stream_builder_destroy(builder);

        return NULL;
    }

    return builder;
#else
    (void) output_path;
    (void) temporary_directory;
    (void) chunk_capacity;
    (void) sort_threads;

    return NULL;
#endif
}

void geo_stream_builder_destroy(GeoStreamBuilder *builder)
{
    if (!builder) {
        return;
    }

    for (size_t i = 0; i < builder->run_count; ++i) {
        fclose(builder->runs[i].file);
    }

    free(builder->runs);
    geo_parallel_sorter_destroy(builder->sorter);
    geo_index_destroy(builder->chunk);
    free(builder->temporary_directory);
    free(builder->output_path);
    free(builder);
}

bool geo_stream_builder_add_batch(GeoStreamBuilder *builder,
                                  const uint64_t *ids,
                                  const double *latitudes,
                                  const double *longitudes,
                                  size_t count)
{
    if (!builder || builder->failed || builder->finished || (count && (!ids || !latitudes || !longitudes))) {
        return false;
    }

    if (count > UINT64_MAX - builder->total_records) {
        builder->failed = true;

        return false;
    }

    if (!geo_simd_validate_points(latitudes, longitudes, count)) {
        return false;
    }

    size_t offset = 0;

    while (offset < count) {
        size_t available = builder->chunk->capacity - builder->chunk->count;

        if (!available) {
            if (!stream_flush_chunk(builder)) {
                builder->failed = true;

                return false;
            }

            available = builder->chunk->capacity;
        }

        size_t batch_count = count - offset < available ? count - offset : available;

        if (!geo_index_add_batch_unchecked(builder->chunk,
                                           ids + offset,
                                           latitudes + offset,
                                           longitudes + offset,
                                           batch_count)) {
            builder->failed = true;

            return false;
        }

        offset += batch_count;
    }

    builder->total_records += count;

    return true;
}

bool geo_stream_builder_add_records(GeoStreamBuilder *builder, const GeoRecord *records, size_t count)
{
    if (!builder || builder->failed || builder->finished || (count && !records)) {
        return false;
    }

    if (count > UINT64_MAX - builder->total_records) {
        builder->failed = true;

        return false;
    }

    size_t offset = 0;

    while (offset < count) {
        size_t available = builder->chunk->capacity - builder->chunk->count;

        if (!available) {
            if (!stream_flush_chunk(builder)) {
                builder->failed = true;

                return false;
            }

            available = builder->chunk->capacity;
        }

        size_t batch_count = count - offset < available ? count - offset : available;

        if (!geo_index_add_records(builder->chunk, records + offset, batch_count)) {
            builder->failed = true;

            return false;
        }

        offset += batch_count;
    }

    builder->total_records += count;

    return true;
}

static bool merge_node_less(const GeoMergeNode *first, const GeoMergeNode *second)
{
    if (first->record.z != second->record.z) {
        return first->record.z < second->record.z;
    }

    if (first->record.id != second->record.id) {
        return first->record.id < second->record.id;
    }

    return first->run_index < second->run_index;
}

static void merge_heap_push(GeoMergeNode *heap, size_t *count, GeoMergeNode node)
{
    size_t position = (*count)++;

    while (position) {
        size_t parent = (position - 1U) >> 1U;

        if (merge_node_less(heap + parent, &node)) {
            break;
        }

        heap[position] = heap[parent];
        position = parent;
    }

    heap[position] = node;
}

static GeoMergeNode merge_heap_pop(GeoMergeNode *heap, size_t *count)
{
    GeoMergeNode minimum = heap[0];
    GeoMergeNode replacement = heap[--(*count)];
    size_t position = 0;

    while (position * 2U + 1U < *count) {
        size_t left = position * 2U + 1U;
        size_t right = left + 1U;
        size_t smallest = right < *count && merge_node_less(heap + right, heap + left) ? right : left;

        if (merge_node_less(&replacement, heap + smallest)) {
            break;
        }

        heap[position] = heap[smallest];
        position = smallest;
    }

    if (*count) {
        heap[position] = replacement;
    }

    return minimum;
}

static void merge_heap_replace_min(GeoMergeNode *heap, size_t count, GeoMergeNode replacement)
{
    size_t position = 0;

    while (position * 2U + 1U < count) {
        size_t left = position * 2U + 1U;
        size_t right = left + 1U;
        size_t smallest = right < count && merge_node_less(heap + right, heap + left) ? right : left;

        if (!merge_node_less(heap + smallest, &replacement)) {
            break;
        }

        heap[position] = heap[smallest];
        position = smallest;
    }

    heap[position] = replacement;
}

static bool stream_cursor_has_next(const GeoMergeCursor *cursor)
{
    return cursor->position < cursor->buffered || cursor->run->remaining > 0;
}

static bool stream_cursor_read(GeoMergeCursor *cursor, GeoRecord *record)
{
    if (cursor->position == cursor->buffered) {
        size_t requested = cursor->run->remaining < GEO_STREAM_INPUT_BUFFER_RECORDS
                               ? (size_t) cursor->run->remaining
                               : GEO_STREAM_INPUT_BUFFER_RECORDS;

        if (!requested || fread(cursor->records, sizeof(*cursor->records), requested, cursor->run->file) != requested) {
            return false;
        }

        cursor->run->remaining -= requested;
        cursor->position = 0;
        cursor->buffered = requested;
    }

    *record = cursor->records[cursor->position++];

    return true;
}

static bool stream_merge_output_write_block(GeoMergeOutput *output, const GeoRecord *records, size_t count)
{
    if (count > UINT64_MAX - output->records_written) {
        return false;
    }

    uint64_t first_position = output->records_written;

    if (output->prefix_bits) {
        for (size_t i = 0; i < count; ++i) {
            size_t record_prefix = (size_t) (records[i].z >> (64U - output->prefix_bits));

            while (output->next_prefix <= record_prefix) {
                output->prefix_offsets[output->next_prefix++] = (size_t) (first_position + i);
            }
        }
    }

    if (output->density_writer &&
        !geo_density_writer_add(output->density_writer, records, count, first_position)) {
        return false;
    }

    if (fwrite(records, sizeof(*records), count, output->file) != count) {
        return false;
    }

    if (output->checksum_enabled) {
        output->records_checksum = geo_persisted_checksum_update(output->records_checksum,
                                                                 records,
                                                                 count * sizeof(*records));
    }

    output->records_written += count;

    return true;
}

static bool stream_merge_output_append(GeoMergeOutput *output, GeoRecord record)
{
    output->buffer[output->buffered++] = record;

    if (output->buffered < GEO_STREAM_OUTPUT_BUFFER_RECORDS) {
        return true;
    }

    bool succeeded = stream_merge_output_write_block(output, output->buffer, output->buffered);

    output->buffered = 0;

    return succeeded;
}

static bool stream_merge_output_finish(GeoMergeOutput *output,
                                       uint64_t expected_records,
                                       GeoFileHeader *persisted_header)
{
    if (output->buffered) {
        if (!stream_merge_output_write_block(output, output->buffer, output->buffered)) {
            return false;
        }

        output->buffered = 0;
    }

    if (output->records_written != expected_records) {
        return false;
    }

    uint64_t prefix_checksum = geo_persisted_checksum_initial();

    if (output->prefix_bits) {
        while (output->next_prefix < output->prefix_count) {
            output->prefix_offsets[output->next_prefix++] = (size_t) output->records_written;
        }

        if (!geo_io_write_u64_offsets(output->file,
                                      output->prefix_offsets,
                                      output->prefix_count,
                                      &prefix_checksum)) {
            return false;
        }
    }

    uint64_t density_bytes = 0;
    uint64_t density_checksum = 0;

    if (persisted_header &&
        !geo_density_writer_append_file(output->density_writer,
                                        output->file,
                                        &density_bytes,
                                        &density_checksum)) {
        return false;
    }

    if (persisted_header) {
        persisted_header->records_checksum = output->records_checksum;
        persisted_header->prefix_checksum = prefix_checksum;
        persisted_header->density_checksum = density_checksum;
        persisted_header->density_bytes = density_bytes;
    }

    return true;
}

static bool stream_copy_single_run(GeoStreamRun *run, GeoMergeOutput *output)
{
    run->remaining = run->count;

    if (fseek(run->file, 0, SEEK_SET) != 0) {
        return false;
    }

    stream_advise_sequential(run->file);

    while (run->remaining) {
        size_t requested = run->remaining < GEO_STREAM_OUTPUT_BUFFER_RECORDS
                               ? (size_t) run->remaining
                               : GEO_STREAM_OUTPUT_BUFFER_RECORDS;

        if (fread(output->buffer, sizeof(*output->buffer), requested, run->file) != requested) {
            return false;
        }

        run->remaining -= requested;

        if (!stream_merge_output_write_block(output, output->buffer, requested)) {
            return false;
        }
    }

    return true;
}

static bool stream_merge_runs(GeoStreamRun *runs,
                              size_t run_count,
                              const GeoRecord *direct_records,
                              FILE *output,
                              uint64_t expected_records,
                              uint8_t prefix_bits,
                              GeoFileHeader *persisted_header)
{
    if ((run_count && direct_records) ||
        (direct_records && expected_records > SIZE_MAX) ||
        run_count > SIZE_MAX / sizeof(GeoMergeNode) ||
        run_count > SIZE_MAX / sizeof(GeoMergeCursor) ||
        run_count > SIZE_MAX / GEO_STREAM_INPUT_BUFFER_RECORDS ||
        run_count * GEO_STREAM_INPUT_BUFFER_RECORDS > SIZE_MAX / sizeof(GeoRecord) ||
        (prefix_bits && expected_records > SIZE_MAX)) {
        return false;
    }

    bool direct_source = direct_records != NULL;
    bool needs_k_way_merge = run_count > 1;
    GeoMergeNode *heap = needs_k_way_merge ? malloc(run_count * sizeof(*heap)) : NULL;
    GeoMergeCursor *cursors = needs_k_way_merge ? calloc(run_count, sizeof(*cursors)) : NULL;
    size_t input_record_count = needs_k_way_merge ? run_count * GEO_STREAM_INPUT_BUFFER_RECORDS : 0;
    GeoRecord *input_buffers = input_record_count ? malloc(input_record_count * sizeof(*input_buffers)) : NULL;
    GeoRecord *output_buffer = direct_source ? NULL : malloc(GEO_STREAM_OUTPUT_BUFFER_RECORDS * sizeof(*output_buffer));
    size_t prefix_count = prefix_bits ? ((size_t) 1 << prefix_bits) + 1 : 0;
    size_t *prefix_offsets = prefix_count ? malloc(prefix_count * sizeof(*prefix_offsets)) : NULL;
    GeoDensityWriter *density_writer = persisted_header
                                           ? geo_density_writer_create(expected_records, prefix_bits)
                                           : NULL;

    if ((needs_k_way_merge && (!heap || !cursors || !input_buffers)) || (!direct_source && !output_buffer) ||
        (prefix_count && !prefix_offsets) || (persisted_header && !density_writer)) {
        geo_density_writer_destroy(density_writer);
        free(prefix_offsets);
        free(output_buffer);
        free(input_buffers);
        free(cursors);
        free(heap);

        return false;
    }

    GeoMergeOutput merge_output = {
        .file = output,
        .buffer = output_buffer,
        .prefix_offsets = prefix_offsets,
        .prefix_count = prefix_count,
        .prefix_bits = prefix_bits,
        .density_writer = density_writer,
        .records_checksum = geo_persisted_checksum_initial(),
        .checksum_enabled = persisted_header != NULL,
    };
    bool succeeded = true;

    if (direct_source) {
        size_t offset = 0;
        size_t direct_count = (size_t) expected_records;

        while (succeeded && offset < direct_count) {
            size_t remaining = direct_count - offset;
            size_t block_count = remaining < GEO_STREAM_OUTPUT_BUFFER_RECORDS
                                     ? remaining
                                     : GEO_STREAM_OUTPUT_BUFFER_RECORDS;

            succeeded = stream_merge_output_write_block(&merge_output, direct_records + offset, block_count);
            offset += block_count;
        }
    }

    if (run_count == 1) {
        succeeded = stream_copy_single_run(runs, &merge_output);
    }

    size_t heap_count = 0;

    for (size_t run_index = 0; succeeded && needs_k_way_merge && run_index < run_count; ++run_index) {
        GeoStreamRun *run = runs + run_index;
        GeoMergeCursor *cursor = cursors + run_index;
        GeoMergeNode node;

        run->remaining = run->count;
        cursor->run = run;
        cursor->records = input_buffers + run_index * GEO_STREAM_INPUT_BUFFER_RECORDS;

        if (fseek(run->file, 0, SEEK_SET) != 0) {
            succeeded = false;
            break;
        }

        stream_advise_sequential(run->file);

        if (run->remaining && stream_cursor_read(cursor, &node.record)) {
            node.run_index = run_index;
            merge_heap_push(heap, &heap_count, node);
        } else if (run->remaining) {
            succeeded = false;
            break;
        }
    }

    uint64_t merged_records = 0;

    while (succeeded && heap_count) {
        if (merged_records == expected_records) {
            succeeded = false;
            break;
        }

        GeoMergeNode node = heap[0];

        succeeded = stream_merge_output_append(&merge_output, node.record);
        merged_records++;

        GeoMergeCursor *cursor = cursors + node.run_index;

        if (succeeded && stream_cursor_has_next(cursor)) {
            if (!stream_cursor_read(cursor, &node.record)) {
                succeeded = false;
                break;
            }

            merge_heap_replace_min(heap, heap_count, node);
        } else if (succeeded) {
            (void) merge_heap_pop(heap, &heap_count);
        }
    }

    if (succeeded) {
        succeeded = stream_merge_output_finish(&merge_output, expected_records, persisted_header);
    }

    geo_density_writer_destroy(density_writer);
    free(prefix_offsets);
    free(output_buffer);
    free(input_buffers);
    free(cursors);
    free(heap);

    return succeeded;
}

static bool stream_compact_tail(GeoStreamBuilder *builder)
{
    while (builder->run_count >= GEO_STREAM_MERGE_FAN_IN) {
        size_t first = builder->run_count - GEO_STREAM_MERGE_FAN_IN;
        uint32_t level = builder->runs[first].level;

        for (size_t i = first + 1; i < builder->run_count; ++i) {
            if (builder->runs[i].level != level) {
                return true;
            }
        }

        if (level == UINT32_MAX) {
            return false;
        }

        uint64_t merged_count = 0;

        for (size_t i = first; i < builder->run_count; ++i) {
            if (builder->runs[i].count > UINT64_MAX - merged_count) {
                return false;
            }

            merged_count += builder->runs[i].count;
        }

        FILE *merged_file = stream_create_temporary_run(builder->temporary_directory);

        if (!merged_file) {
            return false;
        }

        stream_advise_sequential(merged_file);
        stream_preallocate(merged_file, merged_count, 0);

        if (builder->run_count + 1 > builder->peak_open_runs) {
            builder->peak_open_runs = builder->run_count + 1;
        }

        double merge_start = geo_get_time_ms();
        bool succeeded = stream_merge_runs(builder->runs + first,
                                           GEO_STREAM_MERGE_FAN_IN,
                                           NULL,
                                           merged_file,
                                           merged_count,
                                           0,
                                           NULL) &&
                         fflush(merged_file) == 0;

        builder->intermediate_merge_time_ms += geo_get_time_ms() - merge_start;

        if (!succeeded) {
            fclose(merged_file);

            return false;
        }

        for (size_t i = first; i < builder->run_count; ++i) {
            fclose(builder->runs[i].file);
        }

        builder->runs[first] = (GeoStreamRun) {
            .file = merged_file,
            .count = merged_count,
            .remaining = 0,
            .level = level + 1,
        };
        builder->run_count = first + 1;
        builder->intermediate_merges++;
    }

    return true;
}

bool geo_stream_builder_finish(GeoStreamBuilder *builder, GeoStreamBuildStats *stats)
{
#if GEO_STREAM_SUPPORTED
    if (stats) {
        memset(stats, 0, sizeof(*stats));
    }

    if (!builder || builder->failed || builder->finished) {
        return false;
    }

    bool direct_chunk = builder->run_count == 0;

    if (direct_chunk) {
        double sort_start = geo_get_time_ms();

        if (!stream_sort_chunk(builder)) {
            builder->failed = true;

            return false;
        }

        builder->chunk_build_time_ms += geo_get_time_ms() - sort_start;
    } else if (!stream_flush_chunk(builder)) {
        builder->failed = true;

        return false;
    }

    char *temporary_path = NULL;
    FILE *output = geo_io_create_atomic_file(builder->output_path, &temporary_path);

    if (!output) {
        builder->failed = true;

        return false;
    }

    uint8_t prefix_bits = geo_persisted_prefix_bits(builder->total_records);
    GeoFileHeader header = {
        .magic = { 0 },
        .version = GEO_FILE_VERSION,
        .record_size = sizeof(GeoRecord),
        .endian_marker = GEO_FILE_ENDIAN_MARKER,
        .prefix_bits = prefix_bits,
        .count = builder->total_records,
    };

    memcpy(header.magic, GEO_FILE_MAGIC, sizeof(header.magic));

    size_t prefix_count = prefix_bits ? ((size_t) 1 << prefix_bits) + 1 : 0;
    size_t metadata_bytes = sizeof(header) + prefix_count * sizeof(uint64_t);

    stream_advise_sequential(output);
    stream_preallocate(output, builder->total_records, metadata_bytes);

    double merge_start = geo_get_time_ms();
    bool succeeded = fwrite(&header, sizeof(header), 1, output) == 1 &&
                     stream_merge_runs(builder->runs,
                                       direct_chunk ? 0 : builder->run_count,
                                       direct_chunk ? builder->chunk->records : NULL,
                                       output,
                                       builder->total_records,
                                       prefix_bits,
                                       &header);

    if (succeeded) {
        succeeded = fseeko(output, 0, SEEK_SET) == 0 && fwrite(&header, sizeof(header), 1, output) == 1;
    }

    if (succeeded) {
        succeeded = geo_io_publish_atomic_file(output, temporary_path, builder->output_path);
    } else {
        geo_io_discard_atomic_file(output, temporary_path);
    }

    if (!succeeded) {
        builder->failed = true;
    }

    double merge_time = geo_get_time_ms() - merge_start;

    free(temporary_path);

    if (!succeeded) {
        return false;
    }

    builder->finished = true;

    if (stats) {
        stats->records_written = builder->total_records;
        stats->runs_created = builder->runs_created;
        stats->intermediate_merges = builder->intermediate_merges;
        stats->peak_open_runs = builder->peak_open_runs;
        stats->chunk_build_time_ms = builder->chunk_build_time_ms;
        stats->intermediate_merge_time_ms = builder->intermediate_merge_time_ms;
        stats->merge_time_ms = merge_time;
    }

    return true;
#else
    (void) builder;
    (void) stats;

    return false;
#endif
}
