#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "geobolt/geo_index.h"
#include "geo_index_internal.h"
#include "geo_index_io.h"
#include "geo_index_persistence.h"
#include "geo_index_private.h"
#include "geo_thread_pool.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <pthread.h>
#include <sched.h>
#include <sys/stat.h>
#include <unistd.h>
#define GEO_SEGMENTS_SUPPORTED 1
#else
#define GEO_SEGMENTS_SUPPORTED 0
#endif

#define GEO_MANIFEST_VERSION 5U
#define GEO_MANIFEST_ENDIAN_MARKER UINT32_C(0x01020304)
#define GEO_MANIFEST_MAX_SEGMENTS UINT64_C(1048576)
#define GEO_MANIFEST_MAX_PATH_BYTES UINT32_C(1048576)
#define GEO_SEGMENT_OUTPUT_BUFFER_RECORDS 4096
#define GEO_COMPACTION_PARTITION_BITS 16U
#define GEO_COMPACTION_MIN_RECORDS_PER_WORKER UINT64_C(500000)
#define GEO_COMPACTION_MAX_WORKERS 16U
#define GEO_COMPACTION_TASKS_PER_WORKER 4U
#define GEO_BACKGROUND_COMPACTION_FANOUT 4
#define GEO_READER_SHARD_COUNT 128U
#define GEO_CACHE_LINE_SIZE 64U
#define GEO_MUTATION_CHECKPOINT_MIN_RECORDS UINT64_C(4096)
#define GEO_MUTATION_CHECKPOINT_AMPLIFICATION 4U
#define GEO_MUTATION_STATE_BITS 2U
#define GEO_MUTATION_STATE_MASK ((UINT64_C(1) << GEO_MUTATION_STATE_BITS) - 1U)
#define GEO_MUTATION_MAX_GENERATION (UINT64_MAX >> GEO_MUTATION_STATE_BITS)
#define GEO_MUTATION_CONTROL_GROUP_WIDTH 8U
#define GEO_MUTATION_CONTROL_EMPTY 0U
#define GEO_MUTATION_SMALL_ENTRY_CAPACITY 8U
#define GEO_BACKGROUND_DEFAULT_MUTATION_MIN_ENTRIES 4096U
#define GEO_BACKGROUND_DEFAULT_MUTATION_RATIO_NUMERATOR 1U
#define GEO_BACKGROUND_DEFAULT_MUTATION_RATIO_DENOMINATOR 8U
#define GEO_BACKGROUND_DEFAULT_MUTATION_MAX_PASSES 2U

_Static_assert(GEO_MUTATION_CONTROL_GROUP_WIDTH == sizeof(uint64_t), "Mutation control probes operate on one machine word");
_Static_assert((GEO_MUTATION_CONTROL_GROUP_WIDTH & (GEO_MUTATION_CONTROL_GROUP_WIDTH - 1U)) == 0,
               "Mutation control group width must be a power of two");

typedef struct {
    char magic[8];
    uint32_t version;
    uint32_t endian_marker;
    uint64_t generation;
    uint64_t segment_count;
    uint64_t mutation_count;
    uint64_t mutation_checksum;
    uint64_t mutation_slot;
    uint64_t durable_watermark;
    uint64_t checksum;
} GeoManifestHeader;

typedef struct {
    uint32_t path_length;
    uint32_t flags;
    uint64_t record_count;
    uint64_t generation;
    uint64_t content_checksum;
} GeoManifestEntry;

_Static_assert(sizeof(GeoManifestHeader) == 72, "GeoManifestHeader layout is persisted");
_Static_assert(offsetof(GeoManifestHeader, durable_watermark) == 56,
               "GeoManifestHeader.durable_watermark offset is persisted");
_Static_assert(offsetof(GeoManifestHeader, checksum) == 64, "GeoManifestHeader.checksum offset is persisted");
_Static_assert(sizeof(GeoManifestEntry) == 32, "GeoManifestEntry layout is persisted");
_Static_assert(offsetof(GeoManifestEntry, record_count) == 8, "GeoManifestEntry.record_count offset is persisted");
_Static_assert(offsetof(GeoManifestEntry, content_checksum) == 24, "GeoManifestEntry.content_checksum offset is persisted");

typedef enum {
    GEO_MUTATION_RESET = 0,
    GEO_MUTATION_DELETE = 1,
    GEO_MUTATION_LIVE = 2,
} GeoMutationState;

_Static_assert(GEO_MUTATION_LIVE <= GEO_MUTATION_STATE_MASK, "Mutation state must fit packed metadata");

typedef struct {
    uint64_t id;
    uint64_t generation;
    uint64_t state;
} GeoMutationRecord;

typedef struct {
    uint64_t id;
    // Zero marks an unused hash slot. Occupied slots store generation in the high 62 bits and state in the low two bits.
    uint64_t metadata;
} GeoMutationEntry;

_Static_assert(sizeof(GeoMutationEntry) == 16, "Mutation hash entries must pack four per cache line");

typedef struct {
    GeoMutationEntry *entries;
    uint8_t *controls;
    size_t capacity;
    size_t count;
} GeoMutationTable;

typedef struct {
    GeoRecord record;
    size_t segment_index;
    size_t position;
} GeoSegmentMergeNode;

typedef enum {
    GEO_BACKGROUND_COMPACTION_NONE = 0,
    GEO_BACKGROUND_COMPACTION_SIZE_TIER,
    GEO_BACKGROUND_COMPACTION_MUTATION_REWRITE,
} GeoBackgroundCompactionMode;

typedef struct {
    size_t index;
    size_t record_count;
} GeoSegmentCandidate;

typedef struct {
    _Alignas(GEO_CACHE_LINE_SIZE) atomic_size_t active;
} GeoReaderShard;

_Static_assert((GEO_READER_SHARD_COUNT & (GEO_READER_SHARD_COUNT - 1U)) == 0,
               "GEO_READER_SHARD_COUNT must be a power of two");
_Static_assert(sizeof(GeoReaderShard) >= GEO_CACHE_LINE_SIZE, "Reader shards must not share a cache line");
_Static_assert(sizeof(GeoReaderShard) % GEO_CACHE_LINE_SIZE == 0, "Reader shard stride must preserve cache-line isolation");

struct GeoSegmentSet {
    char *manifest_path;
    char **paths;
    GeoIndex **segments;
    uint64_t *segment_generations;
    GeoMutationTable mutations;
    char *mutation_path;
    size_t count;
    size_t capacity;
    uint64_t generation;
    uint64_t record_count;
    uint64_t mutation_count;
    uint64_t mutation_checksum;
    uint64_t mutation_slot;
    uint64_t durable_watermark;
#if GEO_SEGMENTS_SUPPORTED
    // Readers enter a QSBR grace period with atomics only. Writers close the gate,
    // wait for the current epoch to quiesce, publish, and then reclaim old arrays.
    char *background_directory;
    pthread_mutex_t mutation_lock;
    pthread_mutex_t background_lock;
    pthread_cond_t background_condition;
    pthread_t background_thread;
    GeoThreadPool *compaction_pool;
    atomic_bool writer_pending;
    GeoReaderShard *reader_shards;
    GeoSegmentCompactionPolicy background_policy;
    GeoSegmentCompactionStats background_last_stats;
    GeoSegmentCompactionStats background_total_stats;
    uint64_t background_completed_runs;
    uint64_t background_failed_runs;
    bool background_enabled;
    bool background_running;
    bool background_requested;
    bool background_thread_started;
    bool background_stop;
    bool background_last_succeeded;
    bool compaction_running;
    bool locks_initialized;
#endif
};

static const char GEO_MANIFEST_MAGIC[8] = { 'G', 'B', 'M', 'A', 'N', 'I', '1', '\0' };

#if GEO_SEGMENTS_SUPPORTED
static atomic_uint_fast64_t GEO_COMPACTION_SEQUENCE = ATOMIC_VAR_INIT(0);

static void segment_schedule_background_compaction(GeoSegmentSet *set);
static bool segment_set_compact_size_tier(GeoSegmentSet *set,
                                          const char *output_path,
                                          GeoSegmentCompactionStats *stats);
static bool mutation_checkpoint_locked(GeoSegmentSet *set);
static void mutation_checkpoint_if_amplified(GeoSegmentSet *set);
#endif

// =============================================================================
// Manifest encoding and atomic publication
// =============================================================================

static char *segment_copy_string(const char *source)
{
    size_t length = strlen(source);
    char *copy = malloc(length + 1);

    if (copy) {
        memcpy(copy, source, length + 1);
    }

    return copy;
}

static char *segment_mutation_path(const char *manifest_path, uint64_t slot)
{
    static const char primary_suffix[] = ".mutations";
    static const char checkpoint_suffix[] = ".mutations.checkpoint";
    const char *suffix = slot ? checkpoint_suffix : primary_suffix;
    size_t length = strlen(manifest_path);
    size_t suffix_length = strlen(suffix) + 1U;

    if (slot > 1U || length > SIZE_MAX - suffix_length) {
        return NULL;
    }

    char *path = malloc(length + suffix_length);

    if (path) {
        memcpy(path, manifest_path, length);
        memcpy(path + length, suffix, suffix_length);
    }

    return path;
}

static uint64_t segment_hash_id(uint64_t id)
{
    id ^= id >> 30U;
    id *= UINT64_C(0xbf58476d1ce4e5b9);
    id ^= id >> 27U;
    id *= UINT64_C(0x94d049bb133111eb);

    return id ^ (id >> 31U);
}

static inline bool mutation_entry_is_occupied(const GeoMutationEntry *entry)
{
    return entry->metadata != 0;
}

static inline uint64_t mutation_entry_generation(const GeoMutationEntry *entry)
{
    return entry->metadata >> GEO_MUTATION_STATE_BITS;
}

static inline GeoMutationState mutation_entry_state(const GeoMutationEntry *entry)
{
    return (GeoMutationState) (entry->metadata & GEO_MUTATION_STATE_MASK);
}

static inline void mutation_entry_set(GeoMutationEntry *entry, uint64_t id, uint64_t generation, GeoMutationState state)
{
    entry->id = id;
    entry->metadata = (generation << GEO_MUTATION_STATE_BITS) | (uint64_t) state;
}

static inline uint8_t mutation_hash_fingerprint(uint64_t hash)
{
    uint8_t fingerprint = (uint8_t) (hash >> 56U);

    return fingerprint ? fingerprint : 1U;
}

static inline void mutation_control_set(GeoMutationTable *table, size_t position, uint8_t control)
{
    table->controls[position] = control;

    // The mirrored prefix makes every eight-byte probe contiguous, including a probe that wraps at the table boundary.
    if (position < GEO_MUTATION_CONTROL_GROUP_WIDTH - 1U) {
        table->controls[table->capacity + position] = control;
    }
}

static inline uint64_t mutation_zero_byte_candidates(uint64_t value)
{
    const uint64_t low_bits = UINT64_C(0x0101010101010101);
    const uint64_t high_bits = UINT64_C(0x8080808080808080);

    // Borrow propagation can add conservative candidates. Callers verify the selected control byte before using it.
    return (value - low_bits) & ~value & high_bits;
}

static inline unsigned mutation_candidate_lane(uint64_t candidates)
{
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return (unsigned) (__builtin_clzll(candidates) >> 3U);
#else
    return (unsigned) (__builtin_ctzll(candidates) >> 3U);
#endif
}

static inline uint64_t mutation_remove_candidate_lane(uint64_t candidates, unsigned lane)
{
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return candidates & ~(UINT64_C(0x80) << ((7U - lane) * 8U));
#else
    return candidates & ~(UINT64_C(0x80) << (lane * 8U));
#endif
}

static GeoMutationEntry *mutation_table_find_slot(const GeoMutationTable *table, uint64_t id, bool *found)
{
    *found = false;

    if (!table->capacity) {
        return NULL;
    }

    uint64_t hash = segment_hash_id(id);
    uint8_t fingerprint = mutation_hash_fingerprint(hash);
    uint64_t repeated_fingerprint = UINT64_C(0x0101010101010101) * fingerprint;
    size_t mask = table->capacity - 1U;
    size_t group_position = (size_t) hash & mask;

    for (;;) {
        uint64_t control_word;

        memcpy(&control_word, table->controls + group_position, sizeof(control_word));

        uint64_t matching = mutation_zero_byte_candidates(control_word ^ repeated_fingerprint);
        uint64_t empty = mutation_zero_byte_candidates(control_word);
        unsigned first_empty_lane = GEO_MUTATION_CONTROL_GROUP_WIDTH;

        while (empty) {
            unsigned lane = mutation_candidate_lane(empty);

            if (table->controls[group_position + lane] == GEO_MUTATION_CONTROL_EMPTY) {
                first_empty_lane = lane;
                break;
            }

            empty = mutation_remove_candidate_lane(empty, lane);
        }

        while (matching) {
            unsigned lane = mutation_candidate_lane(matching);

            if (lane >= first_empty_lane) {
                break;
            }

            size_t position = (group_position + lane) & mask;

            if (table->controls[group_position + lane] == fingerprint && table->entries[position].id == id) {
                *found = true;

                return table->entries + position;
            }

            matching = mutation_remove_candidate_lane(matching, lane);
        }

        if (first_empty_lane < GEO_MUTATION_CONTROL_GROUP_WIDTH) {
            size_t position = (group_position + first_empty_lane) & mask;

            return table->entries + position;
        }

        group_position = (group_position + GEO_MUTATION_CONTROL_GROUP_WIDTH) & mask;
    }
}

static const GeoMutationEntry *mutation_table_lookup(const GeoMutationTable *table, uint64_t id)
{
    bool found;
    GeoMutationEntry *entry = mutation_table_find_slot(table, id, &found);

    return found ? entry : NULL;
}

static void mutation_table_release_storage(GeoMutationTable *table)
{
    free(table->controls);
    free(table->entries);
    memset(table, 0, sizeof(*table));
}

static bool mutation_table_copy_entry(GeoMutationTable *table, const GeoMutationEntry *source)
{
    bool found;
    GeoMutationEntry *destination = mutation_table_find_slot(table, source->id, &found);

    if (!destination || !table->entries) {
        return false;
    }

    if (!found) {
        size_t position = (size_t) (destination - table->entries);

        mutation_control_set(table, position, mutation_hash_fingerprint(segment_hash_id(source->id)));
        table->count++;
    }

    *destination = *source;

    return true;
}

static bool mutation_table_rehash(GeoMutationTable *table, size_t capacity)
{
    if (capacity > SIZE_MAX / sizeof(GeoMutationEntry) ||
        capacity > SIZE_MAX - (GEO_MUTATION_CONTROL_GROUP_WIDTH - 1U)) {
        return false;
    }

    GeoMutationEntry *entries = calloc(capacity, sizeof(*entries));
    uint8_t *controls = calloc(capacity + GEO_MUTATION_CONTROL_GROUP_WIDTH - 1U, sizeof(*controls));

    if (!entries || !controls) {
        free(controls);
        free(entries);

        return false;
    }

    GeoMutationTable replacement = {
        .entries = entries,
        .controls = controls,
        .capacity = capacity,
    };

    for (size_t i = 0; i < table->capacity; ++i) {
        if (!mutation_entry_is_occupied(table->entries + i)) {
            continue;
        }

        if (!mutation_table_copy_entry(&replacement, table->entries + i)) {
            mutation_table_release_storage(&replacement);

            return false;
        }
    }

    mutation_table_release_storage(table);
    *table = replacement;

    return true;
}

static bool mutation_table_reserve(GeoMutationTable *table, size_t required)
{
    if (required <= table->capacity - table->capacity / 4U) {
        return true;
    }

    size_t capacity = table->capacity ? table->capacity : 16;

    while (required > capacity - capacity / 4U) {
        if (capacity > SIZE_MAX / 2U) {
            return false;
        }

        capacity *= 2U;
    }

    return mutation_table_rehash(table, capacity);
}

static bool mutation_table_apply(GeoMutationTable *table, const GeoMutationRecord *mutation)
{
    bool found;
    GeoMutationEntry *entry = mutation_table_find_slot(table, mutation->id, &found);

    if (!entry) {
        if (!mutation_table_reserve(table, table->count + 1U)) {
            return false;
        }

        entry = mutation_table_find_slot(table, mutation->id, &found);
    }

    if (!found) {
        size_t position = (size_t) (entry - table->entries);

        mutation_entry_set(entry, mutation->id, mutation->generation, (GeoMutationState) mutation->state);
        mutation_control_set(table, position, mutation_hash_fingerprint(segment_hash_id(mutation->id)));
        table->count++;
    }

    if (mutation->generation >= mutation_entry_generation(entry)) {
        mutation_entry_set(entry, mutation->id, mutation->generation, (GeoMutationState) mutation->state);
    }

    return true;
}

static bool mutation_table_clone(const GeoMutationTable *source, GeoMutationTable *destination, size_t additional)
{
    memset(destination, 0, sizeof(*destination));

    if (source->count > SIZE_MAX - additional || !mutation_table_reserve(destination, source->count + additional)) {
        return false;
    }

    for (size_t i = 0; i < source->capacity; ++i) {
        if (mutation_entry_is_occupied(source->entries + i)) {
            if (!mutation_table_copy_entry(destination, source->entries + i)) {
                mutation_table_release_storage(destination);

                return false;
            }
        }
    }

    return true;
}

static void mutation_table_destroy(GeoMutationTable *table)
{
    mutation_table_release_storage(table);
}

static bool mutation_table_discard_resets(GeoMutationTable *table)
{
    size_t active_count = 0;

    for (size_t i = 0; i < table->capacity; ++i) {
        const GeoMutationEntry *entry = table->entries + i;

        active_count += mutation_entry_is_occupied(entry) && mutation_entry_state(entry) != GEO_MUTATION_RESET;
    }

    GeoMutationTable compacted = { 0 };

    if (!mutation_table_reserve(&compacted, active_count)) {
        return false;
    }

    for (size_t i = 0; i < table->capacity; ++i) {
        const GeoMutationEntry *entry = table->entries + i;

        if (!mutation_entry_is_occupied(entry) || mutation_entry_state(entry) == GEO_MUTATION_RESET) {
            continue;
        }

        GeoMutationRecord mutation = {
            .id = entry->id,
            .generation = mutation_entry_generation(entry),
            .state = mutation_entry_state(entry),
        };

        if (!mutation_table_apply(&compacted, &mutation)) {
            mutation_table_destroy(&compacted);

            return false;
        }
    }

    mutation_table_destroy(table);
    *table = compacted;

    return true;
}

static uint64_t manifest_checksum_update(uint64_t checksum, const void *data, size_t size)
{
    const uint8_t *bytes = data;

    for (size_t i = 0; i < size; ++i) {
        checksum ^= bytes[i];
        checksum *= UINT64_C(1099511628211);
    }

    return checksum;
}

static bool mutation_log_create(const char *path)
{
    int descriptor = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0644);

    if (descriptor < 0) {
        return false;
    }

    bool succeeded = fsync(descriptor) == 0;

    if (close(descriptor) != 0) {
        succeeded = false;
    }

    if (!succeeded) {
        unlink(path);
    }

    return succeeded;
}

static bool mutation_log_append(GeoSegmentSet *set,
                                const GeoMutationRecord *mutations,
                                size_t count,
                                uint64_t *checksum)
{
    if (!count) {
        *checksum = set->mutation_checksum;

        return true;
    }

    if (set->mutation_count > (uint64_t) (INT64_MAX / (int64_t) sizeof(*mutations)) ||
        count > UINT64_MAX - set->mutation_count) {
        return false;
    }

    FILE *file = fopen(set->mutation_path, "r+b");

    if (!file) {
        return false;
    }

    off_t committed_size = (off_t) (set->mutation_count * sizeof(*mutations));
    bool succeeded = ftruncate(fileno(file), committed_size) == 0 &&
                     fseeko(file, committed_size, SEEK_SET) == 0 &&
                     fwrite(mutations, sizeof(*mutations), count, file) == count &&
                     fflush(file) == 0 &&
                     fsync(fileno(file)) == 0;

    if (fclose(file) != 0) {
        succeeded = false;
    }

    if (!succeeded) {
        return false;
    }

    uint64_t updated = set->mutation_checksum;

    for (size_t i = 0; i < count; ++i) {
        updated = manifest_checksum_update(updated, mutations + i, sizeof(*mutations));
    }

    *checksum = updated;

    return true;
}

static bool mutation_log_load(GeoSegmentSet *set)
{
    FILE *file = fopen(set->mutation_path, "rb");

    if (!file) {
        return false;
    }

    uint64_t checksum = UINT64_C(1469598103934665603);
    bool succeeded = set->mutation_count <= SIZE_MAX &&
                     mutation_table_reserve(&set->mutations, (size_t) set->mutation_count);

    for (uint64_t i = 0; succeeded && i < set->mutation_count; ++i) {
        GeoMutationRecord mutation;

        succeeded = fread(&mutation, sizeof(mutation), 1, file) == 1 &&
                    mutation.generation > 0 &&
                    mutation.generation <= set->generation &&
                    mutation.state <= GEO_MUTATION_LIVE;

        if (succeeded) {
            checksum = manifest_checksum_update(checksum, &mutation, sizeof(mutation));
            succeeded = mutation_table_apply(&set->mutations, &mutation);
        }
    }

    if (succeeded) {
        GeoMutationRecord ignored_tail;

        // A failed manifest publication may leave a valid uncommitted tail. It is intentionally ignored.
        (void) fread(&ignored_tail, sizeof(ignored_tail), 1, file);
        succeeded = !ferror(file) && checksum == set->mutation_checksum;
    }

    fclose(file);

    if (succeeded) {
        succeeded = mutation_table_discard_resets(&set->mutations);
    }

    return succeeded;
}

static uint64_t manifest_checksum(char *const *paths,
                                  GeoIndex *const *segments,
                                  const uint64_t *segment_generations,
                                  size_t count,
                                  uint64_t generation,
                                  uint64_t mutation_count,
                                  uint64_t mutation_checksum,
                                  uint64_t mutation_slot,
                                  uint64_t durable_watermark)
{
    uint64_t checksum = UINT64_C(1469598103934665603);
    uint64_t segment_count = count;

    checksum = manifest_checksum_update(checksum, &generation, sizeof(generation));
    checksum = manifest_checksum_update(checksum, &segment_count, sizeof(segment_count));
    checksum = manifest_checksum_update(checksum, &mutation_count, sizeof(mutation_count));
    checksum = manifest_checksum_update(checksum, &mutation_checksum, sizeof(mutation_checksum));
    checksum = manifest_checksum_update(checksum, &mutation_slot, sizeof(mutation_slot));
    checksum = manifest_checksum_update(checksum, &durable_watermark, sizeof(durable_watermark));

    for (size_t i = 0; i < count; ++i) {
        GeoManifestEntry entry = {
            .path_length = (uint32_t) strlen(paths[i]),
            .flags = 0,
            .record_count = segments[i]->count,
            .generation = segment_generations[i],
            .content_checksum = segments[i]->content_checksum,
        };

        checksum = manifest_checksum_update(checksum, &entry, sizeof(entry));
        checksum = manifest_checksum_update(checksum, paths[i], entry.path_length);
    }

    return checksum;
}

static bool manifest_write(const char *manifest_path,
                           char *const *paths,
                           GeoIndex *const *segments,
                           const uint64_t *segment_generations,
                           size_t count,
                           uint64_t generation,
                           uint64_t mutation_count,
                           uint64_t mutation_checksum,
                           uint64_t mutation_slot,
                           uint64_t durable_watermark)
{
    if (count > GEO_MANIFEST_MAX_SEGMENTS || mutation_slot > 1U) {
        return false;
    }

    for (size_t i = 0; i < count; ++i) {
        size_t path_length = strlen(paths[i]);

        if (!path_length || path_length > GEO_MANIFEST_MAX_PATH_BYTES || path_length > UINT32_MAX) {
            return false;
        }
    }

    GeoManifestHeader header = {
        .magic = { 0 },
        .version = GEO_MANIFEST_VERSION,
        .endian_marker = GEO_MANIFEST_ENDIAN_MARKER,
        .generation = generation,
        .segment_count = count,
        .mutation_count = mutation_count,
        .mutation_checksum = mutation_checksum,
        .mutation_slot = mutation_slot,
        .durable_watermark = durable_watermark,
        .checksum = manifest_checksum(paths,
                                      segments,
                                      segment_generations,
                                      count,
                                      generation,
                                      mutation_count,
                                      mutation_checksum,
                                      mutation_slot,
                                      durable_watermark),
    };

    memcpy(header.magic, GEO_MANIFEST_MAGIC, sizeof(header.magic));

    char *temporary_path = NULL;
    FILE *file = geo_io_create_atomic_file(manifest_path, &temporary_path);

    if (!file) {
        return false;
    }

    bool succeeded = fwrite(&header, sizeof(header), 1, file) == 1;

    for (size_t i = 0; succeeded && i < count; ++i) {
        GeoManifestEntry entry = {
            .path_length = (uint32_t) strlen(paths[i]),
            .flags = 0,
            .record_count = segments[i]->count,
            .generation = segment_generations[i],
            .content_checksum = segments[i]->content_checksum,
        };

        succeeded = fwrite(&entry, sizeof(entry), 1, file) == 1 &&
                    fwrite(paths[i], 1, entry.path_length, file) == entry.path_length;
    }

    if (succeeded) {
        succeeded = geo_io_publish_atomic_file(file, temporary_path, manifest_path);
    } else {
        geo_io_discard_atomic_file(file, temporary_path);
    }

    free(temporary_path);

    return succeeded;
}

static GeoSegmentSet *segment_set_allocate(const char *manifest_path)
{
    GeoSegmentSet *set = calloc(1, sizeof(*set));

    if (!set) {
        return NULL;
    }

    set->manifest_path = segment_copy_string(manifest_path);
    set->mutation_path = segment_mutation_path(manifest_path, 0);

    if (!set->manifest_path || !set->mutation_path) {
        free(set->mutation_path);
        free(set->manifest_path);
        free(set);

        return NULL;
    }

#if GEO_SEGMENTS_SUPPORTED
    if (pthread_mutex_init(&set->mutation_lock, NULL) != 0) {
        free(set->manifest_path);
        free(set->mutation_path);
        free(set);

        return NULL;
    }

    if (pthread_mutex_init(&set->background_lock, NULL) != 0) {
        pthread_mutex_destroy(&set->mutation_lock);
        free(set->manifest_path);
        free(set->mutation_path);
        free(set);

        return NULL;
    }

    if (pthread_cond_init(&set->background_condition, NULL) != 0) {
        pthread_mutex_destroy(&set->background_lock);
        pthread_mutex_destroy(&set->mutation_lock);
        free(set->manifest_path);
        free(set->mutation_path);
        free(set);

        return NULL;
    }

    set->reader_shards = aligned_alloc(GEO_CACHE_LINE_SIZE,
                                       GEO_READER_SHARD_COUNT * sizeof(*set->reader_shards));

    if (!set->reader_shards) {
        pthread_cond_destroy(&set->background_condition);
        pthread_mutex_destroy(&set->background_lock);
        pthread_mutex_destroy(&set->mutation_lock);
        free(set->manifest_path);
        free(set->mutation_path);
        free(set);

        return NULL;
    }

    atomic_init(&set->writer_pending, false);

    for (size_t shard = 0; shard < GEO_READER_SHARD_COUNT; ++shard) {
        atomic_init(&set->reader_shards[shard].active, 0);
    }

    set->locks_initialized = true;
#endif

    return set;
}

static size_t segment_reader_shard_index(void)
{
    static _Thread_local unsigned char thread_marker;
    uintptr_t value = (uintptr_t) &thread_marker;

    value ^= value >> 17U;
    value *= UINT64_C(0xed5ad4bb);
    value ^= value >> 11U;

    return (size_t) value & (GEO_READER_SHARD_COUNT - 1U);
}

static bool segment_set_read_lock(const GeoSegmentSet *set)
{
#if GEO_SEGMENTS_SUPPORTED
    GeoSegmentSet *mutable_set = (GeoSegmentSet *) set;
    GeoReaderShard *shard = mutable_set->reader_shards + segment_reader_shard_index();

    for (;;) {
        while (atomic_load_explicit(&mutable_set->writer_pending, memory_order_acquire)) {
            sched_yield();
        }

        atomic_fetch_add_explicit(&shard->active, 1, memory_order_acquire);

        if (!atomic_load_explicit(&mutable_set->writer_pending, memory_order_acquire)) {
            return true;
        }

        atomic_fetch_sub_explicit(&shard->active, 1, memory_order_release);
    }
#else
    (void) set;

    return true;
#endif
}

static bool segment_set_write_lock(GeoSegmentSet *set)
{
#if GEO_SEGMENTS_SUPPORTED
    atomic_store_explicit(&set->writer_pending, true, memory_order_release);

    for (;;) {
        bool readers_active = false;

        for (size_t shard = 0; shard < GEO_READER_SHARD_COUNT; ++shard) {
            if (atomic_load_explicit(&set->reader_shards[shard].active, memory_order_acquire) != 0) {
                readers_active = true;
                break;
            }
        }

        if (!readers_active) {
            break;
        }

        sched_yield();
    }

    atomic_thread_fence(memory_order_acquire);

    return true;
#else
    (void) set;

    return true;
#endif
}

static void segment_set_write_unlock(GeoSegmentSet *set)
{
#if GEO_SEGMENTS_SUPPORTED
    atomic_thread_fence(memory_order_release);
    atomic_store_explicit(&set->writer_pending, false, memory_order_release);
#else
    (void) set;
#endif
}

static void segment_set_read_unlock(const GeoSegmentSet *set)
{
#if GEO_SEGMENTS_SUPPORTED
    GeoSegmentSet *mutable_set = (GeoSegmentSet *) set;
    GeoReaderShard *shard = mutable_set->reader_shards + segment_reader_shard_index();

    atomic_fetch_sub_explicit(&shard->active, 1, memory_order_release);
#else
    (void) set;
#endif
}

bool geo_segment_set_read_snapshot_acquire(const GeoSegmentSet *set)
{
    return set && segment_set_read_lock(set);
}

void geo_segment_set_read_snapshot_release(const GeoSegmentSet *set)
{
    if (set) {
        segment_set_read_unlock(set);
    }
}

static bool segment_set_reserve(GeoSegmentSet *set, size_t capacity)
{
    if (capacity <= set->capacity) {
        return true;
    }

    if (capacity > SIZE_MAX / sizeof(*set->paths) ||
        capacity > SIZE_MAX / sizeof(*set->segments) ||
        capacity > SIZE_MAX / sizeof(*set->segment_generations)) {
        return false;
    }

    char **paths = malloc(capacity * sizeof(*paths));
    GeoIndex **segments = malloc(capacity * sizeof(*segments));
    uint64_t *segment_generations = malloc(capacity * sizeof(*segment_generations));

    if (!paths || !segments || !segment_generations) {
        free(segment_generations);
        free(segments);
        free(paths);

        return false;
    }

    if (set->count) {
        memcpy(paths, set->paths, set->count * sizeof(*paths));
        memcpy(segments, set->segments, set->count * sizeof(*segments));
        memcpy(segment_generations,
               set->segment_generations,
               set->count * sizeof(*segment_generations));
    }

    free(set->segments);
    free(set->paths);
    free(set->segment_generations);
    set->paths = paths;
    set->segments = segments;
    set->segment_generations = segment_generations;
    set->capacity = capacity;

    return true;
}

// =============================================================================
// Segment-set lifecycle
// =============================================================================

GeoSegmentSet *geo_segment_set_create(const char *manifest_path)
{
#if GEO_SEGMENTS_SUPPORTED
    if (!manifest_path || !manifest_path[0] || access(manifest_path, F_OK) == 0) {
        return NULL;
    }

    GeoSegmentSet *set = segment_set_allocate(manifest_path);

    if (!set || !mutation_log_create(set->mutation_path) ||
        !manifest_write(manifest_path,
                        NULL,
                        NULL,
                        NULL,
                        0,
                        1,
                        0,
                        UINT64_C(1469598103934665603),
                        0,
                        0)) {
        if (set) {
            unlink(set->mutation_path);
        }
        geo_segment_set_destroy(set);

        return NULL;
    }

    set->generation = 1;
    set->mutation_checksum = UINT64_C(1469598103934665603);

    return set;
#else
    (void) manifest_path;

    return NULL;
#endif
}

GeoSegmentSet *geo_segment_set_open(const char *manifest_path)
{
#if GEO_SEGMENTS_SUPPORTED
    if (!manifest_path || !manifest_path[0]) {
        return NULL;
    }

    FILE *file = fopen(manifest_path, "rb");

    if (!file) {
        return NULL;
    }

    GeoManifestHeader header;
    bool succeeded = fread(&header, sizeof(header), 1, file) == 1 &&
                     memcmp(header.magic, GEO_MANIFEST_MAGIC, sizeof(header.magic)) == 0 &&
                     header.version == GEO_MANIFEST_VERSION &&
                     header.endian_marker == GEO_MANIFEST_ENDIAN_MARKER &&
                     header.generation > 0 &&
                     header.generation <= GEO_MUTATION_MAX_GENERATION &&
                     header.mutation_slot <= 1U &&
                     header.segment_count <= GEO_MANIFEST_MAX_SEGMENTS &&
                     header.segment_count <= SIZE_MAX;
    GeoSegmentSet *set = succeeded ? segment_set_allocate(manifest_path) : NULL;

    if (succeeded && !set) {
        succeeded = false;
    }

    if (succeeded && header.segment_count) {
        succeeded = segment_set_reserve(set, (size_t) header.segment_count);
    }

    for (size_t i = 0; succeeded && i < (size_t) header.segment_count; ++i) {
        GeoManifestEntry entry;

        succeeded = fread(&entry, sizeof(entry), 1, file) == 1 &&
                    entry.path_length > 0 &&
                    entry.path_length <= GEO_MANIFEST_MAX_PATH_BYTES &&
                    entry.flags == 0 &&
                    entry.generation > 0 &&
                    entry.generation <= header.generation;

        if (!succeeded) {
            break;
        }

        char *path = malloc((size_t) entry.path_length + 1);

        if (!path || fread(path, 1, entry.path_length, file) != entry.path_length) {
            free(path);
            succeeded = false;
            break;
        }

        path[entry.path_length] = '\0';

        if (strlen(path) != entry.path_length) {
            free(path);
            succeeded = false;
            break;
        }

        GeoIndex *segment = geo_index_open_mmap(path);

        if (!segment || segment->count != entry.record_count || segment->content_checksum != entry.content_checksum ||
            segment->count > UINT64_MAX - set->record_count) {
            geo_index_destroy(segment);
            free(path);
            succeeded = false;
            break;
        }

        set->paths[set->count] = path;
        set->segments[set->count] = segment;
        set->segment_generations[set->count] = entry.generation;
        set->count++;
        set->record_count += segment->count;
    }

    if (succeeded) {
        succeeded = fgetc(file) == EOF &&
                    !ferror(file) &&
                    manifest_checksum(set->paths,
                                      set->segments,
                                      set->segment_generations,
                                      set->count,
                                      header.generation,
                                      header.mutation_count,
                                      header.mutation_checksum,
                                      header.mutation_slot,
                                      header.durable_watermark) == header.checksum;
    }

    fclose(file);

    if (!succeeded) {
        geo_segment_set_destroy(set);

        return NULL;
    }

    set->generation = header.generation;
    set->mutation_count = header.mutation_count;
    set->mutation_checksum = header.mutation_checksum;
    set->mutation_slot = header.mutation_slot;
    set->durable_watermark = header.durable_watermark;

    if (set->mutation_slot) {
        char *mutation_path = segment_mutation_path(manifest_path, set->mutation_slot);

        if (!mutation_path) {
            geo_segment_set_destroy(set);

            return NULL;
        }

        free(set->mutation_path);
        set->mutation_path = mutation_path;
    }

    if (!mutation_log_load(set)) {
        geo_segment_set_destroy(set);

        return NULL;
    }

    return set;
#else
    (void) manifest_path;

    return NULL;
#endif
}

void geo_segment_set_destroy(GeoSegmentSet *set)
{
    if (!set) {
        return;
    }

#if GEO_SEGMENTS_SUPPORTED
    if (set->locks_initialized) {
        pthread_mutex_lock(&set->background_lock);
        set->background_enabled = false;
        set->background_requested = false;
        set->background_stop = true;
        pthread_cond_broadcast(&set->background_condition);
        pthread_mutex_unlock(&set->background_lock);

        if (set->background_thread_started) {
            (void) pthread_join(set->background_thread, NULL);
        }

        geo_thread_pool_destroy(set->compaction_pool);
    }
#endif

    for (size_t i = 0; i < set->count; ++i) {
        geo_index_destroy(set->segments[i]);
        free(set->paths[i]);
    }

    free(set->segments);
    free(set->paths);
    free(set->segment_generations);
    mutation_table_destroy(&set->mutations);
    free(set->manifest_path);
    free(set->mutation_path);
#if GEO_SEGMENTS_SUPPORTED
    if (set->locks_initialized) {
        free(set->background_directory);
        free(set->reader_shards);
        pthread_cond_destroy(&set->background_condition);
        pthread_mutex_destroy(&set->background_lock);
        pthread_mutex_destroy(&set->mutation_lock);
    }
#endif
    free(set);
}

static bool segment_set_publish_file(GeoSegmentSet *set,
                                     const char *segment_path,
                                     bool upsert_all,
                                     bool advance_watermark,
                                     uint64_t durable_watermark)
{
#if GEO_SEGMENTS_SUPPORTED
    if (!set || !segment_path || !segment_path[0] || pthread_mutex_lock(&set->mutation_lock) != 0) {
        return false;
    }

    if (set->generation >= GEO_MUTATION_MAX_GENERATION ||
        (advance_watermark && durable_watermark < set->durable_watermark)) {
        pthread_mutex_unlock(&set->mutation_lock);

        return false;
    }

    uint64_t published_watermark = advance_watermark ? durable_watermark : set->durable_watermark;

    char *canonical_path = realpath(segment_path, NULL);

    if (!canonical_path || strlen(canonical_path) > GEO_MANIFEST_MAX_PATH_BYTES) {
        free(canonical_path);
        pthread_mutex_unlock(&set->mutation_lock);

        return false;
    }

    for (size_t i = 0; i < set->count; ++i) {
        if (strcmp(set->paths[i], canonical_path) == 0) {
            free(canonical_path);
            pthread_mutex_unlock(&set->mutation_lock);

            return false;
        }
    }

    GeoIndex *segment = geo_index_open_mmap(canonical_path);

    if (!segment || segment->count > UINT64_MAX - set->record_count) {
        geo_index_destroy(segment);
        free(canonical_path);
        pthread_mutex_unlock(&set->mutation_lock);

        return false;
    }

    size_t capacity = set->capacity;

    if (set->count == capacity) {
        capacity = capacity ? capacity * 2 : 4;
    }

    bool capacity_valid = capacity >= set->count + 1 &&
                          capacity <= SIZE_MAX / sizeof(*set->paths) &&
                          capacity <= SIZE_MAX / sizeof(*set->segments) &&
                          capacity <= SIZE_MAX / sizeof(*set->segment_generations);
    char **paths = capacity_valid ? malloc(capacity * sizeof(*paths)) : NULL;
    GeoIndex **segments = capacity_valid ? malloc(capacity * sizeof(*segments)) : NULL;
    uint64_t *segment_generations = capacity_valid ? malloc(capacity * sizeof(*segment_generations)) : NULL;

    if (!paths || !segments || !segment_generations) {
        free(segment_generations);
        free(segments);
        free(paths);
        geo_index_destroy(segment);
        free(canonical_path);
        pthread_mutex_unlock(&set->mutation_lock);

        return false;
    }

    if (set->count) {
        memcpy(paths, set->paths, set->count * sizeof(*paths));
        memcpy(segments, set->segments, set->count * sizeof(*segments));
        memcpy(segment_generations,
               set->segment_generations,
               set->count * sizeof(*segment_generations));
    }

    uint64_t new_generation = set->generation + 1U;
    paths[set->count] = canonical_path;
    segments[set->count] = segment;
    segment_generations[set->count] = new_generation;
    size_t new_count = set->count + 1;

    size_t mutation_capacity = 0;
    size_t mutation_count = 0;
    GeoMutationRecord *mutations = NULL;

    for (size_t i = 0; (upsert_all || set->mutations.count) && i < segment->count; ++i) {
        const GeoMutationEntry *entry = mutation_table_lookup(&set->mutations, segment->records[i].id);

        if (!upsert_all && (!entry || mutation_entry_state(entry) == GEO_MUTATION_RESET)) {
            continue;
        }

        if (mutation_count == mutation_capacity) {
            size_t expanded_capacity = mutation_capacity ? mutation_capacity * 2U : 16U;
            bool expansion_valid = expanded_capacity >= mutation_capacity &&
                                   expanded_capacity <= SIZE_MAX / sizeof(*mutations);
            GeoMutationRecord *expanded = expansion_valid
                                              ? realloc(mutations, expanded_capacity * sizeof(*expanded))
                                              : NULL;

            if (!expanded) {
                free(mutations);
                free(segment_generations);
                free(segments);
                free(paths);
                geo_index_destroy(segment);
                free(canonical_path);
                pthread_mutex_unlock(&set->mutation_lock);

                return false;
            }

            mutations = expanded;
            mutation_capacity = expanded_capacity;
        }

        mutations[mutation_count++] = (GeoMutationRecord) {
            .id = segment->records[i].id,
            .generation = new_generation,
            .state = GEO_MUTATION_LIVE,
        };
    }

    GeoMutationTable replacement = { 0 };
    bool succeeded = mutation_count <= UINT64_MAX - set->mutation_count;

    if (succeeded && upsert_all) {
        succeeded = mutation_table_clone(&set->mutations, &replacement, mutation_count);
    }

    for (size_t i = 0; succeeded && upsert_all && i < mutation_count; ++i) {
        const GeoMutationEntry *entry = mutation_table_lookup(&replacement, mutations[i].id);

        if (entry && mutation_entry_generation(entry) == new_generation) {
            succeeded = false;
            break;
        }

        succeeded = mutation_table_apply(&replacement, mutations + i);
    }

    uint64_t updated_checksum = set->mutation_checksum;

    if (succeeded) {
        succeeded = mutation_log_append(set, mutations, mutation_count, &updated_checksum);
    }

    if (succeeded) {
        succeeded = manifest_write(set->manifest_path,
                                   paths,
                                   segments,
                                   segment_generations,
                                   new_count,
                                   new_generation,
                                   set->mutation_count + mutation_count,
                                   updated_checksum,
                                   set->mutation_slot,
                                   published_watermark);
    }

    if (succeeded) {
        succeeded = segment_set_write_lock(set);
    }

    if (succeeded) {
        if (upsert_all) {
            GeoMutationTable old_mutations = set->mutations;

            set->mutations = replacement;
            replacement = old_mutations;
        } else {
            for (size_t i = 0; i < mutation_count; ++i) {
                // These IDs were selected from this table while mutation_lock remained held, so apply cannot allocate or fail.
                (void) mutation_table_apply(&set->mutations, mutations + i);
            }
        }

        free(set->segments);
        free(set->paths);
        free(set->segment_generations);
        set->paths = paths;
        set->segments = segments;
        set->segment_generations = segment_generations;
        set->count = new_count;
        set->capacity = capacity;
        set->generation = new_generation;
        set->record_count += segment->count;
        set->mutation_count += mutation_count;
        set->mutation_checksum = updated_checksum;
        set->durable_watermark = published_watermark;
        segment_set_write_unlock(set);
    } else {
        free(segment_generations);
        free(segments);
        free(paths);
        geo_index_destroy(segment);
        free(canonical_path);
    }

    mutation_table_destroy(&replacement);
    free(mutations);

    if (succeeded) {
        mutation_checkpoint_if_amplified(set);
    }

    pthread_mutex_unlock(&set->mutation_lock);

    if (succeeded) {
        segment_schedule_background_compaction(set);
    }

    return succeeded;
#else
    (void) set;
    (void) segment_path;
    (void) upsert_all;
    (void) advance_watermark;
    (void) durable_watermark;

    return false;
#endif
}

bool geo_segment_set_add_file(GeoSegmentSet *set, const char *segment_path)
{
    return segment_set_publish_file(set, segment_path, false, false, 0U);
}

bool geo_segment_set_upsert_file(GeoSegmentSet *set, const char *segment_path)
{
    return segment_set_publish_file(set, segment_path, true, false, 0U);
}

bool geo_segment_set_upsert_file_at_watermark(GeoSegmentSet *set,
                                              const char *segment_path,
                                              uint64_t durable_watermark)
{
    return segment_set_publish_file(set, segment_path, true, true, durable_watermark);
}

size_t geo_segment_set_count(const GeoSegmentSet *set)
{
    if (!set || !segment_set_read_lock(set)) {
        return 0;
    }

    size_t count = set->count;

    segment_set_read_unlock(set);

    return count;
}

uint64_t geo_segment_set_record_count(const GeoSegmentSet *set)
{
    if (!set || !segment_set_read_lock(set)) {
        return 0;
    }

    uint64_t record_count = set->record_count;

    segment_set_read_unlock(set);

    return record_count;
}

uint64_t geo_segment_set_durable_watermark(const GeoSegmentSet *set)
{
#if GEO_SEGMENTS_SUPPORTED
    if (!set) {
        return 0U;
    }

    GeoSegmentSet *mutable_set = (GeoSegmentSet *) set;

    if (pthread_mutex_lock(&mutable_set->mutation_lock) != 0) {
        return 0U;
    }

    uint64_t durable_watermark = set->durable_watermark;

    pthread_mutex_unlock(&mutable_set->mutation_lock);
    return durable_watermark;
#else
    (void) set;

    return 0U;
#endif
}

static bool segment_set_remove_ids(GeoSegmentSet *set,
                                   const uint64_t *ids,
                                   size_t count,
                                   bool advance_watermark,
                                   uint64_t durable_watermark)
{
#if GEO_SEGMENTS_SUPPORTED
    if (!set || (!ids && count) || pthread_mutex_lock(&set->mutation_lock) != 0) {
        return false;
    }

    if (advance_watermark && durable_watermark < set->durable_watermark) {
        pthread_mutex_unlock(&set->mutation_lock);

        return false;
    }

    uint64_t published_watermark = advance_watermark ? durable_watermark : set->durable_watermark;

    if (!count) {
        bool succeeded = !advance_watermark || durable_watermark == set->durable_watermark ||
                         manifest_write(set->manifest_path,
                                        set->paths,
                                        set->segments,
                                        set->segment_generations,
                                        set->count,
                                        set->generation,
                                        set->mutation_count,
                                        set->mutation_checksum,
                                        set->mutation_slot,
                                        published_watermark);

        if (succeeded) {
            set->durable_watermark = published_watermark;
        }

        pthread_mutex_unlock(&set->mutation_lock);

        return succeeded;
    }

    if (set->generation >= GEO_MUTATION_MAX_GENERATION || count > SIZE_MAX / sizeof(GeoMutationRecord) ||
        count > UINT64_MAX - set->mutation_count) {
        pthread_mutex_unlock(&set->mutation_lock);

        return false;
    }

    uint64_t new_generation = set->generation + 1U;
    GeoMutationRecord *mutations = malloc(count * sizeof(*mutations));
    bool needs_rehash = count > SIZE_MAX - set->mutations.count ||
                        set->mutations.count + count > set->mutations.capacity - set->mutations.capacity / 4U;
    GeoMutationTable replacement = { 0 };
    bool succeeded = mutations != NULL;

    if (succeeded && needs_rehash) {
        succeeded = mutation_table_clone(&set->mutations, &replacement, count);
    }

    for (size_t i = 0; succeeded && i < count; ++i) {
        mutations[i] = (GeoMutationRecord) {
            .id = ids[i],
            .generation = new_generation,
            .state = GEO_MUTATION_DELETE,
        };
        if (needs_rehash) {
            succeeded = mutation_table_apply(&replacement, mutations + i);
        }
    }

    uint64_t updated_checksum = set->mutation_checksum;

    if (succeeded) {
        succeeded = mutation_log_append(set, mutations, count, &updated_checksum);
    }

    if (succeeded) {
        succeeded = manifest_write(set->manifest_path,
                                   set->paths,
                                   set->segments,
                                   set->segment_generations,
                                   set->count,
                                   new_generation,
                                   set->mutation_count + count,
                                   updated_checksum,
                                   set->mutation_slot,
                                   published_watermark);
    }

    if (succeeded) {
        succeeded = segment_set_write_lock(set);
    }

    if (succeeded) {
        if (needs_rehash) {
            GeoMutationTable old_mutations = set->mutations;

            set->mutations = replacement;
            replacement = old_mutations;
        } else {
            for (size_t i = 0; i < count; ++i) {
                succeeded = mutation_table_apply(&set->mutations, mutations + i);
            }
        }

        set->generation = new_generation;
        set->mutation_count += count;
        set->mutation_checksum = updated_checksum;
        set->durable_watermark = published_watermark;
        segment_set_write_unlock(set);
    }

    mutation_table_destroy(&replacement);
    free(mutations);

    if (succeeded) {
        mutation_checkpoint_if_amplified(set);
    }

    pthread_mutex_unlock(&set->mutation_lock);

    if (succeeded) {
        segment_schedule_background_compaction(set);
    }

    return succeeded;
#else
    (void) set;
    (void) ids;
    (void) count;
    (void) advance_watermark;
    (void) durable_watermark;

    return false;
#endif
}

bool geo_segment_set_remove_ids(GeoSegmentSet *set, const uint64_t *ids, size_t count)
{
    return segment_set_remove_ids(set, ids, count, false, 0U);
}

bool geo_segment_set_remove_ids_at_watermark(GeoSegmentSet *set,
                                             const uint64_t *ids,
                                             size_t count,
                                             uint64_t durable_watermark)
{
    return segment_set_remove_ids(set, ids, count, true, durable_watermark);
}

bool geo_segment_set_remove(GeoSegmentSet *set, uint64_t id)
{
    return geo_segment_set_remove_ids(set, &id, 1);
}

static bool mutation_log_write_checkpoint(const char *path,
                                          const GeoMutationTable *mutations,
                                          uint64_t *record_count,
                                          uint64_t *checksum)
{
    char *temporary_path = NULL;
    FILE *file = geo_io_create_atomic_file(path, &temporary_path);

    if (!file) {
        return false;
    }

    uint64_t written = 0;
    uint64_t value_checksum = UINT64_C(1469598103934665603);
    bool succeeded = true;

    for (size_t slot = 0; succeeded && slot < mutations->capacity; ++slot) {
        const GeoMutationEntry *entry = mutations->entries + slot;

        if (!mutation_entry_is_occupied(entry) || mutation_entry_state(entry) == GEO_MUTATION_RESET) {
            continue;
        }

        GeoMutationRecord record = {
            .id = entry->id,
            .generation = mutation_entry_generation(entry),
            .state = mutation_entry_state(entry),
        };

        succeeded = fwrite(&record, sizeof(record), 1, file) == 1;

        if (succeeded) {
            value_checksum = manifest_checksum_update(value_checksum, &record, sizeof(record));
            written++;
        }
    }

    if (succeeded) {
        succeeded = geo_io_publish_atomic_file(file, temporary_path, path);
    } else {
        geo_io_discard_atomic_file(file, temporary_path);
    }

    free(temporary_path);

    if (succeeded) {
        *record_count = written;
        *checksum = value_checksum;
    }

    return succeeded;
}

static bool mutation_checkpoint_locked(GeoSegmentSet *set)
{
    uint64_t checkpoint_slot = set->mutation_slot ^ 1U;
    char *checkpoint_path = segment_mutation_path(set->manifest_path, checkpoint_slot);
    uint64_t checkpoint_count = 0;
    uint64_t checkpoint_checksum = 0;
    bool succeeded = checkpoint_path &&
                     mutation_log_write_checkpoint(checkpoint_path,
                                                   &set->mutations,
                                                   &checkpoint_count,
                                                   &checkpoint_checksum) &&
                     manifest_write(set->manifest_path,
                                    set->paths,
                                    set->segments,
                                    set->segment_generations,
                                    set->count,
                                    set->generation,
                                    checkpoint_count,
                                    checkpoint_checksum,
                                    checkpoint_slot,
                                    set->durable_watermark);

    if (succeeded) {
        char *obsolete_path = set->mutation_path;

        set->mutation_path = checkpoint_path;
        set->mutation_slot = checkpoint_slot;
        set->mutation_count = checkpoint_count;
        set->mutation_checksum = checkpoint_checksum;
        checkpoint_path = NULL;
        (void) unlink(obsolete_path);
        (void) geo_io_sync_parent_directory(obsolete_path);
        free(obsolete_path);
    } else if (checkpoint_path) {
        (void) unlink(checkpoint_path);
    }

    free(checkpoint_path);

    return succeeded;
}

static void mutation_checkpoint_if_amplified(GeoSegmentSet *set)
{
    uint64_t active_count = set->mutations.count;
    uint64_t amplified_limit = active_count > UINT64_MAX / GEO_MUTATION_CHECKPOINT_AMPLIFICATION
                                   ? UINT64_MAX
                                   : active_count * GEO_MUTATION_CHECKPOINT_AMPLIFICATION;

    if (set->mutation_count >= GEO_MUTATION_CHECKPOINT_MIN_RECORDS && set->mutation_count > amplified_limit) {
        (void) mutation_checkpoint_locked(set);
    }
}

bool geo_segment_set_checkpoint_mutations(GeoSegmentSet *set)
{
#if GEO_SEGMENTS_SUPPORTED
    if (!set || pthread_mutex_lock(&set->mutation_lock) != 0) {
        return false;
    }

    bool succeeded = mutation_checkpoint_locked(set);

    pthread_mutex_unlock(&set->mutation_lock);

    return succeeded;
#else
    (void) set;

    return false;
#endif
}

typedef struct {
    const GeoMutationTable *mutations;
    uint64_t segment_generation;
    GeoMutationEntry small_entries[GEO_MUTATION_SMALL_ENTRY_CAPACITY];
    size_t small_count;
    bool uses_small_entries;
} GeoSegmentVisibility;

static void segment_visibility_initialize(GeoSegmentVisibility *visibility,
                                          const GeoMutationTable *mutations,
                                          uint64_t segment_generation)
{
    memset(visibility, 0, sizeof(*visibility));
    visibility->mutations = mutations;
    visibility->segment_generation = segment_generation;

    if (mutations->count > sizeof(visibility->small_entries) / sizeof(visibility->small_entries[0])) {
        return;
    }

    visibility->uses_small_entries = true;

    for (size_t slot = 0; slot < mutations->capacity; ++slot) {
        if (mutation_entry_is_occupied(mutations->entries + slot)) {
            visibility->small_entries[visibility->small_count++] = mutations->entries[slot];
        }
    }
}

static bool segment_mutation_entry_is_visible(const GeoMutationEntry *entry, uint64_t segment_generation)
{
    return !entry || mutation_entry_state(entry) == GEO_MUTATION_RESET ||
           (mutation_entry_state(entry) == GEO_MUTATION_LIVE && mutation_entry_generation(entry) == segment_generation);
}

static bool segment_record_is_visible_for(const GeoRecord *record, const GeoSegmentVisibility *visibility)
{
    const GeoMutationEntry *entry = NULL;

    if (visibility->uses_small_entries) {
        for (size_t i = 0; i < visibility->small_count; ++i) {
            if (visibility->small_entries[i].id == record->id) {
                entry = visibility->small_entries + i;
                break;
            }
        }
    } else {
        entry = mutation_table_lookup(visibility->mutations, record->id);
    }

    return segment_mutation_entry_is_visible(entry, visibility->segment_generation);
}

static bool segment_record_is_visible(const GeoRecord *record, void *context)
{
    return segment_record_is_visible_for(record, context);
}

static size_t segment_scan_visible_records(const GeoIndex *segment,
                                           const GeoSegmentVisibility *visibility,
                                           uint64_t *visibility_bits)
{
    size_t word_count = geo_internal_bit_word_count(segment->count);
    size_t visible_count = 0;

    for (size_t word = 0; word < word_count; ++word) {
        size_t first_record = word * 64U;
        size_t records_in_word = segment->count - first_record;

        if (records_in_word > 64U) {
            records_in_word = 64U;
        }

        uint64_t visible_word = 0;

        for (size_t lane = 0; lane < records_in_word; ++lane) {
            bool visible = segment_record_is_visible_for(segment->records + first_record + lane, visibility);

            visible_word |= (uint64_t) visible << lane;
        }

        if (visibility_bits) {
            visibility_bits[word] = visible_word;
        }

        visible_count += (size_t) __builtin_popcountll(visible_word);
    }

    return visible_count;
}

static size_t segment_filter_mutation_table_bits(const GeoRecord *records,
                                                 size_t count,
                                                 uint64_t *candidate_bits,
                                                 const GeoSegmentVisibility *visibility)
{
    const GeoMutationTable *mutations = visibility->mutations;
    size_t word_count = geo_internal_bit_word_count(count);
    size_t matched = 0;

    for (size_t word = 0; word < word_count; ++word) {
        uint64_t candidates = candidate_bits[word];

        while (candidates) {
            unsigned lane = (unsigned) __builtin_ctzll(candidates);
            size_t record_index = word * 64U + lane;
            const GeoMutationEntry *entry = mutation_table_lookup(mutations, records[record_index].id);

            if (!segment_mutation_entry_is_visible(entry, visibility->segment_generation)) {
                candidate_bits[word] &= ~(UINT64_C(1) << lane);
            }

            candidates &= candidates - 1U;
        }

        matched += (size_t) __builtin_popcountll(candidate_bits[word]);
    }

    return matched;
}

static size_t segment_filter_visible_bits(const GeoRecord *records,
                                          size_t count,
                                          uint64_t *candidate_bits,
                                          void *context)
{
    const GeoSegmentVisibility *visibility = context;
    size_t word_count = geo_internal_bit_word_count(count);
    size_t matched = 0;

    if (visibility->uses_small_entries && visibility->small_count == 1U) {
        const GeoMutationEntry *entry = visibility->small_entries;

        if (!segment_mutation_entry_is_visible(entry, visibility->segment_generation)) {
            return geo_simd_exclude_id_bits(records, count, entry->id, candidate_bits);
        }

        for (size_t word = 0; word < word_count; ++word) {
            matched += (size_t) __builtin_popcountll(candidate_bits[word]);
        }

        return matched;
    }

    if (!visibility->uses_small_entries) {
        return segment_filter_mutation_table_bits(records, count, candidate_bits, visibility);
    }

    for (size_t word = 0; word < word_count; ++word) {
        uint64_t candidates = candidate_bits[word];

        while (candidates) {
            unsigned lane = (unsigned) __builtin_ctzll(candidates);
            size_t record_index = word * 64U + lane;
            const GeoMutationEntry *entry = NULL;

            for (size_t i = 0; i < visibility->small_count; ++i) {
                if (visibility->small_entries[i].id == records[record_index].id) {
                    entry = visibility->small_entries + i;
                    break;
                }
            }

            if (!segment_mutation_entry_is_visible(entry, visibility->segment_generation)) {
                candidate_bits[word] &= ~(UINT64_C(1) << lane);
            }

            candidates &= candidates - 1U;
        }

        matched += (size_t) __builtin_popcountll(candidate_bits[word]);
    }

    return matched;
}

typedef struct {
    GeoSegmentVisibility visibility;
    GeoRecordFilter external_filter;
    void *external_context;
    bool has_mutations;
} GeoCombinedRecordFilter;

static bool segment_combined_record_is_visible(const GeoRecord *record, void *context)
{
    const GeoCombinedRecordFilter *combined = context;

    return (!combined->has_mutations || segment_record_is_visible_for(record, &combined->visibility)) &&
           combined->external_filter(record, combined->external_context);
}

static size_t segment_filter_combined_bits(const GeoRecord *records,
                                           size_t count,
                                           uint64_t *candidate_bits,
                                           void *context)
{
    GeoCombinedRecordFilter *combined = context;

    if (combined->has_mutations) {
        (void) segment_filter_visible_bits(records, count, candidate_bits, &combined->visibility);
    }

    size_t matched = 0U;
    size_t word_count = geo_internal_bit_word_count(count);

    for (size_t word = 0U; word < word_count; ++word) {
        uint64_t candidates = candidate_bits[word];

        while (candidates) {
            unsigned lane = (unsigned) __builtin_ctzll(candidates);
            size_t record_index = word * 64U + lane;

            if (!combined->external_filter(records + record_index, combined->external_context)) {
                candidate_bits[word] &= ~(UINT64_C(1) << lane);
            }

            candidates &= candidates - 1U;
        }

        matched += (size_t) __builtin_popcountll(candidate_bits[word]);
    }

    return matched;
}

bool geo_segment_set_copy_path(const GeoSegmentSet *set,
                               size_t segment_index,
                               char *buffer,
                               size_t buffer_size)
{
    if (!set || !buffer || !buffer_size || !segment_set_read_lock(set)) {
        return false;
    }

    bool succeeded = false;

    if (segment_index < set->count) {
        size_t path_size = strlen(set->paths[segment_index]) + 1;

        if (path_size <= buffer_size) {
            memcpy(buffer, set->paths[segment_index], path_size);
            succeeded = true;
        }
    }

    segment_set_read_unlock(set);

    return succeeded;
}

// =============================================================================
// Automatic size-tiered background compaction
// =============================================================================

#if GEO_SEGMENTS_SUPPORTED
static char *segment_create_background_output_path(const GeoSegmentSet *set, const char *directory)
{
    if (!segment_set_read_lock(set)) {
        return NULL;
    }

    uint64_t generation = set->generation;

    segment_set_read_unlock(set);

    const char *separator = directory[strlen(directory) - 1] == '/' ? "" : "/";

    for (size_t attempt = 0; attempt < 16; ++attempt) {
        uint64_t sequence = atomic_fetch_add_explicit(&GEO_COMPACTION_SEQUENCE, 1, memory_order_relaxed) + 1;
        int path_length = snprintf(NULL,
                                   0,
                                   "%s%sgeobolt-compacted-%ld-%" PRIu64 "-%" PRIu64 ".geobolt",
                                   directory,
                                   separator,
                                   (long) getpid(),
                                   generation,
                                   sequence);

        if (path_length < 0) {
            return NULL;
        }

        char *path = malloc((size_t) path_length + 1);

        if (!path) {
            return NULL;
        }

        snprintf(path,
                 (size_t) path_length + 1,
                 "%s%sgeobolt-compacted-%ld-%" PRIu64 "-%" PRIu64 ".geobolt",
                 directory,
                 separator,
                 (long) getpid(),
                 generation,
                 sequence);

        if (access(path, F_OK) != 0) {
            return path;
        }

        free(path);
    }

    return NULL;
}

static uint64_t segment_mutation_pressure_threshold(uint64_t record_count, const GeoSegmentCompactionPolicy *policy)
{
    uint64_t denominator = policy->mutation_ratio_denominator;
    uint64_t numerator = policy->mutation_ratio_numerator;
    uint64_t quotient = record_count / denominator;
    uint64_t remainder_product = (record_count % denominator) * numerator;
    uint64_t threshold = quotient * numerator + remainder_product / denominator;

    return threshold + (remainder_product % denominator != 0);
}

static GeoBackgroundCompactionMode segment_background_compaction_mode(const GeoSegmentSet *set,
                                                                      const GeoSegmentCompactionPolicy *policy)
{
    if (!segment_set_read_lock(set)) {
        return GEO_BACKGROUND_COMPACTION_NONE;
    }

    uint64_t mutation_threshold = segment_mutation_pressure_threshold(set->record_count, policy);

    bool mutation_pressure = set->count &&
                             set->record_count &&
                             set->mutations.count >= policy->minimum_mutations &&
                             mutation_threshold <= SIZE_MAX &&
                             set->mutations.count >= (size_t) mutation_threshold;
    bool size_pressure = set->count > policy->max_active_segments;

    segment_set_read_unlock(set);

    if (mutation_pressure) {
        return GEO_BACKGROUND_COMPACTION_MUTATION_REWRITE;
    }

    return size_pressure ? GEO_BACKGROUND_COMPACTION_SIZE_TIER : GEO_BACKGROUND_COMPACTION_NONE;
}

static void *segment_background_compaction_worker(void *argument)
{
    GeoSegmentSet *set = argument;

    pthread_mutex_lock(&set->background_lock);

    while (!set->background_stop) {
        while (!set->background_requested && !set->background_stop) {
            pthread_cond_wait(&set->background_condition, &set->background_lock);
        }

        if (set->background_stop) {
            break;
        }

        set->background_requested = false;
        set->background_running = true;
        memset(&set->background_last_stats, 0, sizeof(set->background_last_stats));
        set->background_last_succeeded = false;

        GeoSegmentCompactionStats total_stats = { 0 };
        size_t mutation_rewrite_passes = 0;
        bool succeeded = true;
        GeoBackgroundCompactionMode mode = set->background_enabled
                                               ? segment_background_compaction_mode(set, &set->background_policy)
                                               : GEO_BACKGROUND_COMPACTION_NONE;
        char *output_path = mode != GEO_BACKGROUND_COMPACTION_NONE
                                ? segment_create_background_output_path(set, set->background_directory)
                                : NULL;

        if (mode != GEO_BACKGROUND_COMPACTION_NONE && !output_path) {
            succeeded = false;
        }

        pthread_mutex_unlock(&set->background_lock);

        while (output_path) {
            GeoSegmentCompactionStats partial_stats;

            if (mode == GEO_BACKGROUND_COMPACTION_MUTATION_REWRITE) {
                succeeded = geo_segment_set_compact(set, output_path, &partial_stats);
                mutation_rewrite_passes++;
            } else {
                succeeded = segment_set_compact_size_tier(set, output_path, &partial_stats);
            }

            free(output_path);
            output_path = NULL;

            if (!succeeded ||
                partial_stats.records_written > UINT64_MAX - total_stats.records_written ||
                partial_stats.input_segments > SIZE_MAX - total_stats.input_segments ||
                partial_stats.partition_count > SIZE_MAX - total_stats.partition_count) {
                succeeded = false;
                break;
            }

            total_stats.records_written += partial_stats.records_written;
            total_stats.input_segments += partial_stats.input_segments;
            total_stats.partition_count += partial_stats.partition_count;
            total_stats.merge_time_ms += partial_stats.merge_time_ms;

            if (partial_stats.worker_count > total_stats.worker_count) {
                total_stats.worker_count = partial_stats.worker_count;
            }

            pthread_mutex_lock(&set->background_lock);

            mode = set->background_enabled
                       ? segment_background_compaction_mode(set, &set->background_policy)
                       : GEO_BACKGROUND_COMPACTION_NONE;

            if (mode == GEO_BACKGROUND_COMPACTION_MUTATION_REWRITE &&
                mutation_rewrite_passes >= set->background_policy.maximum_mutation_rewrite_passes) {
                mode = GEO_BACKGROUND_COMPACTION_NONE;
            }

            if (mode != GEO_BACKGROUND_COMPACTION_NONE) {
                output_path = segment_create_background_output_path(set, set->background_directory);

                if (!output_path) {
                    succeeded = false;
                }
            }

            pthread_mutex_unlock(&set->background_lock);
        }

        pthread_mutex_lock(&set->background_lock);
        set->background_last_stats = total_stats;
        set->background_last_succeeded = succeeded;

        if (total_stats.records_written > UINT64_MAX - set->background_total_stats.records_written) {
            set->background_total_stats.records_written = UINT64_MAX;
        } else {
            set->background_total_stats.records_written += total_stats.records_written;
        }

        if (total_stats.input_segments > SIZE_MAX - set->background_total_stats.input_segments) {
            set->background_total_stats.input_segments = SIZE_MAX;
        } else {
            set->background_total_stats.input_segments += total_stats.input_segments;
        }

        if (total_stats.partition_count > SIZE_MAX - set->background_total_stats.partition_count) {
            set->background_total_stats.partition_count = SIZE_MAX;
        } else {
            set->background_total_stats.partition_count += total_stats.partition_count;
        }

        if (total_stats.worker_count > set->background_total_stats.worker_count) {
            set->background_total_stats.worker_count = total_stats.worker_count;
        }

        set->background_total_stats.merge_time_ms += total_stats.merge_time_ms;

        if (set->background_completed_runs < UINT64_MAX) {
            set->background_completed_runs++;
        }

        if (!succeeded && set->background_failed_runs < UINT64_MAX) {
            set->background_failed_runs++;
        }

        set->background_running = false;
        pthread_cond_broadcast(&set->background_condition);
    }

    set->background_running = false;
    pthread_cond_broadcast(&set->background_condition);
    pthread_mutex_unlock(&set->background_lock);

    return NULL;
}

static void segment_schedule_background_compaction(GeoSegmentSet *set)
{
    if (pthread_mutex_lock(&set->background_lock) != 0) {
        return;
    }

    GeoBackgroundCompactionMode mode = set->background_enabled && !set->background_requested && !set->background_running
                                           ? segment_background_compaction_mode(set, &set->background_policy)
                                           : GEO_BACKGROUND_COMPACTION_NONE;

    if (mode != GEO_BACKGROUND_COMPACTION_NONE) {
        set->background_requested = true;
        pthread_cond_signal(&set->background_condition);
    }

    pthread_mutex_unlock(&set->background_lock);
}
#endif

bool geo_segment_set_configure_background_compaction(GeoSegmentSet *set,
                                                     const char *output_directory,
                                                     const GeoSegmentCompactionPolicy *policy)
{
#if GEO_SEGMENTS_SUPPORTED
    if (!set ||
        !output_directory ||
        !output_directory[0] ||
        !policy ||
        !policy->max_active_segments ||
        !policy->minimum_mutations ||
        !policy->maximum_mutation_rewrite_passes ||
        !policy->mutation_ratio_numerator ||
        !policy->mutation_ratio_denominator ||
        policy->mutation_ratio_numerator > policy->mutation_ratio_denominator) {
        return false;
    }

    char *canonical_directory = realpath(output_directory, NULL);
    struct stat directory_status;

    if (!canonical_directory ||
        stat(canonical_directory, &directory_status) != 0 ||
        !S_ISDIR(directory_status.st_mode)) {
        free(canonical_directory);

        return false;
    }

    if (pthread_mutex_lock(&set->background_lock) != 0) {
        free(canonical_directory);

        return false;
    }

    if (set->background_running || set->background_requested) {
        pthread_mutex_unlock(&set->background_lock);
        free(canonical_directory);

        return false;
    }

    free(set->background_directory);
    set->background_directory = canonical_directory;
    set->background_policy = *policy;
    set->background_enabled = true;

    if (!set->background_thread_started) {
        if (pthread_create(&set->background_thread, NULL, segment_background_compaction_worker, set) != 0) {
            set->background_enabled = false;
            pthread_mutex_unlock(&set->background_lock);

            return false;
        }

        set->background_thread_started = true;
    }

    pthread_mutex_unlock(&set->background_lock);
    segment_schedule_background_compaction(set);

    return true;
#else
    (void) set;
    (void) output_directory;
    (void) policy;

    return false;
#endif
}

bool geo_segment_set_enable_background_compaction(GeoSegmentSet *set,
                                                  const char *output_directory,
                                                  size_t max_active_segments)
{
    GeoSegmentCompactionPolicy policy = {
        .max_active_segments = max_active_segments,
        .minimum_mutations = GEO_BACKGROUND_DEFAULT_MUTATION_MIN_ENTRIES,
        .maximum_mutation_rewrite_passes = GEO_BACKGROUND_DEFAULT_MUTATION_MAX_PASSES,
        .mutation_ratio_numerator = GEO_BACKGROUND_DEFAULT_MUTATION_RATIO_NUMERATOR,
        .mutation_ratio_denominator = GEO_BACKGROUND_DEFAULT_MUTATION_RATIO_DENOMINATOR,
    };

    return geo_segment_set_configure_background_compaction(set, output_directory, &policy);
}

bool geo_segment_set_wait_for_background_compaction(GeoSegmentSet *set,
                                                    GeoSegmentCompactionStats *stats)
{
    if (stats) {
        memset(stats, 0, sizeof(*stats));
    }

#if GEO_SEGMENTS_SUPPORTED
    if (!set || pthread_mutex_lock(&set->background_lock) != 0) {
        return false;
    }

    while (set->background_running || set->background_requested) {
        if (pthread_cond_wait(&set->background_condition, &set->background_lock) != 0) {
            pthread_mutex_unlock(&set->background_lock);

            return false;
        }
    }

    bool succeeded = !set->background_completed_runs || set->background_last_succeeded;

    if (stats) {
        *stats = set->background_last_stats;
    }

    pthread_mutex_unlock(&set->background_lock);

    return succeeded;
#else
    (void) set;

    return false;
#endif
}

bool geo_segment_set_background_compaction_totals(const GeoSegmentSet *set,
                                                  GeoSegmentCompactionStats *stats,
                                                  uint64_t *completed_runs,
                                                  uint64_t *failed_runs)
{
#if GEO_SEGMENTS_SUPPORTED
    if (!set || !stats || !completed_runs || !failed_runs) {
        return false;
    }

    GeoSegmentSet *mutable_set = (GeoSegmentSet *) set;

    if (pthread_mutex_lock(&mutable_set->background_lock) != 0) {
        return false;
    }

    *stats = mutable_set->background_total_stats;
    *completed_runs = mutable_set->background_completed_runs;
    *failed_runs = mutable_set->background_failed_runs;
    pthread_mutex_unlock(&mutable_set->background_lock);

    return true;
#else
    (void) set;
    (void) stats;
    (void) completed_runs;
    (void) failed_runs;

    return false;
#endif
}

bool geo_segment_set_background_compaction_active(const GeoSegmentSet *set)
{
#if GEO_SEGMENTS_SUPPORTED
    if (!set) {
        return false;
    }

    GeoSegmentSet *mutable_set = (GeoSegmentSet *) set;

    if (pthread_mutex_lock(&mutable_set->background_lock) != 0) {
        return false;
    }

    bool active = mutable_set->background_running;

    pthread_mutex_unlock(&mutable_set->background_lock);

    return active;
#else
    (void) set;

    return false;
#endif
}

bool geo_segment_set_compaction_active(const GeoSegmentSet *set)
{
#if GEO_SEGMENTS_SUPPORTED
    if (!set) {
        return false;
    }

    GeoSegmentSet *mutable_set = (GeoSegmentSet *) set;

    if (pthread_mutex_lock(&mutable_set->mutation_lock) != 0) {
        return false;
    }

    bool active = mutable_set->compaction_running;

    pthread_mutex_unlock(&mutable_set->mutation_lock);

    return active;
#else
    (void) set;

    return false;
#endif
}

// =============================================================================
// Queries across active immutable segments
// =============================================================================

static void segment_stats_reset(GeoSearchStats *stats)
{
    if (stats) {
        memset(stats, 0, sizeof(*stats));
    }
}

static void segment_stats_add(GeoSearchStats *total, const GeoSearchStats *partial)
{
    if (!total) {
        return;
    }

    total->records_scanned += partial->records_scanned;
    total->records_matched += partial->records_matched;
    total->ranges_checked += partial->ranges_checked;
}

static unsigned segment_radius_density_refinement(const GeoSegmentSet *set, double latitude, double longitude)
{
    unsigned maximum_refinement = 0;

    for (size_t i = 0; i < set->count; ++i) {
        unsigned refinement = geo_index_density_refinement(set->segments[i], latitude, longitude);

        if (refinement > maximum_refinement) {
            maximum_refinement = refinement;
        }
    }

    return maximum_refinement;
}

static uint64_t segment_saturating_multiply(uint64_t value, uint64_t multiplier)
{
    return value && multiplier > UINT64_MAX / value ? UINT64_MAX : value * multiplier;
}

static uint64_t segment_saturating_add(uint64_t first, uint64_t second)
{
    return second > UINT64_MAX - first ? UINT64_MAX : first + second;
}

static bool segment_prepare_radius_query(const GeoSegmentSet *set,
                                         double latitude,
                                         double longitude,
                                         double radius_km,
                                         size_t output_record_size,
                                         GeoRadiusQueryPlan *plan)
{
    if (set->count == 1U) {
        return geo_index_prepare_radius_query_for_index(set->segments[0],
                                                        latitude,
                                                        longitude,
                                                        radius_km,
                                                        output_record_size,
                                                        plan);
    }

    unsigned maximum_refinement = segment_radius_density_refinement(set, latitude, longitude);
    bool selected = false;

    for (unsigned refinement = 0; refinement <= maximum_refinement; ++refinement) {
        GeoRadiusQueryPlan candidate;

        if (!geo_index_prepare_radius_query(latitude,
                                            longitude,
                                            radius_km,
                                            set->record_count,
                                            refinement,
                                            &candidate)) {
            continue;
        }

        uint64_t candidate_records = 0;

        for (size_t segment = 0; segment < set->count; ++segment) {
            uint64_t partial = geo_index_radius_plan_candidate_count(set->segments[segment], &candidate);

            candidate_records = segment_saturating_add(candidate_records, partial);
        }

        candidate.candidate_records = candidate_records;
        uint64_t output_bytes = segment_saturating_multiply(candidate_records, output_record_size);
        uint64_t range_checks = segment_saturating_multiply((uint64_t) candidate.range_count, set->count);

        candidate.estimated_output_bytes = output_bytes - output_bytes / 4U;
        candidate.estimated_cost = segment_saturating_add(segment_saturating_multiply(candidate_records, 10U),
                                                          segment_saturating_multiply(range_checks, 96U));
        candidate.estimated_cost = segment_saturating_add(candidate.estimated_cost,
                                                          candidate.estimated_output_bytes / 16U);

        if (!selected || candidate.estimated_cost < plan->estimated_cost) {
            *plan = candidate;
            selected = true;
        }

        if (candidate.strategy == GEO_QUERY_STRATEGY_COPY) {
            break;
        }
    }

    return selected;
}

static bool segment_set_search_radius_plan_filtered_snapshot(const GeoSegmentSet *set,
                                                             const GeoRadiusQueryPlan *plan,
                                                             GeoSearchResult *result,
                                                             size_t *count,
                                                             GeoRecordFilter external_filter,
                                                             void *external_context,
                                                             GeoSearchStats *stats)
{
    segment_stats_reset(stats);

    if (!set || !plan || plan->range_count <= 0 || (!result && !count) || (result && count)) {
        return false;
    }

    double start = stats ? geo_get_time_ms() : 0.0;
    size_t total_count = 0U;
    bool succeeded = true;

    if (result) {
        geo_result_clear(result);
    }

    for (size_t i = 0; succeeded && i < set->count; ++i) {
        GeoSearchStats partial_stats;
        size_t partial_count = 0U;
        GeoSegmentVisibility visibility;
        GeoCombinedRecordFilter combined;
        GeoRecordBitFilter bit_filter;
        GeoRecordFilter record_filter;
        void *filter_context;

        segment_visibility_initialize(&visibility, &set->mutations, set->segment_generations[i]);

        if (external_filter) {
            combined = (GeoCombinedRecordFilter) {
                .visibility = visibility,
                .external_filter = external_filter,
                .external_context = external_context,
                .has_mutations = set->mutations.count != 0U,
            };
            bit_filter = segment_filter_combined_bits;
            record_filter = segment_combined_record_is_visible;
            filter_context = &combined;
        } else {
            bit_filter = set->mutations.count ? segment_filter_visible_bits : NULL;
            record_filter = set->mutations.count ? segment_record_is_visible : NULL;
            filter_context = &visibility;
        }

        succeeded = geo_index_search_radius_plan_append_filtered(set->segments[i],
                                                                 plan,
                                                                 result,
                                                                 &partial_count,
                                                                 &partial_stats,
                                                                 bit_filter,
                                                                 record_filter,
                                                                 filter_context) &&
                    partial_count <= SIZE_MAX - total_count;

        if (succeeded) {
            total_count += partial_count;
            segment_stats_add(stats, &partial_stats);
        }
    }

    if (succeeded && count) {
        *count = total_count;
    } else if (!succeeded && result) {
        geo_result_clear(result);
    }

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start;
    }

    return succeeded;
}

static bool segment_set_search_radius_filtered_snapshot(const GeoSegmentSet *set,
                                                        double latitude,
                                                        double longitude,
                                                        double radius_km,
                                                        GeoSearchResult *result,
                                                        size_t *count,
                                                        GeoRecordFilter external_filter,
                                                        void *external_context,
                                                        GeoSearchStats *stats)
{
    GeoRadiusQueryPlan plan;
    size_t output_record_size = result ? sizeof(GeoRecord) : 0U;

    if (!segment_prepare_radius_query(set, latitude, longitude, radius_km, output_record_size, &plan)) {
        return false;
    }

    return segment_set_search_radius_plan_filtered_snapshot(set,
                                                            &plan,
                                                            result,
                                                            count,
                                                            external_filter,
                                                            external_context,
                                                            stats);
}

bool geo_segment_set_search_radius_filtered(const GeoSegmentSet *set,
                                            double latitude,
                                            double longitude,
                                            double radius_km,
                                            GeoSearchResult *result,
                                            size_t *count,
                                            GeoRecordFilter filter,
                                            void *filter_context,
                                            GeoSearchStats *stats)
{
    if (!set || !filter || !segment_set_read_lock(set)) {
        return false;
    }

    bool succeeded = segment_set_search_radius_filtered_snapshot(set,
                                                                 latitude,
                                                                 longitude,
                                                                 radius_km,
                                                                 result,
                                                                 count,
                                                                 filter,
                                                                 filter_context,
                                                                 stats);

    segment_set_read_unlock(set);
    return succeeded;
}

bool geo_segment_set_prepare_radius_query(const GeoSegmentSet *set,
                                          double latitude,
                                          double longitude,
                                          double radius_km,
                                          size_t output_record_size,
                                          GeoRadiusQueryPlan *plan)
{
    if (!set || !plan || !segment_set_read_lock(set)) {
        return false;
    }

    bool succeeded = segment_prepare_radius_query(set,
                                                  latitude,
                                                  longitude,
                                                  radius_km,
                                                  output_record_size,
                                                  plan);

    segment_set_read_unlock(set);
    return succeeded;
}

bool geo_segment_set_search_radius_plan_reuse(const GeoSegmentSet *set,
                                              const GeoRadiusQueryPlan *plan,
                                              GeoSearchResult *result,
                                              GeoSearchStats *stats)
{
    if (!set || !plan || !result || !segment_set_read_lock(set)) {
        return false;
    }

    bool succeeded = segment_set_search_radius_plan_filtered_snapshot(set,
                                                                      plan,
                                                                      result,
                                                                      NULL,
                                                                      NULL,
                                                                      NULL,
                                                                      stats);

    segment_set_read_unlock(set);
    return succeeded;
}

bool geo_segment_set_search_radius_reuse(const GeoSegmentSet *set,
                                         double lat,
                                         double lng,
                                         double radius_km,
                                         GeoSearchResult *result,
                                         GeoSearchStats *stats)
{
    if (!set || !segment_set_read_lock(set)) {
        return false;
    }

    bool succeeded = segment_set_search_radius_filtered_snapshot(set,
                                                                 lat,
                                                                 lng,
                                                                 radius_km,
                                                                 result,
                                                                 NULL,
                                                                 NULL,
                                                                 NULL,
                                                                 stats);

    segment_set_read_unlock(set);
    return succeeded;
}

GeoSearchResult *geo_segment_set_search_radius(const GeoSegmentSet *set,
                                               double lat,
                                               double lng,
                                               double radius_km,
                                               GeoSearchStats *stats)
{
    GeoSearchResult *result = geo_result_create(64);

    if (!result) {
        return NULL;
    }

    if (!geo_segment_set_search_radius_reuse(set, lat, lng, radius_km, result, stats)) {
        geo_result_destroy(result);

        return NULL;
    }

    return result;
}

bool geo_segment_set_search_radius_count_snapshot(const GeoSegmentSet *set,
                                                  double latitude,
                                                  double longitude,
                                                  double radius_km,
                                                  size_t *count)
{
    return segment_set_search_radius_filtered_snapshot(set,
                                                       latitude,
                                                       longitude,
                                                       radius_km,
                                                       NULL,
                                                       count,
                                                       NULL,
                                                       NULL,
                                                       NULL);
}

bool geo_segment_set_search_radius_count(const GeoSegmentSet *set,
                                         double lat,
                                         double lng,
                                         double radius_km,
                                         size_t *count,
                                         GeoSearchStats *stats)
{
    if (!set || !segment_set_read_lock(set)) {
        return false;
    }

    bool succeeded = segment_set_search_radius_filtered_snapshot(set,
                                                                 lat,
                                                                 lng,
                                                                 radius_km,
                                                                 NULL,
                                                                 count,
                                                                 NULL,
                                                                 NULL,
                                                                 stats);

    segment_set_read_unlock(set);

    return succeeded;
}

static bool segment_set_search_radius_ids_snapshot(const GeoSegmentSet *set,
                                                   double latitude,
                                                   double longitude,
                                                   double radius_km,
                                                   GeoIdResult *result,
                                                   bool allow_growth,
                                                   GeoSearchStats *stats)
{
    segment_stats_reset(stats);

    if (!set || !result || (result->capacity && !result->ids)) {
        return false;
    }

    GeoRadiusQueryPlan plan;

    if (!segment_prepare_radius_query(set, latitude, longitude, radius_km, sizeof(uint64_t), &plan)) {
        return false;
    }

    double start = stats ? geo_get_time_ms() : 0.0;
    result->count = 0;
    size_t total_count = 0;
    bool succeeded = true;

    for (size_t i = 0; succeeded && i < set->count; ++i) {
        GeoSearchStats partial_stats;
        size_t partial_count = 0;
        GeoSegmentVisibility visibility;

        segment_visibility_initialize(&visibility, &set->mutations, set->segment_generations[i]);

        succeeded = geo_index_search_radius_plan_ids_append_filtered(
                        set->segments[i],
                        &plan,
                        result,
                        allow_growth,
                        &partial_count,
                        &partial_stats,
                        set->mutations.count ? segment_filter_visible_bits : NULL,
                        set->mutations.count ? segment_record_is_visible : NULL,
                        &visibility) &&
                    partial_count <= SIZE_MAX - total_count;

        if (succeeded) {
            total_count += partial_count;
            segment_stats_add(stats, &partial_stats);
        }
    }

    if (!succeeded || result->count != total_count) {
        result->count = 0;
        succeeded = false;
    }

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start;
    }

    return succeeded;
}

bool geo_segment_set_search_radius_ids_reuse(const GeoSegmentSet *set,
                                             double lat,
                                             double lng,
                                             double radius_km,
                                             GeoIdResult *result,
                                             GeoSearchStats *stats)
{
    if (!set || !segment_set_read_lock(set)) {
        return false;
    }

    bool succeeded = segment_set_search_radius_ids_snapshot(set, lat, lng, radius_km, result, true, stats);

    segment_set_read_unlock(set);

    return succeeded;
}

bool geo_segment_set_search_radius_ids_into_snapshot(const GeoSegmentSet *set,
                                                     double latitude,
                                                     double longitude,
                                                     double radius_km,
                                                     uint64_t *ids,
                                                     size_t capacity,
                                                     size_t *count)
{
    if (!set || !count || (capacity && !ids)) {
        return false;
    }

    GeoIdResult result = {
        .ids = ids,
        .count = 0,
        .capacity = capacity,
    };

    if (!segment_set_search_radius_ids_snapshot(set, latitude, longitude, radius_km, &result, false, NULL)) {
        return false;
    }

    *count = result.count;

    return true;
}

bool geo_segment_set_search_radius_ids_into(const GeoSegmentSet *set,
                                            double lat,
                                            double lng,
                                            double radius_km,
                                            uint64_t *ids,
                                            size_t capacity,
                                            size_t *count,
                                            GeoSearchStats *stats)
{
    if (!set || !count || (capacity && !ids) || !segment_set_read_lock(set)) {
        return false;
    }

    GeoIdResult result = {
        .ids = ids,
        .count = 0,
        .capacity = capacity,
    };

    bool succeeded = segment_set_search_radius_ids_snapshot(set, lat, lng, radius_km, &result, false, stats);

    if (succeeded) {
        *count = result.count;
    }

    segment_set_read_unlock(set);

    return succeeded;
}

bool geo_segment_set_search_bbox_reuse(const GeoSegmentSet *set,
                                       double min_lat,
                                       double max_lat,
                                       double min_lng,
                                       double max_lng,
                                       GeoSearchResult *result,
                                       GeoSearchStats *stats)
{
    segment_stats_reset(stats);

    GeoBboxQueryPlan plan;

    if (!set || !result || !geo_index_prepare_bbox_query(min_lat, max_lat, min_lng, max_lng, &plan)) {
        return false;
    }

    if (!segment_set_read_lock(set)) {
        return false;
    }

    double start = stats ? geo_get_time_ms() : 0.0;
    geo_result_clear(result);
    bool succeeded = true;

    for (size_t i = 0; succeeded && i < set->count; ++i) {
        GeoSearchStats partial_stats;
        GeoSegmentVisibility visibility;

        segment_visibility_initialize(&visibility, &set->mutations, set->segment_generations[i]);

        succeeded = geo_index_search_bbox_plan_append_filtered(set->segments[i],
                                                               &plan,
                                                               result,
                                                               NULL,
                                                               &partial_stats,
                                                               set->mutations.count ? segment_filter_visible_bits : NULL,
                                                               set->mutations.count ? segment_record_is_visible : NULL,
                                                               &visibility);

        if (succeeded) {
            segment_stats_add(stats, &partial_stats);
        }
    }

    if (!succeeded) {
        geo_result_clear(result);
    }

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start;
    }

    segment_set_read_unlock(set);

    return succeeded;
}

GeoSearchResult *geo_segment_set_search_bbox(const GeoSegmentSet *set,
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

    if (!geo_segment_set_search_bbox_reuse(set, min_lat, max_lat, min_lng, max_lng, result, stats)) {
        geo_result_destroy(result);

        return NULL;
    }

    return result;
}

bool geo_segment_set_search_bbox_count(const GeoSegmentSet *set,
                                       double min_lat,
                                       double max_lat,
                                       double min_lng,
                                       double max_lng,
                                       size_t *count,
                                       GeoSearchStats *stats)
{
    segment_stats_reset(stats);

    GeoBboxQueryPlan plan;

    if (!set || !count || !geo_index_prepare_bbox_query(min_lat, max_lat, min_lng, max_lng, &plan)) {
        return false;
    }

    if (!segment_set_read_lock(set)) {
        return false;
    }

    double start = stats ? geo_get_time_ms() : 0.0;
    size_t total_count = 0;
    bool succeeded = true;

    for (size_t i = 0; succeeded && i < set->count; ++i) {
        GeoSearchStats partial_stats;
        size_t partial_count;
        GeoSegmentVisibility visibility;

        segment_visibility_initialize(&visibility, &set->mutations, set->segment_generations[i]);

        succeeded = geo_index_search_bbox_plan_append_filtered(set->segments[i],
                                                               &plan,
                                                               NULL,
                                                               &partial_count,
                                                               &partial_stats,
                                                               set->mutations.count ? segment_filter_visible_bits : NULL,
                                                               set->mutations.count ? segment_record_is_visible : NULL,
                                                               &visibility) &&
                    partial_count <= SIZE_MAX - total_count;

        if (succeeded) {
            total_count += partial_count;
            segment_stats_add(stats, &partial_stats);
        }
    }

    if (succeeded) {
        *count = total_count;
    }

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start;
    }

    segment_set_read_unlock(set);

    return succeeded;
}

GeoSearchResult *geo_segment_set_search_knn(const GeoSegmentSet *set,
                                            double lat,
                                            double lng,
                                            size_t k,
                                            double max_radius_km,
                                            GeoSearchStats *stats)
{
    segment_stats_reset(stats);

    if (!set || !geo_is_valid_point(lat, lng) || !k || !isfinite(max_radius_km) || max_radius_km < 0.0) {
        return NULL;
    }

    if (!segment_set_read_lock(set)) {
        return NULL;
    }

    GeoSearchResult *result = geo_result_create(k);
    GeoKnnWorkspace *workspace = geo_knn_workspace_create(k, 64);
    GeoKnnSource *sources = set->count ? malloc(set->count * sizeof(*sources)) : NULL;
    GeoSegmentVisibility *visibilities = set->count ? malloc(set->count * sizeof(*visibilities)) : NULL;

    if (!result || !workspace || (set->count && (!sources || !visibilities))) {
        free(visibilities);
        free(sources);
        geo_knn_workspace_destroy(workspace);
        geo_result_destroy(result);
        segment_set_read_unlock(set);

        return NULL;
    }

    for (size_t i = 0; i < set->count; ++i) {
        segment_visibility_initialize(visibilities + i, &set->mutations, set->segment_generations[i]);
        sources[i] = (GeoKnnSource) {
            .index = set->segments[i],
            .filter = set->mutations.count ? segment_record_is_visible : NULL,
            .filter_context = visibilities + i,
        };
    }

    bool succeeded = geo_index_search_knn_sources(sources,
                                                  set->count,
                                                  lat,
                                                  lng,
                                                  k,
                                                  max_radius_km,
                                                  result,
                                                  workspace,
                                                  stats);

    free(visibilities);
    free(sources);
    geo_knn_workspace_destroy(workspace);

    if (!succeeded) {
        geo_result_destroy(result);
        segment_set_read_unlock(set);

        return NULL;
    }

    segment_set_read_unlock(set);

    return result;
}

// =============================================================================
// Direct Morton compaction
// =============================================================================

#if GEO_SEGMENTS_SUPPORTED
static bool segment_merge_node_less(const GeoSegmentMergeNode *first, const GeoSegmentMergeNode *second)
{
    if (first->record.z != second->record.z) {
        return first->record.z < second->record.z;
    }

    if (first->record.id != second->record.id) {
        return first->record.id < second->record.id;
    }

    return first->segment_index < second->segment_index;
}

static void segment_merge_heap_push(GeoSegmentMergeNode *heap, size_t *count, GeoSegmentMergeNode node)
{
    size_t position = (*count)++;

    while (position) {
        size_t parent = (position - 1) >> 1;

        if (!segment_merge_node_less(&node, heap + parent)) {
            break;
        }

        heap[position] = heap[parent];
        position = parent;
    }

    heap[position] = node;
}

static void segment_merge_heap_replace_min(GeoSegmentMergeNode *heap, size_t count, GeoSegmentMergeNode replacement)
{
    size_t position = 0;

    while (position * 2 + 1 < count) {
        size_t left = position * 2 + 1;
        size_t right = left + 1;
        size_t smallest = right < count && segment_merge_node_less(heap + right, heap + left) ? right : left;

        if (!segment_merge_node_less(heap + smallest, &replacement)) {
            break;
        }

        heap[position] = heap[smallest];
        position = smallest;
    }

    heap[position] = replacement;
}

static void segment_merge_heap_pop(GeoSegmentMergeNode *heap, size_t *count)
{
    GeoSegmentMergeNode replacement = heap[--(*count)];

    if (*count) {
        segment_merge_heap_replace_min(heap, *count, replacement);
    }
}

static void segment_preallocate_output(FILE *file, uint64_t record_count, size_t metadata_bytes)
{
#if defined(__linux__)
    if (metadata_bytes > (size_t) INT64_MAX) {
        return;
    }

    uint64_t maximum_records = ((uint64_t) INT64_MAX - metadata_bytes) / sizeof(GeoRecord);

    if (record_count <= maximum_records) {
        off_t file_size = (off_t) (record_count * sizeof(GeoRecord) + metadata_bytes);

        (void) posix_fallocate(fileno(file), 0, file_size);
    }
#else
    (void) file;
    (void) record_count;
    (void) metadata_bytes;
#endif
}

typedef struct {
    uint64_t range_min;
    uint64_t range_end;
    size_t output_begin;
    size_t record_count;
    uint64_t zero_seed_checksum;
    GeoDensityWriter *density_writer;
    bool reaches_morton_end;
    bool succeeded;
} GeoCompactionPartition;

typedef struct {
    GeoIndex *const *segments;
    const uint64_t *segment_generations;
    const GeoMutationTable *mutations;
    const uint64_t *const *visibility_bits;
    GeoCompactionPartition *partitions;
    size_t partition_count;
    size_t segment_count;
    size_t *prefix_offsets;
    GeoDensityWriter *single_pass_density;
    int descriptor;
    uint8_t prefix_bits;
    atomic_size_t next_partition;
    atomic_bool failed;
} GeoCompactionMerge;

typedef struct {
    GeoSegmentMergeNode *heap;
    GeoSegmentVisibility *visibilities;
    size_t *segment_ends;
    GeoRecord *output_buffer;
} GeoCompactionScratch;

static bool segment_pwrite_all(int descriptor, const void *data, size_t size, off_t offset)
{
    const unsigned char *bytes = data;

    while (size) {
        ssize_t written = pwrite(descriptor, bytes, size, offset);

        if (written < 0 && errno == EINTR) {
            continue;
        }

        if (written <= 0) {
            return false;
        }

        bytes += (size_t) written;
        size -= (size_t) written;
        offset += written;
    }

    return true;
}

static bool segment_compaction_scratch_initialize(const GeoCompactionMerge *merge, GeoCompactionScratch *scratch)
{
    memset(scratch, 0, sizeof(*scratch));

    scratch->heap = malloc(merge->segment_count * sizeof(*scratch->heap));
    scratch->segment_ends = malloc(merge->segment_count * sizeof(*scratch->segment_ends));

    if (merge->mutations->count && !merge->visibility_bits && merge->segment_count) {
        scratch->visibilities = malloc(merge->segment_count * sizeof(*scratch->visibilities));
    }

    scratch->output_buffer = malloc(GEO_SEGMENT_OUTPUT_BUFFER_RECORDS * sizeof(*scratch->output_buffer));

    if (!scratch->heap ||
        !scratch->segment_ends ||
        (merge->mutations->count && !merge->visibility_bits && merge->segment_count && !scratch->visibilities) ||
        !scratch->output_buffer) {
        return false;
    }

    for (size_t segment = 0; scratch->visibilities && segment < merge->segment_count; ++segment) {
        segment_visibility_initialize(scratch->visibilities + segment,
                                      merge->mutations,
                                      merge->segment_generations[segment]);
    }

    return true;
}

static void segment_compaction_scratch_destroy(GeoCompactionScratch *scratch)
{
    free(scratch->output_buffer);
    free(scratch->segment_ends);
    free(scratch->visibilities);
    free(scratch->heap);
}

static bool segment_compaction_flush_records(const GeoCompactionMerge *merge,
                                             GeoCompactionPartition *partition,
                                             GeoCompactionScratch *scratch,
                                             size_t first_local_record,
                                             size_t count)
{
    size_t first_output_record = partition->output_begin + first_local_record;
    size_t bytes = count * sizeof(*scratch->output_buffer);
    uint64_t byte_offset = sizeof(GeoFileHeader) + (uint64_t) first_output_record * sizeof(GeoRecord);

    GeoDensityWriter *density_writer = merge->single_pass_density
                                           ? merge->single_pass_density
                                           : partition->density_writer;

    if (density_writer &&
        !geo_density_writer_add(density_writer, scratch->output_buffer, count, first_output_record)) {
        return false;
    }

    partition->zero_seed_checksum = geo_persisted_checksum_update(partition->zero_seed_checksum,
                                                                   scratch->output_buffer,
                                                                   bytes);

    return segment_pwrite_all(merge->descriptor, scratch->output_buffer, bytes, (off_t) byte_offset);
}

static bool segment_compaction_record_visible(const GeoCompactionMerge *merge,
                                              const GeoCompactionScratch *scratch,
                                              const GeoSegmentMergeNode *node)
{
    size_t record_position = node->position - 1U;

    if (merge->visibility_bits) {
        return (merge->visibility_bits[node->segment_index][record_position >> 6U] &
                (UINT64_C(1) << (record_position & 63U))) != 0;
    }

    if (scratch->visibilities) {
        return segment_record_is_visible_for(&node->record, scratch->visibilities + node->segment_index);
    }

    return true;
}

static bool segment_merge_partition(GeoCompactionMerge *merge,
                                    GeoCompactionPartition *partition,
                                    GeoCompactionScratch *scratch)
{
    size_t heap_count = 0;

    for (size_t segment_index = 0; segment_index < merge->segment_count; ++segment_index) {
        const GeoIndex *segment = merge->segments[segment_index];
        size_t begin = partition->range_min
                           ? geo_lower_bound(segment->records, segment->count, partition->range_min)
                           : 0;
        size_t end = partition->reaches_morton_end
                         ? segment->count
                         : geo_lower_bound(segment->records, segment->count, partition->range_end);

        scratch->segment_ends[segment_index] = end;

        if (begin < end) {
            segment_merge_heap_push(scratch->heap,
                                    &heap_count,
                                    (GeoSegmentMergeNode) {
                                        .record = segment->records[begin],
                                        .segment_index = segment_index,
                                        .position = begin + 1U,
                                    });
        }
    }

    size_t buffered = 0;
    size_t records_written = 0;
    size_t next_prefix = 0;
    bool succeeded = true;

    while (succeeded && heap_count) {
        GeoSegmentMergeNode node = scratch->heap[0];
        bool visible = segment_compaction_record_visible(merge, scratch, &node);

        if (visible && merge->prefix_offsets) {
            size_t record_prefix = (size_t) (node.record.z >> (64U - merge->prefix_bits));

            while (next_prefix <= record_prefix) {
                merge->prefix_offsets[next_prefix++] = partition->output_begin + records_written;
            }
        }

        if (visible) {
            scratch->output_buffer[buffered++] = node.record;
            records_written++;
        }

        if (buffered == GEO_SEGMENT_OUTPUT_BUFFER_RECORDS) {
            succeeded = segment_compaction_flush_records(merge,
                                                         partition,
                                                         scratch,
                                                         records_written - buffered,
                                                         buffered);
            buffered = 0;
        }

        if (succeeded && node.position < scratch->segment_ends[node.segment_index]) {
            node.record = merge->segments[node.segment_index]->records[node.position++];
            segment_merge_heap_replace_min(scratch->heap, heap_count, node);
        } else if (succeeded) {
            segment_merge_heap_pop(scratch->heap, &heap_count);
        }
    }

    if (succeeded && buffered) {
        succeeded = segment_compaction_flush_records(merge,
                                                     partition,
                                                     scratch,
                                                     records_written - buffered,
                                                     buffered);
    }

    if (succeeded && merge->prefix_offsets) {
        size_t prefix_count = ((size_t) 1 << merge->prefix_bits) + 1U;

        while (next_prefix < prefix_count) {
            merge->prefix_offsets[next_prefix++] = partition->output_begin + records_written;
        }
    }

    partition->succeeded = succeeded && records_written == partition->record_count;

    return partition->succeeded;
}

static bool segment_compaction_merge_worker(void *argument, size_t worker_index, size_t worker_count)
{
    GeoCompactionMerge *merge = argument;
    GeoCompactionScratch scratch;

    (void) worker_index;
    (void) worker_count;

    if (!segment_compaction_scratch_initialize(merge, &scratch)) {
        segment_compaction_scratch_destroy(&scratch);
        atomic_store_explicit(&merge->failed, true, memory_order_release);

        return false;
    }

    while (!atomic_load_explicit(&merge->failed, memory_order_acquire)) {
        size_t partition_index = atomic_fetch_add_explicit(&merge->next_partition, 1U, memory_order_relaxed);

        if (partition_index >= merge->partition_count) {
            break;
        }

        if (!segment_merge_partition(merge, merge->partitions + partition_index, &scratch)) {
            atomic_store_explicit(&merge->failed, true, memory_order_release);
            break;
        }
    }

    segment_compaction_scratch_destroy(&scratch);

    return !atomic_load_explicit(&merge->failed, memory_order_acquire);
}

static size_t segment_compaction_worker_count(uint64_t record_count, size_t requested_workers)
{
    uint64_t capacity_workers = record_count / GEO_COMPACTION_MIN_RECORDS_PER_WORKER;
    size_t workers = requested_workers ? requested_workers : capacity_workers > SIZE_MAX ? SIZE_MAX : (size_t) capacity_workers;
    long online_processors = sysconf(_SC_NPROCESSORS_ONLN);
    size_t online_workers = online_processors > 0 ? (size_t) online_processors : 1U;

    if (!workers) {
        workers = 1;
    }

    if (workers > online_workers) {
        workers = online_workers;
    }

    if (workers > GEO_COMPACTION_MAX_WORKERS) {
        workers = GEO_COMPACTION_MAX_WORKERS;
    }

    return workers;
}

static bool segment_compaction_count_ranges(const GeoCompactionMerge *merge,
                                            size_t expected_records,
                                            size_t *prefix_counts)
{
    GeoSegmentVisibility *visibilities = NULL;

    if (merge->mutations->count && !merge->visibility_bits && merge->segment_count) {
        visibilities = malloc(merge->segment_count * sizeof(*visibilities));

        if (!visibilities) {
            return false;
        }

        for (size_t segment = 0; segment < merge->segment_count; ++segment) {
            segment_visibility_initialize(visibilities + segment,
                                          merge->mutations,
                                          merge->segment_generations[segment]);
        }
    }

    size_t records_counted = 0;

    for (size_t segment_index = 0; segment_index < merge->segment_count; ++segment_index) {
        const GeoIndex *segment = merge->segments[segment_index];
        size_t partition_index = 0;

        for (size_t position = 0; position < segment->count; ++position) {
            bool visible = true;

            while (!merge->partitions[partition_index].reaches_morton_end &&
                   segment->records[position].z >= merge->partitions[partition_index].range_end) {
                partition_index++;
            }

            if (merge->visibility_bits) {
                visible = (merge->visibility_bits[segment_index][position >> 6U] &
                           (UINT64_C(1) << (position & 63U))) != 0;
            } else if (visibilities) {
                visible = segment_record_is_visible_for(segment->records + position, visibilities + segment_index);
            }

            if (!visible) {
                continue;
            }

            size_t prefix = (size_t) (segment->records[position].z >> (64U - GEO_COMPACTION_PARTITION_BITS));

            prefix_counts[prefix]++;
            merge->partitions[partition_index].record_count++;
            records_counted++;
        }
    }

    free(visibilities);

    return records_counted == expected_records;
}

static void segment_compaction_build_prefix_offsets(const size_t *fine_prefix_offsets,
                                                    uint8_t prefix_bits,
                                                    size_t *prefix_offsets)
{
    size_t prefix_count = (size_t) 1 << prefix_bits;
    unsigned fine_bits = GEO_COMPACTION_PARTITION_BITS - prefix_bits;

    for (size_t prefix = 0; prefix <= prefix_count; ++prefix) {
        prefix_offsets[prefix] = fine_prefix_offsets[prefix << fine_bits];
    }
}

static int segment_compaction_sample_compare(const void *first_value, const void *second_value)
{
    uint64_t first = *(const uint64_t *) first_value;
    uint64_t second = *(const uint64_t *) second_value;

    return (first > second) - (first < second);
}

static size_t segment_compaction_build_partitions(const GeoCompactionMerge *merge,
                                                  size_t worker_count,
                                                  GeoCompactionPartition *partitions)
{
    const size_t samples_per_partition = 64U;
    size_t desired_partitions = worker_count * GEO_COMPACTION_TASKS_PER_WORKER;
    size_t maximum_samples = desired_partitions * samples_per_partition;
    size_t physical_records = 0;

    for (size_t segment = 0; segment < merge->segment_count; ++segment) {
        if (merge->segments[segment]->count > SIZE_MAX - physical_records) {
            return 0;
        }

        physical_records += merge->segments[segment]->count;
    }

    size_t sample_count = physical_records < maximum_samples ? physical_records : maximum_samples;
    uint64_t *samples = sample_count ? malloc(sample_count * sizeof(*samples)) : NULL;

    if (sample_count && !samples) {
        return 0;
    }

    size_t segment_index = 0;
    size_t segment_begin = 0;

    for (size_t sample = 0; sample < sample_count; ++sample) {
        size_t ordinal = (physical_records / sample_count) * sample +
                         ((physical_records % sample_count) * sample) / sample_count;

        while (ordinal - segment_begin >= merge->segments[segment_index]->count) {
            segment_begin += merge->segments[segment_index]->count;
            segment_index++;
        }

        samples[sample] = merge->segments[segment_index]->records[ordinal - segment_begin].z;
    }

    if (sample_count) {
        qsort(samples, sample_count, sizeof(*samples), segment_compaction_sample_compare);
    }

    size_t partition_count = 0;
    uint64_t range_begin = 0;

    for (size_t partition = 1; sample_count && partition < desired_partitions; ++partition) {
        size_t sample_index = (sample_count / desired_partitions) * partition +
                              ((sample_count % desired_partitions) * partition) / desired_partitions;
        uint64_t range_end = samples[sample_index];

        if (range_end <= range_begin) {
            continue;
        }

        partitions[partition_count++] = (GeoCompactionPartition) {
            .range_min = range_begin,
            .range_end = range_end,
        };
        range_begin = range_end;
    }

    partitions[partition_count++] = (GeoCompactionPartition) {
        .range_min = range_begin,
        .reaches_morton_end = true,
    };

    free(samples);

    return partition_count;
}

static bool segment_compaction_execute(GeoCompactionMerge *merge,
                                       GeoThreadPool *pool,
                                       size_t worker_limit,
                                       size_t *workers_executed)
{
    atomic_init(&merge->next_partition, 0);
    atomic_init(&merge->failed, false);

    if (worker_limit == 1U) {
        *workers_executed = 1U;

        return segment_compaction_merge_worker(merge, 0U, 1U);
    }

    return pool &&
           geo_thread_pool_run(pool,
                               worker_limit,
                               segment_compaction_merge_worker,
                               merge,
                               workers_executed) &&
           !atomic_load_explicit(&merge->failed, memory_order_acquire);
}

static void segment_compaction_destroy_density_partitions(GeoDensityWriter **density_partitions, size_t partition_count)
{
    if (!density_partitions) {
        return;
    }

    for (size_t partition = 0; partition < partition_count; ++partition) {
        geo_density_writer_destroy(density_partitions[partition]);
    }

    free(density_partitions);
}

static bool segment_write_compacted_file(GeoIndex *const *segments,
                                         const uint64_t *segment_generations,
                                         const GeoMutationTable *mutations,
                                         const uint64_t *const *visibility_bits,
                                         size_t segment_count,
                                         uint64_t record_count,
                                         const char *output_path,
                                         GeoThreadPool *compaction_pool,
                                         size_t requested_workers,
                                         size_t *workers_used,
                                         size_t *partitions_used)
{
    if (!segments || !segment_generations || !mutations || !segment_count || !output_path || !workers_used || !partitions_used ||
        record_count > SIZE_MAX ||
        record_count > ((uint64_t) INT64_MAX - sizeof(GeoFileHeader)) / sizeof(GeoRecord) ||
        segment_count > SIZE_MAX / sizeof(GeoSegmentMergeNode) ||
        segment_count > SIZE_MAX / sizeof(size_t) ||
        (mutations->count && !visibility_bits && segment_count > SIZE_MAX / sizeof(GeoSegmentVisibility))) {
        return false;
    }

    uint8_t prefix_bits = geo_persisted_prefix_bits(record_count);
    size_t prefix_count = prefix_bits ? ((size_t) 1 << prefix_bits) + 1U : 0;
    size_t metadata_bytes = sizeof(GeoFileHeader) + prefix_count * sizeof(uint64_t);
    size_t *prefix_offsets = prefix_count ? malloc(prefix_count * sizeof(*prefix_offsets)) : NULL;
    GeoDensityWriter *density_writer = geo_density_writer_create(record_count, prefix_bits);
    size_t worker_count = segment_compaction_worker_count(record_count, requested_workers);
    bool uses_partitioned_metadata = worker_count > 1U;
    size_t maximum_partitions = worker_count * GEO_COMPACTION_TASKS_PER_WORKER;
    GeoCompactionPartition *partitions = calloc(maximum_partitions, sizeof(*partitions));
    size_t *fine_prefix_counts = NULL;
    size_t *fine_prefix_offsets = NULL;
    GeoDensityWriter **density_partitions = NULL;

    if ((prefix_count && !prefix_offsets) || !density_writer || !partitions) {
        free(partitions);
        geo_density_writer_destroy(density_writer);
        free(prefix_offsets);

        return false;
    }

    GeoCompactionMerge merge = {
        .segments = segments,
        .segment_generations = segment_generations,
        .mutations = mutations,
        .visibility_bits = visibility_bits,
        .partitions = partitions,
        .partition_count = 1,
        .segment_count = segment_count,
        .prefix_offsets = uses_partitioned_metadata ? NULL : prefix_offsets,
        .single_pass_density = uses_partitioned_metadata ? NULL : density_writer,
        .prefix_bits = prefix_bits,
    };

    partitions[0] = (GeoCompactionPartition) {
        .record_count = (size_t) record_count,
        .reaches_morton_end = true,
    };

    if (uses_partitioned_metadata) {
        size_t fine_prefix_count = (size_t) 1 << GEO_COMPACTION_PARTITION_BITS;

        fine_prefix_counts = calloc(fine_prefix_count, sizeof(*fine_prefix_counts));
        fine_prefix_offsets = malloc((fine_prefix_count + 1U) * sizeof(*fine_prefix_offsets));

        if (fine_prefix_counts && fine_prefix_offsets) {
            merge.partition_count = segment_compaction_build_partitions(&merge, worker_count, partitions);
        }

        bool counted = fine_prefix_counts &&
                       fine_prefix_offsets &&
                       merge.partition_count &&
                       segment_compaction_count_ranges(&merge, (size_t) record_count, fine_prefix_counts);

        size_t output_begin = 0;

        for (size_t partition = 0; counted && partition < merge.partition_count; ++partition) {
            partitions[partition].output_begin = output_begin;

            if (partitions[partition].record_count > (size_t) record_count - output_begin) {
                counted = false;
            } else {
                output_begin += partitions[partition].record_count;
            }
        }

        counted = counted && output_begin == (size_t) record_count;

        if (!counted) {
            free(fine_prefix_offsets);
            free(fine_prefix_counts);
            free(partitions);
            geo_density_writer_destroy(density_writer);
            free(prefix_offsets);

            return false;
        }

        fine_prefix_offsets[0] = 0;

        for (size_t prefix = 0; prefix < fine_prefix_count; ++prefix) {
            fine_prefix_offsets[prefix + 1U] = fine_prefix_offsets[prefix] + fine_prefix_counts[prefix];
        }

        segment_compaction_build_prefix_offsets(fine_prefix_offsets, prefix_bits, prefix_offsets);

        density_partitions = calloc(merge.partition_count, sizeof(*density_partitions));

        for (size_t partition = 0; density_partitions && partition < merge.partition_count; ++partition) {
            density_partitions[partition] = geo_density_writer_create_partition(record_count,
                                                                                 prefix_bits,
                                                                                 partitions[partition].output_begin,
                                                                                 partitions[partition].record_count,
                                                                                 prefix_offsets);
            partitions[partition].density_writer = density_partitions[partition];

            if (!density_partitions[partition]) {
                segment_compaction_destroy_density_partitions(density_partitions, merge.partition_count);
                density_partitions = NULL;
                break;
            }
        }

        if (!density_partitions) {
            free(fine_prefix_offsets);
            free(fine_prefix_counts);
            free(partitions);
            geo_density_writer_destroy(density_writer);
            free(prefix_offsets);

            return false;
        }

        if (worker_count > merge.partition_count) {
            worker_count = merge.partition_count;
        }
    }

    char *temporary_path = NULL;
    FILE *output = geo_io_create_atomic_file(output_path, &temporary_path);

    if (!output) {
        segment_compaction_destroy_density_partitions(density_partitions, merge.partition_count);
        free(fine_prefix_offsets);
        free(fine_prefix_counts);
        free(partitions);
        geo_density_writer_destroy(density_writer);
        free(prefix_offsets);

        return false;
    }

#if defined(POSIX_FADV_SEQUENTIAL)
    (void) posix_fadvise(fileno(output), 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
    segment_preallocate_output(output, record_count, metadata_bytes);

    GeoFileHeader header = {
        .magic = { 0 },
        .version = GEO_FILE_VERSION,
        .record_size = sizeof(GeoRecord),
        .endian_marker = GEO_FILE_ENDIAN_MARKER,
        .prefix_bits = prefix_bits,
        .count = record_count,
    };

    memcpy(header.magic, GEO_FILE_MAGIC, sizeof(header.magic));
    merge.descriptor = fileno(output);
    bool succeeded = fwrite(&header, sizeof(header), 1, output) == 1 && fflush(output) == 0;

    if (succeeded) {
        succeeded = segment_compaction_execute(&merge, compaction_pool, worker_count, &worker_count);
    }

    uint64_t records_checksum = geo_persisted_checksum_initial();

    for (size_t partition = 0; succeeded && partition < merge.partition_count; ++partition) {
        size_t suffix_bytes = partitions[partition].record_count * sizeof(GeoRecord);

        succeeded = partitions[partition].succeeded;
        records_checksum = geo_persisted_checksum_combine_aligned(records_checksum,
                                                                   partitions[partition].zero_seed_checksum,
                                                                   suffix_bytes);
    }

    if (succeeded && uses_partitioned_metadata) {
        succeeded = geo_density_writer_merge_partitions(density_writer, density_partitions, merge.partition_count);
    }

    uint64_t records_bytes = record_count * sizeof(GeoRecord);
    uint64_t metadata_offset = sizeof(GeoFileHeader) + records_bytes;

    if (succeeded) {
        succeeded = fseeko(output, (off_t) metadata_offset, SEEK_SET) == 0;
    }

    uint64_t prefix_checksum = geo_persisted_checksum_initial();

    if (succeeded && prefix_offsets) {
        succeeded = geo_io_write_u64_offsets(output, prefix_offsets, prefix_count, &prefix_checksum);
    }

    uint64_t density_bytes = 0;
    uint64_t density_checksum = 0;

    if (succeeded) {
        succeeded = geo_density_writer_append_file(density_writer,
                                                   output,
                                                   &density_bytes,
                                                   &density_checksum);
    }

    if (succeeded) {
        header.records_checksum = records_checksum;
        header.prefix_checksum = prefix_checksum;
        header.density_checksum = density_checksum;
        header.density_bytes = density_bytes;
        succeeded = fseeko(output, 0, SEEK_SET) == 0 && fwrite(&header, sizeof(header), 1, output) == 1;
    }

    if (succeeded) {
        succeeded = geo_io_publish_atomic_file(output, temporary_path, output_path);
    } else {
        geo_io_discard_atomic_file(output, temporary_path);
    }

    free(temporary_path);
    segment_compaction_destroy_density_partitions(density_partitions, merge.partition_count);
    free(fine_prefix_offsets);
    free(fine_prefix_counts);
    free(partitions);
    geo_density_writer_destroy(density_writer);
    free(prefix_offsets);

    if (succeeded) {
        *workers_used = worker_count;
        *partitions_used = merge.partition_count;
    }

    return succeeded;
}

static int segment_candidate_compare(const void *first_value, const void *second_value)
{
    const GeoSegmentCandidate *first = first_value;
    const GeoSegmentCandidate *second = second_value;

    if (first->record_count < second->record_count) {
        return -1;
    }

    if (first->record_count > second->record_count) {
        return 1;
    }

    return (first->index > second->index) - (first->index < second->index);
}

static int segment_index_compare(const void *first_value, const void *second_value)
{
    size_t first = *(const size_t *) first_value;
    size_t second = *(const size_t *) second_value;

    return (first > second) - (first < second);
}

static bool segment_set_compact_selection_locked(GeoSegmentSet *set,
                                                 const char *output_path,
                                                 const size_t *selected_indices,
                                                 size_t selected_count,
                                                 size_t requested_workers,
                                                 GeoSegmentCompactionStats *stats)
{
    double start = geo_get_time_ms();
    uint64_t selected_physical_count = 0;

    if (!selected_count ||
        selected_count > SIZE_MAX / sizeof(GeoIndex *) ||
        selected_count > SIZE_MAX / sizeof(uint64_t)) {
        return false;
    }

    for (size_t i = 0; i < selected_count; ++i) {
        size_t segment_index = selected_indices[i];

        if (segment_index >= set->count ||
            (i && selected_indices[i - 1] >= segment_index) ||
            set->segments[segment_index]->count > UINT64_MAX - selected_physical_count) {
            return false;
        }

        selected_physical_count += set->segments[segment_index]->count;
    }

    if (selected_physical_count > SIZE_MAX) {
        return false;
    }

    size_t desired_workers = segment_compaction_worker_count(selected_physical_count, requested_workers);

    if (desired_workers > 1U && !set->compaction_pool) {
        set->compaction_pool = geo_thread_pool_create(GEO_COMPACTION_MAX_WORKERS);

        if (!set->compaction_pool) {
            return false;
        }
    }

    GeoIndex **merge_segments = malloc(selected_count * sizeof(*merge_segments));
    uint64_t *merge_generations = malloc(selected_count * sizeof(*merge_generations));
    GeoMutationTable snapshot_mutations = { 0 };

    if (!merge_segments || !merge_generations || !mutation_table_clone(&set->mutations, &snapshot_mutations, 0)) {
        mutation_table_destroy(&snapshot_mutations);
        free(merge_generations);
        free(merge_segments);

        return false;
    }

    for (size_t i = 0; i < selected_count; ++i) {
        merge_segments[i] = set->segments[selected_indices[i]];
        merge_generations[i] = set->segment_generations[selected_indices[i]];
    }

    uint64_t **visibility_bits = NULL;
    uint64_t *visibility_storage = NULL;

    bool materialize_visibility = snapshot_mutations.count > GEO_MUTATION_SMALL_ENTRY_CAPACITY;

    if (materialize_visibility) {
        size_t total_visibility_words = 0;
        bool visibility_size_valid = selected_count <= SIZE_MAX / sizeof(*visibility_bits);

        for (size_t i = 0; visibility_size_valid && i < selected_count; ++i) {
            size_t segment_words = geo_internal_bit_word_count(merge_segments[i]->count);

            visibility_size_valid = segment_words <= SIZE_MAX - total_visibility_words;
            total_visibility_words += visibility_size_valid ? segment_words : 0;
        }

        visibility_size_valid = visibility_size_valid &&
                                total_visibility_words <= SIZE_MAX / sizeof(*visibility_storage);

        visibility_bits = visibility_size_valid ? malloc(selected_count * sizeof(*visibility_bits)) : NULL;

        if (visibility_size_valid && total_visibility_words) {
            visibility_storage = malloc(total_visibility_words * sizeof(*visibility_storage));
        }

        if (!visibility_size_valid || !visibility_bits || (total_visibility_words && !visibility_storage)) {
            free(visibility_storage);
            free(visibility_bits);
            mutation_table_destroy(&snapshot_mutations);
            free(merge_generations);
            free(merge_segments);

            return false;
        }

        size_t visibility_offset = 0;

        for (size_t i = 0; i < selected_count; ++i) {
            size_t segment_words = geo_internal_bit_word_count(merge_segments[i]->count);

            visibility_bits[i] = visibility_storage ? visibility_storage + visibility_offset : NULL;
            visibility_offset += segment_words;
        }
    }

    uint64_t snapshot_generation = set->generation;
    pthread_mutex_unlock(&set->mutation_lock);

    uint64_t selected_record_count = selected_physical_count;

    if (snapshot_mutations.count) {
        selected_record_count = 0;

        for (size_t i = 0; i < selected_count; ++i) {
            GeoSegmentVisibility visibility;

            segment_visibility_initialize(&visibility, &snapshot_mutations, merge_generations[i]);

            selected_record_count += segment_scan_visible_records(merge_segments[i],
                                                                  &visibility,
                                                                  materialize_visibility ? visibility_bits[i] : NULL);
        }
    }

    size_t workers_used = 0;
    size_t partitions_used = 0;
    bool succeeded = segment_write_compacted_file(merge_segments,
                                                  merge_generations,
                                                  &snapshot_mutations,
                                                  (const uint64_t *const *) visibility_bits,
                                                  selected_count,
                                                  selected_record_count,
                                                  output_path,
                                                  set->compaction_pool,
                                                  requested_workers,
                                                  &workers_used,
                                                  &partitions_used);
    free(visibility_storage);
    free(visibility_bits);
    char *canonical_path = succeeded ? realpath(output_path, NULL) : NULL;
    GeoIndex *compacted = canonical_path ? geo_index_open_mmap(canonical_path) : NULL;

    succeeded = canonical_path && compacted && compacted->count == selected_record_count;

    bool lock_reacquired = pthread_mutex_lock(&set->mutation_lock) == 0;

    if (!succeeded || !lock_reacquired) {
        geo_index_destroy(compacted);

        if (canonical_path) {
            (void) unlink(canonical_path);
        }

        free(canonical_path);
        mutation_table_destroy(&snapshot_mutations);
        free(merge_generations);
        free(merge_segments);

        return false;
    }

    succeeded = succeeded && set->generation < GEO_MUTATION_MAX_GENERATION;

    for (size_t i = 0; succeeded && i < selected_count; ++i) {
        size_t segment_index = selected_indices[i];

        succeeded = segment_index < set->count &&
                    set->segments[segment_index] == merge_segments[i] &&
                    set->segment_generations[segment_index] == merge_generations[i];
    }

    mutation_table_destroy(&snapshot_mutations);
    free(merge_generations);
    free(merge_segments);

    if (!succeeded) {
        geo_index_destroy(compacted);
        (void) unlink(canonical_path);
        free(canonical_path);

        return false;
    }

    size_t new_count = set->count - selected_count + 1;
    size_t new_capacity = set->capacity > new_count ? set->capacity : new_count;
    char **new_paths = malloc(new_capacity * sizeof(*new_paths));
    GeoIndex **new_segments = malloc(new_capacity * sizeof(*new_segments));
    uint64_t *new_segment_generations = malloc(new_capacity * sizeof(*new_segment_generations));

    if (!new_paths || !new_segments || !new_segment_generations) {
        free(new_segment_generations);
        free(new_segments);
        free(new_paths);
        geo_index_destroy(compacted);
        (void) unlink(canonical_path);
        free(canonical_path);

        return false;
    }

    size_t selected_position = 0;
    size_t output_position = 0;

    for (size_t i = 0; i < set->count; ++i) {
        if (selected_position < selected_count && selected_indices[selected_position] == i) {
            selected_position++;
        } else {
            new_paths[output_position] = set->paths[i];
            new_segments[output_position] = set->segments[i];
            new_segment_generations[output_position] = set->segment_generations[i];
            output_position++;
        }
    }

    new_paths[output_position] = canonical_path;
    new_segments[output_position] = compacted;
    uint64_t new_generation = set->generation + 1U;
    new_segment_generations[output_position] = new_generation;
    bool compacts_all_segments = selected_count == set->count && snapshot_generation == set->generation;

    size_t carried_live_count = 0;

    for (size_t entry_index = 0; !compacts_all_segments && entry_index < set->mutations.capacity; ++entry_index) {
        const GeoMutationEntry *entry = set->mutations.entries + entry_index;

        if (!mutation_entry_is_occupied(entry) || mutation_entry_state(entry) != GEO_MUTATION_LIVE) {
            continue;
        }

        for (size_t selected = 0; selected < selected_count; ++selected) {
            if (mutation_entry_generation(entry) == set->segment_generations[selected_indices[selected]]) {
                carried_live_count++;
                break;
            }
        }
    }

    size_t carry_allocation_count = carried_live_count ? carried_live_count : 1U;
    GeoMutationRecord *carried_live = malloc(carry_allocation_count * sizeof(*carried_live));
    succeeded = carried_live != NULL;
    size_t carry_position = 0;

    for (size_t entry_index = 0;
         succeeded && !compacts_all_segments && entry_index < set->mutations.capacity;
         ++entry_index) {
        const GeoMutationEntry *entry = set->mutations.entries + entry_index;

        if (!mutation_entry_is_occupied(entry) || mutation_entry_state(entry) != GEO_MUTATION_LIVE) {
            continue;
        }

        for (size_t selected = 0; selected < selected_count; ++selected) {
            if (mutation_entry_generation(entry) != set->segment_generations[selected_indices[selected]]) {
                continue;
            }

            /*
             * The selected live record moves into the compacted segment, but an older copy can still exist in an
             * unselected segment. Carrying LIVE forward to the output generation keeps every older generation hidden.
             * Only a full, generation-stable compaction may discard the mutation table entirely.
             */
            carried_live[carry_position] = (GeoMutationRecord) {
                .id = entry->id,
                .generation = new_generation,
                .state = GEO_MUTATION_LIVE,
            };
            carry_position++;
            break;
        }
    }

    succeeded = succeeded && carry_position == carried_live_count;

    GeoMutationTable compacted_mutations = { 0 };

    if (succeeded && !compacts_all_segments) {
        succeeded = mutation_table_clone(&set->mutations, &compacted_mutations, 0);
    }

    for (size_t i = 0; succeeded && i < carried_live_count; ++i) {
        succeeded = mutation_table_apply(&compacted_mutations, carried_live + i);
    }

    uint64_t updated_checksum = compacts_all_segments
                                    ? UINT64_C(1469598103934665603)
                                    : set->mutation_checksum;

    if (succeeded && carried_live_count > UINT64_MAX - set->mutation_count) {
        succeeded = false;
    }

    if (succeeded && !compacts_all_segments) {
        succeeded = mutation_log_append(set, carried_live, carried_live_count, &updated_checksum);
    }

    uint64_t updated_mutation_count = compacts_all_segments ? 0 : set->mutation_count + carried_live_count;

    if (succeeded) {
        succeeded = manifest_write(set->manifest_path,
                                   new_paths,
                                   new_segments,
                                   new_segment_generations,
                                   new_count,
                                   new_generation,
                                   updated_mutation_count,
                                   updated_checksum,
                                   set->mutation_slot,
                                   set->durable_watermark);
    }

    if (succeeded) {
        succeeded = segment_set_write_lock(set);
    }

    if (!succeeded) {
        mutation_table_destroy(&compacted_mutations);
        free(carried_live);
        free(new_segment_generations);
        free(new_segments);
        free(new_paths);
        geo_index_destroy(compacted);
        (void) unlink(canonical_path);
        free(canonical_path);

        return false;
    }

    char **old_paths = set->paths;
    GeoIndex **old_segments = set->segments;
    uint64_t *old_segment_generations = set->segment_generations;

    GeoMutationTable old_mutations = set->mutations;

    set->paths = new_paths;
    set->segments = new_segments;
    set->segment_generations = new_segment_generations;
    set->mutations = compacted_mutations;
    set->count = new_count;
    set->capacity = new_capacity;
    set->generation = new_generation;
    set->mutation_count = updated_mutation_count;
    set->mutation_checksum = updated_checksum;
    set->record_count = set->record_count - selected_physical_count + selected_record_count;

    segment_set_write_unlock(set);

    if (compacts_all_segments) {
        // The durable manifest already commits an empty mutation prefix. Truncating afterwards is crash-safe:
        // an old tail is ignored if truncation fails or the process exits before this best-effort reclamation.
        (void) mutation_log_create(set->mutation_path);
        (void) geo_io_sync_parent_directory(set->mutation_path);
    }

    for (size_t i = 0; i < selected_count; ++i) {
        size_t segment_index = selected_indices[i];

        geo_index_destroy(old_segments[segment_index]);
        free(old_paths[segment_index]);
    }

    free(old_segments);
    free(old_paths);
    free(old_segment_generations);
    mutation_table_destroy(&old_mutations);
    free(carried_live);

    if (stats) {
        stats->records_written = selected_record_count;
        stats->input_segments = selected_count;
        stats->worker_count = workers_used;
        stats->partition_count = partitions_used;
        stats->merge_time_ms = geo_get_time_ms() - start;
    }

    return true;
}

static bool segment_set_compact_size_tier(GeoSegmentSet *set,
                                          const char *output_path,
                                          GeoSegmentCompactionStats *stats)
{
    if (stats) {
        memset(stats, 0, sizeof(*stats));
    }

    if (!set || !output_path || !output_path[0] || access(output_path, F_OK) == 0 ||
        pthread_mutex_lock(&set->mutation_lock) != 0) {
        return false;
    }

    if (set->count < 2 || set->generation >= GEO_MUTATION_MAX_GENERATION || set->compaction_running) {
        pthread_mutex_unlock(&set->mutation_lock);

        return false;
    }

    GeoSegmentCandidate *candidates = malloc(set->count * sizeof(*candidates));

    if (!candidates) {
        pthread_mutex_unlock(&set->mutation_lock);

        return false;
    }

    for (size_t i = 0; i < set->count; ++i) {
        candidates[i] = (GeoSegmentCandidate) {
            .index = i,
            .record_count = set->segments[i]->count,
        };
    }

    qsort(candidates, set->count, sizeof(*candidates), segment_candidate_compare);

    size_t maximum_selection = set->count < GEO_BACKGROUND_COMPACTION_FANOUT
                                   ? set->count
                                   : GEO_BACKGROUND_COMPACTION_FANOUT;
    size_t selection_count = 1;
    size_t tier_limit = candidates[0].record_count > SIZE_MAX / GEO_BACKGROUND_COMPACTION_FANOUT
                            ? SIZE_MAX
                            : candidates[0].record_count * GEO_BACKGROUND_COMPACTION_FANOUT;

    while (selection_count < maximum_selection && candidates[selection_count].record_count <= tier_limit) {
        selection_count++;
    }

    if (selection_count < 2) {
        selection_count = 2;
    }

    size_t *selected_indices = malloc(selection_count * sizeof(*selected_indices));

    if (!selected_indices) {
        free(candidates);
        pthread_mutex_unlock(&set->mutation_lock);

        return false;
    }

    for (size_t i = 0; i < selection_count; ++i) {
        selected_indices[i] = candidates[i].index;
    }

    qsort(selected_indices, selection_count, sizeof(*selected_indices), segment_index_compare);

    set->compaction_running = true;
    bool succeeded = segment_set_compact_selection_locked(set,
                                                          output_path,
                                                          selected_indices,
                                                          selection_count,
                                                          0,
                                                          stats);
    set->compaction_running = false;

    free(selected_indices);
    free(candidates);
    pthread_mutex_unlock(&set->mutation_lock);

    return succeeded;
}
#endif

static bool segment_set_compact_with_worker_limit(GeoSegmentSet *set,
                                                  const char *output_path,
                                                  size_t requested_workers,
                                                  GeoSegmentCompactionStats *stats)
{
    if (stats) {
        memset(stats, 0, sizeof(*stats));
    }

#if GEO_SEGMENTS_SUPPORTED
    if (!set || !output_path || !output_path[0] || access(output_path, F_OK) == 0 ||
        pthread_mutex_lock(&set->mutation_lock) != 0) {
        return false;
    }

    if (!set->count || set->generation >= GEO_MUTATION_MAX_GENERATION || set->count > SIZE_MAX / sizeof(size_t) ||
        set->compaction_running) {
        pthread_mutex_unlock(&set->mutation_lock);

        return false;
    }

    size_t *selected_indices = malloc(set->count * sizeof(*selected_indices));

    if (!selected_indices) {
        pthread_mutex_unlock(&set->mutation_lock);

        return false;
    }

    for (size_t i = 0; i < set->count; ++i) {
        selected_indices[i] = i;
    }

    set->compaction_running = true;
    bool succeeded = segment_set_compact_selection_locked(set,
                                                          output_path,
                                                          selected_indices,
                                                          set->count,
                                                          requested_workers,
                                                          stats);
    set->compaction_running = false;

    free(selected_indices);
    pthread_mutex_unlock(&set->mutation_lock);

    return succeeded;
#else
    (void) set;
    (void) output_path;
    (void) requested_workers;

    return false;
#endif
}

bool geo_segment_set_compact(GeoSegmentSet *set,
                             const char *output_path,
                             GeoSegmentCompactionStats *stats)
{
    return segment_set_compact_with_worker_limit(set, output_path, 0, stats);
}

bool geo_segment_set_compact_with_workers(GeoSegmentSet *set,
                                          const char *output_path,
                                          size_t worker_count,
                                          GeoSegmentCompactionStats *stats)
{
    if (!worker_count || worker_count > GEO_COMPACTION_MAX_WORKERS) {
        if (stats) {
            memset(stats, 0, sizeof(*stats));
        }

        return false;
    }

    return segment_set_compact_with_worker_limit(set, output_path, worker_count, stats);
}
