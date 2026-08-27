#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "geobolt/geobolt.h"

#include "geo_index_io.h"
#include "geo_index_persistence.h"
#include "geo_wal.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#define GEO_DATABASE_STATE_VERSION 1U
#define GEO_DATABASE_WAL_BATCH_VERSION 1U
#define GEO_DATABASE_DEFAULT_WAL_SEGMENT_SIZE (64U * 1024U * 1024U)
#define GEO_DATABASE_DEFAULT_CHECKPOINT_OPERATIONS (256U * 1024U)
#define GEO_DATABASE_DEFAULT_MAX_SEGMENTS 16U
#define GEO_DATABASE_MAX_WAL_PAYLOAD_SIZE (64U * 1024U * 1024U)

typedef struct {
    char magic[8];
    uint32_t version;
    uint32_t endian_marker;
    uint32_t header_size;
    uint32_t reserved32;
    uint64_t applied_next_sequence;
    uint64_t committed_operations;
    uint64_t checksum;
    uint8_t padding[16];
} GeoDatabaseState;

typedef struct {
    uint32_t version;
    uint32_t header_size;
    uint32_t operation_size;
    uint32_t operation_count;
    uint64_t first_sequence;
    uint64_t reserved;
} GeoDatabaseWalBatch;

typedef struct {
    uint64_t object_id;
    uint64_t morton_code;
    uint32_t operation;
    uint32_t reserved32;
    uint64_t reserved64;
} GeoDatabaseWalOperation;

_Static_assert(sizeof(GeoDatabaseState) == 64U, "Database state layout is persisted");
_Static_assert(offsetof(GeoDatabaseState, applied_next_sequence) == 24U, "Database state sequence offset is persisted");
_Static_assert(offsetof(GeoDatabaseState, checksum) == 40U, "Database state checksum offset is persisted");
_Static_assert(sizeof(GeoDatabaseWalBatch) == 32U, "Database WAL batch layout is persisted");
_Static_assert(sizeof(GeoDatabaseWalOperation) == 32U, "Database WAL operation layout is persisted");
_Static_assert(offsetof(GeoDatabaseWalOperation, operation) == 16U, "Database WAL operation type offset is persisted");

struct GeoDatabase {
    char *directory;
    char *state_path;
    char *manifest_path;
    char *segments_directory;
    char *wal_directory;
    GeoSegmentSet *segments;
    GeoWal *wal;
    pthread_rwlock_t visibility_lock;
    pthread_mutex_t writer_lock;
    GeoDatabaseConfig config;
    uint64_t applied_next_sequence;
    uint64_t committed_operations;
    uint64_t recovered_operations;
    uint64_t operations_since_checkpoint;
    uint64_t checkpoints;
    uint64_t maintenance_failures;
    bool visibility_lock_initialized;
    bool writer_lock_initialized;
    bool failed;
};

static const char GEO_DATABASE_STATE_MAGIC[8] = { 'G', 'B', 'D', 'B', 'S', '1', '\0', '\0' };

static uint64_t database_state_checksum(const GeoDatabaseState *state)
{
    GeoDatabaseState copy = *state;

    copy.checksum = 0;

    return geo_persisted_checksum_update(geo_persisted_checksum_initial(), &copy, sizeof(copy));
}

static char *database_join_path(const char *directory, const char *name)
{
    size_t directory_length = strlen(directory);
    size_t name_length = strlen(name);
    bool needs_separator = directory_length && directory[directory_length - 1U] != '/';
    size_t separator_length = needs_separator ? 1U : 0U;

    if (name_length == SIZE_MAX ||
        directory_length > SIZE_MAX - separator_length ||
        directory_length + separator_length > SIZE_MAX - name_length - 1U) {
        return NULL;
    }

    size_t path_length = directory_length + separator_length + name_length;
    char *path = malloc(path_length + 1U);

    if (!path) {
        return NULL;
    }

    memcpy(path, directory, directory_length);

    if (needs_separator) {
        path[directory_length] = '/';
    }

    memcpy(path + directory_length + separator_length, name, name_length + 1U);

    return path;
}

static bool database_ensure_directory(const char *path, bool create_if_missing, bool *created)
{
    struct stat status;

    *created = false;

    if (stat(path, &status) == 0) {
        return S_ISDIR(status.st_mode);
    }

    if (errno != ENOENT || !create_if_missing || mkdir(path, S_IRWXU) != 0) {
        return false;
    }

    *created = true;

    return geo_io_sync_parent_directory(path);
}

static bool database_write_state(GeoDatabase *database, uint64_t applied_next_sequence, uint64_t committed_operations)
{
    GeoDatabaseState state = {
        .version = GEO_DATABASE_STATE_VERSION,
        .endian_marker = GEO_FILE_ENDIAN_MARKER,
        .header_size = sizeof(GeoDatabaseState),
        .applied_next_sequence = applied_next_sequence,
        .committed_operations = committed_operations,
    };

    memcpy(state.magic, GEO_DATABASE_STATE_MAGIC, sizeof(state.magic));
    state.checksum = database_state_checksum(&state);

    char *temporary_path = NULL;
    FILE *file = geo_io_create_atomic_file(database->state_path, &temporary_path);
    bool succeeded = file && fwrite(&state, sizeof(state), 1U, file) == 1U;

    if (succeeded) {
        succeeded = geo_io_publish_atomic_file(file, temporary_path, database->state_path);
        file = NULL;
    }

    if (file) {
        geo_io_discard_atomic_file(file, temporary_path);
    }

    free(temporary_path);

    return succeeded;
}

static bool database_read_state(GeoDatabase *database)
{
    FILE *file = fopen(database->state_path, "rb");

    if (!file) {
        return false;
    }

    GeoDatabaseState state;
    bool succeeded = fread(&state, sizeof(state), 1U, file) == 1U &&
                     fgetc(file) == EOF &&
                     !ferror(file) &&
                     memcmp(state.magic, GEO_DATABASE_STATE_MAGIC, sizeof(state.magic)) == 0 &&
                     state.version == GEO_DATABASE_STATE_VERSION &&
                     state.endian_marker == GEO_FILE_ENDIAN_MARKER &&
                     state.header_size == sizeof(state) &&
                     state.reserved32 == 0 &&
                     state.applied_next_sequence > 0 &&
                     state.committed_operations <= state.applied_next_sequence - 1U &&
                     state.checksum == database_state_checksum(&state);

    if (fclose(file) != 0) {
        succeeded = false;
    }

    if (succeeded) {
        database->applied_next_sequence = state.applied_next_sequence;
        database->committed_operations = state.committed_operations;
    }

    return succeeded;
}

static bool database_segment_is_published(const GeoDatabase *database, const char *path)
{
    char *canonical_path = realpath(path, NULL);

    if (!canonical_path) {
        return false;
    }

    size_t path_size = strlen(canonical_path) + 1U;
    char *active_path = malloc(path_size);

    if (!active_path) {
        free(canonical_path);

        return false;
    }

    size_t segment_count = geo_segment_set_count(database->segments);
    bool published = false;

    for (size_t segment = 0; !published && segment < segment_count; ++segment) {
        published = geo_segment_set_copy_path(database->segments, segment, active_path, path_size) &&
                    strcmp(active_path, canonical_path) == 0;
    }

    free(active_path);
    free(canonical_path);

    return published;
}

static char *database_segment_path(const GeoDatabase *database, uint64_t first_sequence)
{
    char name[64];
    int length = snprintf(name, sizeof(name), "write-%020" PRIu64 ".gbi", first_sequence);

    if (length < 0 || (size_t) length >= sizeof(name)) {
        return NULL;
    }

    return database_join_path(database->segments_directory, name);
}

static bool database_build_segment(const char *path,
                                   const GeoDatabaseWalOperation *operations,
                                   uint32_t operation_count,
                                   size_t upsert_count)
{
    if (!upsert_count) {
        return true;
    }

    GeoIndex *index = geo_index_create(upsert_count);
    GeoRecord *records = malloc(upsert_count * sizeof(*records));

    if (!index || !records) {
        free(records);
        geo_index_destroy(index);

        return false;
    }

    size_t record_count = 0;

    for (uint32_t operation = 0; operation < operation_count; ++operation) {
        if (operations[operation].operation != GEO_DATABASE_UPSERT) {
            continue;
        }

        records[record_count++] = (GeoRecord) {
            .id = operations[operation].object_id,
            .z = operations[operation].morton_code,
        };
    }

    bool succeeded = record_count == upsert_count &&
                     geo_index_add_records(index, records, record_count) &&
                     geo_index_build(index) &&
                     geo_index_save(index, path);

    free(records);
    geo_index_destroy(index);

    return succeeded;
}

static bool database_apply_operations(GeoDatabase *database,
                                      uint64_t first_sequence,
                                      const GeoDatabaseWalOperation *operations,
                                      uint32_t operation_count,
                                      bool segment_prebuilt)
{
    size_t upsert_count = 0;
    size_t delete_count = 0;

    for (uint32_t operation = 0; operation < operation_count; ++operation) {
        if (operations[operation].operation == GEO_DATABASE_UPSERT) {
            upsert_count++;
        } else {
            delete_count++;
        }
    }

    char *segment_path = upsert_count ? database_segment_path(database, first_sequence) : NULL;
    uint64_t *delete_ids = delete_count ? malloc(delete_count * sizeof(*delete_ids)) : NULL;
    bool succeeded = (!upsert_count || segment_path) && (!delete_count || delete_ids);

    for (uint32_t operation = 0, delete_index = 0; succeeded && operation < operation_count; ++operation) {
        if (operations[operation].operation == GEO_DATABASE_DELETE) {
            delete_ids[delete_index++] = operations[operation].object_id;
        }
    }

    if (succeeded && upsert_count && !segment_prebuilt) {
        succeeded = database_build_segment(segment_path, operations, operation_count, upsert_count);
    }

    bool segment_published = succeeded && upsert_count && database_segment_is_published(database, segment_path);

    if (succeeded && pthread_rwlock_wrlock(&database->visibility_lock) != 0) {
        succeeded = false;
    } else if (succeeded) {
        if (delete_count) {
            succeeded = geo_segment_set_remove_ids(database->segments, delete_ids, delete_count);
        }

        if (succeeded && upsert_count && !segment_published) {
            succeeded = geo_segment_set_upsert_file(database->segments, segment_path);
        }

        if (pthread_rwlock_unlock(&database->visibility_lock) != 0) {
            succeeded = false;
        }
    }

    free(delete_ids);
    free(segment_path);

    return succeeded;
}

static uint64_t database_hash_id(uint64_t id)
{
    id ^= id >> 30U;
    id *= UINT64_C(0xbf58476d1ce4e5b9);
    id ^= id >> 27U;
    id *= UINT64_C(0x94d049bb133111eb);

    return id ^ (id >> 31U);
}

static GeoDatabaseStatus database_validate_mutations(const GeoDatabaseMutation *mutations,
                                                     size_t mutation_count)
{
    size_t maximum_operations = (GEO_DATABASE_MAX_WAL_PAYLOAD_SIZE - sizeof(GeoDatabaseWalBatch)) /
                                sizeof(GeoDatabaseWalOperation);

    if (!mutations || !mutation_count || mutation_count > maximum_operations) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    size_t capacity = 2U;

    while (capacity < mutation_count + mutation_count / 2U) {
        if (capacity > SIZE_MAX / 2U) {
            return GEO_DATABASE_INVALID_ARGUMENT;
        }

        capacity *= 2U;
    }

    uint64_t *ids = calloc(capacity, sizeof(*ids));

    if (!ids) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    size_t mask = capacity - 1U;

    for (size_t index = 0; index < mutation_count; ++index) {
        bool operation_valid = mutations[index].operation == GEO_DATABASE_UPSERT ||
                               mutations[index].operation == GEO_DATABASE_DELETE;

        if (!mutations[index].object_id || !operation_valid || mutations[index].reserved != 0) {
            free(ids);

            return GEO_DATABASE_INVALID_ARGUMENT;
        }

        size_t slot = (size_t) database_hash_id(mutations[index].object_id) & mask;

        while (ids[slot] && ids[slot] != mutations[index].object_id) {
            slot = (slot + 1U) & mask;
        }

        if (ids[slot]) {
            free(ids);

            return GEO_DATABASE_INVALID_ARGUMENT;
        }

        ids[slot] = mutations[index].object_id;
    }

    free(ids);

    return GEO_DATABASE_OK;
}

static bool database_decode_wal_batch(uint64_t first_sequence,
                                      uint32_t entry_count,
                                      const void *payload,
                                      uint32_t payload_size,
                                      const GeoDatabaseWalOperation **operations)
{
    if (!payload || payload_size < sizeof(GeoDatabaseWalBatch)) {
        return false;
    }

    const GeoDatabaseWalBatch *batch = payload;
    size_t operations_size = (size_t) entry_count * sizeof(GeoDatabaseWalOperation);
    size_t expected_size = sizeof(*batch) + operations_size;

    if (batch->version != GEO_DATABASE_WAL_BATCH_VERSION ||
        batch->header_size != sizeof(*batch) ||
        batch->operation_size != sizeof(GeoDatabaseWalOperation) ||
        batch->operation_count != entry_count ||
        batch->first_sequence != first_sequence ||
        batch->reserved != 0 ||
        expected_size != payload_size) {
        return false;
    }

    const GeoDatabaseWalOperation *decoded = (const GeoDatabaseWalOperation *) (batch + 1);

    for (uint32_t operation = 0; operation < entry_count; ++operation) {
        if (!decoded[operation].object_id ||
            (decoded[operation].operation != GEO_DATABASE_UPSERT && decoded[operation].operation != GEO_DATABASE_DELETE) ||
            decoded[operation].reserved32 != 0 ||
            decoded[operation].reserved64 != 0) {
            return false;
        }
    }

    *operations = decoded;

    return true;
}

static bool database_replay_wal(void *context,
                                uint64_t first_sequence,
                                uint32_t entry_count,
                                const void *payload,
                                uint32_t payload_size)
{
    GeoDatabase *database = context;
    const GeoDatabaseWalOperation *operations;

    if (entry_count > UINT64_MAX - first_sequence ||
        !database_decode_wal_batch(first_sequence, entry_count, payload, payload_size, &operations)) {
        return false;
    }

    uint64_t next_sequence = first_sequence + entry_count;

    if (next_sequence <= database->applied_next_sequence) {
        return true;
    }

    if (first_sequence < database->applied_next_sequence ||
        !database_apply_operations(database, first_sequence, operations, entry_count, false)) {
        return false;
    }

    database->applied_next_sequence = next_sequence;
    database->committed_operations += entry_count;
    database->recovered_operations += entry_count;

    return true;
}

GeoDatabaseConfig geo_database_default_config(void)
{
    return (GeoDatabaseConfig) {
        .wal_segment_size = GEO_DATABASE_DEFAULT_WAL_SEGMENT_SIZE,
        .checkpoint_interval_operations = GEO_DATABASE_DEFAULT_CHECKPOINT_OPERATIONS,
        .max_active_segments = GEO_DATABASE_DEFAULT_MAX_SEGMENTS,
        .create_if_missing = true,
    };
}

static GeoDatabaseStatus database_initialize_paths(GeoDatabase *database, const char *directory)
{
    bool root_created;

    if (!database_ensure_directory(directory, database->config.create_if_missing, &root_created)) {
        return GEO_DATABASE_IO_ERROR;
    }

    database->directory = realpath(directory, NULL);

    if (!database->directory) {
        return GEO_DATABASE_IO_ERROR;
    }

    database->state_path = database_join_path(database->directory, "state.gbs");
    database->manifest_path = database_join_path(database->directory, "segments.gbm");
    database->segments_directory = database_join_path(database->directory, "segments");
    database->wal_directory = database_join_path(database->directory, "wal");

    if (!database->state_path || !database->manifest_path || !database->segments_directory || !database->wal_directory) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    bool segments_created;
    bool wal_created;

    if (!database_ensure_directory(database->segments_directory, database->config.create_if_missing, &segments_created) ||
        !database_ensure_directory(database->wal_directory, database->config.create_if_missing, &wal_created)) {
        return GEO_DATABASE_IO_ERROR;
    }

    (void) root_created;

    return GEO_DATABASE_OK;
}

GeoDatabase *geo_database_open(const char *directory,
                               const GeoDatabaseConfig *config,
                               GeoDatabaseStatus *status)
{
    GeoDatabaseStatus open_status = GEO_DATABASE_INVALID_ARGUMENT;

    if (!directory || !directory[0]) {
        if (status) {
            *status = open_status;
        }

        return NULL;
    }

    GeoDatabaseConfig selected_config = config ? *config : geo_database_default_config();

    if (selected_config.wal_segment_size < 64U * 1024U ||
        !selected_config.checkpoint_interval_operations ||
        !selected_config.max_active_segments) {
        if (status) {
            *status = open_status;
        }

        return NULL;
    }

    GeoDatabase *database = calloc(1, sizeof(*database));

    if (!database) {
        open_status = GEO_DATABASE_OUT_OF_MEMORY;
        goto complete;
    }

    database->config = selected_config;
    open_status = database_initialize_paths(database, directory);

    if (open_status != GEO_DATABASE_OK) {
        goto complete;
    }

    bool state_exists = access(database->state_path, F_OK) == 0;
    bool manifest_exists = access(database->manifest_path, F_OK) == 0;

    if (state_exists != manifest_exists || (!state_exists && !database->config.create_if_missing)) {
        open_status = state_exists || manifest_exists ? GEO_DATABASE_CORRUPTION : GEO_DATABASE_IO_ERROR;
        goto complete;
    }

    if (pthread_rwlock_init(&database->visibility_lock, NULL) != 0) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }

    database->visibility_lock_initialized = true;

    if (pthread_mutex_init(&database->writer_lock, NULL) != 0) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }

    database->writer_lock_initialized = true;

    if (state_exists) {
        if (!database_read_state(database)) {
            open_status = GEO_DATABASE_CORRUPTION;
            goto complete;
        }

        database->segments = geo_segment_set_open(database->manifest_path);
    } else {
        database->applied_next_sequence = 1U;
        database->segments = geo_segment_set_create(database->manifest_path);

        if (database->segments && !database_write_state(database, 1U, 0U)) {
            geo_segment_set_destroy(database->segments);
            database->segments = NULL;
        }
    }

    if (!database->segments) {
        open_status = GEO_DATABASE_CORRUPTION;
        goto complete;
    }

    uint64_t wal_next_sequence = 0;
    database->wal = geo_wal_open(database->wal_directory,
                                 database->config.wal_segment_size,
                                 database_replay_wal,
                                 database,
                                 &wal_next_sequence);

    if (!database->wal || wal_next_sequence < database->applied_next_sequence) {
        open_status = GEO_DATABASE_CORRUPTION;
        goto complete;
    }

    if (database->recovered_operations) {
        if (!database_write_state(database, database->applied_next_sequence, database->committed_operations) ||
            !geo_wal_checkpoint(database->wal, wal_next_sequence)) {
            open_status = GEO_DATABASE_IO_ERROR;
            goto complete;
        }

        database->checkpoints++;
    }

    if (!geo_segment_set_enable_background_compaction(database->segments,
                                                       database->segments_directory,
                                                       database->config.max_active_segments)) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }

    open_status = GEO_DATABASE_OK;

complete:
    if (open_status != GEO_DATABASE_OK) {
        geo_database_close(database);
        database = NULL;
    }

    if (status) {
        *status = open_status;
    }

    return database;
}

void geo_database_close(GeoDatabase *database)
{
    if (!database) {
        return;
    }

    geo_wal_close(database->wal);
    geo_segment_set_destroy(database->segments);

    if (database->writer_lock_initialized) {
        (void) pthread_mutex_destroy(&database->writer_lock);
    }

    if (database->visibility_lock_initialized) {
        (void) pthread_rwlock_destroy(&database->visibility_lock);
    }

    free(database->wal_directory);
    free(database->segments_directory);
    free(database->manifest_path);
    free(database->state_path);
    free(database->directory);
    free(database);
}

static GeoDatabaseStatus database_checkpoint_locked(GeoDatabase *database)
{
    if (!geo_segment_set_checkpoint_mutations(database->segments) ||
        !database_write_state(database, database->applied_next_sequence, database->committed_operations) ||
        !geo_wal_checkpoint(database->wal, database->applied_next_sequence)) {
        return GEO_DATABASE_IO_ERROR;
    }

    database->operations_since_checkpoint = 0;
    database->checkpoints++;

    return GEO_DATABASE_OK;
}

GeoDatabaseStatus geo_database_write(GeoDatabase *database,
                                     const GeoDatabaseMutation *mutations,
                                     size_t mutation_count)
{
    if (!database) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoDatabaseStatus status = database_validate_mutations(mutations, mutation_count);

    if (status != GEO_DATABASE_OK) {
        return status;
    }

    if (mutation_count > (UINT32_MAX - sizeof(GeoDatabaseWalBatch)) / sizeof(GeoDatabaseWalOperation)) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    size_t payload_size = sizeof(GeoDatabaseWalBatch) + mutation_count * sizeof(GeoDatabaseWalOperation);
    GeoDatabaseWalBatch *batch = malloc(payload_size);

    if (!batch) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    if (pthread_mutex_lock(&database->writer_lock) != 0) {
        free(batch);

        return GEO_DATABASE_FAILED_STATE;
    }

    if (database->failed || mutation_count > UINT64_MAX - database->applied_next_sequence) {
        status = GEO_DATABASE_FAILED_STATE;
        goto complete;
    }

    uint64_t first_sequence = database->applied_next_sequence;
    *batch = (GeoDatabaseWalBatch) {
        .version = GEO_DATABASE_WAL_BATCH_VERSION,
        .header_size = sizeof(*batch),
        .operation_size = sizeof(GeoDatabaseWalOperation),
        .operation_count = (uint32_t) mutation_count,
        .first_sequence = first_sequence,
    };

    GeoDatabaseWalOperation *operations = (GeoDatabaseWalOperation *) (batch + 1);
    size_t upsert_count = 0;

    for (size_t operation = 0; operation < mutation_count; ++operation) {
        operations[operation] = (GeoDatabaseWalOperation) {
            .object_id = mutations[operation].object_id,
            .morton_code = mutations[operation].morton_code,
            .operation = (uint32_t) mutations[operation].operation,
        };
        upsert_count += mutations[operation].operation == GEO_DATABASE_UPSERT;
    }

    char *segment_path = upsert_count ? database_segment_path(database, first_sequence) : NULL;

    if (upsert_count && (!segment_path || !database_build_segment(segment_path,
                                                                  operations,
                                                                  (uint32_t) mutation_count,
                                                                  upsert_count))) {
        free(segment_path);
        status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }

    free(segment_path);

    bool wal_committed = geo_wal_append(database->wal,
                                        first_sequence,
                                        (uint32_t) mutation_count,
                                        batch,
                                        (uint32_t) payload_size) &&
                         geo_wal_sync(database->wal);

    if (!wal_committed) {
        database->failed = true;
        status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }

    uint64_t next_sequence = first_sequence + mutation_count;

    if (!database_apply_operations(database,
                                   first_sequence,
                                   operations,
                                   (uint32_t) mutation_count,
                                   true) ||
        !database_write_state(database,
                              next_sequence,
                              database->committed_operations + mutation_count)) {
        database->failed = true;
        status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }

    database->applied_next_sequence = next_sequence;
    database->committed_operations += mutation_count;
    database->operations_since_checkpoint += mutation_count;
    status = GEO_DATABASE_OK;

    if (database->operations_since_checkpoint >= database->config.checkpoint_interval_operations) {
        GeoDatabaseStatus checkpoint_status = database_checkpoint_locked(database);

        if (checkpoint_status != GEO_DATABASE_OK) {
            // The commit and its state watermark are already durable. Retaining the WAL is safe and recovery will retry reclamation.
            database->maintenance_failures++;
        }
    }

complete:
    pthread_mutex_unlock(&database->writer_lock);
    free(batch);

    return status;
}

GeoDatabaseStatus geo_database_upsert(GeoDatabase *database,
                                      uint64_t object_id,
                                      double latitude,
                                      double longitude)
{
    if (!isfinite(latitude) || !isfinite(longitude) ||
        latitude < GEO_MIN_LAT || latitude > GEO_MAX_LAT ||
        longitude < GEO_MIN_LNG || longitude > GEO_MAX_LNG) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoDatabaseMutation mutation = {
        .object_id = object_id,
        .morton_code = geo_encode(latitude, longitude),
        .operation = GEO_DATABASE_UPSERT,
    };

    return geo_database_write(database, &mutation, 1U);
}

GeoDatabaseStatus geo_database_remove(GeoDatabase *database, uint64_t object_id)
{
    GeoDatabaseMutation mutation = {
        .object_id = object_id,
        .operation = GEO_DATABASE_DELETE,
    };

    return geo_database_write(database, &mutation, 1U);
}

bool geo_database_search_radius_reuse(const GeoDatabase *database,
                                      double latitude,
                                      double longitude,
                                      double radius_km,
                                      GeoSearchResult *result,
                                      GeoSearchStats *stats)
{
    if (!database || !result) {
        return false;
    }

    GeoDatabase *mutable_database = (GeoDatabase *) database;

    if (pthread_rwlock_rdlock(&mutable_database->visibility_lock) != 0) {
        return false;
    }

    bool succeeded = geo_segment_set_search_radius_reuse(database->segments,
                                                          latitude,
                                                          longitude,
                                                          radius_km,
                                                          result,
                                                          stats);

    if (pthread_rwlock_unlock(&mutable_database->visibility_lock) != 0) {
        succeeded = false;
    }

    return succeeded;
}

bool geo_database_search_radius_count(const GeoDatabase *database,
                                      double latitude,
                                      double longitude,
                                      double radius_km,
                                      size_t *count,
                                      GeoSearchStats *stats)
{
    if (!database || !count) {
        return false;
    }

    GeoDatabase *mutable_database = (GeoDatabase *) database;

    if (pthread_rwlock_rdlock(&mutable_database->visibility_lock) != 0) {
        return false;
    }

    bool succeeded = geo_segment_set_search_radius_count(database->segments,
                                                          latitude,
                                                          longitude,
                                                          radius_km,
                                                          count,
                                                          stats);

    if (pthread_rwlock_unlock(&mutable_database->visibility_lock) != 0) {
        succeeded = false;
    }

    return succeeded;
}

GeoDatabaseStatus geo_database_checkpoint(GeoDatabase *database)
{
    if (!database) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    if (pthread_mutex_lock(&database->writer_lock) != 0) {
        return GEO_DATABASE_FAILED_STATE;
    }

    GeoDatabaseStatus status = database->failed ? GEO_DATABASE_FAILED_STATE : database_checkpoint_locked(database);

    pthread_mutex_unlock(&database->writer_lock);

    return status;
}

bool geo_database_get_stats(const GeoDatabase *database, GeoDatabaseStats *stats)
{
    if (!database || !stats) {
        return false;
    }

    GeoDatabase *mutable_database = (GeoDatabase *) database;

    if (pthread_mutex_lock(&mutable_database->writer_lock) != 0) {
        return false;
    }

    *stats = (GeoDatabaseStats) {
        .committed_operations = database->committed_operations,
        .recovered_operations = database->recovered_operations,
        .checkpoints = database->checkpoints,
        .maintenance_failures = database->maintenance_failures,
        .next_sequence = database->applied_next_sequence,
        .physical_records = geo_segment_set_record_count(database->segments),
        .active_segments = geo_segment_set_count(database->segments),
    };

    pthread_mutex_unlock(&mutable_database->writer_lock);

    return true;
}
