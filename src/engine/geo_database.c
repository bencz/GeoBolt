#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "geo_database_internal.h"

#include "geobolt/geodoc.h"
#include "geo_db_format.h"
#include "geo_index_private.h"
#include "geo_index_io.h"
#include "geo_rocks_bridge.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#define GEO_DATABASE_DEFAULT_BLOCK_CACHE_SIZE (256U * 1024U * 1024U)
#define GEO_DATABASE_DEFAULT_OBJECT_CACHE_SIZE (128U * 1024U * 1024U)
#define GEO_DATABASE_DEFAULT_WRITE_BUFFER_SIZE (64U * 1024U * 1024U)
#define GEO_DATABASE_DEFAULT_MAX_SEGMENTS 16U
#define GEO_DATABASE_DEFAULT_GROUP_COMMIT_OPERATIONS (64U * 1024U)
#define GEO_DATABASE_DEFAULT_SPATIAL_MEMTABLE_OPERATIONS (256U * 1024U)
#define GEO_DATABASE_DEFAULT_GROUP_COMMIT_DELAY_US 25U
#define GEO_DATABASE_REPLAY_BATCH_SIZE (256U * 1024U)

static char *database_join_path(const char *directory, const char *name)
{
    size_t directory_length = strlen(directory);
    size_t name_length = strlen(name);
    bool needs_separator = directory_length && directory[directory_length - 1U] != '/';
    size_t separator_length = needs_separator ? 1U : 0U;

    if (name_length == SIZE_MAX || directory_length > SIZE_MAX - separator_length ||
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

static bool database_ensure_directory(const char *path, bool create_if_missing)
{
    struct stat status;

    if (stat(path, &status) == 0) {
        return S_ISDIR(status.st_mode);
    }

    if (errno != ENOENT || !create_if_missing || mkdir(path, S_IRWXU) != 0) {
        return false;
    }

    return geo_io_sync_parent_directory(path);
}

static GeoDatabaseStatus database_initialize_paths(GeoDatabase *database, const char *directory)
{
    if (!database_ensure_directory(directory, database->config.create_if_missing)) {
        return GEO_DATABASE_IO_ERROR;
    }

    database->directory = realpath(directory, NULL);

    if (!database->directory) {
        return GEO_DATABASE_IO_ERROR;
    }

    database->rocks_directory = database_join_path(database->directory, "rocksdb");
    database->manifest_path = database_join_path(database->directory, "spatial-manifest.gbm");
    database->segments_directory = database_join_path(database->directory, "spatial-segments");

    if (!database->rocks_directory || !database->manifest_path || !database->segments_directory) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    return database_ensure_directory(database->segments_directory, true) ? GEO_DATABASE_OK : GEO_DATABASE_IO_ERROR;
}

static uint64_t database_hash_id(uint64_t id)
{
    id ^= id >> 30U;
    id *= UINT64_C(0xbf58476d1ce4e5b9);
    id ^= id >> 27U;
    id *= UINT64_C(0x94d049bb133111eb);
    return id ^ (id >> 31U);
}

static void database_fingerprint_bytes(uint64_t *low, uint64_t *high, const void *data, size_t size)
{
    const unsigned char *bytes = data;

    for (size_t index = 0; index < size; ++index) {
        *low = (*low ^ bytes[index]) * UINT64_C(0x100000001b3);
        *high = (*high + bytes[index] + UINT64_C(0x9e3779b97f4a7c15)) * UINT64_C(0xbf58476d1ce4e5b9);
        *high = *high << 27U | *high >> 37U;
    }
}

static void database_fingerprint_u64(uint64_t *low, uint64_t *high, uint64_t value)
{
    unsigned char encoded[8];

    for (size_t index = 0; index < sizeof(encoded); ++index) {
        encoded[index] = (unsigned char) (value >> (index * 8U));
    }

    database_fingerprint_bytes(low, high, encoded, sizeof(encoded));
}

static GeoDbIdempotencyRecord database_idempotency_record(const GeoDatabaseMutationView *view,
                                                          size_t mutation_count,
                                                          uint64_t first_sequence,
                                                          uint64_t expires_at_seconds)
{
    uint64_t low = UINT64_C(0xcbf29ce484222325);
    uint64_t high = UINT64_C(0x6eed0e9da4d94a4f);

    for (size_t index = 0; index < mutation_count; ++index) {
        GeoDatabaseObjectMutation mutation = geo_database_mutation_at(view, index);
        database_fingerprint_u64(&low, &high, mutation.object_id);
        database_fingerprint_u64(&low, &high, mutation.morton_code);
        database_fingerprint_u64(&low, &high, mutation.document_size);
        database_fingerprint_u64(&low, &high, mutation.operation);

        if (mutation.document_size) {
            database_fingerprint_bytes(&low, &high, mutation.document, mutation.document_size);
        }
    }

    return (GeoDbIdempotencyRecord) {
        .fingerprint_low = low ^ (high >> 29U),
        .fingerprint_high = high ^ (low << 31U),
        .expires_at_seconds = expires_at_seconds,
        .first_sequence = first_sequence,
        .operation_count = (uint32_t) mutation_count,
    };
}

static bool database_realtime_seconds(uint64_t *seconds)
{
    struct timespec now;

    if (clock_gettime(CLOCK_REALTIME, &now) != 0 || now.tv_sec < 0) {
        return false;
    }

    *seconds = (uint64_t) now.tv_sec;
    return true;
}

GeoDatabaseObjectMutation geo_database_mutation_at(const GeoDatabaseMutationView *view, size_t index)
{
    if (view->layout == GEO_DATABASE_MUTATIONS_OBJECT) {
        return ((const GeoDatabaseObjectMutation *) view->mutations)[index];
    }

    const GeoDatabaseMutation *mutation = (const GeoDatabaseMutation *) view->mutations + index;

    return (GeoDatabaseObjectMutation) {
        .object_id = mutation->object_id,
        .morton_code = mutation->morton_code,
        .operation = mutation->operation,
        .reserved = mutation->reserved,
    };
}

GeoDatabaseStatus geo_database_validate_mutations(const GeoDatabaseMutationView *view, size_t mutation_count)
{
    if (!view || !view->mutations || !mutation_count || mutation_count > UINT32_MAX) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    bool has_idempotency = view->idempotency_key_size != 0U;

    if ((has_idempotency && (!view->idempotency_key ||
                             view->idempotency_key_size > GEO_DATABASE_MAX_IDEMPOTENCY_KEY_SIZE ||
                             !view->idempotency_retention_seconds ||
                             view->idempotency_retention_seconds > GEO_DATABASE_MAX_IDEMPOTENCY_RETENTION_SECONDS)) ||
        (!has_idempotency && (view->idempotency_key || view->idempotency_retention_seconds))) {
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
    size_t generated_count = 0U;
    GeoDatabaseStatus status = GEO_DATABASE_OK;

    for (size_t index = 0; index < mutation_count; ++index) {
        GeoDatabaseObjectMutation mutation = geo_database_mutation_at(view, index);
        bool generated = mutation.operation == GEO_DATABASE_INTERNAL_INSERT_GENERATED;
        bool operation_valid = (mutation.operation >= GEO_DATABASE_UPSERT && mutation.operation <= GEO_DATABASE_UPDATE) || generated;
        bool document_valid = mutation.operation != GEO_DATABASE_DELETE
                                  ? (!mutation.document_size || mutation.document)
                                  : (!mutation.document && mutation.document_size == 0U && mutation.morton_code == 0U);

        generated_count += generated;

        if ((generated ? mutation.object_id != 0U : mutation.object_id == 0U) ||
            (generated && view->layout != GEO_DATABASE_MUTATIONS_OBJECT) || generated_count > 1U || !operation_valid ||
            mutation.reserved != 0U || !document_valid ||
            mutation.document_size > UINT32_MAX) {
            status = GEO_DATABASE_INVALID_ARGUMENT;
            break;
        }

        if (mutation.document_size) {
            GeoDocView document;

            if (geo_doc_open(mutation.document, mutation.document_size, &document) != GEO_DOC_OK) {
                status = GEO_DATABASE_INVALID_ARGUMENT;
                break;
            }
        }

        if (generated) {
            continue;
        }

        size_t slot = (size_t) database_hash_id(mutation.object_id) & mask;

        while (ids[slot] && ids[slot] != mutation.object_id) {
            slot = (slot + 1U) & mask;
        }

        if (ids[slot]) {
            status = GEO_DATABASE_INVALID_ARGUMENT;
            break;
        }

        ids[slot] = mutation.object_id;
    }

    free(ids);
    return status;
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

static bool database_build_segment(const char *path,
                                   const GeoDatabaseMutation *operations,
                                   size_t operation_count,
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

    size_t record_count = 0U;

    for (size_t operation = 0; operation < operation_count; ++operation) {
        if (operations[operation].operation == GEO_DATABASE_UPSERT) {
            records[record_count++] = (GeoRecord) {
                .id = operations[operation].object_id,
                .z = operations[operation].morton_code,
            };
        }
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
                                      uint64_t next_sequence,
                                      const GeoDatabaseMutation *operations,
                                      size_t operation_count,
                                      bool segment_prebuilt,
                                      bool visibility_locked)
{
    if (next_sequence <= first_sequence || operation_count > next_sequence - first_sequence) {
        return false;
    }

    if (geo_segment_set_durable_watermark(database->segments) >= next_sequence) {
        return true;
    }

    size_t upsert_count = 0U;
    size_t delete_count = 0U;

    for (size_t operation = 0; operation < operation_count; ++operation) {
        upsert_count += operations[operation].operation == GEO_DATABASE_UPSERT;
        delete_count += operations[operation].operation == GEO_DATABASE_DELETE;
    }

    char *segment_path = upsert_count ? database_segment_path(database, first_sequence) : NULL;
    uint64_t *delete_ids = delete_count ? malloc(delete_count * sizeof(*delete_ids)) : NULL;
    bool succeeded = (!upsert_count || segment_path) && (!delete_count || delete_ids);

    for (size_t operation = 0U, delete_index = 0U; succeeded && operation < operation_count; ++operation) {
        if (operations[operation].operation == GEO_DATABASE_DELETE) {
            delete_ids[delete_index++] = operations[operation].object_id;
        }
    }

    if (succeeded && upsert_count && !segment_prebuilt) {
        succeeded = database_build_segment(segment_path, operations, operation_count, upsert_count);
    }

    bool segment_published = succeeded && upsert_count && database_segment_is_published(database, segment_path);

    if (succeeded && !visibility_locked && !geo_visibility_gate_write_lock(&database->visibility_gate)) {
        succeeded = false;
    } else if (succeeded) {
        if (delete_count) {
            uint64_t delete_watermark = upsert_count ? geo_segment_set_durable_watermark(database->segments) : next_sequence;

            succeeded = geo_segment_set_remove_ids_at_watermark(database->segments,
                                                                 delete_ids,
                                                                 delete_count,
                                                                 delete_watermark);
        }

        if (succeeded && upsert_count && !segment_published) {
            succeeded = geo_segment_set_upsert_file_at_watermark(database->segments, segment_path, next_sequence);
        }

        if (!visibility_locked && !geo_visibility_gate_write_unlock(&database->visibility_gate)) {
            succeeded = false;
        }
    }

    free(delete_ids);
    free(segment_path);
    return succeeded;
}

static bool database_store_catalog(GeoDatabase *database, bool synchronize)
{
    GeoRocksStatus rocks_status;
    GeoRocksBatch *batch = geo_rocks_batch_create(database->rocks, 128U, &rocks_status);
    GeoDatabaseStatus status = GEO_DATABASE_OK;
    unsigned char first_delta[8] = { 0 };
    unsigned char first_unapplied_delta[8];

    geo_db_store_u64_be(first_unapplied_delta, database->catalog.spatial_applied_sequence);

    bool succeeded = batch &&
                     geo_db_catalog_put(batch, &database->catalog, &status) &&
                     geo_rocks_batch_delete_range(batch,
                                                  GEO_ROCKS_CF_SPATIAL_DELTA,
                                                  first_delta,
                                                  sizeof(first_delta),
                                                  first_unapplied_delta,
                                                  sizeof(first_unapplied_delta),
                                                  &rocks_status) &&
                     geo_rocks_write(database->rocks, batch, synchronize, &rocks_status);

    geo_rocks_batch_destroy(batch);
    return succeeded;
}

GeoDatabaseStatus geo_database_object_load(GeoDatabase *database,
                                           uint64_t object_id,
                                           GeoRocksBuffer *value,
                                           bool *exists)
{
    unsigned char key[8];
    GeoRocksStatus rocks_status;

    *value = (GeoRocksBuffer) { 0 };
    geo_db_store_u64_be(key, object_id);

    if (geo_rocks_get(database->rocks, GEO_ROCKS_CF_OBJECTS, NULL, key, sizeof(key), value, &rocks_status)) {
        *exists = true;
        return GEO_DATABASE_OK;
    }

    if (rocks_status.code == GEO_ROCKS_NOT_FOUND) {
        *exists = false;
        return GEO_DATABASE_OK;
    }

    return geo_db_status_from_rocks(&rocks_status);
}

static GeoDatabaseStatus database_object_exists(GeoDatabase *database, uint64_t object_id, bool *exists)
{
    GeoRocksBuffer value;
    GeoDatabaseStatus status = geo_database_object_load(database, object_id, &value, exists);

    geo_rocks_buffer_release(&value);
    return status;
}

static bool database_replay_spatial_deltas(GeoDatabase *database)
{
    if (database->catalog.spatial_applied_sequence == database->catalog.next_sequence) {
        return true;
    }

    GeoDatabaseMutation *operations = malloc(GEO_DATABASE_REPLAY_BATCH_SIZE * sizeof(*operations));
    GeoRocksStatus rocks_status;
    GeoRocksIterator *iterator = geo_rocks_iterator_create(database->rocks, GEO_ROCKS_CF_SPATIAL_DELTA, NULL, &rocks_status);

    if (!operations || !iterator) {
        free(operations);
        geo_rocks_iterator_destroy(iterator);
        return false;
    }

    unsigned char seek_key[8];
    geo_db_store_u64_be(seek_key, database->catalog.spatial_applied_sequence);
    geo_rocks_iterator_seek(iterator, seek_key, sizeof(seek_key));
    bool succeeded = true;

    while (succeeded && database->catalog.spatial_applied_sequence < database->catalog.next_sequence) {
        uint64_t first_sequence = database->catalog.spatial_applied_sequence;
        size_t count = 0U;

        while (count < GEO_DATABASE_REPLAY_BATCH_SIZE && geo_rocks_iterator_valid(iterator)) {
            size_t key_size = 0U;
            size_t value_size = 0U;
            const void *key = geo_rocks_iterator_key(iterator, &key_size);
            const void *value = geo_rocks_iterator_value(iterator, &value_size);
            GeoDbSpatialDelta delta;

            if (!geo_db_spatial_delta_decode(key, key_size, value, value_size, &delta) ||
                delta.sequence != first_sequence + count || delta.sequence >= database->catalog.next_sequence) {
                succeeded = false;
                break;
            }

            operations[count++] = (GeoDatabaseMutation) {
                .object_id = delta.object_id,
                .morton_code = delta.morton_code,
                .operation = delta.operation,
            };
            geo_rocks_iterator_next(iterator);
        }

        if (!succeeded || !count ||
            !database_apply_operations(database,
                                       first_sequence,
                                       first_sequence + count,
                                       operations,
                                       count,
                                       false,
                                       false)) {
            succeeded = false;
            break;
        }

        database->catalog.spatial_applied_sequence += count;
        database->recovered_operations += count;
        succeeded = database_store_catalog(database, true);
    }

    if (succeeded) {
        succeeded = geo_rocks_iterator_status(iterator, &rocks_status);
    }

    free(operations);
    geo_rocks_iterator_destroy(iterator);
    return succeeded;
}

static bool database_rebuild_spatial_index(GeoDatabase *database)
{
    GeoDatabaseMutation *operations = malloc(GEO_DATABASE_REPLAY_BATCH_SIZE * sizeof(*operations));
    GeoRocksStatus rocks_status;
    GeoRocksIterator *iterator = geo_rocks_iterator_create(database->rocks, GEO_ROCKS_CF_OBJECTS, NULL, &rocks_status);

    if (!operations || !iterator) {
        free(operations);
        geo_rocks_iterator_destroy(iterator);
        return false;
    }

    geo_rocks_iterator_seek_first(iterator);
    bool succeeded = true;
    uint64_t synthetic_sequence = 1U;

    while (succeeded && geo_rocks_iterator_valid(iterator)) {
        size_t count = 0U;

        while (count < GEO_DATABASE_REPLAY_BATCH_SIZE && geo_rocks_iterator_valid(iterator)) {
            size_t key_size = 0U;
            size_t value_size = 0U;
            const void *key = geo_rocks_iterator_key(iterator, &key_size);
            const void *value = geo_rocks_iterator_value(iterator, &value_size);
            uint64_t sequence;
            uint64_t morton_code;
            const void *document;
            size_t document_size;

            if (!key || key_size != sizeof(uint64_t) ||
                !geo_db_object_decode(value, value_size, &sequence, &morton_code, &document, &document_size)) {
                succeeded = false;
                break;
            }

            operations[count++] = (GeoDatabaseMutation) {
                .object_id = geo_db_load_u64_be(key),
                .morton_code = morton_code,
                .operation = GEO_DATABASE_UPSERT,
            };
            geo_rocks_iterator_next(iterator);
        }

        if (count && !database_apply_operations(database,
                                                synthetic_sequence,
                                                synthetic_sequence + count,
                                                operations,
                                                count,
                                                false,
                                                false)) {
            succeeded = false;
        }

        synthetic_sequence += count;
    }

    if (succeeded) {
        succeeded = geo_rocks_iterator_status(iterator, &rocks_status);
    }

    if (succeeded) {
        succeeded = geo_segment_set_remove_ids_at_watermark(database->segments,
                                                             NULL,
                                                             0U,
                                                             database->catalog.next_sequence);
    }

    if (succeeded) {
        database->catalog.spatial_applied_sequence = database->catalog.next_sequence;
        succeeded = database_store_catalog(database, true);
    }

    free(operations);
    geo_rocks_iterator_destroy(iterator);
    return succeeded;
}

static bool database_remove_derived_manifest(const GeoDatabase *database)
{
    static const char *const suffixes[] = { "", ".mutations", ".mutations.checkpoint" };

    for (size_t index = 0; index < sizeof(suffixes) / sizeof(suffixes[0]); ++index) {
        size_t path_size = strlen(database->manifest_path) + strlen(suffixes[index]) + 1U;
        char *path = malloc(path_size);

        if (!path) {
            return false;
        }

        int length = snprintf(path, path_size, "%s%s", database->manifest_path, suffixes[index]);
        bool removed = length > 0 && (size_t) length < path_size && (unlink(path) == 0 || errno == ENOENT);
        free(path);

        if (!removed) {
            return false;
        }
    }

    return geo_io_sync_parent_directory(database->manifest_path);
}

static bool database_rotate_spatial_memtable_locked(GeoDatabase *database)
{
    if (pthread_mutex_lock(&database->spatial_flush_lock) != 0) {
        return false;
    }

    bool succeeded = !geo_spatial_memtable_has_frozen(database->spatial_memtable) &&
                     geo_spatial_memtable_rotate(database->spatial_memtable);

    if (succeeded) {
        (void) pthread_cond_signal(&database->spatial_flush_ready);
    }

    (void) pthread_mutex_unlock(&database->spatial_flush_lock);
    return succeeded;
}

static bool database_wait_for_frozen_flush(GeoDatabase *database)
{
    if (pthread_mutex_lock(&database->spatial_flush_lock) != 0) {
        return false;
    }

    while (geo_spatial_memtable_has_frozen(database->spatial_memtable) && !database->spatial_flush_failed) {
        if (pthread_cond_wait(&database->spatial_flush_ready, &database->spatial_flush_lock) != 0) {
            database->spatial_flush_failed = true;
        }
    }

    bool succeeded = !database->spatial_flush_failed;
    (void) pthread_mutex_unlock(&database->spatial_flush_lock);
    return succeeded;
}

static bool database_prepare_spatial_memtable(GeoDatabase *database, size_t operation_count)
{
    while (!geo_spatial_memtable_can_append(database->spatial_memtable, operation_count)) {
        if (!geo_visibility_gate_write_lock(&database->visibility_gate)) {
            return false;
        }

        bool has_frozen = geo_spatial_memtable_has_frozen(database->spatial_memtable);
        bool rotated = has_frozen || database_rotate_spatial_memtable_locked(database);
        (void) geo_visibility_gate_write_unlock(&database->visibility_gate);

        if (!rotated || !database_wait_for_frozen_flush(database)) {
            return false;
        }
    }

    return true;
}

static bool database_publish_frozen_memtable(GeoDatabase *database)
{
    GeoDatabaseMutation *operations = NULL;
    size_t operation_count = 0U;
    uint64_t first_sequence = 0U;
    uint64_t next_sequence = 0U;
    bool succeeded = geo_spatial_memtable_export_frozen(database->spatial_memtable,
                                                        &operations,
                                                        &operation_count,
                                                        &first_sequence,
                                                        &next_sequence);
    size_t upsert_count = 0U;

    for (size_t index = 0U; succeeded && index < operation_count; ++index) {
        upsert_count += operations[index].operation == GEO_DATABASE_UPSERT;
    }

    char *segment_path = succeeded && upsert_count ? database_segment_path(database, first_sequence) : NULL;

    if (succeeded && upsert_count) {
        succeeded = segment_path && database_build_segment(segment_path, operations, operation_count, upsert_count);
    }

    if (succeeded && geo_visibility_gate_write_lock(&database->visibility_gate)) {
        succeeded = database_apply_operations(database,
                                              first_sequence,
                                              next_sequence,
                                              operations,
                                              operation_count,
                                              true,
                                              true);

        if (pthread_mutex_lock(&database->spatial_flush_lock) == 0) {
            if (succeeded) {
                geo_spatial_memtable_release_frozen(database->spatial_memtable);
            } else {
                database->spatial_flush_failed = true;
            }

            database->spatial_flush_busy = false;
            (void) pthread_cond_broadcast(&database->spatial_flush_ready);
            (void) pthread_mutex_unlock(&database->spatial_flush_lock);
        } else {
            succeeded = false;
        }

        (void) geo_visibility_gate_write_unlock(&database->visibility_gate);
    } else if (succeeded) {
        succeeded = false;
    }

    free(segment_path);
    free(operations);
    return succeeded;
}

static void *database_spatial_flush_thread_main(void *argument)
{
    GeoDatabase *database = argument;

    for (;;) {
        if (pthread_mutex_lock(&database->spatial_flush_lock) != 0) {
            atomic_store_explicit(&database->failed, true, memory_order_release);
            return NULL;
        }

        while (!geo_spatial_memtable_has_frozen(database->spatial_memtable) && !database->spatial_flush_stop) {
            if (pthread_cond_wait(&database->spatial_flush_ready, &database->spatial_flush_lock) != 0) {
                database->spatial_flush_failed = true;
                database->spatial_flush_stop = true;
            }
        }

        bool stopped = database->spatial_flush_stop && !geo_spatial_memtable_has_frozen(database->spatial_memtable);
        database->spatial_flush_busy = !stopped;
        (void) pthread_mutex_unlock(&database->spatial_flush_lock);

        if (stopped) {
            return NULL;
        }

        if (!database_publish_frozen_memtable(database)) {
            if (pthread_mutex_lock(&database->spatial_flush_lock) == 0) {
                database->spatial_flush_failed = true;
                database->spatial_flush_busy = false;
                (void) pthread_cond_broadcast(&database->spatial_flush_ready);
                (void) pthread_mutex_unlock(&database->spatial_flush_lock);
            }

            atomic_store_explicit(&database->failed, true, memory_order_release);
            return NULL;
        }
    }
}

static bool database_flush_all_spatial_memtables(GeoDatabase *database)
{
    for (;;) {
        if (!geo_visibility_gate_write_lock(&database->visibility_gate)) {
            return false;
        }

        bool has_active = geo_spatial_memtable_has_active(database->spatial_memtable);
        bool has_frozen = geo_spatial_memtable_has_frozen(database->spatial_memtable);
        bool succeeded = !has_active || has_frozen || database_rotate_spatial_memtable_locked(database);
        (void) geo_visibility_gate_write_unlock(&database->visibility_gate);

        if (!succeeded || (has_frozen && !database_wait_for_frozen_flush(database))) {
            return false;
        }

        if (!has_active) {
            return !has_frozen || database_wait_for_frozen_flush(database);
        }

        if (!database_wait_for_frozen_flush(database)) {
            return false;
        }
    }
}

GeoDatabaseConfig geo_database_default_config(void)
{
    GeoRocksConfig rocks_config = geo_rocks_default_config();

    return (GeoDatabaseConfig) {
        .block_cache_bytes = GEO_DATABASE_DEFAULT_BLOCK_CACHE_SIZE,
        .object_cache_bytes = GEO_DATABASE_DEFAULT_OBJECT_CACHE_SIZE,
        .write_buffer_bytes = GEO_DATABASE_DEFAULT_WRITE_BUFFER_SIZE,
        .max_active_segments = GEO_DATABASE_DEFAULT_MAX_SEGMENTS,
        .group_commit_max_operations = GEO_DATABASE_DEFAULT_GROUP_COMMIT_OPERATIONS,
        .spatial_memtable_max_operations = GEO_DATABASE_DEFAULT_SPATIAL_MEMTABLE_OPERATIONS,
        .group_commit_delay_us = GEO_DATABASE_DEFAULT_GROUP_COMMIT_DELAY_US,
        .background_jobs = rocks_config.background_jobs,
        .create_if_missing = true,
    };
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

    if (selected_config.block_cache_bytes < 1024U * 1024U ||
        (selected_config.object_cache_bytes && selected_config.object_cache_bytes < 1024U * 1024U) ||
        selected_config.write_buffer_bytes < 1024U * 1024U ||
        !selected_config.max_active_segments || !selected_config.group_commit_max_operations ||
        !selected_config.spatial_memtable_max_operations ||
        selected_config.group_commit_max_operations > selected_config.spatial_memtable_max_operations ||
        selected_config.spatial_memtable_max_operations > UINT32_MAX ||
        selected_config.group_commit_delay_us > 1000000U ||
        selected_config.background_jobs <= 0) {
        if (status) {
            *status = open_status;
        }
        return NULL;
    }

    GeoDatabase *database = calloc(1U, sizeof(*database));

    if (!database) {
        open_status = GEO_DATABASE_OUT_OF_MEMORY;
        goto complete;
    }

    atomic_init(&database->failed, false);

    database->config = selected_config;

    if (selected_config.object_cache_bytes) {
        database->object_cache = geo_object_cache_create(selected_config.object_cache_bytes);

        if (!database->object_cache) {
            open_status = GEO_DATABASE_OUT_OF_MEMORY;
            goto complete;
        }
    }

    open_status = database_initialize_paths(database, directory);

    if (open_status != GEO_DATABASE_OK) {
        goto complete;
    }

    if (!geo_visibility_gate_initialize(&database->visibility_gate)) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }
    database->visibility_gate_initialized = true;

    if (pthread_mutex_init(&database->writer_lock, NULL) != 0) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }
    database->writer_lock_initialized = true;

    if (pthread_mutex_init(&database->commit_queue_lock, NULL) != 0) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }
    database->commit_queue_lock_initialized = true;

    if (!geo_database_commit_queue_initialize(database)) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }
    database->commit_queue_ready_initialized = true;

    GeoRocksConfig rocks_config = {
        .block_cache_bytes = selected_config.block_cache_bytes,
        .write_buffer_bytes = selected_config.write_buffer_bytes,
        .background_jobs = selected_config.background_jobs,
        .create_if_missing = selected_config.create_if_missing,
    };
    GeoRocksStatus rocks_status;

    database->rocks = geo_rocks_open(database->rocks_directory, &rocks_config, &rocks_status);

    if (!database->rocks) {
        open_status = geo_db_status_from_rocks(&rocks_status);
        goto complete;
    }

    if (!geo_db_catalog_load(database->rocks, selected_config.create_if_missing, &database->catalog, &open_status)) {
        goto complete;
    }

    database->secondary_indexes = geo_secondary_catalog_open(database->rocks,
                                                              database->catalog.next_secondary_index_id,
                                                              &open_status);

    if (!database->secondary_indexes) {
        goto complete;
    }

    open_status = geo_secondary_catalog_recover(database->rocks, database->secondary_indexes);

    if (open_status != GEO_DATABASE_OK) {
        goto complete;
    }

    bool manifest_exists = access(database->manifest_path, F_OK) == 0;
    database->segments = manifest_exists ? geo_segment_set_open(database->manifest_path)
                                         : geo_segment_set_create(database->manifest_path);

    if (!database->segments && manifest_exists && database_remove_derived_manifest(database)) {
        database->segments = geo_segment_set_create(database->manifest_path);
        manifest_exists = false;
    }

    if (!database->segments) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }

    bool rebuild_spatial_index = !manifest_exists;

    if (manifest_exists) {
        uint64_t durable_watermark = geo_segment_set_durable_watermark(database->segments);

        if (durable_watermark > database->catalog.next_sequence) {
            open_status = GEO_DATABASE_CORRUPTION;
            goto complete;
        }

        if (durable_watermark < database->catalog.spatial_applied_sequence) {
            geo_segment_set_destroy(database->segments);
            database->segments = NULL;

            if (!database_remove_derived_manifest(database)) {
                open_status = GEO_DATABASE_IO_ERROR;
                goto complete;
            }

            database->segments = geo_segment_set_create(database->manifest_path);

            if (!database->segments) {
                open_status = GEO_DATABASE_IO_ERROR;
                goto complete;
            }

            rebuild_spatial_index = true;
        } else if (durable_watermark > database->catalog.spatial_applied_sequence) {
            database->catalog.spatial_applied_sequence = durable_watermark;

            if (!database_store_catalog(database, true)) {
                open_status = GEO_DATABASE_IO_ERROR;
                goto complete;
            }
        }
    }

    if (rebuild_spatial_index) {
        if (!database_rebuild_spatial_index(database)) {
            open_status = GEO_DATABASE_CORRUPTION;
            goto complete;
        }
    } else if (!database_replay_spatial_deltas(database)) {
        open_status = GEO_DATABASE_CORRUPTION;
        goto complete;
    }

    if (!geo_segment_set_enable_background_compaction(database->segments,
                                                       database->segments_directory,
                                                       database->config.max_active_segments)) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }

    database->spatial_memtable = geo_spatial_memtable_create(database->config.spatial_memtable_max_operations);

    if (!database->spatial_memtable) {
        open_status = GEO_DATABASE_OUT_OF_MEMORY;
        goto complete;
    }

    if (pthread_mutex_init(&database->spatial_flush_lock, NULL) != 0) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }
    database->spatial_flush_lock_initialized = true;

    if (pthread_cond_init(&database->spatial_flush_ready, NULL) != 0) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }
    database->spatial_flush_ready_initialized = true;

    if (pthread_create(&database->spatial_flush_thread, NULL, database_spatial_flush_thread_main, database) != 0) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }
    database->spatial_flush_thread_started = true;

    if (pthread_create(&database->commit_thread, NULL, geo_database_commit_thread_main, database) != 0) {
        open_status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }
    database->commit_thread_started = true;

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

    if (database->commit_thread_started) {
        if (pthread_mutex_lock(&database->commit_queue_lock) == 0) {
            database->commit_stop = true;
            (void) pthread_cond_signal(&database->commit_queue_ready);
            (void) pthread_mutex_unlock(&database->commit_queue_lock);
        }
        (void) pthread_join(database->commit_thread, NULL);
    }

    if (database->spatial_flush_thread_started) {
        bool flushed = database->writer_lock_initialized && pthread_mutex_lock(&database->writer_lock) == 0;

        if (flushed) {
            flushed = database_flush_all_spatial_memtables(database);

            if (flushed) {
                database->catalog.spatial_applied_sequence = geo_segment_set_durable_watermark(database->segments);
                flushed = database_store_catalog(database, true);
            }

            (void) pthread_mutex_unlock(&database->writer_lock);
        }

        if (!flushed) {
            atomic_store_explicit(&database->failed, true, memory_order_release);
        }

        if (pthread_mutex_lock(&database->spatial_flush_lock) == 0) {
            database->spatial_flush_stop = true;
            (void) pthread_cond_broadcast(&database->spatial_flush_ready);
            (void) pthread_mutex_unlock(&database->spatial_flush_lock);
        }

        (void) pthread_join(database->spatial_flush_thread, NULL);
    }

    geo_spatial_memtable_destroy(database->spatial_memtable);
    geo_segment_set_destroy(database->segments);
    geo_secondary_catalog_destroy(database->secondary_indexes);
    geo_rocks_close(database->rocks);
    geo_object_cache_destroy(database->object_cache);

    if (database->writer_lock_initialized) {
        (void) pthread_mutex_destroy(&database->writer_lock);
    }

    if (database->commit_queue_ready_initialized) {
        (void) pthread_cond_destroy(&database->commit_queue_ready);
    }

    if (database->commit_queue_lock_initialized) {
        (void) pthread_mutex_destroy(&database->commit_queue_lock);
    }

    if (database->spatial_flush_ready_initialized) {
        (void) pthread_cond_destroy(&database->spatial_flush_ready);
    }

    if (database->spatial_flush_lock_initialized) {
        (void) pthread_mutex_destroy(&database->spatial_flush_lock);
    }

    if (database->visibility_gate_initialized) {
        geo_visibility_gate_destroy(&database->visibility_gate);
    }

    free(database->segments_directory);
    free(database->manifest_path);
    free(database->rocks_directory);
    free(database->directory);
    free(database);
}

GeoDatabaseStatus geo_database_commit_view(GeoDatabase *database,
                                           const GeoDatabaseMutationView *view,
                                           size_t mutation_count,
                                           uint64_t *generated_object_id,
                                           bool *replayed)
{
    if (!database) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoDatabaseStatus status = geo_database_validate_mutations(view, mutation_count);

    if (status != GEO_DATABASE_OK) {
        return status;
    }

    if (mutation_count > SIZE_MAX / sizeof(GeoDatabaseMutation) ||
        mutation_count > SIZE_MAX / sizeof(GeoDatabaseObjectMutation) ||
        mutation_count > SIZE_MAX / sizeof(GeoObjectCacheEntry *)) {
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    GeoDatabaseMutation *spatial_operations = malloc(mutation_count * sizeof(*spatial_operations));
    GeoDatabaseObjectMutation *normalized_mutations = malloc(mutation_count * sizeof(*normalized_mutations));
    GeoObjectCacheEntry **cache_entries = database->object_cache
                                             ? calloc(mutation_count, sizeof(*cache_entries))
                                             : NULL;
    GeoRocksBuffer *old_objects = NULL;
    GeoSecondaryWriteState secondary_write = { 0 };
    bool visibility_locked = false;

    if (!spatial_operations || !normalized_mutations || (database->object_cache && !cache_entries)) {
        free(cache_entries);
        free(normalized_mutations);
        free(spatial_operations);
        return GEO_DATABASE_OUT_OF_MEMORY;
    }

    if (pthread_mutex_lock(&database->writer_lock) != 0) {
        free(cache_entries);
        free(old_objects);
        free(normalized_mutations);
        free(spatial_operations);
        return GEO_DATABASE_FAILED_STATE;
    }

    bool has_secondary_indexes = geo_secondary_catalog_has_indexes(database->secondary_indexes);

    if (has_secondary_indexes) {
        old_objects = calloc(mutation_count, sizeof(*old_objects));

        if (!old_objects) {
            status = GEO_DATABASE_OUT_OF_MEMORY;
            goto complete;
        }
    }

    if (replayed) {
        *replayed = false;
    }

    if (atomic_load_explicit(&database->failed, memory_order_acquire) ||
        mutation_count > UINT64_MAX - database->catalog.next_sequence ||
        mutation_count > UINT64_MAX - database->catalog.committed_operations) {
        status = GEO_DATABASE_FAILED_STATE;
        goto complete;
    }

    uint64_t first_sequence = database->catalog.next_sequence;
    uint64_t next_server_id = database->catalog.next_server_id;
    uint64_t allocated_object_id = 0U;
    GeoDbIdempotencyRecord idempotency_record = { 0 };
    bool has_idempotency = view->idempotency_key_size != 0U;
    size_t reserved_bytes = 256U;

    if (has_idempotency) {
        uint64_t now_seconds;
        bool found;

        if (!database_realtime_seconds(&now_seconds) ||
            view->idempotency_retention_seconds > UINT64_MAX - now_seconds) {
            status = GEO_DATABASE_FAILED_STATE;
            goto complete;
        }

        idempotency_record = database_idempotency_record(view,
                                                         mutation_count,
                                                         first_sequence,
                                                         now_seconds + view->idempotency_retention_seconds);
        GeoDbIdempotencyRecord stored_record;

        if (!geo_db_idempotency_get(database->rocks,
                                    view->idempotency_key,
                                    view->idempotency_key_size,
                                    &stored_record,
                                    &found,
                                    &status)) {
            goto complete;
        }

        if (found && stored_record.expires_at_seconds > now_seconds) {
            bool matches = stored_record.operation_count == idempotency_record.operation_count &&
                           stored_record.fingerprint_low == idempotency_record.fingerprint_low &&
                           stored_record.fingerprint_high == idempotency_record.fingerprint_high;

            if (!matches) {
                status = GEO_DATABASE_IDEMPOTENCY_CONFLICT;
                goto complete;
            }

            if (replayed) {
                *replayed = true;
            }

            status = GEO_DATABASE_OK;
            goto complete;
        }
    }

    for (size_t index = 0; index < mutation_count; ++index) {
        GeoDatabaseObjectMutation mutation = geo_database_mutation_at(view, index);
        bool exists = false;
        bool generated = mutation.operation == GEO_DATABASE_INTERNAL_INSERT_GENERATED;

        if (generated) {
            do {
                if (next_server_id == UINT64_MAX) {
                    status = GEO_DATABASE_FAILED_STATE;
                    goto complete;
                }

                mutation.object_id = next_server_id++;
                status = database_object_exists(database, mutation.object_id, &exists);

                if (status != GEO_DATABASE_OK) {
                    goto complete;
                }
            } while (exists);

            mutation.operation = GEO_DATABASE_INSERT;
            allocated_object_id = mutation.object_id;
        }

        bool conditional = mutation.operation == GEO_DATABASE_INSERT || mutation.operation == GEO_DATABASE_UPDATE;

        if (!generated && (has_secondary_indexes || conditional)) {
            GeoRocksBuffer transient = { 0 };
            GeoRocksBuffer *loaded = has_secondary_indexes ? old_objects + index : &transient;

            status = geo_database_object_load(database, mutation.object_id, loaded, &exists);

            if (!has_secondary_indexes) {
                geo_rocks_buffer_release(&transient);
            }

            if (status != GEO_DATABASE_OK) {
                goto complete;
            }
        }

        if (conditional) {
            if (mutation.operation == GEO_DATABASE_INSERT && exists) {
                status = GEO_DATABASE_ALREADY_EXISTS;
                goto complete;
            }

            if (mutation.operation == GEO_DATABASE_UPDATE && !exists) {
                status = GEO_DATABASE_NOT_FOUND;
                goto complete;
            }
        }

        spatial_operations[index] = (GeoDatabaseMutation) {
            .object_id = mutation.object_id,
            .morton_code = mutation.morton_code,
            .operation = mutation.operation == GEO_DATABASE_DELETE ? GEO_DATABASE_DELETE : GEO_DATABASE_UPSERT,
        };
        normalized_mutations[index] = mutation;

        if (cache_entries && mutation.operation != GEO_DATABASE_DELETE) {
            GeoDocView document = {
                .data = mutation.document,
                .size = mutation.document_size,
            };

            cache_entries[index] = geo_object_cache_entry_create(mutation.object_id,
                                                                  first_sequence + index,
                                                                  mutation.morton_code,
                                                                  document);
        }

        size_t mutation_bytes = mutation.document_size + GEO_DB_OBJECT_HEADER_SIZE + GEO_DB_SPATIAL_DELTA_SIZE + 16U;

        if (mutation_bytes <= SIZE_MAX - reserved_bytes) {
            reserved_bytes += mutation_bytes;
        }
    }

    if (!database_prepare_spatial_memtable(database, mutation_count)) {
        status = GEO_DATABASE_FAILED_STATE;
        goto complete;
    }

    GeoRocksStatus rocks_status;
    GeoRocksBatch *batch = geo_rocks_batch_create(database->rocks, reserved_bytes, &rocks_status);

    if (!batch) {
        status = geo_db_status_from_rocks(&rocks_status);
        goto complete;
    }

    status = geo_secondary_write_begin(database->secondary_indexes, &secondary_write);

    if (status != GEO_DATABASE_OK) {
        geo_rocks_batch_destroy(batch);
        goto complete;
    }

    GeoDbCatalog committed_catalog = database->catalog;
    committed_catalog.next_sequence += mutation_count;
    committed_catalog.committed_operations += mutation_count;
    committed_catalog.next_server_id = next_server_id;

    for (size_t index = 0; status == GEO_DATABASE_OK && index < mutation_count; ++index) {
        GeoDatabaseObjectMutation mutation = normalized_mutations[index];
        uint64_t sequence = first_sequence + index;
        GeoDocView old_document = { 0 };

        if (old_objects && old_objects[index].data) {
            uint64_t old_sequence;
            uint64_t old_morton_code;

            if (!geo_db_object_decode_view(old_objects[index].data,
                                           old_objects[index].size,
                                           &old_sequence,
                                           &old_morton_code,
                                           &old_document)) {
                status = GEO_DATABASE_CORRUPTION;
                break;
            }
        }

        /* validate_mutations established this view before any transaction preparation.
         * Reuse it across indexes and serialization; persisted old objects are validated on decode. */
        GeoDocView new_document = {
            .data = mutation.document,
            .size = mutation.document_size,
        };

        status = geo_secondary_write_object(database->secondary_indexes,
                                            batch,
                                            &secondary_write,
                                            mutation.object_id,
                                            old_document,
                                            new_document);

        if (status != GEO_DATABASE_OK) {
            break;
        }

        bool object_stored = mutation.operation != GEO_DATABASE_DELETE
                                 ? geo_db_object_put(batch,
                                                     mutation.object_id,
                                                     sequence,
                                                     mutation.morton_code,
                                                     new_document,
                                                     &status)
                                 : geo_db_object_delete(batch, mutation.object_id, &status);
        GeoDbSpatialDelta delta = {
            .object_id = mutation.object_id,
            .sequence = sequence,
            .morton_code = mutation.morton_code,
            .operation = mutation.operation == GEO_DATABASE_DELETE ? GEO_DATABASE_DELETE : GEO_DATABASE_UPSERT,
        };

        if (object_stored && !geo_db_spatial_delta_put(batch, &delta, &status)) {
            object_stored = false;
        }

        if (!object_stored && status == GEO_DATABASE_OK) {
            status = GEO_DATABASE_IO_ERROR;
        }
    }

    if (status == GEO_DATABASE_OK) {
        status = geo_secondary_write_statistics(database->secondary_indexes, batch, &secondary_write);
    }

    if (status == GEO_DATABASE_OK && !geo_db_catalog_put(batch, &committed_catalog, &status)) {
        status = GEO_DATABASE_IO_ERROR;
    }

    if (status == GEO_DATABASE_OK && has_idempotency &&
        !geo_db_idempotency_put(batch,
                                view->idempotency_key,
                                view->idempotency_key_size,
                                &idempotency_record,
                                &status)) {
        status = GEO_DATABASE_IO_ERROR;
    }

    if (status == GEO_DATABASE_OK && !geo_visibility_gate_write_lock(&database->visibility_gate)) {
        status = GEO_DATABASE_FAILED_STATE;
    } else if (status == GEO_DATABASE_OK) {
        visibility_locked = true;
    }

    if (status == GEO_DATABASE_OK && !geo_rocks_write(database->rocks, batch, true, &rocks_status)) {
        status = geo_db_status_from_rocks(&rocks_status);
    }
    geo_rocks_batch_destroy(batch);

    if (status != GEO_DATABASE_OK) {
        goto complete;
    }

    database->catalog = committed_catalog;

    if (!geo_secondary_write_commit(database->secondary_indexes, &secondary_write)) {
        atomic_store_explicit(&database->failed, true, memory_order_release);
        status = GEO_DATABASE_FAILED_STATE;
        goto complete;
    }

    if (!geo_spatial_memtable_append(database->spatial_memtable,
                                     spatial_operations,
                                     mutation_count,
                                     first_sequence)) {
        atomic_store_explicit(&database->failed, true, memory_order_release);
        status = GEO_DATABASE_IO_ERROR;
        goto complete;
    }

    if (cache_entries) {
        for (size_t index = 0U; index < mutation_count; ++index) {
            uint64_t object_id = normalized_mutations[index].object_id;

            if (cache_entries[index]) {
                geo_object_cache_publish(database->object_cache, cache_entries[index]);
                cache_entries[index] = NULL;
            } else {
                geo_object_cache_remove(database->object_cache, object_id);
            }
        }
    }

    if (geo_spatial_memtable_should_rotate(database->spatial_memtable) &&
        !geo_spatial_memtable_has_frozen(database->spatial_memtable) &&
        !database_rotate_spatial_memtable_locked(database)) {
        database->maintenance_failures++;
    }

    if (generated_object_id) {
        *generated_object_id = allocated_object_id;
    }

complete:
    if (visibility_locked) {
        (void) geo_visibility_gate_write_unlock(&database->visibility_gate);
    }
    (void) pthread_mutex_unlock(&database->writer_lock);
    geo_secondary_write_destroy(&secondary_write);
    if (old_objects) {
        for (size_t index = 0U; index < mutation_count; ++index) {
            geo_rocks_buffer_release(old_objects + index);
        }
    }
    free(old_objects);

    if (cache_entries) {
        for (size_t index = 0U; index < mutation_count; ++index) {
            geo_object_cache_entry_destroy(cache_entries[index]);
        }
    }

    free(cache_entries);
    free(normalized_mutations);
    free(spatial_operations);
    return status;
}

GeoDatabaseStatus geo_database_write_objects(GeoDatabase *database,
                                             const GeoDatabaseObjectMutation *mutations,
                                             size_t mutation_count)
{
    GeoDatabaseMutationView view = {
        .mutations = mutations,
        .layout = GEO_DATABASE_MUTATIONS_OBJECT,
    };

    return geo_database_write_view(database, &view, mutation_count, NULL, NULL);
}

GeoDatabaseStatus geo_database_write_objects_idempotent(GeoDatabase *database,
                                                        const GeoDatabaseObjectMutation *mutations,
                                                        size_t mutation_count,
                                                        const void *idempotency_key,
                                                        size_t idempotency_key_size,
                                                        uint32_t retention_seconds,
                                                        bool *replayed)
{
    if (!replayed) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoDatabaseMutationView view = {
        .mutations = mutations,
        .idempotency_key = idempotency_key,
        .idempotency_key_size = idempotency_key_size,
        .idempotency_retention_seconds = retention_seconds,
        .layout = GEO_DATABASE_MUTATIONS_OBJECT,
    };

    *replayed = false;
    return geo_database_write_view(database, &view, mutation_count, NULL, replayed);
}

GeoDatabaseStatus geo_database_write(GeoDatabase *database,
                                     const GeoDatabaseMutation *mutations,
                                     size_t mutation_count)
{
    GeoDatabaseMutationView view = {
        .mutations = mutations,
        .layout = GEO_DATABASE_MUTATIONS_COMPACT,
    };

    return geo_database_write_view(database, &view, mutation_count, NULL, NULL);
}

static GeoDatabaseStatus database_conditional_write_document(GeoDatabase *database,
                                                             uint64_t object_id,
                                                             double latitude,
                                                             double longitude,
                                                             const void *document,
                                                             size_t document_size,
                                                             GeoDatabaseOperation operation)
{
    if (!isfinite(latitude) || !isfinite(longitude) || latitude < GEO_MIN_LAT || latitude > GEO_MAX_LAT ||
        longitude < GEO_MIN_LNG || longitude > GEO_MAX_LNG) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoDatabaseObjectMutation mutation = {
        .object_id = object_id,
        .morton_code = geo_encode(latitude, longitude),
        .document = document,
        .document_size = document_size,
        .operation = operation,
    };

    return geo_database_write_objects(database, &mutation, 1U);
}

GeoDatabaseStatus geo_database_upsert_document(GeoDatabase *database,
                                               uint64_t object_id,
                                               double latitude,
                                               double longitude,
                                               const void *document,
                                               size_t document_size)
{
    return database_conditional_write_document(database,
                                               object_id,
                                               latitude,
                                               longitude,
                                               document,
                                               document_size,
                                               GEO_DATABASE_UPSERT);
}

GeoDatabaseStatus geo_database_insert_document(GeoDatabase *database,
                                               uint64_t object_id,
                                               double latitude,
                                               double longitude,
                                               const void *document,
                                               size_t document_size)
{
    return database_conditional_write_document(database,
                                               object_id,
                                               latitude,
                                               longitude,
                                               document,
                                               document_size,
                                               GEO_DATABASE_INSERT);
}

GeoDatabaseStatus geo_database_insert_generated_document(GeoDatabase *database,
                                                         double latitude,
                                                         double longitude,
                                                         const void *document,
                                                         size_t document_size,
                                                         uint64_t *object_id)
{
    if (!object_id || !isfinite(latitude) || !isfinite(longitude) || latitude < GEO_MIN_LAT || latitude > GEO_MAX_LAT ||
        longitude < GEO_MIN_LNG || longitude > GEO_MAX_LNG) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    return geo_database_insert_generated_object(database,
                                                geo_encode(latitude, longitude),
                                                document,
                                                document_size,
                                                object_id);
}

GeoDatabaseStatus geo_database_insert_generated_object(GeoDatabase *database,
                                                       uint64_t morton_code,
                                                       const void *document,
                                                       size_t document_size,
                                                       uint64_t *object_id)
{
    if (!object_id) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoDatabaseObjectMutation mutation = {
        .morton_code = morton_code,
        .document = document,
        .document_size = document_size,
        .operation = GEO_DATABASE_INTERNAL_INSERT_GENERATED,
    };
    GeoDatabaseMutationView view = {
        .mutations = &mutation,
        .layout = GEO_DATABASE_MUTATIONS_OBJECT,
    };

    *object_id = 0U;
    return geo_database_write_view(database, &view, 1U, object_id, NULL);
}

GeoDatabaseStatus geo_database_update_document(GeoDatabase *database,
                                               uint64_t object_id,
                                               double latitude,
                                               double longitude,
                                               const void *document,
                                               size_t document_size)
{
    return database_conditional_write_document(database,
                                               object_id,
                                               latitude,
                                               longitude,
                                               document,
                                               document_size,
                                               GEO_DATABASE_UPDATE);
}

GeoDatabaseStatus geo_database_upsert(GeoDatabase *database,
                                      uint64_t object_id,
                                      double latitude,
                                      double longitude)
{
    return geo_database_upsert_document(database, object_id, latitude, longitude, NULL, 0U);
}

GeoDatabaseStatus geo_database_remove(GeoDatabase *database, uint64_t object_id)
{
    GeoDatabaseMutation mutation = {
        .object_id = object_id,
        .operation = GEO_DATABASE_DELETE,
    };

    return geo_database_write(database, &mutation, 1U);
}

GeoDatabaseStatus geo_database_get(const GeoDatabase *database, uint64_t object_id, GeoDatabaseObject *object)
{
    if (!database || !object || !object_id) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    *object = (GeoDatabaseObject) { 0 };
    GeoDatabase *mutable_database = (GeoDatabase *) database;

    if (!geo_visibility_gate_read_lock(&mutable_database->visibility_gate)) {
        return GEO_DATABASE_FAILED_STATE;
    }

    unsigned char key[8];
    GeoRocksBuffer value = { 0 };
    GeoRocksStatus rocks_status;
    GeoDatabaseStatus status = GEO_DATABASE_OK;

    geo_db_store_u64_be(key, object_id);

    if (!geo_rocks_get(database->rocks, GEO_ROCKS_CF_OBJECTS, NULL, key, sizeof(key), &value, &rocks_status)) {
        status = rocks_status.code == GEO_ROCKS_NOT_FOUND ? GEO_DATABASE_NOT_FOUND : geo_db_status_from_rocks(&rocks_status);
        goto complete;
    }

    uint64_t sequence;
    uint64_t morton_code;
    const void *document;
    size_t document_size;

    if (!geo_db_object_decode(value.data, value.size, &sequence, &morton_code, &document, &document_size)) {
        status = GEO_DATABASE_CORRUPTION;
        goto complete;
    }

    void *document_copy = document_size ? malloc(document_size) : NULL;

    if (document_size && !document_copy) {
        status = GEO_DATABASE_OUT_OF_MEMORY;
        goto complete;
    }

    if (document_size) {
        memcpy(document_copy, document, document_size);
    }

    *object = (GeoDatabaseObject) {
        .object_id = object_id,
        .sequence = sequence,
        .morton_code = morton_code,
        .document = document_copy,
        .document_size = document_size,
    };

complete:
    geo_rocks_buffer_release(&value);

    if (!geo_visibility_gate_read_unlock(&mutable_database->visibility_gate)) {
        geo_database_object_release(object);
        return GEO_DATABASE_FAILED_STATE;
    }

    return status;
}

void geo_database_object_release(GeoDatabaseObject *object)
{
    if (!object) {
        return;
    }

    free(object->document);
    *object = (GeoDatabaseObject) { 0 };
}

GeoDatabaseStatus geo_database_create_index(GeoDatabase *database,
                                            const char *name,
                                            const char *json_pointer,
                                            GeoDatabaseIndexType type)
{
    if (!database) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    if (pthread_mutex_lock(&database->writer_lock) != 0) {
        return GEO_DATABASE_FAILED_STATE;
    }

    GeoDatabaseStatus status = atomic_load_explicit(&database->failed, memory_order_acquire)
                                   ? GEO_DATABASE_FAILED_STATE
                                   : geo_secondary_catalog_create(database->rocks,
                                                                  database->secondary_indexes,
                                                                  &database->catalog,
                                                                  name,
                                                                  json_pointer,
                                                                  type);

    if (status == GEO_DATABASE_FAILED_STATE) {
        atomic_store_explicit(&database->failed, true, memory_order_release);
    }

    (void) pthread_mutex_unlock(&database->writer_lock);
    return status;
}

GeoDatabaseStatus geo_database_drop_index(GeoDatabase *database, const char *name)
{
    if (!database) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    if (pthread_mutex_lock(&database->writer_lock) != 0) {
        return GEO_DATABASE_FAILED_STATE;
    }

    GeoDatabaseStatus status = GEO_DATABASE_OK;

    if (!geo_visibility_gate_write_lock(&database->visibility_gate)) {
        status = GEO_DATABASE_FAILED_STATE;
    } else {
        status = atomic_load_explicit(&database->failed, memory_order_acquire)
                     ? GEO_DATABASE_FAILED_STATE
                     : geo_secondary_catalog_drop(database->rocks, database->secondary_indexes, name);
        (void) geo_visibility_gate_write_unlock(&database->visibility_gate);
    }

    (void) pthread_mutex_unlock(&database->writer_lock);
    return status;
}

GeoDatabaseStatus geo_database_list_indexes(const GeoDatabase *database,
                                            GeoDatabaseIndexInfo *indexes,
                                            size_t capacity,
                                            size_t *count)
{
    if (!database) {
        return GEO_DATABASE_INVALID_ARGUMENT;
    }

    return geo_secondary_catalog_list(database->secondary_indexes, indexes, capacity, count);
}

GeoDatabaseStatus geo_database_query_index_reuse(const GeoDatabase *database,
                                                 const GeoDatabaseIndexPredicate *predicate,
                                                 GeoIdResult *result,
                                                 GeoDatabaseIndexQueryStats *stats)
{
    if (!database || atomic_load_explicit(&database->failed, memory_order_acquire)) {
        return database ? GEO_DATABASE_FAILED_STATE : GEO_DATABASE_INVALID_ARGUMENT;
    }

    GeoDatabase *mutable_database = (GeoDatabase *) database;

    if (!geo_visibility_gate_read_lock(&mutable_database->visibility_gate)) {
        return GEO_DATABASE_FAILED_STATE;
    }

    GeoDatabaseStatus status = geo_secondary_query(database->secondary_indexes,
                                                   database->rocks,
                                                   predicate,
                                                   result,
                                                   stats);

    if (!geo_visibility_gate_read_unlock(&mutable_database->visibility_gate)) {
        return GEO_DATABASE_FAILED_STATE;
    }

    return status;
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

    if (!geo_visibility_gate_read_lock(&mutable_database->visibility_gate)) {
        return false;
    }

    double start = stats ? geo_get_time_ms() : 0.0;
    bool succeeded = geo_segment_set_search_radius_filtered(database->segments,
                                                            latitude,
                                                            longitude,
                                                            radius_km,
                                                            result,
                                                            NULL,
                                                            geo_spatial_memtable_allows_persisted_record,
                                                            database->spatial_memtable,
                                                            stats) &&
                     geo_spatial_memtable_search_radius(database->spatial_memtable,
                                                        latitude,
                                                        longitude,
                                                        radius_km,
                                                        result,
                                                        NULL,
                                                        NULL,
                                                        NULL,
                                                        stats);

    if (!succeeded) {
        geo_result_clear(result);
    }

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start;
    }

    if (!geo_visibility_gate_read_unlock(&mutable_database->visibility_gate)) {
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

    if (!geo_visibility_gate_read_lock(&mutable_database->visibility_gate)) {
        return false;
    }

    double start = stats ? geo_get_time_ms() : 0.0;
    bool succeeded = geo_segment_set_search_radius_filtered(database->segments,
                                                            latitude,
                                                            longitude,
                                                            radius_km,
                                                            NULL,
                                                            count,
                                                            geo_spatial_memtable_allows_persisted_record,
                                                            database->spatial_memtable,
                                                            stats) &&
                     geo_spatial_memtable_search_radius(database->spatial_memtable,
                                                        latitude,
                                                        longitude,
                                                        radius_km,
                                                        NULL,
                                                        count,
                                                        NULL,
                                                        NULL,
                                                        stats);

    if (stats) {
        stats->search_time_ms = geo_get_time_ms() - start;
    }

    if (!geo_visibility_gate_read_unlock(&mutable_database->visibility_gate)) {
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

    GeoRocksStatus rocks_status;
    bool succeeded = !atomic_load_explicit(&database->failed, memory_order_acquire) &&
                     database_flush_all_spatial_memtables(database);

    if (succeeded) {
        database->catalog.spatial_applied_sequence = geo_segment_set_durable_watermark(database->segments);
        succeeded = database_store_catalog(database, true) &&
                     geo_segment_set_checkpoint_mutations(database->segments) &&
                     geo_rocks_flush(database->rocks, true, &rocks_status);
    }

    if (succeeded) {
        database->checkpoints++;
    } else {
        database->maintenance_failures++;
    }

    (void) pthread_mutex_unlock(&database->writer_lock);
    return succeeded ? GEO_DATABASE_OK : GEO_DATABASE_IO_ERROR;
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

    if (!geo_visibility_gate_read_lock(&mutable_database->visibility_gate)) {
        (void) pthread_mutex_unlock(&mutable_database->writer_lock);
        return false;
    }

    *stats = (GeoDatabaseStats) {
        .committed_operations = database->catalog.committed_operations,
        .recovered_operations = database->recovered_operations,
        .checkpoints = database->checkpoints,
        .maintenance_failures = database->maintenance_failures,
        .next_sequence = database->catalog.next_sequence,
        .physical_records = geo_segment_set_record_count(database->segments),
        .active_segments = geo_segment_set_count(database->segments),
        .active_spatial_operations = geo_spatial_memtable_active_operations(database->spatial_memtable),
        .frozen_spatial_operations = geo_spatial_memtable_frozen_operations(database->spatial_memtable),
    };

    (void) geo_visibility_gate_read_unlock(&mutable_database->visibility_gate);
    (void) pthread_mutex_unlock(&mutable_database->writer_lock);
    return true;
}
