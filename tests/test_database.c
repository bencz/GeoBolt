#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "geobolt/geobolt.h"
#include "test_support.h"

#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdatomic.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#define TEST_ASSERT(condition, message)                                                                                   \
    do {                                                                                                                  \
        if (!(condition)) {                                                                                               \
            fprintf(stderr, "FAIL: %s (%s:%d)\n", message, __FILE__, __LINE__);                                         \
            succeeded = false;                                                                                            \
            goto cleanup;                                                                                                 \
        }                                                                                                                 \
    } while (0)

static bool test_write_all(int descriptor, const void *data, size_t size)
{
    const unsigned char *bytes = data;

    while (size) {
        ssize_t written = write(descriptor, bytes, size);

        if (written < 0 && errno == EINTR) {
            continue;
        }

        if (written <= 0) {
            return false;
        }

        bytes += (size_t) written;
        size -= (size_t) written;
    }

    return true;
}

static bool test_read_file(const char *path, unsigned char **contents, size_t *size)
{
    struct stat status;
    int descriptor = open(path, O_RDONLY);

    if (descriptor < 0 || fstat(descriptor, &status) != 0 || status.st_size < 0 || (uintmax_t) status.st_size > SIZE_MAX) {
        if (descriptor >= 0) {
            (void) close(descriptor);
        }

        return false;
    }

    unsigned char *buffer = status.st_size ? malloc((size_t) status.st_size) : NULL;
    size_t offset = 0;

    while (offset < (size_t) status.st_size) {
        ssize_t received = read(descriptor, buffer + offset, (size_t) status.st_size - offset);

        if (received < 0 && errno == EINTR) {
            continue;
        }

        if (received <= 0) {
            free(buffer);
            (void) close(descriptor);

            return false;
        }

        offset += (size_t) received;
    }

    bool succeeded = close(descriptor) == 0;

    if (!succeeded) {
        free(buffer);

        return false;
    }

    *contents = buffer;
    *size = (size_t) status.st_size;

    return true;
}

static bool test_replace_file(const char *path, const void *contents, size_t size)
{
    int descriptor = open(path, O_WRONLY | O_TRUNC);
    bool succeeded = descriptor >= 0 && test_write_all(descriptor, contents, size) && fsync(descriptor) == 0;

    if (descriptor >= 0 && close(descriptor) != 0) {
        succeeded = false;
    }

    return succeeded;
}

static bool result_contains(const GeoSearchResult *result, uint64_t object_id, uint64_t morton_code)
{
    for (size_t index = 0; index < result->count; ++index) {
        if (result->results[index].id == object_id && result->results[index].z == morton_code) {
            return true;
        }
    }

    return false;
}

static bool test_database_recovery_and_last_write_wins(void)
{
    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-database-test-XXXXXX";
    char *directory = mkdtemp(directory_template);
    GeoDatabase *database = NULL;
    GeoSearchResult *result = NULL;
    unsigned char *initial_state = NULL;
    size_t initial_state_size = 0;
    char state_path[256];

    TEST_ASSERT(directory != NULL, "temporary database directory must be created");
    TEST_ASSERT(snprintf(state_path, sizeof(state_path), "%s/state.gbs", directory) > 0, "state path must fit");

    GeoDatabaseConfig config = geo_database_default_config();

    config.wal_segment_size = 64U * 1024U;
    config.checkpoint_interval_operations = 1024U;
    config.max_active_segments = 4U;

    GeoDatabaseStatus status;
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "new database must open");

    geo_database_close(database);
    database = NULL;
    TEST_ASSERT(test_read_file(state_path, &initial_state, &initial_state_size), "initial durable state must be readable");

    config.create_if_missing = false;
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "existing database must reopen");

    GeoDatabaseMutation initial[] = {
        { .object_id = 101U, .morton_code = geo_encode(-23.5505, -46.6333), .operation = GEO_DATABASE_UPSERT },
        { .object_id = 202U, .morton_code = geo_encode(-22.9068, -43.1729), .operation = GEO_DATABASE_UPSERT },
        { .object_id = 303U, .morton_code = geo_encode(40.7128, -74.0060), .operation = GEO_DATABASE_UPSERT },
    };

    TEST_ASSERT(geo_database_write(database, initial, 3U) == GEO_DATABASE_OK, "initial batch must commit");
    geo_database_close(database);
    database = NULL;

    // Models a crash after segment publication and before the state watermark rename. Recovery must replay idempotently.
    TEST_ASSERT(test_replace_file(state_path, initial_state, initial_state_size), "pre-commit state watermark must be restored");
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "database must recover a fully durable unapplied WAL frame");

    GeoDatabaseStats stats;
    TEST_ASSERT(geo_database_get_stats(database, &stats), "database stats must be available");
    TEST_ASSERT(stats.recovered_operations == 3U, "recovery must account for the replayed operations");
    TEST_ASSERT(stats.committed_operations == 3U && stats.next_sequence == 4U, "recovery sequence must advance exactly once");

    result = geo_result_create(8U);
    TEST_ASSERT(result != NULL, "reusable result must be allocated");
    TEST_ASSERT(geo_database_search_radius_reuse(database, 0.0, 0.0, 25000.0, result, NULL), "global query must succeed");
    TEST_ASSERT(result->count == 3U, "idempotent recovery must not expose duplicate objects");

    GeoDatabaseMutation update_and_delete[] = {
        { .object_id = 101U, .morton_code = geo_encode(51.5074, -0.1278), .operation = GEO_DATABASE_UPSERT },
        { .object_id = 202U, .operation = GEO_DATABASE_DELETE },
    };

    TEST_ASSERT(geo_database_write(database, update_and_delete, 2U) == GEO_DATABASE_OK, "mixed update/delete batch must commit");
    TEST_ASSERT(geo_database_search_radius_reuse(database, 0.0, 0.0, 25000.0, result, NULL), "post-update query must succeed");
    TEST_ASSERT(result->count == 2U, "delete and replacement must have last-write-wins visibility");
    TEST_ASSERT(result_contains(result, 101U, update_and_delete[0].morton_code), "updated Morton coordinate must be visible");
    TEST_ASSERT(!result_contains(result, 101U, initial[0].morton_code), "superseded Morton coordinate must be invisible");

    GeoDatabaseMutation duplicate_ids[] = {
        { .object_id = 404U, .morton_code = geo_encode(1.0, 1.0), .operation = GEO_DATABASE_UPSERT },
        { .object_id = 404U, .operation = GEO_DATABASE_DELETE },
    };

    TEST_ASSERT(geo_database_write(database, duplicate_ids, 2U) == GEO_DATABASE_INVALID_ARGUMENT,
                "ambiguous duplicate IDs inside one generation must be rejected");
    TEST_ASSERT(geo_database_checkpoint(database) == GEO_DATABASE_OK, "explicit checkpoint must durably reclaim the WAL prefix");

    geo_database_close(database);
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "checkpointed database must reopen");
    TEST_ASSERT(geo_database_search_radius_reuse(database, 0.0, 0.0, 25000.0, result, NULL), "reopened query must succeed");
    TEST_ASSERT(result->count == 2U && result_contains(result, 101U, update_and_delete[0].morton_code),
                "checkpoint must preserve exact logical state");

cleanup:
    geo_result_destroy(result);
    geo_database_close(database);
    free(initial_state);

    if (directory && !geobolt_test_remove_tree(directory)) {
        fprintf(stderr, "FAIL: temporary database tree could not be removed\n");
        succeeded = false;
    }

    return succeeded;
}

typedef struct {
    const GeoDatabase *database;
    atomic_bool *start;
    atomic_bool *stop;
    atomic_bool *failed;
    size_t expected_count;
    size_t queries;
} DatabaseReaderContext;

static void *database_reader_main(void *argument)
{
    DatabaseReaderContext *context = argument;

    while (!atomic_load_explicit(context->start, memory_order_acquire)) {
    }

    while (!atomic_load_explicit(context->stop, memory_order_acquire)) {
        size_t count = 0;

        if (!geo_database_search_radius_count(context->database, 0.0, 0.0, 25000.0, &count, NULL) ||
            count != context->expected_count) {
            atomic_store_explicit(context->failed, true, memory_order_release);
            break;
        }

        context->queries++;
    }

    return NULL;
}

static bool test_database_concurrent_batches_and_queries(void)
{
    enum {
        OBJECT_COUNT = 256,
        READER_COUNT = 3,
        WRITE_ROUNDS = 12,
    };

    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-database-concurrent-XXXXXX";
    char *directory = mkdtemp(directory_template);
    GeoDatabase *database = NULL;
    GeoDatabaseMutation *mutations = NULL;
    pthread_t readers[READER_COUNT];
    DatabaseReaderContext contexts[READER_COUNT];
    size_t readers_started = 0;
    atomic_bool start = ATOMIC_VAR_INIT(false);
    atomic_bool stop = ATOMIC_VAR_INIT(false);
    atomic_bool failed = ATOMIC_VAR_INIT(false);

    TEST_ASSERT(directory != NULL, "concurrent test directory must be created");

    GeoDatabaseConfig config = geo_database_default_config();

    config.wal_segment_size = 64U * 1024U;
    config.checkpoint_interval_operations = 1024U;
    config.max_active_segments = 4U;

    GeoDatabaseStatus status;
    database = geo_database_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status == GEO_DATABASE_OK, "concurrent database must open");

    mutations = calloc(OBJECT_COUNT, sizeof(*mutations));
    TEST_ASSERT(mutations != NULL, "concurrent mutation batch must allocate");

    for (size_t object = 0; object < OBJECT_COUNT; ++object) {
        mutations[object] = (GeoDatabaseMutation) {
            .object_id = object + 1U,
            .morton_code = geo_encode(-40.0 + (double) object * 0.25, -120.0 + (double) object * 0.5),
            .operation = GEO_DATABASE_UPSERT,
        };
    }

    TEST_ASSERT(geo_database_write(database, mutations, OBJECT_COUNT) == GEO_DATABASE_OK, "initial concurrent batch must commit");

    for (size_t reader = 0; reader < READER_COUNT; ++reader) {
        contexts[reader] = (DatabaseReaderContext) {
            .database = database,
            .start = &start,
            .stop = &stop,
            .failed = &failed,
            .expected_count = OBJECT_COUNT,
        };

        TEST_ASSERT(pthread_create(readers + reader, NULL, database_reader_main, contexts + reader) == 0,
                    "database query reader must start");
        readers_started++;
    }

    atomic_store_explicit(&start, true, memory_order_release);

    for (size_t round = 0; round < WRITE_ROUNDS; ++round) {
        for (size_t object = 0; object < OBJECT_COUNT; ++object) {
            double latitude = -60.0 + (double) ((object + round * 17U) % 480U) * 0.25;
            double longitude = -170.0 + (double) ((object * 7U + round * 29U) % 680U) * 0.5;

            mutations[object].morton_code = geo_encode(latitude, longitude);
        }

        TEST_ASSERT(geo_database_write(database, mutations, OBJECT_COUNT) == GEO_DATABASE_OK,
                    "concurrent replacement batch must commit");
    }

    atomic_store_explicit(&stop, true, memory_order_release);

    for (size_t reader = 0; reader < readers_started; ++reader) {
        TEST_ASSERT(pthread_join(readers[reader], NULL) == 0, "database query reader must join");
    }

    readers_started = 0;
    TEST_ASSERT(!atomic_load_explicit(&failed, memory_order_acquire), "readers must never observe a partial replacement generation");

    size_t total_queries = 0;

    for (size_t reader = 0; reader < READER_COUNT; ++reader) {
        total_queries += contexts[reader].queries;
    }

    TEST_ASSERT(total_queries > 0, "concurrent readers must execute queries while batches publish");

cleanup:
    atomic_store_explicit(&stop, true, memory_order_release);
    atomic_store_explicit(&start, true, memory_order_release);

    for (size_t reader = 0; reader < readers_started; ++reader) {
        (void) pthread_join(readers[reader], NULL);
    }

    free(mutations);
    geo_database_close(database);

    if (directory && !geobolt_test_remove_tree(directory)) {
        fprintf(stderr, "FAIL: concurrent database tree could not be removed\n");
        succeeded = false;
    }

    return succeeded;
}

int main(void)
{
    printf("========================================\n");
    printf("GEOBOLT DATABASE TEST SUITE\n");
    printf("========================================\n\n");

    if (!test_database_recovery_and_last_write_wins()) {
        return EXIT_FAILURE;
    }

    printf("[PASS] recovery, checkpoint, update, delete, and last-write-wins\n");

    if (!test_database_concurrent_batches_and_queries()) {
        return EXIT_FAILURE;
    }

    printf("[PASS] concurrent batch publication and snapshot queries\n");

    return EXIT_SUCCESS;
}
