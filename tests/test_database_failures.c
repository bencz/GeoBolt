#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "geobolt/geobolt.h"
#include "geobolt/geodoc.h"
#include "geo_rocks_bridge.h"
#include "geo_database_internal.h"
#include "test_support.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>

/* Linker interposition stays in this executable; the production library has no fault hooks.
 * Only the calling thread is armed, so maintenance workers never consume a test fault. */
static _Thread_local bool fail_iterator;
static _Thread_local bool interrupt_build;
static _Thread_local bool fault_observed;
static _Thread_local bool fail_rollback;
static _Thread_local bool rollback_fault_observed;
static _Atomic(pthread_cond_t *) watched_condition = NULL;
static atomic_size_t timed_wait_count = 0U;

int __real_pthread_cond_timedwait(pthread_cond_t *condition, pthread_mutex_t *mutex, const struct timespec *deadline);

int __wrap_pthread_cond_timedwait(pthread_cond_t *condition, pthread_mutex_t *mutex, const struct timespec *deadline)
{
    if (condition == atomic_load_explicit(&watched_condition, memory_order_acquire)) {
        atomic_fetch_add_explicit(&timed_wait_count, 1U, memory_order_relaxed);
    }

    return __real_pthread_cond_timedwait(condition, mutex, deadline);
}

bool __real_geo_rocks_write(GeoRocksDatabase *database, GeoRocksBatch *batch, bool synchronize, GeoRocksStatus *status);

bool __wrap_geo_rocks_write(GeoRocksDatabase *database, GeoRocksBatch *batch, bool synchronize, GeoRocksStatus *status)
{
    if (fail_rollback && fault_observed) {
        fail_rollback = false;
        rollback_fault_observed = true;
        *status = (GeoRocksStatus) { .code = GEO_ROCKS_IO_ERROR };
        return false;
    }

    return __real_geo_rocks_write(database, batch, synchronize, status);
}

GeoRocksIterator *__real_geo_rocks_iterator_create(GeoRocksDatabase *database,
                                                  GeoRocksColumnFamily column_family,
                                                  const GeoRocksSnapshot *snapshot,
                                                  GeoRocksStatus *status);

GeoRocksIterator *__wrap_geo_rocks_iterator_create(GeoRocksDatabase *database,
                                                  GeoRocksColumnFamily column_family,
                                                  const GeoRocksSnapshot *snapshot,
                                                  GeoRocksStatus *status)
{
    if (interrupt_build) {
        _exit(77);
    }

    if (fail_iterator) {
        fail_iterator = false;
        fault_observed = true;
        *status = (GeoRocksStatus) { .code = GEO_ROCKS_OUT_OF_MEMORY };
        return NULL;
    }

    return __real_geo_rocks_iterator_create(database, column_family, snapshot, status);
}

static bool verify_object(GeoDatabase *database)
{
    GeoDatabaseObject object = { 0 };
    bool succeeded = geo_database_get(database, 1U, &object) == GEO_DATABASE_OK;
    geo_database_object_release(&object);
    return succeeded;
}

static bool verify_index(GeoDatabase *database, bool exists)
{
    GeoDatabaseIndexPredicate predicate = {
        .index_name = "value_idx",
        .operation = GEO_DATABASE_INDEX_EQUAL,
        .lower = { .type = GEO_DATABASE_INDEX_BOOL, .as.boolean = true },
    };
    GeoIdResult result = { 0 };
    GeoDatabaseStatus status = geo_database_query_index_reuse(database, &predicate, &result, NULL);
    bool succeeded = exists ? status == GEO_DATABASE_OK && result.count == 1U && result.ids[0] == 1U
                            : status == GEO_DATABASE_NOT_FOUND && result.count == 0U;
    free(result.ids);
    return succeeded;
}

static bool run_build_failures(void)
{
    char path[] = "/tmp/geobolt-index-failures-XXXXXX";
    char *directory = mkdtemp(path);
    GeoDatabaseConfig config = geo_database_default_config();
    config.block_cache_bytes = 8U * 1024U * 1024U;
    config.object_cache_bytes = 0U;
    config.write_buffer_bytes = 4U * 1024U * 1024U;
    config.group_commit_max_operations = 1U;
    GeoDatabaseStatus status;
    GeoDatabase *database = directory ? geo_database_open(directory, &config, &status) : NULL;

    if (database) {
        atomic_store_explicit(&watched_condition, &database->commit_queue_ready, memory_order_release);
    }

    GeoDocBuilder *builder = geo_doc_builder_create();
    GeoDocBuffer document = { 0 };
    bool succeeded = database && builder &&
                     geo_doc_builder_add_bool(builder, 0U, "value", true, NULL) == GEO_DOC_OK &&
                     geo_doc_builder_finish(builder, &document) == GEO_DOC_OK &&
                     geo_database_upsert_document(database, 1U, 0.0, 0.0, document.data, document.size) == GEO_DATABASE_OK;

    /* Assert scheduling behavior, not elapsed time: a full one-operation group must never
     * enter the collection timer, even on a contended test host. */
    succeeded = succeeded && atomic_load_explicit(&timed_wait_count, memory_order_relaxed) == 0U;
    atomic_store_explicit(&watched_condition, NULL, memory_order_release);

    if (succeeded) {
        fail_iterator = true;
        status = geo_database_create_index(database, "value_idx", "/value", GEO_DATABASE_INDEX_BOOL);
        succeeded = fault_observed && status == GEO_DATABASE_OUT_OF_MEMORY && verify_index(database, false);
    }

    if (succeeded) {
        status = geo_database_create_index(database, "value_idx", "/value", GEO_DATABASE_INDEX_STRING);
        succeeded = status == GEO_DATABASE_INVALID_ARGUMENT && verify_index(database, false);
    }

    geo_database_close(database);
    database = NULL;

    if (succeeded) {
        database = geo_database_open(directory, &config, &status);
        succeeded = database && verify_object(database) && verify_index(database, false);
    }

    geo_database_close(database);
    database = NULL;

    if (succeeded) {
        /* No database threads exist at fork. The child persists BUILDING and dies before
         * scanning an incompatible canonical field, modeling interruption before rollback. */
        pid_t child = fork();

        if (child == 0) {
            database = geo_database_open(directory, &config, &status);

            if (!database) {
                _exit(78);
            }

            interrupt_build = true;
            (void) geo_database_create_index(database, "value_idx", "/value", GEO_DATABASE_INDEX_STRING);
            _exit(79);
        }

        int child_status = 0;
        pid_t waited;

        do {
            waited = child > 0 ? waitpid(child, &child_status, 0) : -1;
        } while (child > 0 && waited < 0 && errno == EINTR);

        succeeded = waited == child && child > 0 && WIFEXITED(child_status) && WEXITSTATUS(child_status) == 77;
    }

    if (succeeded) {
        database = geo_database_open(directory, &config, &status);
        succeeded = database && verify_object(database) && verify_index(database, false) &&
                    geo_database_create_index(database, "value_idx", "/value", GEO_DATABASE_INDEX_BOOL) == GEO_DATABASE_OK &&
                    verify_index(database, true);
    }

    geo_database_close(database);
    database = NULL;

    if (succeeded) {
        database = geo_database_open(directory, &config, &status);
        succeeded = database && verify_object(database) && verify_index(database, true);
    }

    if (succeeded) {
        succeeded = geo_database_drop_index(database, "value_idx") == GEO_DATABASE_OK;
    }

    if (succeeded) {
        fault_observed = false;
        fail_iterator = true;
        fail_rollback = true;
        status = geo_database_create_index(database, "value_idx", "/value", GEO_DATABASE_INDEX_BOOL);
        succeeded = fault_observed && rollback_fault_observed && status == GEO_DATABASE_FAILED_STATE &&
                    geo_database_upsert(database, 2U, 0.0, 0.0) == GEO_DATABASE_FAILED_STATE;
    }

    geo_database_close(database);
    database = NULL;

    if (succeeded) {
        database = geo_database_open(directory, &config, &status);
        succeeded = database && verify_object(database) && verify_index(database, true);
    }

    geo_database_close(database);
    geo_doc_buffer_release(&document);
    geo_doc_builder_destroy(builder);

    if (directory && !geobolt_test_remove_tree(directory)) {
        succeeded = false;
    }

    return succeeded;
}

int main(void)
{
    if (!run_build_failures()) {
        fprintf(stderr, "FAIL: index allocation failure, rejection, interrupted build, or durable retry\n");
        return EXIT_FAILURE;
    }

    printf("[PASS] index OOM, rejected type, interrupted recovery, failed rollback quarantine, and durable rebuild\n");
    return EXIT_SUCCESS;
}
