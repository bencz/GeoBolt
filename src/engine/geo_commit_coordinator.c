#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "geo_database_internal.h"

#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

bool geo_database_commit_queue_initialize(GeoDatabase *database)
{
    if (!database) {
        return false;
    }

    pthread_condattr_t attributes;

    if (pthread_condattr_init(&attributes) != 0) {
        return false;
    }

#if defined(__APPLE__)
    /* Darwin uses a relative timed wait below rather than a selectable condition clock. */
    bool supported = true;
#elif defined(_POSIX_CLOCK_SELECTION) && _POSIX_CLOCK_SELECTION >= 0
    bool supported = pthread_condattr_setclock(&attributes, CLOCK_MONOTONIC) == 0;
#else
    bool supported = false;
#endif

    bool initialized = supported && pthread_cond_init(&database->commit_queue_ready, &attributes) == 0;
    int destroy_status = pthread_condattr_destroy(&attributes);

    if (initialized && destroy_status != 0) {
        (void) pthread_cond_destroy(&database->commit_queue_ready);
        return false;
    }

    return initialized;
}

static int commit_wait_until(GeoDatabase *database, const struct timespec *deadline)
{
#if defined(__APPLE__)
    struct timespec now;

    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) {
        return errno;
    }
    if (now.tv_sec > deadline->tv_sec || (now.tv_sec == deadline->tv_sec && now.tv_nsec >= deadline->tv_nsec)) {
        return ETIMEDOUT;
    }

    struct timespec remaining = {
        .tv_sec = deadline->tv_sec - now.tv_sec,
        .tv_nsec = deadline->tv_nsec - now.tv_nsec,
    };

    if (remaining.tv_nsec < 0) {
        remaining.tv_sec--;
        remaining.tv_nsec += 1000000000L;
    }

    return pthread_cond_timedwait_relative_np(&database->commit_queue_ready, &database->commit_queue_lock, &remaining);
#else
    return pthread_cond_timedwait(&database->commit_queue_ready, &database->commit_queue_lock, deadline);
#endif
}

static bool commit_view_is_groupable(const GeoDatabaseMutationView *view, size_t mutation_count)
{
    if (view->idempotency_key_size) {
        return false;
    }

    for (size_t index = 0; index < mutation_count; ++index) {
        GeoDatabaseObjectMutation mutation = geo_database_mutation_at(view, index);

        if (mutation.operation != GEO_DATABASE_UPSERT && mutation.operation != GEO_DATABASE_DELETE) {
            return false;
        }
    }

    return true;
}

static bool commit_deadline(uint32_t delay_us, struct timespec *deadline)
{
    if (clock_gettime(CLOCK_MONOTONIC, deadline) != 0) {
        return false;
    }

    deadline->tv_nsec += (long) delay_us * 1000L;

    if (deadline->tv_nsec >= 1000000000L) {
        deadline->tv_sec++;
        deadline->tv_nsec -= 1000000000L;
    }

    return true;
}

static GeoDatabaseWriteRequest *commit_dequeue_group(GeoDatabase *database,
                                                     size_t *request_count,
                                                     size_t *mutation_count)
{
    while (!database->commit_head && !database->commit_stop) {
        int wait_status = pthread_cond_wait(&database->commit_queue_ready, &database->commit_queue_lock);

        if (wait_status != 0) {
            atomic_store_explicit(&database->failed, true, memory_order_release);
            database->commit_stop = true;
        }
    }

    if (!database->commit_head) {
        return NULL;
    }

    bool groupable = commit_view_is_groupable(&database->commit_head->view, database->commit_head->mutation_count);

    struct timespec deadline;
    bool may_wait = groupable && database->config.group_commit_delay_us &&
                    database->commit_head->mutation_count < database->config.group_commit_max_operations;

    if (may_wait && !commit_deadline(database->config.group_commit_delay_us, &deadline)) {
        atomic_store_explicit(&database->failed, true, memory_order_release);
        may_wait = false;
    }

    GeoDatabaseWriteRequest *first = NULL;
    GeoDatabaseWriteRequest *last = NULL;
    size_t requests = 0U;
    size_t mutations = 0U;

    for (;;) {
        if (!database->commit_head) {
            if (!may_wait || database->commit_stop || mutations >= database->config.group_commit_max_operations) {
                break;
            }

            int wait_status = commit_wait_until(database, &deadline);

            if (wait_status != 0) {
                if (wait_status != ETIMEDOUT) {
                    atomic_store_explicit(&database->failed, true, memory_order_release);
                }

                break;
            }

            continue;
        }

        GeoDatabaseWriteRequest *candidate = database->commit_head;
        bool candidate_groupable = commit_view_is_groupable(&candidate->view, candidate->mutation_count);
        bool capacity_exhausted = mutations >= database->config.group_commit_max_operations;
        bool candidate_exceeds_capacity = !capacity_exhausted &&
                                          candidate->mutation_count > database->config.group_commit_max_operations - mutations;

        if (requests && (!groupable || !candidate_groupable || capacity_exhausted || candidate_exceeds_capacity)) {
            break;
        }

        database->commit_head = candidate->next;

        /* Producers may enqueue while the collector waits. Never leave the queue tail
         * pointing into the detached group, which has separate completion ownership. */
        if (!database->commit_head) {
            database->commit_tail = NULL;
        }

        candidate->next = NULL;

        if (last) {
            last->next = candidate;
        } else {
            first = candidate;
        }

        last = candidate;
        requests++;
        mutations += candidate->mutation_count;

        if (!groupable || mutations == database->config.group_commit_max_operations) {
            break;
        }
    }

    *request_count = requests;
    *mutation_count = mutations;
    return first;
}

static void commit_complete_group_locked(GeoDatabaseWriteRequest *requests, GeoDatabaseStatus status)
{
    for (GeoDatabaseWriteRequest *request = requests; request; request = request->next) {
        request->status = status;
        request->completed = true;
        (void) pthread_cond_signal(&request->completion);
    }
}

static void commit_complete_group(GeoDatabase *database,
                                  GeoDatabaseWriteRequest *requests,
                                  GeoDatabaseStatus status)
{
    if (pthread_mutex_lock(&database->commit_queue_lock) != 0) {
        atomic_store_explicit(&database->failed, true, memory_order_release);
        return;
    }

    commit_complete_group_locked(requests, status);
    (void) pthread_mutex_unlock(&database->commit_queue_lock);
}

static void commit_process_group(GeoDatabase *database,
                                 GeoDatabaseWriteRequest *requests,
                                 size_t request_count,
                                 size_t mutation_count)
{
    if (request_count == 1U) {
        GeoDatabaseStatus status = geo_database_commit_view(database,
                                                            &requests->view,
                                                            requests->mutation_count,
                                                            &requests->generated_object_id,
                                                            &requests->replayed);

        commit_complete_group(database, requests, status);
        return;
    }

    if (mutation_count > SIZE_MAX / sizeof(GeoDatabaseObjectMutation)) {
        commit_complete_group(database, requests, GEO_DATABASE_OUT_OF_MEMORY);
        return;
    }

    GeoDatabaseObjectMutation *combined = malloc(mutation_count * sizeof(*combined));

    if (!combined) {
        commit_complete_group(database, requests, GEO_DATABASE_OUT_OF_MEMORY);
        return;
    }

    size_t output = 0U;

    for (GeoDatabaseWriteRequest *request = requests; request; request = request->next) {
        for (size_t index = 0; index < request->mutation_count; ++index) {
            combined[output++] = geo_database_mutation_at(&request->view, index);
        }
    }

    GeoDatabaseMutationView combined_view = {
        .mutations = combined,
        .layout = GEO_DATABASE_MUTATIONS_OBJECT,
    };
    GeoDatabaseStatus status = geo_database_commit_view(database, &combined_view, mutation_count, NULL, NULL);
    free(combined);

    // Duplicate IDs across independent calls make the combined batch invalid, but do not make either call invalid on its own.
    if (status == GEO_DATABASE_INVALID_ARGUMENT) {
        for (GeoDatabaseWriteRequest *request = requests; request;) {
            GeoDatabaseWriteRequest *next = request->next;
            request->next = NULL;
            GeoDatabaseStatus request_status = geo_database_commit_view(database,
                                                                        &request->view,
                                                                        request->mutation_count,
                                                                        &request->generated_object_id,
                                                                        &request->replayed);

            commit_complete_group(database, request, request_status);
            request = next;
        }

        return;
    }

    commit_complete_group(database, requests, status);
}

void *geo_database_commit_thread_main(void *argument)
{
    GeoDatabase *database = argument;

    for (;;) {
        if (pthread_mutex_lock(&database->commit_queue_lock) != 0) {
            atomic_store_explicit(&database->failed, true, memory_order_release);
            return NULL;
        }

        size_t request_count = 0U;
        size_t mutation_count = 0U;
        GeoDatabaseWriteRequest *requests = commit_dequeue_group(database, &request_count, &mutation_count);
        bool stopped = !requests && database->commit_stop;
        (void) pthread_mutex_unlock(&database->commit_queue_lock);

        if (stopped) {
            return NULL;
        }

        if (requests) {
            commit_process_group(database, requests, request_count, mutation_count);
        }
    }
}

GeoDatabaseStatus geo_database_write_view(GeoDatabase *database,
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

    GeoDatabaseWriteRequest request = {
        .view = *view,
        .mutation_count = mutation_count,
        .status = GEO_DATABASE_FAILED_STATE,
    };

    if (pthread_cond_init(&request.completion, NULL) != 0) {
        return GEO_DATABASE_FAILED_STATE;
    }

    if (pthread_mutex_lock(&database->commit_queue_lock) != 0) {
        (void) pthread_cond_destroy(&request.completion);
        return GEO_DATABASE_FAILED_STATE;
    }

    if (database->commit_stop || !database->commit_thread_started) {
        (void) pthread_mutex_unlock(&database->commit_queue_lock);
        (void) pthread_cond_destroy(&request.completion);
        return GEO_DATABASE_FAILED_STATE;
    }

    if (database->commit_tail) {
        database->commit_tail->next = &request;
    } else {
        database->commit_head = &request;
    }

    database->commit_tail = &request;
    (void) pthread_cond_signal(&database->commit_queue_ready);

    while (!request.completed) {
        int wait_status = pthread_cond_wait(&request.completion, &database->commit_queue_lock);

        if (wait_status != 0) {
            atomic_store_explicit(&database->failed, true, memory_order_release);
            status = GEO_DATABASE_FAILED_STATE;
            break;
        }
    }

    if (request.completed) {
        status = request.status;

        if (status == GEO_DATABASE_OK && generated_object_id) {
            *generated_object_id = request.generated_object_id;
        }

        if (status == GEO_DATABASE_OK && replayed) {
            *replayed = request.replayed;
        }
    }

    (void) pthread_mutex_unlock(&database->commit_queue_lock);
    (void) pthread_cond_destroy(&request.completion);
    return status;
}
