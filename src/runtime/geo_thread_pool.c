#include "geo_thread_pool.h"

#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
    struct GeoThreadPool *pool;
    size_t worker_index;
} GeoThreadPoolSlot;

struct GeoThreadPool {
    pthread_mutex_t run_lock;
    pthread_mutex_t state_lock;
    pthread_cond_t work_available;
    pthread_cond_t work_complete;
    pthread_t *threads;
    GeoThreadPoolSlot *slots;
    GeoThreadPoolWorker worker;
    void *context;
    uint64_t generation;
    size_t capacity;
    size_t active_workers;
    size_t completed_workers;
    bool failed;
    bool stopping;
};

static void *geo_thread_pool_worker_main(void *argument)
{
    GeoThreadPoolSlot *slot = argument;
    GeoThreadPool *pool = slot->pool;
    uint64_t observed_generation = 0;

    pthread_mutex_lock(&pool->state_lock);

    while (!pool->stopping) {
        while (observed_generation == pool->generation && !pool->stopping) {
            if (pthread_cond_wait(&pool->work_available, &pool->state_lock) != 0) {
                pool->failed = true;
                pool->stopping = true;
                pthread_cond_broadcast(&pool->work_complete);
                break;
            }
        }

        if (pool->stopping) {
            break;
        }

        observed_generation = pool->generation;

        if (slot->worker_index >= pool->active_workers) {
            continue;
        }

        GeoThreadPoolWorker worker = pool->worker;
        void *context = pool->context;
        size_t worker_count = pool->active_workers;

        pthread_mutex_unlock(&pool->state_lock);
        bool succeeded = worker(context, slot->worker_index, worker_count);
        pthread_mutex_lock(&pool->state_lock);

        if (!succeeded) {
            pool->failed = true;
        }

        pool->completed_workers++;
        pthread_cond_signal(&pool->work_complete);
    }

    pthread_mutex_unlock(&pool->state_lock);

    return NULL;
}

static void geo_thread_pool_destroy_initialized(GeoThreadPool *pool, size_t thread_count)
{
    pthread_mutex_lock(&pool->state_lock);
    pool->stopping = true;
    pthread_cond_broadcast(&pool->work_available);
    pthread_mutex_unlock(&pool->state_lock);

    for (size_t thread = 0; thread < thread_count; ++thread) {
        (void) pthread_join(pool->threads[thread], NULL);
    }

    pthread_cond_destroy(&pool->work_complete);
    pthread_cond_destroy(&pool->work_available);
    pthread_mutex_destroy(&pool->state_lock);
    pthread_mutex_destroy(&pool->run_lock);
    free(pool->slots);
    free(pool->threads);
    free(pool);
}

GeoThreadPool *geo_thread_pool_create(size_t worker_capacity)
{
    if (!worker_capacity || worker_capacity > SIZE_MAX / sizeof(pthread_t) + 1U ||
        worker_capacity > SIZE_MAX / sizeof(GeoThreadPoolSlot) + 1U) {
        return NULL;
    }

    GeoThreadPool *pool = calloc(1, sizeof(*pool));

    if (!pool) {
        return NULL;
    }

    size_t background_capacity = worker_capacity - 1U;
    pool->threads = background_capacity ? malloc(background_capacity * sizeof(*pool->threads)) : NULL;
    pool->slots = background_capacity ? malloc(background_capacity * sizeof(*pool->slots)) : NULL;

    if ((background_capacity && (!pool->threads || !pool->slots)) ||
        pthread_mutex_init(&pool->run_lock, NULL) != 0) {
        free(pool->slots);
        free(pool->threads);
        free(pool);

        return NULL;
    }

    if (pthread_mutex_init(&pool->state_lock, NULL) != 0) {
        pthread_mutex_destroy(&pool->run_lock);
        free(pool->slots);
        free(pool->threads);
        free(pool);

        return NULL;
    }

    if (pthread_cond_init(&pool->work_available, NULL) != 0) {
        pthread_mutex_destroy(&pool->state_lock);
        pthread_mutex_destroy(&pool->run_lock);
        free(pool->slots);
        free(pool->threads);
        free(pool);

        return NULL;
    }

    if (pthread_cond_init(&pool->work_complete, NULL) != 0) {
        pthread_cond_destroy(&pool->work_available);
        pthread_mutex_destroy(&pool->state_lock);
        pthread_mutex_destroy(&pool->run_lock);
        free(pool->slots);
        free(pool->threads);
        free(pool);

        return NULL;
    }

    size_t threads_started = 0;

    while (threads_started < background_capacity) {
        pool->slots[threads_started] = (GeoThreadPoolSlot) {
            .pool = pool,
            .worker_index = threads_started + 1U,
        };

        if (pthread_create(pool->threads + threads_started,
                           NULL,
                           geo_thread_pool_worker_main,
                           pool->slots + threads_started) != 0) {
            break;
        }

        threads_started++;
    }

    pool->capacity = threads_started + 1U;

    if (worker_capacity > 1U && pool->capacity == 1U) {
        geo_thread_pool_destroy_initialized(pool, 0);

        return NULL;
    }

    return pool;
}

void geo_thread_pool_destroy(GeoThreadPool *pool)
{
    if (!pool) {
        return;
    }

    pthread_mutex_lock(&pool->run_lock);
    size_t thread_count = pool->capacity - 1U;
    pthread_mutex_lock(&pool->state_lock);
    pool->stopping = true;
    pthread_cond_broadcast(&pool->work_available);
    pthread_mutex_unlock(&pool->state_lock);

    for (size_t thread = 0; thread < thread_count; ++thread) {
        (void) pthread_join(pool->threads[thread], NULL);
    }

    pthread_mutex_unlock(&pool->run_lock);
    pthread_cond_destroy(&pool->work_complete);
    pthread_cond_destroy(&pool->work_available);
    pthread_mutex_destroy(&pool->state_lock);
    pthread_mutex_destroy(&pool->run_lock);
    free(pool->slots);
    free(pool->threads);
    free(pool);
}

size_t geo_thread_pool_capacity(const GeoThreadPool *pool)
{
    return pool ? pool->capacity : 0;
}

bool geo_thread_pool_run(GeoThreadPool *pool,
                         size_t requested_workers,
                         GeoThreadPoolWorker worker,
                         void *context,
                         size_t *workers_used)
{
    if (!pool || !requested_workers || !worker || !workers_used || pthread_mutex_lock(&pool->run_lock) != 0) {
        return false;
    }

    size_t active_workers = requested_workers < pool->capacity ? requested_workers : pool->capacity;

    if (pthread_mutex_lock(&pool->state_lock) != 0) {
        pthread_mutex_unlock(&pool->run_lock);

        return false;
    }

    if (pool->stopping) {
        pthread_mutex_unlock(&pool->state_lock);
        pthread_mutex_unlock(&pool->run_lock);

        return false;
    }

    pool->worker = worker;
    pool->context = context;
    pool->active_workers = active_workers;
    pool->completed_workers = 0;
    pool->failed = false;
    pool->generation++;
    pthread_cond_broadcast(&pool->work_available);
    pthread_mutex_unlock(&pool->state_lock);

    bool succeeded = worker(context, 0, active_workers);

    pthread_mutex_lock(&pool->state_lock);

    while (pool->completed_workers + 1U < active_workers && !pool->stopping) {
        if (pthread_cond_wait(&pool->work_complete, &pool->state_lock) != 0) {
            pool->failed = true;
            break;
        }
    }

    succeeded = succeeded && !pool->failed && !pool->stopping;
    *workers_used = active_workers;
    pthread_mutex_unlock(&pool->state_lock);
    pthread_mutex_unlock(&pool->run_lock);

    return succeeded;
}
