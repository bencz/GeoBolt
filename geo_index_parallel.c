#include "geo_index.h"
#include "geo_index_private.h"

#include <stdlib.h>
#include <string.h>

#if defined(__unix__) || defined(__APPLE__)
#include <pthread.h>
#include <unistd.h>
#define GEO_PARALLEL_SORT_SUPPORTED 1
#else
#define GEO_PARALLEL_SORT_SUPPORTED 0
#endif

#define GEO_PARALLEL_RADIX_BITS 11
#define GEO_PARALLEL_RADIX_BUCKETS (1U << GEO_PARALLEL_RADIX_BITS)
#define GEO_PARALLEL_MIN_RECORDS_PER_THREAD 131072
#define GEO_PARALLEL_MAX_THREADS 32

#if GEO_PARALLEL_SORT_SUPPORTED

typedef struct GeoRadixPool GeoRadixPool;

typedef struct {
    GeoRadixPool *pool;
    size_t index;
    const GeoRecord *source;
    GeoRecord *destination;
    size_t begin;
    size_t end;
    unsigned shift;
    uint64_t bucket_mask;
    size_t *buckets;
} GeoRadixWorker;

typedef void (*GeoWorkerFunction)(GeoRadixWorker *worker);

struct GeoRadixPool {
    GeoRadixWorker *workers;
    pthread_t *threads;
    size_t thread_count;
    size_t active_thread_count;
    size_t created_threads;
    size_t completed_threads;
    uint64_t generation;
    GeoWorkerFunction function;
    pthread_mutex_t lock;
    pthread_cond_t work_available;
    pthread_cond_t work_complete;
    bool shutdown;
    bool lock_initialized;
    bool work_available_initialized;
    bool work_complete_initialized;
};

struct GeoParallelSorter {
    GeoRadixPool pool;
    GeoRadixWorker *workers;
    pthread_t *threads;
    size_t *worker_buckets;
    GeoRecord *scratch;
    size_t record_capacity;
    size_t thread_count;
};

static void radix_histogram_worker(GeoRadixWorker *worker)
{
    memset(worker->buckets, 0, GEO_PARALLEL_RADIX_BUCKETS * sizeof(*worker->buckets));

    for (size_t i = worker->begin; i < worker->end; ++i) {
        size_t bucket = (size_t) ((worker->source[i].z >> worker->shift) & worker->bucket_mask);

        worker->buckets[bucket]++;
    }
}

static void radix_scatter_worker(GeoRadixWorker *worker)
{
    for (size_t i = worker->begin; i < worker->end; ++i) {
        size_t bucket = (size_t) ((worker->source[i].z >> worker->shift) & worker->bucket_mask);

        worker->destination[worker->buckets[bucket]++] = worker->source[i];
    }
}

static void *radix_pool_worker_main(void *argument)
{
    GeoRadixWorker *worker = argument;
    GeoRadixPool *pool = worker->pool;
    uint64_t observed_generation = 0;

    for (;;) {
        pthread_mutex_lock(&pool->lock);

        while (!pool->shutdown && observed_generation == pool->generation) {
            pthread_cond_wait(&pool->work_available, &pool->lock);
        }

        if (pool->shutdown) {
            pthread_mutex_unlock(&pool->lock);

            return NULL;
        }

        observed_generation = pool->generation;
        GeoWorkerFunction function = pool->function;
        bool active = worker->index < pool->active_thread_count;

        pthread_mutex_unlock(&pool->lock);

        if (!active) {
            continue;
        }

        function(worker);
        pthread_mutex_lock(&pool->lock);
        pool->completed_threads++;

        if (pool->completed_threads == pool->active_thread_count) {
            pthread_cond_signal(&pool->work_complete);
        }

        pthread_mutex_unlock(&pool->lock);
    }
}

static void radix_pool_destroy(GeoRadixPool *pool)
{
    if (pool->created_threads && pool->lock_initialized) {
        pthread_mutex_lock(&pool->lock);
        pool->shutdown = true;

        if (pool->work_available_initialized) {
            pthread_cond_broadcast(&pool->work_available);
        }

        pthread_mutex_unlock(&pool->lock);

        for (size_t thread = 0; thread < pool->created_threads; ++thread) {
            (void) pthread_join(pool->threads[thread], NULL);
        }

        pool->created_threads = 0;
    }

    if (pool->work_complete_initialized) {
        pthread_cond_destroy(&pool->work_complete);
        pool->work_complete_initialized = false;
    }

    if (pool->work_available_initialized) {
        pthread_cond_destroy(&pool->work_available);
        pool->work_available_initialized = false;
    }

    if (pool->lock_initialized) {
        pthread_mutex_destroy(&pool->lock);
        pool->lock_initialized = false;
    }
}

static bool radix_pool_create(GeoRadixPool *pool,
                              GeoRadixWorker *workers,
                              pthread_t *threads,
                              size_t thread_count)
{
    memset(pool, 0, sizeof(*pool));
    pool->workers = workers;
    pool->threads = threads;
    pool->thread_count = thread_count;

    if (pthread_mutex_init(&pool->lock, NULL) != 0) {
        return false;
    }

    pool->lock_initialized = true;

    if (pthread_cond_init(&pool->work_available, NULL) != 0) {
        radix_pool_destroy(pool);

        return false;
    }

    pool->work_available_initialized = true;

    if (pthread_cond_init(&pool->work_complete, NULL) != 0) {
        radix_pool_destroy(pool);

        return false;
    }

    pool->work_complete_initialized = true;

    for (size_t thread = 0; thread < thread_count; ++thread) {
        workers[thread].pool = pool;
        workers[thread].index = thread;

        if (pthread_create(threads + thread, NULL, radix_pool_worker_main, workers + thread) != 0) {
            radix_pool_destroy(pool);

            return false;
        }

        pool->created_threads++;
    }

    return true;
}

static void radix_pool_run(GeoRadixPool *pool, GeoWorkerFunction function, size_t active_thread_count)
{
    pthread_mutex_lock(&pool->lock);
    pool->function = function;
    pool->active_thread_count = active_thread_count;
    pool->completed_threads = 0;
    pool->generation++;
    pthread_cond_broadcast(&pool->work_available);

    while (pool->completed_threads != active_thread_count) {
        pthread_cond_wait(&pool->work_complete, &pool->lock);
    }

    pthread_mutex_unlock(&pool->lock);
}

static size_t radix_effective_thread_count(size_t record_count, size_t requested_threads)
{
    bool automatic = requested_threads == 0;

    if (!requested_threads) {
        long online_cpus = sysconf(_SC_NPROCESSORS_ONLN);

        requested_threads = online_cpus > 0 ? (size_t) online_cpus : 1;
    }

    if (automatic && requested_threads > 8) {
        requested_threads = 8;
    }

    if (requested_threads > GEO_PARALLEL_MAX_THREADS) {
        requested_threads = GEO_PARALLEL_MAX_THREADS;
    }

    size_t work_limited_threads = record_count / GEO_PARALLEL_MIN_RECORDS_PER_THREAD;

    if (requested_threads > work_limited_threads) {
        requested_threads = work_limited_threads;
    }

    return requested_threads > 1 ? requested_threads : 1;
}

GeoParallelSorter *geo_parallel_sorter_create(size_t record_capacity, size_t requested_threads)
{
    if (!record_capacity || record_capacity > SIZE_MAX / sizeof(GeoRecord)) {
        return NULL;
    }

    GeoParallelSorter *sorter = calloc(1, sizeof(*sorter));

    if (!sorter) {
        return NULL;
    }

    sorter->record_capacity = record_capacity;
    sorter->thread_count = radix_effective_thread_count(record_capacity, requested_threads);

    if (sorter->thread_count == 1) {
        return sorter;
    }

    size_t bucket_slots = sorter->thread_count * GEO_PARALLEL_RADIX_BUCKETS;
    size_t bucket_bytes = bucket_slots * sizeof(*sorter->worker_buckets);

    sorter->scratch = malloc(record_capacity * sizeof(*sorter->scratch));
    sorter->worker_buckets = aligned_alloc(64U, bucket_bytes);
    sorter->workers = calloc(sorter->thread_count, sizeof(*sorter->workers));
    sorter->threads = malloc(sorter->thread_count * sizeof(*sorter->threads));

    if (!sorter->scratch || !sorter->worker_buckets || !sorter->workers || !sorter->threads ||
        !radix_pool_create(&sorter->pool, sorter->workers, sorter->threads, sorter->thread_count)) {
        geo_parallel_sorter_destroy(sorter);

        return NULL;
    }

    return sorter;
}

void geo_parallel_sorter_destroy(GeoParallelSorter *sorter)
{
    if (!sorter) {
        return;
    }

    if (sorter->thread_count > 1) {
        radix_pool_destroy(&sorter->pool);
    }

    free(sorter->threads);
    free(sorter->workers);
    free(sorter->worker_buckets);
    free(sorter->scratch);
    free(sorter);
}

static bool radix_sort_records_with_sorter(GeoRecord *records, size_t count, GeoParallelSorter *sorter)
{
    size_t active_thread_count = radix_effective_thread_count(count, sorter->thread_count);

    if (active_thread_count == 1) {
        return false;
    }

    size_t record_bytes = count * sizeof(*records);
    uint64_t varying_bits = 0;
    uint64_t first_z = records[0].z;

    for (size_t i = 1; i < count; ++i) {
        varying_bits |= first_z ^ records[i].z;
    }

    GeoRecord *source = records;
    GeoRecord *destination = sorter->scratch;
    unsigned passes = 0;

    for (unsigned shift = 0; shift < 64; shift += GEO_PARALLEL_RADIX_BITS) {
        unsigned pass_bits = 64U - shift < GEO_PARALLEL_RADIX_BITS ? 64U - shift : GEO_PARALLEL_RADIX_BITS;
        size_t bucket_count = (size_t) 1 << pass_bits;
        uint64_t bucket_mask = bucket_count - 1;

        if (!((varying_bits >> shift) & bucket_mask)) {
            continue;
        }

        for (size_t thread = 0; thread < active_thread_count; ++thread) {
            GeoRadixWorker *worker = sorter->workers + thread;

            worker->source = source;
            worker->destination = destination;
            worker->begin = count * thread / active_thread_count;
            worker->end = count * (thread + 1) / active_thread_count;
            worker->shift = shift;
            worker->bucket_mask = bucket_mask;
            worker->buckets = sorter->worker_buckets + thread * GEO_PARALLEL_RADIX_BUCKETS;
        }

        radix_pool_run(&sorter->pool, radix_histogram_worker, active_thread_count);

        size_t offset = 0;

        for (size_t bucket = 0; bucket < bucket_count; ++bucket) {
            for (size_t thread = 0; thread < active_thread_count; ++thread) {
                size_t slot = thread * GEO_PARALLEL_RADIX_BUCKETS + bucket;
                size_t records_in_bucket = sorter->worker_buckets[slot];

                sorter->worker_buckets[slot] = offset;
                offset += records_in_bucket;
            }
        }

        radix_pool_run(&sorter->pool, radix_scatter_worker, active_thread_count);

        GeoRecord *temporary = source;

        source = destination;
        destination = temporary;
        passes++;
    }

    if (passes & 1U) {
        memcpy(records, source, record_bytes);
    }

    return true;
}

#else

struct GeoParallelSorter {
    size_t record_capacity;
};

GeoParallelSorter *geo_parallel_sorter_create(size_t record_capacity, size_t requested_threads)
{
    (void) requested_threads;

    if (!record_capacity) {
        return NULL;
    }

    GeoParallelSorter *sorter = malloc(sizeof(*sorter));

    if (sorter) {
        sorter->record_capacity = record_capacity;
    }

    return sorter;
}

void geo_parallel_sorter_destroy(GeoParallelSorter *sorter)
{
    free(sorter);
}

#endif

bool geo_index_sort_transient_with_sorter(GeoIndex *index, GeoParallelSorter *sorter)
{
    if (!index || !sorter || index->read_only || index->count > sorter->record_capacity) {
        return false;
    }

    if (index->sorted) {
        return true;
    }

#if GEO_PARALLEL_SORT_SUPPORTED
    if (index->count > 1 && sorter->thread_count > 1 && radix_sort_records_with_sorter(index->records, index->count, sorter)) {
        index->sorted = true;

        return true;
    }
#endif

    geo_index_sort_transient(index);

    return index->sorted;
}

bool geo_index_sort_transient_parallel(GeoIndex *index, size_t thread_count)
{
    if (!index || index->read_only) {
        return false;
    }

    if (index->sorted) {
        return true;
    }

#if GEO_PARALLEL_SORT_SUPPORTED
    size_t capacity = index->count ? index->count : 1;
    GeoParallelSorter *sorter = geo_parallel_sorter_create(capacity, thread_count);

    if (sorter) {
        bool succeeded = geo_index_sort_transient_with_sorter(index, sorter);

        geo_parallel_sorter_destroy(sorter);

        if (succeeded) {
            return true;
        }
    }
#else
    (void) thread_count;
#endif

    geo_index_sort_transient(index);

    return index->sorted;
}

bool geo_index_build_parallel(GeoIndex *index, size_t thread_count)
{
    if (!index || index->read_only) {
        return false;
    }

    if (index->sorted) {
        return true;
    }

    if (!geo_index_sort_transient_parallel(index, thread_count)) {
        return false;
    }

    return geo_index_finalize_sorted(index);
}
