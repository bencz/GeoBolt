#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "geo_index.h"
#include "geo_index_private.h"

#include <pthread.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#if defined(__linux__)
#include <sched.h>
#define GEO_BATCH_AFFINITY_SUPPORTED 1
#else
#define GEO_BATCH_AFFINITY_SUPPORTED 0
#endif

#define GEO_BATCH_DEFAULT_CHUNK 8U
#define GEO_BATCH_MAX_THREADS 256U
#define GEO_BATCH_CACHE_LINE 64U

typedef enum {
    GEO_BATCH_PHASE_IDLE,
    GEO_BATCH_PHASE_COUNT,
    GEO_BATCH_PHASE_FILL,
} GeoBatchPhase;

typedef struct {
    _Alignas(GEO_BATCH_CACHE_LINE) atomic_size_t next_query;
} GeoBatchQueue;

typedef struct {
    bool (*count)(const void *source, double latitude, double longitude, double radius_km, size_t *count);
    bool (*fill)(const void *source,
                 double latitude,
                 double longitude,
                 double radius_km,
                 uint64_t *ids,
                 size_t capacity,
                 size_t *count);
    bool (*snapshot_acquire)(const void *source);
    void (*snapshot_release)(const void *source);
} GeoBatchSourceOps;

_Static_assert(sizeof(GeoBatchQueue) >= GEO_BATCH_CACHE_LINE, "GeoBatchQueue must not share a cache line");
_Static_assert(sizeof(GeoBatchQueue) % GEO_BATCH_CACHE_LINE == 0, "GeoBatchQueue size must satisfy aligned_alloc");

typedef struct {
    struct GeoQueryExecutor *executor;
    size_t local_index;
    unsigned cpu;
} GeoQueryWorker;

struct GeoQueryExecutor {
    const GeoIndex **replicas;
    const GeoSegmentSet *segment_set;
    const GeoBatchSourceOps *source_ops;
    pthread_t *threads;
    GeoQueryWorker *workers;
    size_t replica_count;
    size_t thread_count;
    size_t scheduling_chunk;
    size_t created_threads;
    pthread_mutex_t submit_lock;
    pthread_mutex_t state_lock;
    pthread_cond_t work_available;
    pthread_cond_t work_complete;
    GeoBatchQueue *queues;
    atomic_bool failed;
    const double *latitudes;
    const double *longitudes;
    const double *radii_km;
    size_t *counts;
    size_t *partial_offsets;
    size_t count_slot_capacity;
    GeoBatchIdResult *result;
    size_t query_count;
    size_t completed_workers;
    uint64_t generation;
    GeoBatchPhase phase;
    bool pin_workers;
    bool sharded;
    bool shutdown;
    bool submit_lock_initialized;
    bool state_lock_initialized;
    bool work_available_initialized;
    bool work_complete_initialized;
};

static bool batch_index_count(const void *source,
                              double latitude,
                              double longitude,
                              double radius_km,
                              size_t *count)
{
    return geo_search_radius_count(source, latitude, longitude, radius_km, count, NULL);
}

static bool batch_index_fill(const void *source,
                             double latitude,
                             double longitude,
                             double radius_km,
                             uint64_t *ids,
                             size_t capacity,
                             size_t *count)
{
    return geo_search_radius_ids_into(source, latitude, longitude, radius_km, ids, capacity, count, NULL);
}

static bool batch_segment_count(const void *source,
                                double latitude,
                                double longitude,
                                double radius_km,
                                size_t *count)
{
    return geo_segment_set_search_radius_count_snapshot(source, latitude, longitude, radius_km, count);
}

static bool batch_segment_fill(const void *source,
                               double latitude,
                               double longitude,
                               double radius_km,
                               uint64_t *ids,
                               size_t capacity,
                               size_t *count)
{
    return geo_segment_set_search_radius_ids_into_snapshot(source, latitude, longitude, radius_km, ids, capacity, count);
}

static bool batch_segment_snapshot_acquire(const void *source)
{
    return geo_segment_set_read_snapshot_acquire(source);
}

static void batch_segment_snapshot_release(const void *source)
{
    geo_segment_set_read_snapshot_release(source);
}

static const GeoBatchSourceOps g_index_source_ops = {
    .count = batch_index_count,
    .fill = batch_index_fill,
};

static const GeoBatchSourceOps g_segment_source_ops = {
    .count = batch_segment_count,
    .fill = batch_segment_fill,
    .snapshot_acquire = batch_segment_snapshot_acquire,
    .snapshot_release = batch_segment_snapshot_release,
};

static bool batch_size_add(size_t first, size_t second, size_t *result)
{
    if (second > SIZE_MAX - first) {
        return false;
    }

    *result = first + second;

    return true;
}

static bool batch_size_multiply(size_t first, size_t second, size_t *result)
{
    if (first && second > SIZE_MAX / first) {
        return false;
    }

    *result = first * second;

    return true;
}

GeoBatchIdResult *geo_batch_id_result_create(size_t query_capacity, size_t id_capacity)
{
    GeoBatchIdResult *result = calloc(1, sizeof(*result));

    if (!result) {
        return NULL;
    }

    if (query_capacity) {
        size_t offset_count;
        size_t offset_bytes;

        if (!batch_size_add(query_capacity, 1, &offset_count) ||
            !batch_size_multiply(offset_count, sizeof(*result->offsets), &offset_bytes)) {
            free(result);

            return NULL;
        }

        result->offsets = malloc(offset_bytes);

        if (!result->offsets) {
            free(result);

            return NULL;
        }

        result->offset_capacity = offset_count;
    }

    if (id_capacity) {
        size_t id_bytes;

        if (!batch_size_multiply(id_capacity, sizeof(*result->ids), &id_bytes)) {
            geo_batch_id_result_destroy(result);

            return NULL;
        }

        result->ids = malloc(id_bytes);

        if (!result->ids) {
            geo_batch_id_result_destroy(result);

            return NULL;
        }

        result->id_capacity = id_capacity;
    }

    return result;
}

void geo_batch_id_result_destroy(GeoBatchIdResult *result)
{
    if (result) {
        free(result->offsets);
        free(result->ids);
        free(result);
    }
}

void geo_batch_id_result_clear(GeoBatchIdResult *result)
{
    if (result) {
        result->query_count = 0;
        result->id_count = 0;
    }
}

static bool batch_result_reserve_offsets(GeoBatchIdResult *result, size_t query_count)
{
    size_t required;

    if (!batch_size_add(query_count, 1, &required)) {
        return false;
    }

    if (required <= result->offset_capacity) {
        return true;
    }

    size_t bytes;

    if (!batch_size_multiply(required, sizeof(*result->offsets), &bytes)) {
        return false;
    }

    size_t *offsets = realloc(result->offsets, bytes);

    if (!offsets) {
        return false;
    }

    result->offsets = offsets;
    result->offset_capacity = required;

    return true;
}

static bool batch_result_reserve_ids(GeoBatchIdResult *result, size_t id_count)
{
    if (id_count <= result->id_capacity) {
        return true;
    }

    size_t bytes;

    if (!batch_size_multiply(id_count, sizeof(*result->ids), &bytes)) {
        return false;
    }

    uint64_t *ids = realloc(result->ids, bytes);

    if (!ids) {
        return false;
    }

    result->ids = ids;
    result->id_capacity = id_count;

    return true;
}

static bool batch_executor_reserve_count_slots(GeoQueryExecutor *executor, size_t required)
{
    if (required <= executor->count_slot_capacity) {
        return true;
    }

    size_t bytes;

    if (!batch_size_multiply(required, sizeof(*executor->counts), &bytes)) {
        return false;
    }

    size_t *counts = malloc(bytes);
    size_t *partial_offsets = executor->sharded ? malloc(bytes) : NULL;

    if (!counts || (executor->sharded && !partial_offsets)) {
        free(partial_offsets);
        free(counts);

        return false;
    }

    free(executor->partial_offsets);
    free(executor->counts);
    executor->counts = counts;
    executor->partial_offsets = partial_offsets;
    executor->count_slot_capacity = required;

    return true;
}

static size_t batch_default_thread_count(void)
{
    long online_cpus = sysconf(_SC_NPROCESSORS_ONLN);
    size_t thread_count = online_cpus > 0 ? (size_t) online_cpus : 1;

    if (thread_count > GEO_BATCH_MAX_THREADS) {
        thread_count = GEO_BATCH_MAX_THREADS;
    }

    return thread_count;
}

static void batch_pin_worker(unsigned cpu)
{
#if GEO_BATCH_AFFINITY_SUPPORTED
    long online_cpus = sysconf(_SC_NPROCESSORS_ONLN);

    if (online_cpus <= 0) {
        return;
    }

    cpu_set_t affinity;

    CPU_ZERO(&affinity);
    CPU_SET((size_t) cpu % (size_t) online_cpus, &affinity);
    (void) pthread_setaffinity_np(pthread_self(), sizeof(affinity), &affinity);
#else
    (void) cpu;
#endif
}

static void batch_execute_count(GeoQueryWorker *worker)
{
    GeoQueryExecutor *executor = worker->executor;
    const void *source = executor->segment_set ? (const void *) executor->segment_set
                                               : (const void *) executor->replicas[worker->local_index];
    size_t queue = executor->sharded ? worker->local_index : 0;

    for (;;) {
        size_t begin = atomic_fetch_add_explicit(&executor->queues[queue].next_query,
                                                 executor->scheduling_chunk,
                                                 memory_order_relaxed);

        if (begin >= executor->query_count) {
            return;
        }

        size_t end = begin + executor->scheduling_chunk;

        if (end < begin || end > executor->query_count) {
            end = executor->query_count;
        }

        for (size_t query = begin; query < end; ++query) {
            size_t count_slot = executor->sharded
                                    ? worker->local_index * executor->query_count + query
                                    : query;

            if (!executor->source_ops->count(source,
                                             executor->latitudes[query],
                                             executor->longitudes[query],
                                             executor->radii_km[query],
                                             executor->counts + count_slot)) {
                atomic_store_explicit(&executor->failed, true, memory_order_relaxed);
            }
        }
    }
}

static void batch_execute_fill(GeoQueryWorker *worker)
{
    GeoQueryExecutor *executor = worker->executor;
    const void *source = executor->segment_set ? (const void *) executor->segment_set
                                               : (const void *) executor->replicas[worker->local_index];
    size_t queue = executor->sharded ? worker->local_index : 0;

    for (;;) {
        size_t begin = atomic_fetch_add_explicit(&executor->queues[queue].next_query,
                                                 executor->scheduling_chunk,
                                                 memory_order_relaxed);

        if (begin >= executor->query_count) {
            return;
        }

        size_t end = begin + executor->scheduling_chunk;

        if (end < begin || end > executor->query_count) {
            end = executor->query_count;
        }

        for (size_t query = begin; query < end; ++query) {
            size_t count_slot = executor->sharded
                                    ? worker->local_index * executor->query_count + query
                                    : query;
            size_t offset = executor->sharded
                                ? executor->partial_offsets[count_slot]
                                : executor->result->offsets[query];
            size_t capacity = executor->counts[count_slot];
            size_t written = 0;
            uint64_t *destination = executor->result->ids ? executor->result->ids + offset : NULL;

            if (!executor->source_ops->fill(source,
                                            executor->latitudes[query],
                                            executor->longitudes[query],
                                            executor->radii_km[query],
                                            destination,
                                            capacity,
                                            &written) ||
                written != capacity) {
                atomic_store_explicit(&executor->failed, true, memory_order_relaxed);
            }
        }
    }
}

static void *batch_worker_main(void *argument)
{
    GeoQueryWorker *worker = argument;
    GeoQueryExecutor *executor = worker->executor;
    uint64_t observed_generation = 0;

    if (executor->pin_workers) {
        batch_pin_worker(worker->cpu);
    }

    for (;;) {
        pthread_mutex_lock(&executor->state_lock);

        while (!executor->shutdown && observed_generation == executor->generation) {
            pthread_cond_wait(&executor->work_available, &executor->state_lock);
        }

        if (executor->shutdown) {
            pthread_mutex_unlock(&executor->state_lock);

            return NULL;
        }

        observed_generation = executor->generation;
        GeoBatchPhase phase = executor->phase;

        pthread_mutex_unlock(&executor->state_lock);

        if (phase == GEO_BATCH_PHASE_COUNT) {
            batch_execute_count(worker);
        } else if (phase == GEO_BATCH_PHASE_FILL) {
            batch_execute_fill(worker);
        }

        pthread_mutex_lock(&executor->state_lock);
        executor->completed_workers++;

        if (executor->completed_workers == executor->thread_count) {
            pthread_cond_signal(&executor->work_complete);
        }

        pthread_mutex_unlock(&executor->state_lock);
    }
}

static void batch_executor_release(GeoQueryExecutor *executor)
{
    if (!executor) {
        return;
    }

    if (executor->created_threads) {
        pthread_mutex_lock(&executor->state_lock);
        executor->shutdown = true;
        pthread_cond_broadcast(&executor->work_available);
        pthread_mutex_unlock(&executor->state_lock);

        for (size_t thread = 0; thread < executor->created_threads; ++thread) {
            (void) pthread_join(executor->threads[thread], NULL);
        }
    }

    if (executor->work_complete_initialized) {
        pthread_cond_destroy(&executor->work_complete);
    }

    if (executor->work_available_initialized) {
        pthread_cond_destroy(&executor->work_available);
    }

    if (executor->state_lock_initialized) {
        pthread_mutex_destroy(&executor->state_lock);
    }

    if (executor->submit_lock_initialized) {
        pthread_mutex_destroy(&executor->submit_lock);
    }

    free(executor->workers);
    free(executor->threads);
    free(executor->queues);
    free(executor->replicas);
    free(executor->partial_offsets);
    free(executor->counts);
    free(executor);
}

static GeoQueryExecutor *batch_executor_create(const GeoIndex *const *replicas,
                                               size_t replica_count,
                                               const GeoSegmentSet *segment_set,
                                               const GeoBatchSourceOps *source_ops,
                                               const GeoQueryExecutorConfig *config,
                                               bool sharded)
{
    if (!source_ops || !replica_count || (segment_set && (replicas || replica_count != 1U)) || (!segment_set && !replicas)) {
        return NULL;
    }

    if (!segment_set) {
        for (size_t replica = 0; replica < replica_count; ++replica) {
            if (!replicas[replica] || !replicas[replica]->sorted) {
                return NULL;
            }
        }
    }

    size_t thread_count = config && config->thread_count ? config->thread_count : batch_default_thread_count();

    if (!thread_count || thread_count > GEO_BATCH_MAX_THREADS) {
        return NULL;
    }

    if (sharded && thread_count < replica_count) {
        return NULL;
    }

    GeoQueryExecutor *executor = calloc(1, sizeof(*executor));

    if (!executor) {
        return NULL;
    }

    executor->replicas = segment_set ? NULL : malloc(replica_count * sizeof(*executor->replicas));
    executor->threads = malloc(thread_count * sizeof(*executor->threads));
    executor->workers = malloc(thread_count * sizeof(*executor->workers));
    executor->queues = aligned_alloc(GEO_BATCH_CACHE_LINE,
                                     (sharded ? replica_count : 1) * sizeof(*executor->queues));

    if ((!segment_set && !executor->replicas) || !executor->threads || !executor->workers || !executor->queues) {
        batch_executor_release(executor);

        return NULL;
    }

    if (replicas) {
        memcpy(executor->replicas, replicas, replica_count * sizeof(*executor->replicas));
    }

    executor->segment_set = segment_set;
    executor->source_ops = source_ops;
    executor->replica_count = replica_count;
    executor->thread_count = thread_count;
    executor->scheduling_chunk = config && config->scheduling_chunk ? config->scheduling_chunk : GEO_BATCH_DEFAULT_CHUNK;
    executor->pin_workers = config && config->pin_workers;
    executor->sharded = sharded;

    if (pthread_mutex_init(&executor->submit_lock, NULL) != 0) {
        batch_executor_release(executor);

        return NULL;
    }

    executor->submit_lock_initialized = true;

    if (pthread_mutex_init(&executor->state_lock, NULL) != 0) {
        batch_executor_release(executor);

        return NULL;
    }

    executor->state_lock_initialized = true;

    if (pthread_cond_init(&executor->work_available, NULL) != 0) {
        batch_executor_release(executor);

        return NULL;
    }

    executor->work_available_initialized = true;

    if (pthread_cond_init(&executor->work_complete, NULL) != 0) {
        batch_executor_release(executor);

        return NULL;
    }

    executor->work_complete_initialized = true;
    size_t queue_count = sharded ? replica_count : 1;

    for (size_t queue = 0; queue < queue_count; ++queue) {
        atomic_init(&executor->queues[queue].next_query, 0);
    }

    atomic_init(&executor->failed, false);

    for (size_t thread = 0; thread < thread_count; ++thread) {
        executor->workers[thread] = (GeoQueryWorker) {
            .executor = executor,
            .local_index = thread % replica_count,
            .cpu = config && config->worker_cpus && thread < config->worker_cpu_count
                       ? config->worker_cpus[thread]
                       : (unsigned) thread,
        };

        if (pthread_create(executor->threads + thread, NULL, batch_worker_main, executor->workers + thread) != 0) {
            batch_executor_release(executor);

            return NULL;
        }

        executor->created_threads++;
    }

    return executor;
}

GeoQueryExecutor *geo_query_executor_create_replicated(const GeoIndex *const *replicas,
                                                       size_t replica_count,
                                                       const GeoQueryExecutorConfig *config)
{
    return batch_executor_create(replicas, replica_count, NULL, &g_index_source_ops, config, false);
}

GeoQueryExecutor *geo_query_executor_create_sharded(const GeoIndex *const *shards,
                                                    size_t shard_count,
                                                    const GeoQueryExecutorConfig *config)
{
    return batch_executor_create(shards, shard_count, NULL, &g_index_source_ops, config, true);
}

GeoQueryExecutor *geo_query_executor_create(const GeoIndex *index, const GeoQueryExecutorConfig *config)
{
    return geo_query_executor_create_replicated(&index, 1, config);
}

GeoQueryExecutor *geo_query_executor_create_segment_set(const GeoSegmentSet *set,
                                                        const GeoQueryExecutorConfig *config)
{
    if (!set) {
        return NULL;
    }

    return batch_executor_create(NULL, 1, set, &g_segment_source_ops, config, false);
}

void geo_query_executor_destroy(GeoQueryExecutor *executor)
{
    batch_executor_release(executor);
}

static bool batch_executor_dispatch(GeoQueryExecutor *executor, GeoBatchPhase phase)
{
    pthread_mutex_lock(&executor->state_lock);
    executor->phase = phase;
    executor->completed_workers = 0;
    size_t queue_count = executor->sharded ? executor->replica_count : 1;

    for (size_t queue = 0; queue < queue_count; ++queue) {
        atomic_store_explicit(&executor->queues[queue].next_query, 0, memory_order_relaxed);
    }
    atomic_store_explicit(&executor->failed, false, memory_order_relaxed);
    executor->generation++;
    pthread_cond_broadcast(&executor->work_available);

    while (executor->completed_workers != executor->thread_count) {
        pthread_cond_wait(&executor->work_complete, &executor->state_lock);
    }

    pthread_mutex_unlock(&executor->state_lock);

    return !atomic_load_explicit(&executor->failed, memory_order_relaxed);
}

bool geo_query_executor_search_radius_ids(GeoQueryExecutor *executor,
                                          const double *latitudes,
                                          const double *longitudes,
                                          const double *radii_km,
                                          size_t query_count,
                                          GeoBatchIdResult *result)
{
    if (!executor || !result ||
        (query_count && (!latitudes || !longitudes || !radii_km)) ||
        pthread_mutex_lock(&executor->submit_lock) != 0) {
        return false;
    }

    geo_batch_id_result_clear(result);

    if (!query_count) {
        pthread_mutex_unlock(&executor->submit_lock);

        return true;
    }

    size_t count_slots;
    size_t slots_per_query = executor->sharded ? executor->replica_count : 1;
    bool succeeded = batch_result_reserve_offsets(result, query_count) &&
                     batch_size_multiply(query_count, slots_per_query, &count_slots) &&
                     batch_executor_reserve_count_slots(executor, count_slots);

    bool snapshot_acquired = false;

    if (succeeded && executor->source_ops->snapshot_acquire) {
        const void *source = executor->segment_set ? (const void *) executor->segment_set
                                                   : (const void *) executor->replicas[0];

        snapshot_acquired = executor->source_ops->snapshot_acquire(source);
        succeeded = snapshot_acquired;
    }

    if (!succeeded) {
        pthread_mutex_unlock(&executor->submit_lock);

        return false;
    }

    executor->latitudes = latitudes;
    executor->longitudes = longitudes;
    executor->radii_km = radii_km;
    executor->result = result;
    executor->query_count = query_count;

    succeeded = batch_executor_dispatch(executor, GEO_BATCH_PHASE_COUNT);
    result->offsets[0] = 0;

    for (size_t query = 0; succeeded && query < query_count; ++query) {
        size_t next_offset = result->offsets[query];

        for (size_t slot = 0; succeeded && slot < slots_per_query; ++slot) {
            size_t count_slot = executor->sharded ? slot * query_count + query : query;

            if (executor->sharded) {
                executor->partial_offsets[count_slot] = next_offset;
            }

            succeeded = batch_size_add(next_offset, executor->counts[count_slot], &next_offset);
        }

        result->offsets[query + 1] = next_offset;
    }

    size_t id_count = succeeded ? result->offsets[query_count] : 0;

    if (succeeded) {
        succeeded = batch_result_reserve_ids(result, id_count);
    }

    if (succeeded) {
        succeeded = batch_executor_dispatch(executor, GEO_BATCH_PHASE_FILL);
    }

    if (succeeded) {
        result->query_count = query_count;
        result->id_count = id_count;
    } else {
        geo_batch_id_result_clear(result);
    }

    if (snapshot_acquired) {
        const void *source = executor->segment_set ? (const void *) executor->segment_set
                                                   : (const void *) executor->replicas[0];

        executor->source_ops->snapshot_release(source);
    }

    executor->latitudes = NULL;
    executor->longitudes = NULL;
    executor->radii_km = NULL;
    executor->result = NULL;
    executor->query_count = 0;

    pthread_mutex_unlock(&executor->submit_lock);

    return succeeded;
}
