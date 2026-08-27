#ifndef GEO_THREAD_POOL_H
#define GEO_THREAD_POOL_H

#include <stdbool.h>
#include <stddef.h>

typedef struct GeoThreadPool GeoThreadPool;
typedef bool (*GeoThreadPoolWorker)(void *context, size_t worker_index, size_t worker_count);

GeoThreadPool *geo_thread_pool_create(size_t worker_capacity);
void geo_thread_pool_destroy(GeoThreadPool *pool);
size_t geo_thread_pool_capacity(const GeoThreadPool *pool);

// Runs one callback on each selected worker. The calling thread is worker zero and participates in the operation.
bool geo_thread_pool_run(GeoThreadPool *pool,
                         size_t requested_workers,
                         GeoThreadPoolWorker worker,
                         void *context,
                         size_t *workers_used);

#endif
