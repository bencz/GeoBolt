#ifndef GEO_JOB_QUEUE_H
#define GEO_JOB_QUEUE_H

#include <stdbool.h>
#include <stddef.h>

typedef struct GeoJobQueue GeoJobQueue;
typedef void (*GeoJobFunction)(void *context);

GeoJobQueue *geo_job_queue_create(size_t thread_count, size_t queue_capacity);
void geo_job_queue_destroy(GeoJobQueue *queue);

// Never blocks a reactor or producer. False means that bounded backpressure must be applied by the caller.
bool geo_job_queue_try_submit(GeoJobQueue *queue, GeoJobFunction function, void *context);
bool geo_job_queue_wait_idle(GeoJobQueue *queue);
size_t geo_job_queue_capacity(const GeoJobQueue *queue);

#endif
