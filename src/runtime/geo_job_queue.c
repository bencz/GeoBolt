#include "geo_job_queue.h"

#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
    GeoJobFunction function;
    void *context;
} GeoJob;

struct GeoJobQueue {
    pthread_mutex_t lock;
    pthread_cond_t work_available;
    pthread_cond_t idle;
    pthread_t *threads;
    GeoJob *jobs;
    size_t thread_count;
    size_t queue_capacity;
    size_t head;
    size_t count;
    size_t active;
    bool stopping;
};

static void *job_queue_worker_main(void *context)
{
    GeoJobQueue *queue = context;

    if (pthread_mutex_lock(&queue->lock) != 0) {
        return NULL;
    }

    for (;;) {
        while (!queue->count && !queue->stopping) {
            if (pthread_cond_wait(&queue->work_available, &queue->lock) != 0) {
                queue->stopping = true;
                break;
            }
        }

        if (!queue->count && queue->stopping) {
            break;
        }

        GeoJob job = queue->jobs[queue->head];

        queue->head = (queue->head + 1U) % queue->queue_capacity;
        queue->count--;
        queue->active++;
        pthread_mutex_unlock(&queue->lock);

        job.function(job.context);

        pthread_mutex_lock(&queue->lock);
        queue->active--;

        if (!queue->count && !queue->active) {
            pthread_cond_broadcast(&queue->idle);
        }
    }

    pthread_mutex_unlock(&queue->lock);

    return NULL;
}

static void job_queue_release(GeoJobQueue *queue, size_t threads_started)
{
    if (!queue) {
        return;
    }

    if (threads_started) {
        pthread_mutex_lock(&queue->lock);
        queue->stopping = true;
        pthread_cond_broadcast(&queue->work_available);
        pthread_mutex_unlock(&queue->lock);
    }

    for (size_t thread = 0; thread < threads_started; ++thread) {
        (void) pthread_join(queue->threads[thread], NULL);
    }

    pthread_cond_destroy(&queue->idle);
    pthread_cond_destroy(&queue->work_available);
    pthread_mutex_destroy(&queue->lock);
    free(queue->jobs);
    free(queue->threads);
    free(queue);
}

GeoJobQueue *geo_job_queue_create(size_t thread_count, size_t queue_capacity)
{
    if (!thread_count || !queue_capacity ||
        thread_count > SIZE_MAX / sizeof(pthread_t) ||
        queue_capacity > SIZE_MAX / sizeof(GeoJob)) {
        return NULL;
    }

    GeoJobQueue *queue = calloc(1, sizeof(*queue));

    if (!queue) {
        return NULL;
    }

    queue->threads = malloc(thread_count * sizeof(*queue->threads));
    queue->jobs = malloc(queue_capacity * sizeof(*queue->jobs));

    if (!queue->threads || !queue->jobs || pthread_mutex_init(&queue->lock, NULL) != 0) {
        free(queue->jobs);
        free(queue->threads);
        free(queue);

        return NULL;
    }

    if (pthread_cond_init(&queue->work_available, NULL) != 0) {
        pthread_mutex_destroy(&queue->lock);
        free(queue->jobs);
        free(queue->threads);
        free(queue);

        return NULL;
    }

    if (pthread_cond_init(&queue->idle, NULL) != 0) {
        pthread_cond_destroy(&queue->work_available);
        pthread_mutex_destroy(&queue->lock);
        free(queue->jobs);
        free(queue->threads);
        free(queue);

        return NULL;
    }

    queue->thread_count = thread_count;
    queue->queue_capacity = queue_capacity;
    size_t threads_started = 0;

    while (threads_started < thread_count &&
           pthread_create(queue->threads + threads_started, NULL, job_queue_worker_main, queue) == 0) {
        threads_started++;
    }

    if (threads_started != thread_count) {
        job_queue_release(queue, threads_started);

        return NULL;
    }

    return queue;
}

void geo_job_queue_destroy(GeoJobQueue *queue)
{
    if (!queue) {
        return;
    }

    if (pthread_mutex_lock(&queue->lock) == 0) {
        queue->stopping = true;
        pthread_cond_broadcast(&queue->work_available);
        pthread_mutex_unlock(&queue->lock);
    }

    for (size_t thread = 0; thread < queue->thread_count; ++thread) {
        (void) pthread_join(queue->threads[thread], NULL);
    }

    pthread_cond_destroy(&queue->idle);
    pthread_cond_destroy(&queue->work_available);
    pthread_mutex_destroy(&queue->lock);
    free(queue->jobs);
    free(queue->threads);
    free(queue);
}

bool geo_job_queue_try_submit(GeoJobQueue *queue, GeoJobFunction function, void *context)
{
    if (!queue || !function || pthread_mutex_lock(&queue->lock) != 0) {
        return false;
    }

    if (queue->stopping || queue->count == queue->queue_capacity) {
        pthread_mutex_unlock(&queue->lock);

        return false;
    }

    size_t tail = (queue->head + queue->count) % queue->queue_capacity;

    queue->jobs[tail] = (GeoJob) {
        .function = function,
        .context = context,
    };
    queue->count++;
    pthread_cond_signal(&queue->work_available);
    pthread_mutex_unlock(&queue->lock);

    return true;
}

bool geo_job_queue_wait_idle(GeoJobQueue *queue)
{
    if (!queue || pthread_mutex_lock(&queue->lock) != 0) {
        return false;
    }

    while ((queue->count || queue->active) && !queue->stopping) {
        if (pthread_cond_wait(&queue->idle, &queue->lock) != 0) {
            pthread_mutex_unlock(&queue->lock);

            return false;
        }
    }

    bool succeeded = !queue->count && !queue->active;

    pthread_mutex_unlock(&queue->lock);

    return succeeded;
}

size_t geo_job_queue_capacity(const GeoJobQueue *queue)
{
    return queue ? queue->queue_capacity : 0U;
}
