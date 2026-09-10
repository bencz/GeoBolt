#include "geo_visibility_gate.h"

#include <limits.h>
#include <string.h>

static bool visibility_gate_writer_arrive(GeoVisibilityGate *gate)
{
    unsigned waiting = atomic_load_explicit(&gate->waiting_writers, memory_order_relaxed);

    for (;;) {
        if (waiting == UINT_MAX) {
            return false;
        }

        if (atomic_compare_exchange_weak_explicit(&gate->waiting_writers,
                                                  &waiting,
                                                  waiting + 1U,
                                                  memory_order_acq_rel,
                                                  memory_order_relaxed)) {
            return true;
        }
    }
}

static bool visibility_gate_writer_depart(GeoVisibilityGate *gate)
{
    if (pthread_mutex_lock(&gate->wait_lock) != 0) {
        return false;
    }

    unsigned previous = atomic_fetch_sub_explicit(&gate->waiting_writers, 1U, memory_order_release);

    if (previous == 1U) {
        (void) pthread_cond_broadcast(&gate->readers_ready);
    }

    return pthread_mutex_unlock(&gate->wait_lock) == 0;
}

bool geo_visibility_gate_initialize(GeoVisibilityGate *gate)
{
    if (!gate) {
        return false;
    }

    memset(gate, 0, sizeof(*gate));

    if (pthread_rwlock_init(&gate->lock, NULL) != 0) {
        return false;
    }
    if (pthread_mutex_init(&gate->wait_lock, NULL) != 0) {
        (void) pthread_rwlock_destroy(&gate->lock);
        return false;
    }
    if (pthread_cond_init(&gate->readers_ready, NULL) != 0) {
        (void) pthread_mutex_destroy(&gate->wait_lock);
        (void) pthread_rwlock_destroy(&gate->lock);
        return false;
    }

    atomic_init(&gate->waiting_writers, 0U);
    return true;
}

void geo_visibility_gate_destroy(GeoVisibilityGate *gate)
{
    if (!gate) {
        return;
    }

    (void) pthread_cond_destroy(&gate->readers_ready);
    (void) pthread_mutex_destroy(&gate->wait_lock);
    (void) pthread_rwlock_destroy(&gate->lock);
}

bool geo_visibility_gate_read_lock(GeoVisibilityGate *gate)
{
    if (!gate) {
        return false;
    }

    for (;;) {
        /* A reader that races the writer announcement may enter, but the number of such
         * readers is bounded by the threads already in this fast path. Subsequent calls wait. */
        if (!atomic_load_explicit(&gate->waiting_writers, memory_order_acquire)) {
            return pthread_rwlock_rdlock(&gate->lock) == 0;
        }

        if (pthread_mutex_lock(&gate->wait_lock) != 0) {
            return false;
        }

        int wait_status = 0;

        while (!wait_status && atomic_load_explicit(&gate->waiting_writers, memory_order_acquire)) {
            wait_status = pthread_cond_wait(&gate->readers_ready, &gate->wait_lock);
        }

        int unlock_status = pthread_mutex_unlock(&gate->wait_lock);

        if (wait_status || unlock_status) {
            return false;
        }
    }
}

bool geo_visibility_gate_read_unlock(GeoVisibilityGate *gate)
{
    return gate && pthread_rwlock_unlock(&gate->lock) == 0;
}

bool geo_visibility_gate_write_lock(GeoVisibilityGate *gate)
{
    if (!gate) {
        return false;
    }

    if (!visibility_gate_writer_arrive(gate)) {
        return false;
    }

    if (pthread_rwlock_wrlock(&gate->lock) == 0) {
        return true;
    }

    (void) visibility_gate_writer_depart(gate);
    return false;
}

bool geo_visibility_gate_write_unlock(GeoVisibilityGate *gate)
{
    if (!gate) {
        return false;
    }

    bool succeeded = pthread_rwlock_unlock(&gate->lock) == 0;

    return visibility_gate_writer_depart(gate) && succeeded;
}
