#ifndef GEO_VISIBILITY_GATE_H
#define GEO_VISIBILITY_GATE_H

#include <pthread.h>
#include <stdbool.h>
#include <stdatomic.h>

typedef struct {
    pthread_rwlock_t lock;
    pthread_mutex_t wait_lock;
    pthread_cond_t readers_ready;
    atomic_uint waiting_writers;
} GeoVisibilityGate;

bool geo_visibility_gate_initialize(GeoVisibilityGate *gate);
void geo_visibility_gate_destroy(GeoVisibilityGate *gate);

bool geo_visibility_gate_read_lock(GeoVisibilityGate *gate);
bool geo_visibility_gate_read_unlock(GeoVisibilityGate *gate);
bool geo_visibility_gate_write_lock(GeoVisibilityGate *gate);
bool geo_visibility_gate_write_unlock(GeoVisibilityGate *gate);

#endif
