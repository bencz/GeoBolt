#ifndef GEO_DATABASE_INTERNAL_H
#define GEO_DATABASE_INTERNAL_H

#include "geobolt/geobolt.h"
#include "geo_db_format.h"
#include "geo_object_cache.h"
#include "geo_rocks_bridge.h"
#include "geo_secondary_index.h"
#include "geo_spatial_memtable.h"
#include "geo_visibility_gate.h"

#include <pthread.h>
#include <stdatomic.h>

typedef enum {
    GEO_DATABASE_MUTATIONS_COMPACT,
    GEO_DATABASE_MUTATIONS_OBJECT,
} GeoDatabaseMutationLayout;

#define GEO_DATABASE_INTERNAL_INSERT_GENERATED 5U

typedef struct {
    const void *mutations;
    const void *idempotency_key;
    size_t idempotency_key_size;
    uint32_t idempotency_retention_seconds;
    GeoDatabaseMutationLayout layout;
} GeoDatabaseMutationView;

typedef struct GeoDatabaseWriteRequest GeoDatabaseWriteRequest;

struct GeoDatabaseWriteRequest {
    GeoDatabaseWriteRequest *next;
    GeoDatabaseMutationView view;
    size_t mutation_count;
    uint64_t generated_object_id;
    bool replayed;
    GeoDatabaseStatus status;
    pthread_cond_t completion;
    bool completed;
};

struct GeoDatabase {
    char *directory;
    char *rocks_directory;
    char *manifest_path;
    char *segments_directory;
    GeoRocksDatabase *rocks;
    GeoSegmentSet *segments;
    GeoVisibilityGate visibility_gate;
    pthread_mutex_t writer_lock;
    pthread_mutex_t commit_queue_lock;
    pthread_cond_t commit_queue_ready;
    pthread_t commit_thread;
    pthread_mutex_t spatial_flush_lock;
    pthread_cond_t spatial_flush_ready;
    pthread_t spatial_flush_thread;
    GeoDatabaseWriteRequest *commit_head;
    GeoDatabaseWriteRequest *commit_tail;
    GeoDatabaseConfig config;
    GeoDbCatalog catalog;
    GeoSecondaryCatalog *secondary_indexes;
    GeoSpatialMemtable *spatial_memtable;
    GeoObjectCache *object_cache;
    uint64_t recovered_operations;
    uint64_t checkpoints;
    uint64_t maintenance_failures;
    bool visibility_gate_initialized;
    bool writer_lock_initialized;
    bool commit_queue_lock_initialized;
    bool commit_queue_ready_initialized;
    bool commit_thread_started;
    bool spatial_flush_lock_initialized;
    bool spatial_flush_ready_initialized;
    bool spatial_flush_thread_started;
    bool spatial_flush_stop;
    bool spatial_flush_busy;
    bool spatial_flush_failed;
    bool commit_stop;
    atomic_bool failed;
};

GeoDatabaseObjectMutation geo_database_mutation_at(const GeoDatabaseMutationView *view, size_t index);
GeoDatabaseStatus geo_database_validate_mutations(const GeoDatabaseMutationView *view, size_t mutation_count);
GeoDatabaseStatus geo_database_object_load(GeoDatabase *database,
                                           uint64_t object_id,
                                           GeoRocksBuffer *value,
                                           bool *exists);
GeoDatabaseStatus geo_database_commit_view(GeoDatabase *database,
                                           const GeoDatabaseMutationView *view,
                                           size_t mutation_count,
                                           uint64_t *generated_object_id,
                                           bool *replayed);
GeoDatabaseStatus geo_database_write_view(GeoDatabase *database,
                                          const GeoDatabaseMutationView *view,
                                          size_t mutation_count,
                                          uint64_t *generated_object_id,
                                          bool *replayed);
void *geo_database_commit_thread_main(void *argument);
bool geo_database_commit_queue_initialize(GeoDatabase *database);

#endif
