#ifndef GEOBOLT_H
#define GEOBOLT_H

#include "geobolt/geo_index.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define GEO_DATABASE_MAX_IDEMPOTENCY_KEY_SIZE 128U
#define GEO_DATABASE_MAX_IDEMPOTENCY_RETENTION_SECONDS (7U * 24U * 60U * 60U)
#define GEO_DATABASE_MAX_INDEX_NAME_SIZE 63U
#define GEO_DATABASE_MAX_INDEX_POINTER_SIZE 511U
#define GEO_DATABASE_MAX_INDEXED_VALUE_SIZE 8192U
#define GEO_DATABASE_MAX_QUERY_PREDICATES 16U

typedef struct GeoDatabase GeoDatabase;
typedef struct GeoDatabaseQueryWorkspace GeoDatabaseQueryWorkspace;

typedef enum {
    GEO_DATABASE_UPSERT = 1,
    GEO_DATABASE_DELETE = 2,
    GEO_DATABASE_INSERT = 3,
    GEO_DATABASE_UPDATE = 4,
} GeoDatabaseOperation;

typedef enum {
    GEO_DATABASE_OK = 0,
    GEO_DATABASE_INVALID_ARGUMENT,
    GEO_DATABASE_OUT_OF_MEMORY,
    GEO_DATABASE_IO_ERROR,
    GEO_DATABASE_CORRUPTION,
    GEO_DATABASE_NOT_FOUND,
    GEO_DATABASE_ALREADY_EXISTS,
    GEO_DATABASE_IDEMPOTENCY_CONFLICT,
    GEO_DATABASE_FAILED_STATE,
} GeoDatabaseStatus;

typedef struct {
    uint64_t object_id;
    uint64_t morton_code;
    uint32_t operation;
    uint32_t reserved;
} GeoDatabaseMutation;

typedef struct {
    uint64_t object_id;
    uint64_t morton_code;
    const void *document;
    size_t document_size;
    uint32_t operation;
    uint32_t reserved;
} GeoDatabaseObjectMutation;

typedef struct {
    uint64_t object_id;
    uint64_t sequence;
    uint64_t morton_code;
    void *document;
    size_t document_size;
} GeoDatabaseObject;

typedef enum {
    GEO_DATABASE_INDEX_BOOL = 1,
    GEO_DATABASE_INDEX_INT64,
    GEO_DATABASE_INDEX_UINT64,
    GEO_DATABASE_INDEX_DOUBLE,
    GEO_DATABASE_INDEX_DATETIME,
    GEO_DATABASE_INDEX_STRING,
    GEO_DATABASE_INDEX_BYTES,
} GeoDatabaseIndexType;

typedef enum {
    GEO_DATABASE_INDEX_EQUAL = 1,
    GEO_DATABASE_INDEX_LESS,
    GEO_DATABASE_INDEX_LESS_EQUAL,
    GEO_DATABASE_INDEX_GREATER,
    GEO_DATABASE_INDEX_GREATER_EQUAL,
    GEO_DATABASE_INDEX_BETWEEN,
} GeoDatabaseIndexOperator;

typedef struct {
    GeoDatabaseIndexType type;
    union {
        bool boolean;
        int64_t signed_integer;
        uint64_t unsigned_integer;
        double floating_point;
        int64_t datetime;
        struct {
            const void *data;
            size_t size;
        } bytes;
    } as;
} GeoDatabaseIndexValue;

typedef struct {
    const char *index_name;
    GeoDatabaseIndexOperator operation;
    GeoDatabaseIndexValue lower;
    GeoDatabaseIndexValue upper;
} GeoDatabaseIndexPredicate;

typedef struct {
    uint64_t index_id;
    uint64_t entry_count;
    GeoDatabaseIndexType type;
    char name[GEO_DATABASE_MAX_INDEX_NAME_SIZE + 1U];
    char json_pointer[GEO_DATABASE_MAX_INDEX_POINTER_SIZE + 1U];
} GeoDatabaseIndexInfo;

typedef struct {
    uint64_t scanned_entries;
    uint64_t matched_entries;
} GeoDatabaseIndexQueryStats;

typedef enum {
    GEO_DATABASE_QUERY_PLAN_SPATIAL = 1,
    GEO_DATABASE_QUERY_PLAN_SECONDARY = 2,
} GeoDatabaseQueryPlan;

typedef struct {
    const GeoDatabaseIndexPredicate *predicates;
    size_t predicate_count;
    double latitude;
    double longitude;
    double radius_km;
} GeoDatabaseRadiusQuery;

typedef struct {
    GeoSearchStats spatial;
    uint64_t secondary_entries_scanned;
    uint64_t metadata_candidates;
    uint64_t estimated_spatial_candidates;
    // Canonical RocksDB keys requested after process-local object-cache hits are removed.
    uint64_t object_lookups;
    GeoDatabaseQueryPlan plan;
    uint32_t predicate_count;
} GeoDatabaseQueryStats;

#if defined(__cplusplus)
static_assert(sizeof(GeoDatabaseMutation) == 24U, "GeoDatabaseMutation must remain compact and naturally aligned");
static_assert(offsetof(GeoDatabaseMutation, morton_code) == 8U, "GeoDatabaseMutation Morton offset is part of the batch ABI");
static_assert(offsetof(GeoDatabaseMutation, operation) == 16U, "GeoDatabaseMutation operation offset is part of the batch ABI");
#else
_Static_assert(sizeof(GeoDatabaseMutation) == 24U, "GeoDatabaseMutation must remain compact and naturally aligned");
_Static_assert(offsetof(GeoDatabaseMutation, morton_code) == 8U, "GeoDatabaseMutation Morton offset is part of the batch ABI");
_Static_assert(offsetof(GeoDatabaseMutation, operation) == 16U, "GeoDatabaseMutation operation offset is part of the batch ABI");
#endif

typedef struct {
    size_t block_cache_bytes;
    // Process-local cache for validated canonical objects. Zero disables it.
    size_t object_cache_bytes;
    size_t write_buffer_bytes;
    size_t max_active_segments;
    size_t group_commit_max_operations;
    // Hard bound for each of the active and frozen spatial generations. Must cover one maximum group commit.
    size_t spatial_memtable_max_operations;
    uint32_t group_commit_delay_us;
    int background_jobs;
    bool create_if_missing;
} GeoDatabaseConfig;

typedef struct {
    uint64_t committed_operations;
    uint64_t recovered_operations;
    uint64_t checkpoints;
    uint64_t maintenance_failures;
    uint64_t next_sequence;
    uint64_t physical_records;
    size_t active_segments;
    size_t active_spatial_operations;
    size_t frozen_spatial_operations;
} GeoDatabaseStats;

// Uses production defaults for the shared RocksDB block cache, memtables, background jobs, and derived spatial segments.
GeoDatabaseConfig geo_database_default_config(void);

// Opens one single-server database rooted at directory. All paths below the root are owned and maintained by GeoBolt.
GeoDatabase *geo_database_open(const char *directory,
                               const GeoDatabaseConfig *config,
                               GeoDatabaseStatus *status);
void geo_database_close(GeoDatabase *database);

// Every batch atomically commits object state, metadata, spatial deltas, and catalog watermarks in the RocksDB WAL.
// IDs must be nonzero and unique inside one batch. GeoDoc payloads are validated before entering the write batch.
// Mutation arrays and referenced GeoDoc bytes must remain immutable until the synchronous call returns.
GeoDatabaseStatus geo_database_write_objects(GeoDatabase *database,
                                             const GeoDatabaseObjectMutation *mutations,
                                             size_t mutation_count);

// Atomically stores the successful batch and its opaque retry key. Matching live-key retries do not consume sequences; a different
// batch under the same key returns GEO_DATABASE_IDEMPOTENCY_CONFLICT. Expired keys are reclaimed during RocksDB compaction.
GeoDatabaseStatus geo_database_write_objects_idempotent(GeoDatabase *database,
                                                        const GeoDatabaseObjectMutation *mutations,
                                                        size_t mutation_count,
                                                        const void *idempotency_key,
                                                        size_t idempotency_key_size,
                                                        uint32_t retention_seconds,
                                                        bool *replayed);

// Compact ingestion path for objects without metadata.
// Coordinates are supplied as packed 64-bit Morton codes, avoiding double conversion in pre-encoded ingestion paths.
GeoDatabaseStatus geo_database_write(GeoDatabase *database,
                                     const GeoDatabaseMutation *mutations,
                                     size_t mutation_count);

// Convenience path for applications that have latitude/longitude rather than pre-encoded Morton coordinates.
GeoDatabaseStatus geo_database_upsert(GeoDatabase *database,
                                      uint64_t object_id,
                                      double latitude,
                                      double longitude);
GeoDatabaseStatus geo_database_upsert_document(GeoDatabase *database,
                                               uint64_t object_id,
                                               double latitude,
                                               double longitude,
                                               const void *document,
                                               size_t document_size);
GeoDatabaseStatus geo_database_insert_document(GeoDatabase *database,
                                               uint64_t object_id,
                                               double latitude,
                                               double longitude,
                                               const void *document,
                                               size_t document_size);
GeoDatabaseStatus geo_database_insert_generated_document(GeoDatabase *database,
                                                         double latitude,
                                                         double longitude,
                                                         const void *document,
                                                         size_t document_size,
                                                         uint64_t *object_id);
GeoDatabaseStatus geo_database_insert_generated_object(GeoDatabase *database,
                                                       uint64_t morton_code,
                                                       const void *document,
                                                       size_t document_size,
                                                       uint64_t *object_id);
GeoDatabaseStatus geo_database_update_document(GeoDatabase *database,
                                               uint64_t object_id,
                                               double latitude,
                                               double longitude,
                                               const void *document,
                                               size_t document_size);
GeoDatabaseStatus geo_database_remove(GeoDatabase *database, uint64_t object_id);

GeoDatabaseStatus geo_database_get(const GeoDatabase *database, uint64_t object_id, GeoDatabaseObject *object);
void geo_database_object_release(GeoDatabaseObject *object);

// Secondary indexes are strict and sparse: a missing or null path emits no key; a present value with the wrong type rejects the write.
// DATETIME indexes consume a GeoDoc INT64 value and preserve signed chronological ordering.
// Failed construction rolls back its private index state when storage permits; allocated index IDs are never reused.
// Recovery discards an interrupted type-incompatible build, but does not ignore canonical corruption or I/O failures.
GeoDatabaseStatus geo_database_create_index(GeoDatabase *database,
                                            const char *name,
                                            const char *json_pointer,
                                            GeoDatabaseIndexType type);
GeoDatabaseStatus geo_database_drop_index(GeoDatabase *database, const char *name);
GeoDatabaseStatus geo_database_list_indexes(const GeoDatabase *database,
                                            GeoDatabaseIndexInfo *indexes,
                                            size_t capacity,
                                            size_t *count);
GeoDatabaseStatus geo_database_query_index_reuse(const GeoDatabase *database,
                                                 const GeoDatabaseIndexPredicate *predicate,
                                                 GeoIdResult *result,
                                                 GeoDatabaseIndexQueryStats *stats);

// A workspace is reusable across calls but must not be shared by concurrent queries. Result and intermediate capacities are retained.
GeoDatabaseQueryWorkspace *geo_database_query_workspace_create(size_t predicate_capacity);
void geo_database_query_workspace_destroy(GeoDatabaseQueryWorkspace *workspace);
GeoDatabaseStatus geo_database_query_radius_reuse(const GeoDatabase *database,
                                                  const GeoDatabaseRadiusQuery *query,
                                                  GeoDatabaseQueryWorkspace *workspace,
                                                  GeoSearchResult *result,
                                                  GeoDatabaseQueryStats *stats);

bool geo_database_search_radius_reuse(const GeoDatabase *database,
                                      double latitude,
                                      double longitude,
                                      double radius_km,
                                      GeoSearchResult *result,
                                      GeoSearchStats *stats);

bool geo_database_search_radius_count(const GeoDatabase *database,
                                      double latitude,
                                      double longitude,
                                      double radius_km,
                                      size_t *count,
                                      GeoSearchStats *stats);

// Flushes RocksDB and checkpoints the derived spatial mutation map. Normal compaction remains autonomous.
GeoDatabaseStatus geo_database_checkpoint(GeoDatabase *database);
bool geo_database_get_stats(const GeoDatabase *database, GeoDatabaseStats *stats);

#ifdef __cplusplus
}
#endif

#endif
