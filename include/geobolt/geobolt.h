#ifndef GEOBOLT_H
#define GEOBOLT_H

#include "geobolt/geo_index.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GeoDatabase GeoDatabase;

typedef enum {
    GEO_DATABASE_UPSERT = 1,
    GEO_DATABASE_DELETE = 2,
} GeoDatabaseOperation;

typedef enum {
    GEO_DATABASE_OK = 0,
    GEO_DATABASE_INVALID_ARGUMENT,
    GEO_DATABASE_OUT_OF_MEMORY,
    GEO_DATABASE_IO_ERROR,
    GEO_DATABASE_CORRUPTION,
    GEO_DATABASE_FAILED_STATE,
} GeoDatabaseStatus;

typedef struct {
    uint64_t object_id;
    uint64_t morton_code;
    uint32_t operation;
    uint32_t reserved;
} GeoDatabaseMutation;

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
    size_t wal_segment_size;
    size_t checkpoint_interval_operations;
    size_t max_active_segments;
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
} GeoDatabaseStats;

// Uses production defaults: 64 MiB WAL segments, checkpointing after 262,144 operations, and at most 16 active data segments.
GeoDatabaseConfig geo_database_default_config(void);

// Opens one single-server database rooted at directory. All paths below the root are owned and maintained by GeoBolt.
GeoDatabase *geo_database_open(const char *directory,
                               const GeoDatabaseConfig *config,
                               GeoDatabaseStatus *status);
void geo_database_close(GeoDatabase *database);

// Every batch is durably ordered in the WAL before becoming visible. IDs must be nonzero and unique inside one batch.
// Coordinates are supplied as packed 64-bit Morton codes, avoiding double conversion in pre-encoded ingestion paths.
GeoDatabaseStatus geo_database_write(GeoDatabase *database,
                                     const GeoDatabaseMutation *mutations,
                                     size_t mutation_count);

// Convenience path for applications that have latitude/longitude rather than pre-encoded Morton coordinates.
GeoDatabaseStatus geo_database_upsert(GeoDatabase *database,
                                      uint64_t object_id,
                                      double latitude,
                                      double longitude);
GeoDatabaseStatus geo_database_remove(GeoDatabase *database, uint64_t object_id);

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

// Forces a durable mutation checkpoint and advances WAL reclamation. Normal operation invokes this autonomously by policy.
GeoDatabaseStatus geo_database_checkpoint(GeoDatabase *database);
bool geo_database_get_stats(const GeoDatabase *database, GeoDatabaseStats *stats);

#ifdef __cplusplus
}
#endif

#endif
