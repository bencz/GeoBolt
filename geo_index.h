#ifndef GEO_INDEX_H
#define GEO_INDEX_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// =========================================================
// Constants
// =========================================================

#define GEO_MIN_LAT -90.0
#define GEO_MAX_LAT  90.0
#define GEO_MIN_LNG -180.0
#define GEO_MAX_LNG  180.0
#define GEO_EARTH_RADIUS_KM 6371.0088

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

// Precision: 32 bits per coordinate, about 4.2e-8 degrees of latitude per step.
#define GEO_COORD_BITS 32
#define GEO_COORD_MAX  4294967295ULL

// =========================================================
// Structures
// =========================================================

typedef struct {
    double lat;
    double lng;
} GeoPoint;

typedef struct GeoRecord {
    uint64_t id;
    uint64_t z;  // Two normalized 32-bit coordinates interleaved into one Morton code.
} GeoRecord;

typedef struct GeoIndex GeoIndex;
typedef struct GeoDensityIndex GeoDensityIndex;

#if defined(__cplusplus)
static_assert(sizeof(GeoRecord) == 16, "GeoRecord must remain a compact 16-byte persisted record");
static_assert(alignof(GeoRecord) >= alignof(uint64_t), "GeoRecord must preserve naturally aligned 64-bit fields");
static_assert(offsetof(GeoRecord, id) == 0, "GeoRecord.id must remain the first field");
static_assert(offsetof(GeoRecord, z) == 8, "GeoRecord.z must remain naturally aligned after id");
#else
_Static_assert(sizeof(GeoRecord) == 16, "GeoRecord must remain a compact 16-byte persisted record");
_Static_assert(_Alignof(GeoRecord) >= _Alignof(uint64_t), "GeoRecord must preserve naturally aligned 64-bit fields");
_Static_assert(offsetof(GeoRecord, id) == 0, "GeoRecord.id must remain the first field");
_Static_assert(offsetof(GeoRecord, z) == 8, "GeoRecord.z must remain naturally aligned after id");
#endif

typedef struct {
    uint64_t min;
    uint64_t max;
} ZRange;

typedef struct {
    GeoRecord *results;
    size_t count;
    size_t capacity;
} GeoSearchResult;

typedef struct {
    uint64_t *ids;
    size_t count;
    size_t capacity;
} GeoIdResult;

typedef struct {
    uint64_t *ids;
    size_t *offsets;
    size_t query_count;
    size_t id_count;
    size_t id_capacity;
    size_t offset_capacity;
} GeoBatchIdResult;

typedef struct {
    uint64_t records_scanned;
    uint64_t records_matched;
    uint64_t ranges_checked;
    double search_time_ms;
} GeoSearchStats;

typedef struct GeoStreamBuilder GeoStreamBuilder;
typedef struct GeoSegmentSet GeoSegmentSet;
typedef struct GeoQueryExecutor GeoQueryExecutor;
typedef struct GeoKnnWorkspace GeoKnnWorkspace;

typedef enum {
    GEO_INDEX_ADVICE_NORMAL,
    GEO_INDEX_ADVICE_RANDOM,
    GEO_INDEX_ADVICE_SEQUENTIAL,
    GEO_INDEX_ADVICE_WILL_NEED,
    GEO_INDEX_ADVICE_HUGE_PAGE,
} GeoIndexMemoryAdvice;

typedef struct {
    bool verify_checksums;
    bool validate_sorted_order;
    GeoIndexMemoryAdvice memory_advice;
} GeoIndexOpenOptions;

typedef struct {
    size_t thread_count;
    size_t scheduling_chunk;
    const unsigned *worker_cpus;
    size_t worker_cpu_count;
    bool pin_workers;
} GeoQueryExecutorConfig;

typedef struct {
    uint64_t records_written;
    size_t runs_created;
    size_t intermediate_merges;
    size_t peak_open_runs;
    double chunk_build_time_ms;
    double intermediate_merge_time_ms;
    double merge_time_ms;
} GeoStreamBuildStats;

typedef struct {
    uint64_t records_written;
    size_t input_segments;
    size_t worker_count;
    size_t partition_count;
    double merge_time_ms;
} GeoSegmentCompactionStats;

typedef struct {
    size_t max_active_segments;
    size_t minimum_mutations;
    size_t maximum_mutation_rewrite_passes;
    uint32_t mutation_ratio_numerator;
    uint32_t mutation_ratio_denominator;
} GeoSegmentCompactionPolicy;

// =========================================================
// Morton Code Functions (Z-order curve)
// =========================================================

uint64_t geo_spread_bits(uint32_t v);
uint32_t geo_compact_bits(uint64_t v);

// =========================================================
// Coordinate Normalization
// =========================================================

uint32_t geo_normalize_lat(double lat);
uint32_t geo_normalize_lng(double lng);
double geo_denormalize_lat(uint32_t v);
double geo_denormalize_lng(uint32_t v);

// =========================================================
// Encode / Decode
// =========================================================

uint64_t geo_encode(double lat, double lng);
GeoPoint geo_decode(uint64_t z);

// =========================================================
// Distance Calculations
// =========================================================

double geo_to_radians(double degrees);
double geo_to_degrees(double radians);
double geo_haversine_km(double lat1, double lng1, double lat2, double lng2);
double geo_haversine_m(double lat1, double lng1, double lat2, double lng2);

// Fast approximate distance (Equirectangular approximation)
double geo_fast_distance_km(double lat1, double lng1, double lat2, double lng2);

// Bounding box calculations
void geo_bounding_box(double lat, double lng, double radius_km,
                      double *min_lat, double *max_lat,
                      double *min_lng, double *max_lng);

// =========================================================
// Index Management
// =========================================================

GeoIndex* geo_index_create(size_t initial_capacity);
void geo_index_destroy(GeoIndex *index);
bool geo_index_reserve(GeoIndex *index, size_t capacity);
bool geo_index_add(GeoIndex *index, uint64_t id, double lat, double lng);
bool geo_index_add_batch(GeoIndex *index, const uint64_t *ids,
                         const double *lats, const double *lngs, size_t count);
bool geo_index_add_records(GeoIndex *index, const GeoRecord *records, size_t count);
bool geo_index_build(GeoIndex *index);  // Sort and prepare for queries.
bool geo_index_build_parallel(GeoIndex *index, size_t thread_count);
void geo_index_clear(GeoIndex *index);
bool geo_index_is_read_only(const GeoIndex *index);

// Persistent read-only indexes. The mapped representation is directly
// queryable and does not copy the record array into heap memory.
bool geo_index_save(const GeoIndex *index, const char *path);
GeoIndex* geo_index_open_mmap(const char *path);
GeoIndex *geo_index_open_mmap_with_options(const char *path, const GeoIndexOpenOptions *options);
bool geo_index_advise(GeoIndex *index, GeoIndexMemoryAdvice advice);

// External-memory construction. At most chunk_capacity records are held in the
// mutable chunk, while sorted temporary runs are merged into output_path.
GeoStreamBuilder *geo_stream_builder_create(const char *output_path,
                                            const char *temporary_directory,
                                            size_t chunk_capacity);

GeoStreamBuilder *geo_stream_builder_create_parallel(const char *output_path,
                                                     const char *temporary_directory,
                                                     size_t chunk_capacity,
                                                     size_t sort_threads);

void geo_stream_builder_destroy(GeoStreamBuilder *builder);

bool geo_stream_builder_add_batch(GeoStreamBuilder *builder,
                                  const uint64_t *ids,
                                  const double *latitudes,
                                  const double *longitudes,
                                  size_t count);

bool geo_stream_builder_add_records(GeoStreamBuilder *builder, const GeoRecord *records, size_t count);

bool geo_stream_builder_finish(GeoStreamBuilder *builder, GeoStreamBuildStats *stats);

// Immutable segment collections for continuously published indexes.
GeoSegmentSet *geo_segment_set_create(const char *manifest_path);
GeoSegmentSet *geo_segment_set_open(const char *manifest_path);
void geo_segment_set_destroy(GeoSegmentSet *set);
bool geo_segment_set_add_file(GeoSegmentSet *set, const char *segment_path);

// Publishes one new location per ID and makes every older physical copy of those IDs invisible at the same linearization point.
// Duplicate IDs inside the input segment are rejected because a single generation cannot order them.
bool geo_segment_set_upsert_file(GeoSegmentSet *set, const char *segment_path);

// Linearizable logical deletion. Removed records remain in immutable segment files until compaction,
// but are excluded from radius, bounding-box, count-only, and kNN queries immediately after publication.
bool geo_segment_set_remove(GeoSegmentSet *set, uint64_t id);
bool geo_segment_set_remove_ids(GeoSegmentSet *set, const uint64_t *ids, size_t count);
bool geo_segment_set_checkpoint_mutations(GeoSegmentSet *set);
size_t geo_segment_set_count(const GeoSegmentSet *set);
uint64_t geo_segment_set_record_count(const GeoSegmentSet *set);
bool geo_segment_set_compact(GeoSegmentSet *set,
                             const char *output_path,
                             GeoSegmentCompactionStats *stats);

// Executes the same durable compaction with an explicit worker limit. The calling thread participates and stats reports actual usage.
bool geo_segment_set_compact_with_workers(GeoSegmentSet *set,
                                          const char *output_path,
                                          size_t worker_count,
                                          GeoSegmentCompactionStats *stats);

// Configures both size-tiered merging and mutation-pressure reclamation. Every field must be nonzero and the mutation ratio must be in
// the inclusive range (0, 1]. Reconfiguration is rejected while a background worker is active.
bool geo_segment_set_configure_background_compaction(GeoSegmentSet *set,
                                                     const char *output_directory,
                                                     const GeoSegmentCompactionPolicy *policy);

// Enables the default policy: at most max_active_segments, at least 4096 mutations, a 1/8 mutation ratio, and two bounded rewrites.
bool geo_segment_set_enable_background_compaction(GeoSegmentSet *set,
                                                  const char *output_directory,
                                                  size_t max_active_segments);

// Waits for the currently scheduled compaction, if any, and joins its worker thread.
bool geo_segment_set_wait_for_background_compaction(GeoSegmentSet *set,
                                                    GeoSegmentCompactionStats *stats);

// Returns cumulative background work since create/open. Values are monitoring counters and do not reset when a worker is joined.
bool geo_segment_set_background_compaction_totals(const GeoSegmentSet *set,
                                                  GeoSegmentCompactionStats *stats,
                                                  uint64_t *completed_runs,
                                                  uint64_t *failed_runs);

bool geo_segment_set_background_compaction_active(const GeoSegmentSet *set);

// Reports both explicit and background merge activity. The result is an instantaneous observation,
// intended for monitoring and coordination rather than as a synchronization primitive.
bool geo_segment_set_compaction_active(const GeoSegmentSet *set);

// Copies one active segment path while holding a consistent segment-set snapshot.
bool geo_segment_set_copy_path(const GeoSegmentSet *set,
                               size_t segment_index,
                               char *buffer,
                               size_t buffer_size);

GeoSearchResult *geo_segment_set_search_radius(const GeoSegmentSet *set,
                                               double lat,
                                               double lng,
                                               double radius_km,
                                               GeoSearchStats *stats);

bool geo_segment_set_search_radius_reuse(const GeoSegmentSet *set,
                                         double lat,
                                         double lng,
                                         double radius_km,
                                         GeoSearchResult *result,
                                         GeoSearchStats *stats);

bool geo_segment_set_search_radius_count(const GeoSegmentSet *set,
                                         double lat,
                                         double lng,
                                         double radius_km,
                                         size_t *count,
                                         GeoSearchStats *stats);

bool geo_segment_set_search_radius_ids_reuse(const GeoSegmentSet *set,
                                             double lat,
                                             double lng,
                                             double radius_km,
                                             GeoIdResult *result,
                                             GeoSearchStats *stats);

bool geo_segment_set_search_radius_ids_into(const GeoSegmentSet *set,
                                            double lat,
                                            double lng,
                                            double radius_km,
                                            uint64_t *ids,
                                            size_t capacity,
                                            size_t *count,
                                            GeoSearchStats *stats);

GeoSearchResult *geo_segment_set_search_bbox(const GeoSegmentSet *set,
                                             double min_lat,
                                             double max_lat,
                                             double min_lng,
                                             double max_lng,
                                             GeoSearchStats *stats);

bool geo_segment_set_search_bbox_reuse(const GeoSegmentSet *set,
                                       double min_lat,
                                       double max_lat,
                                       double min_lng,
                                       double max_lng,
                                       GeoSearchResult *result,
                                       GeoSearchStats *stats);

bool geo_segment_set_search_bbox_count(const GeoSegmentSet *set,
                                       double min_lat,
                                       double max_lat,
                                       double min_lng,
                                       double max_lng,
                                       size_t *count,
                                       GeoSearchStats *stats);

GeoSearchResult *geo_segment_set_search_knn(const GeoSegmentSet *set,
                                            double lat,
                                            double lng,
                                            size_t k,
                                            double max_radius_km,
                                            GeoSearchStats *stats);

// =========================================================
// Search Functions
// =========================================================

// Binary search utilities
size_t geo_lower_bound(const GeoRecord *records, size_t n, uint64_t key);
size_t geo_upper_bound(const GeoRecord *records, size_t n, uint64_t key);

// Build Z-ranges for radius search
int geo_build_ranges(double lat, double lng, double radius_km, 
                     ZRange *out, int max_ranges);

// Optimized range building with adaptive subdivision
int geo_build_ranges_adaptive(double lat, double lng, double radius_km,
                              ZRange *out, int max_ranges, int precision_level);

// Search within radius
GeoSearchResult* geo_search_radius(const GeoIndex *index, 
                                   double lat, double lng, double radius_km,
                                   GeoSearchStats *stats);

// Allocation-reusing variant. Existing result storage is preserved and reused.
bool geo_search_radius_reuse(const GeoIndex *index,
                             double lat, double lng, double radius_km,
                             GeoSearchResult *result, GeoSearchStats *stats);

// Count-only query. It runs the same exact predicate without allocating or copying records.
bool geo_search_radius_count(const GeoIndex *index,
                             double lat, double lng, double radius_km,
                             size_t *count, GeoSearchStats *stats);

// ID-only materialization halves result bandwidth when callers do not need the Morton code.
bool geo_search_radius_ids_reuse(const GeoIndex *index,
                                 double lat,
                                 double lng,
                                 double radius_km,
                                 GeoIdResult *result,
                                 GeoSearchStats *stats);

// Fills caller-owned storage without allocation. capacity must be at least the count returned
// by geo_search_radius_count() for the same immutable index and query.
bool geo_search_radius_ids_into(const GeoIndex *index,
                                double lat,
                                double lng,
                                double radius_km,
                                uint64_t *ids,
                                size_t capacity,
                                size_t *count,
                                GeoSearchStats *stats);

// Search within bounding box
GeoSearchResult* geo_search_bbox(const GeoIndex *index,
                                 double min_lat, double max_lat,
                                 double min_lng, double max_lng,
                                 GeoSearchStats *stats);

bool geo_search_bbox_reuse(const GeoIndex *index,
                           double min_lat, double max_lat,
                           double min_lng, double max_lng,
                           GeoSearchResult *result, GeoSearchStats *stats);

bool geo_search_bbox_count(const GeoIndex *index,
                           double min_lat, double max_lat,
                           double min_lng, double max_lng,
                           size_t *count, GeoSearchStats *stats);

// K-nearest neighbors search
GeoSearchResult* geo_search_knn(const GeoIndex *index,
                                double lat, double lng, size_t k,
                                double max_radius_km,
                                GeoSearchStats *stats);

GeoKnnWorkspace *geo_knn_workspace_create(size_t neighbor_capacity, size_t cell_capacity);
void geo_knn_workspace_destroy(GeoKnnWorkspace *workspace);
bool geo_knn_workspace_reserve(GeoKnnWorkspace *workspace, size_t neighbor_capacity, size_t cell_capacity);

bool geo_search_knn_reuse(const GeoIndex *index,
                          double lat,
                          double lng,
                          size_t k,
                          double max_radius_km,
                          GeoSearchResult *result,
                          GeoKnnWorkspace *workspace,
                          GeoSearchStats *stats);

// =========================================================
// Result Management
// =========================================================

GeoSearchResult* geo_result_create(size_t initial_capacity);
void geo_result_destroy(GeoSearchResult *result);
void geo_result_clear(GeoSearchResult *result);
bool geo_result_reserve(GeoSearchResult *result, size_t capacity);
bool geo_result_add(GeoSearchResult *result, const GeoRecord *record);
void geo_result_sort_by_distance(GeoSearchResult *result, double lat, double lng);

GeoIdResult *geo_id_result_create(size_t initial_capacity);
void geo_id_result_destroy(GeoIdResult *result);
void geo_id_result_clear(GeoIdResult *result);
bool geo_id_result_reserve(GeoIdResult *result, size_t capacity);

GeoBatchIdResult *geo_batch_id_result_create(size_t initial_query_capacity, size_t initial_id_capacity);
void geo_batch_id_result_destroy(GeoBatchIdResult *result);
void geo_batch_id_result_clear(GeoBatchIdResult *result);

// Persistent query workers consume SoA inputs. The same immutable index may be supplied once,
// or one complete replica per locality domain may be supplied for core-local/NUMA reads.
GeoQueryExecutor *geo_query_executor_create(const GeoIndex *index, const GeoQueryExecutorConfig *config);
GeoQueryExecutor *geo_query_executor_create_replicated(const GeoIndex *const *replicas,
                                                       size_t replica_count,
                                                       const GeoQueryExecutorConfig *config);

// Shards are disjoint partitions of one logical index. Each shard owns a local worker queue;
// query results are concatenated in shard order into the same two-pass contiguous output.
GeoQueryExecutor *geo_query_executor_create_sharded(const GeoIndex *const *shards,
                                                    size_t shard_count,
                                                    const GeoQueryExecutorConfig *config);

// Segment executors retain a consistent QSBR snapshot across both materialization passes.
// The segment set must outlive the executor. Concurrent mutations publish after an active submission releases its snapshot.
GeoQueryExecutor *geo_query_executor_create_segment_set(const GeoSegmentSet *set,
                                                        const GeoQueryExecutorConfig *config);

void geo_query_executor_destroy(GeoQueryExecutor *executor);

// Executes count and fill passes into one contiguous ID array. offsets has query_count + 1
// elements, so the IDs for query i are ids[offsets[i] ... offsets[i + 1]).
bool geo_query_executor_search_radius_ids(GeoQueryExecutor *executor,
                                          const double *latitudes,
                                          const double *longitudes,
                                          const double *radii_km,
                                          size_t query_count,
                                          GeoBatchIdResult *result);

// =========================================================
// Comparison Functions
// =========================================================

int geo_compare_records_by_z(const void *a, const void *b);
int geo_compare_records_by_id(const void *a, const void *b);

// =========================================================
// Utility Functions
// =========================================================

// Get current time in milliseconds (for benchmarking)
double geo_get_time_ms(void);

// Validate coordinates
bool geo_is_valid_lat(double lat);
bool geo_is_valid_lng(double lng);
bool geo_is_valid_point(double lat, double lng);

// Clamp coordinates to valid range
double geo_clamp_lat(double lat);
double geo_clamp_lng(double lng);

// Wrap longitude to [-180, 180]
double geo_wrap_lng(double lng);

// =========================================================
// SIMD Optimizations (optional)
// =========================================================

#ifndef GEO_INDEX_NO_SIMD
#include "geo_index_simd.h"
#endif

#ifdef __cplusplus
}
#endif

#endif // GEO_INDEX_H
