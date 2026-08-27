#ifndef GEO_INDEX_PRIVATE_H
#define GEO_INDEX_PRIVATE_H

#include "geobolt/geo_index.h"
#include "geo_index_persistence.h"

#include <stdio.h>

#define GEO_QUERY_MAX_RANGES 64

struct GeoIndex {
    GeoRecord *records;
    size_t *prefix_offsets;
    size_t count;
    union {
        size_t capacity;
        size_t mapping_size;
    };
    uint64_t content_checksum;
    uint8_t prefix_bits;
    uint8_t maximum_density_refinement;
    bool sorted;
    bool read_only;
    GeoDensityIndex *density;
    void *mapping;
};

#if SIZE_MAX == UINT64_MAX
_Static_assert(sizeof(GeoIndex) == 64, "GeoIndex should fit in one 64-byte cache line on supported 64-bit targets");
_Static_assert(offsetof(GeoIndex, records) == 0, "GeoIndex.records is the primary hot field");
_Static_assert(offsetof(GeoIndex, count) == 16, "GeoIndex.count must remain in the hot field group");
_Static_assert(offsetof(GeoIndex, content_checksum) == 32, "GeoIndex persisted identity must remain naturally aligned");
_Static_assert(offsetof(GeoIndex, prefix_bits) == 40, "GeoIndex query metadata must remain in the hot field group");
#endif

typedef enum {
    GEO_QUERY_STRATEGY_FILTER,
    GEO_QUERY_STRATEGY_COPY,
} GeoQueryStrategy;

typedef bool (*GeoRecordFilter)(const GeoRecord *record, void *context);
typedef size_t (*GeoRecordBitFilter)(const GeoRecord *records,
                                     size_t count,
                                     uint64_t *candidate_bits,
                                     void *context);
typedef struct GeoDensityWriter GeoDensityWriter;

typedef struct {
    const GeoIndex *index;
    GeoRecordFilter filter;
    void *filter_context;
} GeoKnnSource;

typedef struct {
    ZRange ranges[GEO_QUERY_MAX_RANGES];
    size_t range_begins[GEO_QUERY_MAX_RANGES];
    size_t range_ends[GEO_QUERY_MAX_RANGES];
    GeoSimdRadiusQuery radius_query;
    const GeoIndex *bounds_index;
    double latitude;
    double longitude;
    double radius_km;
    uint64_t candidate_records;
    uint64_t estimated_output_bytes;
    uint64_t estimated_cost;
    int range_count;
    uint8_t cover_depth;
    GeoQueryStrategy strategy;
    bool bounds_cached;
} GeoRadiusQueryPlan;

typedef struct {
    ZRange ranges[GEO_QUERY_MAX_RANGES];
    uint8_t needs_filter[GEO_QUERY_MAX_RANGES];
    double min_latitude;
    double max_latitude;
    double min_longitude;
    double max_longitude;
    int range_count;
} GeoBboxQueryPlan;

// Sorts a transient mutable index without constructing query-only metadata.
// This is intended for external-memory runs that are written and cleared immediately.
typedef struct GeoParallelSorter GeoParallelSorter;

void geo_index_sort_transient(GeoIndex *index);
bool geo_index_sort_transient_parallel(GeoIndex *index, size_t thread_count);
GeoParallelSorter *geo_parallel_sorter_create(size_t record_capacity, size_t thread_count);
void geo_parallel_sorter_destroy(GeoParallelSorter *sorter);
bool geo_index_sort_transient_with_sorter(GeoIndex *index, GeoParallelSorter *sorter);
bool geo_index_finalize_sorted(GeoIndex *index);

// Appends a batch whose coordinates were already validated by the caller.
bool geo_index_add_batch_unchecked(GeoIndex *index,
                                   const uint64_t *ids,
                                   const double *latitudes,
                                   const double *longitudes,
                                   size_t count);

bool geo_index_search_radius_append(const GeoIndex *index,
                                    double lat,
                                    double lng,
                                    double radius_km,
                                    GeoSearchResult *result,
                                    GeoSearchStats *stats);

bool geo_index_prepare_radius_query(double latitude,
                                    double longitude,
                                    double radius_km,
                                    uint64_t record_count,
                                    unsigned density_refinement,
                                    GeoRadiusQueryPlan *plan);

bool geo_index_prepare_radius_query_for_index(const GeoIndex *index,
                                              double latitude,
                                              double longitude,
                                              double radius_km,
                                              size_t output_record_size,
                                              GeoRadiusQueryPlan *plan);

uint64_t geo_index_radius_plan_candidate_count(const GeoIndex *index, const GeoRadiusQueryPlan *plan);

unsigned geo_index_density_refinement(const GeoIndex *index, double latitude, double longitude);

bool geo_density_index_build(GeoIndex *index);
void geo_density_index_destroy(GeoIndex *index);
bool geo_density_index_attach(GeoIndex *index, const void *data, size_t available, size_t *consumed);
size_t geo_density_index_serialized_size(const GeoIndex *index);
bool geo_density_index_serialize(const GeoIndex *index, void *destination, size_t size);
GeoDensityWriter *geo_density_writer_create(uint64_t record_count, uint8_t prefix_bits);
GeoDensityWriter *geo_density_writer_create_partition(uint64_t record_count,
                                                       uint8_t prefix_bits,
                                                       uint64_t first_position,
                                                       uint64_t partition_records,
                                                       const size_t *root_offsets);
void geo_density_writer_destroy(GeoDensityWriter *writer);
bool geo_density_writer_add(GeoDensityWriter *writer,
                            const GeoRecord *records,
                            size_t count,
                            uint64_t first_position);
bool geo_density_writer_append_file(GeoDensityWriter *writer,
                                    FILE *file,
                                    uint64_t *serialized_bytes,
                                    uint64_t *checksum);
bool geo_density_writer_merge_partitions(GeoDensityWriter *writer,
                                         GeoDensityWriter *const *partitions,
                                         size_t partition_count);
unsigned geo_density_index_query(const GeoIndex *index, uint64_t morton, size_t *local_records);

bool geo_index_search_radius_plan_append(const GeoIndex *index,
                                         const GeoRadiusQueryPlan *plan,
                                         GeoSearchResult *result,
                                         size_t *matched_count,
                                         GeoSearchStats *stats);

bool geo_index_search_radius_plan_append_filtered(const GeoIndex *index,
                                                  const GeoRadiusQueryPlan *plan,
                                                  GeoSearchResult *result,
                                                  size_t *matched_count,
                                                  GeoSearchStats *stats,
                                                  GeoRecordBitFilter bit_filter,
                                                  GeoRecordFilter filter,
                                                  void *filter_context);

bool geo_index_search_radius_plan_ids_append_filtered(const GeoIndex *index,
                                                      const GeoRadiusQueryPlan *plan,
                                                      GeoIdResult *result,
                                                      bool allow_growth,
                                                      size_t *matched_count,
                                                      GeoSearchStats *stats,
                                                      GeoRecordBitFilter bit_filter,
                                                      GeoRecordFilter filter,
                                                      void *filter_context);

bool geo_segment_set_read_snapshot_acquire(const GeoSegmentSet *set);
void geo_segment_set_read_snapshot_release(const GeoSegmentSet *set);
bool geo_segment_set_search_radius_count_snapshot(const GeoSegmentSet *set,
                                                  double latitude,
                                                  double longitude,
                                                  double radius_km,
                                                  size_t *count);
bool geo_segment_set_search_radius_ids_into_snapshot(const GeoSegmentSet *set,
                                                     double latitude,
                                                     double longitude,
                                                     double radius_km,
                                                     uint64_t *ids,
                                                     size_t capacity,
                                                     size_t *count);

bool geo_index_search_bbox_append(const GeoIndex *index,
                                  double min_lat,
                                  double max_lat,
                                  double min_lng,
                                  double max_lng,
                                  GeoSearchResult *result,
                                  GeoSearchStats *stats);

bool geo_index_prepare_bbox_query(double min_latitude,
                                  double max_latitude,
                                  double min_longitude,
                                  double max_longitude,
                                  GeoBboxQueryPlan *plan);

bool geo_index_search_bbox_plan_append(const GeoIndex *index,
                                       const GeoBboxQueryPlan *plan,
                                       GeoSearchResult *result,
                                       size_t *matched_count,
                                       GeoSearchStats *stats);

bool geo_index_search_bbox_plan_append_filtered(const GeoIndex *index,
                                                const GeoBboxQueryPlan *plan,
                                                GeoSearchResult *result,
                                                size_t *matched_count,
                                                GeoSearchStats *stats,
                                                GeoRecordBitFilter bit_filter,
                                                GeoRecordFilter filter,
                                                void *filter_context);

GeoSearchResult *geo_index_search_knn_filtered(const GeoIndex *index,
                                               double latitude,
                                               double longitude,
                                               size_t k,
                                               double max_radius_km,
                                               GeoSearchStats *stats,
                                               GeoRecordFilter filter,
                                               void *filter_context);

bool geo_index_search_knn_reuse_filtered(const GeoIndex *index,
                                         double latitude,
                                         double longitude,
                                         size_t k,
                                         double max_radius_km,
                                         GeoSearchResult *result,
                                         GeoKnnWorkspace *workspace,
                                         GeoSearchStats *stats,
                                         GeoRecordFilter filter,
                                         void *filter_context);

bool geo_index_search_knn_sources(const GeoKnnSource *sources,
                                  size_t source_count,
                                  double latitude,
                                  double longitude,
                                  size_t k,
                                  double max_radius_km,
                                  GeoSearchResult *result,
                                  GeoKnnWorkspace *workspace,
                                  GeoSearchStats *stats);

#endif
