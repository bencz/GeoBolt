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

// Precision: 32-bit per coordinate = ~0.00001 degree (~1.1m at equator)
#define GEO_COORD_BITS 32
#define GEO_COORD_MAX  4294967295ULL

// =========================================================
// Structures
// =========================================================

typedef struct {
    double lat;
    double lng;
} GeoPoint;

typedef struct __attribute__((aligned(16))) {
    uint64_t id;
    uint64_t z;  // Morton/Z-order code
} GeoRecord;

typedef struct {
    uint64_t min;
    uint64_t max;
} ZRange;

typedef struct {
    GeoRecord *records;
    size_t count;
    size_t capacity;
    bool sorted;
} GeoIndex;

typedef struct {
    GeoRecord *results;
    size_t count;
    size_t capacity;
} GeoSearchResult;

typedef struct {
    uint64_t records_scanned;
    uint64_t records_matched;
    uint64_t ranges_checked;
    double search_time_ms;
} GeoSearchStats;

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
bool geo_index_add(GeoIndex *index, uint64_t id, double lat, double lng);
bool geo_index_add_batch(GeoIndex *index, const uint64_t *ids, 
                         const double *lats, const double *lngs, size_t count);
void geo_index_build(GeoIndex *index);  // Sort and prepare for queries
void geo_index_clear(GeoIndex *index);

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

// Search within bounding box
GeoSearchResult* geo_search_bbox(const GeoIndex *index,
                                 double min_lat, double max_lat,
                                 double min_lng, double max_lng,
                                 GeoSearchStats *stats);

// K-nearest neighbors search
GeoSearchResult* geo_search_knn(const GeoIndex *index,
                                double lat, double lng, size_t k,
                                double max_radius_km,
                                GeoSearchStats *stats);

// =========================================================
// Result Management
// =========================================================

GeoSearchResult* geo_result_create(size_t initial_capacity);
void geo_result_destroy(GeoSearchResult *result);
bool geo_result_add(GeoSearchResult *result, const GeoRecord *record);
void geo_result_sort_by_distance(GeoSearchResult *result, double lat, double lng);

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

#ifdef __cplusplus
}
#endif

#endif // GEO_INDEX_H
