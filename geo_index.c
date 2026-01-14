#include "geo_index.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef __APPLE__
#include <mach/mach_time.h>
#endif

// =========================================================
// Internal Macros
// =========================================================

#define GEO_MIN(a, b) ((a) < (b) ? (a) : (b))
#define GEO_MAX(a, b) ((a) > (b) ? (a) : (b))
#define GEO_CLAMP(v, lo, hi) GEO_MIN(GEO_MAX(v, lo), hi)

// Prefetch hint for cache optimization
#ifdef __GNUC__
#define GEO_PREFETCH(addr) __builtin_prefetch(addr, 0, 3)
#else
#define GEO_PREFETCH(addr) ((void)0)
#endif

// Branch prediction hints
#ifdef __GNUC__
#define GEO_LIKELY(x)   __builtin_expect(!!(x), 1)
#define GEO_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#define GEO_LIKELY(x)   (x)
#define GEO_UNLIKELY(x) (x)
#endif

// =========================================================
// Morton Code Functions (Z-order curve)
// Optimized bit interleaving using magic numbers
// =========================================================

uint64_t geo_spread_bits(uint32_t v) {
    uint64_t x = v;
    x = (x | (x << 32)) & 0x00000000FFFFFFFFULL;
    x = (x | (x << 16)) & 0x0000FFFF0000FFFFULL;
    x = (x | (x << 8))  & 0x00FF00FF00FF00FFULL;
    x = (x | (x << 4))  & 0x0F0F0F0F0F0F0F0FULL;
    x = (x | (x << 2))  & 0x3333333333333333ULL;
    x = (x | (x << 1))  & 0x5555555555555555ULL;
    return x;
}

uint32_t geo_compact_bits(uint64_t v) {
    uint64_t x = v & 0x5555555555555555ULL;
    x = (x | (x >> 1))  & 0x3333333333333333ULL;
    x = (x | (x >> 2))  & 0x0F0F0F0F0F0F0F0FULL;
    x = (x | (x >> 4))  & 0x00FF00FF00FF00FFULL;
    x = (x | (x >> 8))  & 0x0000FFFF0000FFFFULL;
    x = (x | (x >> 16)) & 0x00000000FFFFFFFFULL;
    return (uint32_t)x;
}

// =========================================================
// Coordinate Normalization
// Using fixed-point arithmetic for precision
// =========================================================

uint32_t geo_normalize_lat(double lat) {
    lat = GEO_CLAMP(lat, GEO_MIN_LAT, GEO_MAX_LAT);
    double normalized = (lat - GEO_MIN_LAT) / (GEO_MAX_LAT - GEO_MIN_LAT);
    return (uint32_t)(normalized * (double)GEO_COORD_MAX);
}

uint32_t geo_normalize_lng(double lng) {
    lng = GEO_CLAMP(lng, GEO_MIN_LNG, GEO_MAX_LNG);
    double normalized = (lng - GEO_MIN_LNG) / (GEO_MAX_LNG - GEO_MIN_LNG);
    return (uint32_t)(normalized * (double)GEO_COORD_MAX);
}

double geo_denormalize_lat(uint32_t v) {
    return ((double)v / (double)GEO_COORD_MAX) * (GEO_MAX_LAT - GEO_MIN_LAT) + GEO_MIN_LAT;
}

double geo_denormalize_lng(uint32_t v) {
    return ((double)v / (double)GEO_COORD_MAX) * (GEO_MAX_LNG - GEO_MIN_LNG) + GEO_MIN_LNG;
}

// =========================================================
// Encode / Decode
// =========================================================

uint64_t geo_encode(double lat, double lng) {
    uint32_t lat_norm = geo_normalize_lat(lat);
    uint32_t lng_norm = geo_normalize_lng(lng);
    return geo_spread_bits(lat_norm) | (geo_spread_bits(lng_norm) << 1);
}

GeoPoint geo_decode(uint64_t z) {
    GeoPoint p;
    p.lat = geo_denormalize_lat(geo_compact_bits(z));
    p.lng = geo_denormalize_lng(geo_compact_bits(z >> 1));
    return p;
}

// =========================================================
// Distance Calculations
// =========================================================

double geo_to_radians(double degrees) {
    return degrees * (M_PI / 180.0);
}

double geo_to_degrees(double radians) {
    return radians * (180.0 / M_PI);
}

double geo_haversine_km(double lat1, double lng1, double lat2, double lng2) {
    double lat1_rad = geo_to_radians(lat1);
    double lat2_rad = geo_to_radians(lat2);
    double dlat = geo_to_radians(lat2 - lat1);
    double dlng = geo_to_radians(lng2 - lng1);
    
    double sin_dlat_2 = sin(dlat * 0.5);
    double sin_dlng_2 = sin(dlng * 0.5);
    
    double a = sin_dlat_2 * sin_dlat_2 +
               cos(lat1_rad) * cos(lat2_rad) * sin_dlng_2 * sin_dlng_2;
    
    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
    
    return GEO_EARTH_RADIUS_KM * c;
}

double geo_haversine_m(double lat1, double lng1, double lat2, double lng2) {
    return geo_haversine_km(lat1, lng1, lat2, lng2) * 1000.0;
}

double geo_fast_distance_km(double lat1, double lng1, double lat2, double lng2) {
    double lat_mid = geo_to_radians((lat1 + lat2) * 0.5);
    double cos_lat = cos(lat_mid);
    
    double dlat = (lat2 - lat1) * 111.32;  // km per degree latitude
    double dlng = (lng2 - lng1) * 111.32 * cos_lat;  // adjusted for latitude
    
    return sqrt(dlat * dlat + dlng * dlng);
}

void geo_bounding_box(double lat, double lng, double radius_km,
                      double *min_lat, double *max_lat,
                      double *min_lng, double *max_lng) {
    // Latitude: 1 degree ≈ 111.32 km
    double lat_delta = radius_km / 111.32;
    
    // Longitude: depends on latitude
    double cos_lat = cos(geo_to_radians(lat));
    if (fabs(cos_lat) < 1e-10) cos_lat = 1e-10;
    double lng_delta = radius_km / (111.32 * cos_lat);
    
    *min_lat = GEO_MAX(lat - lat_delta, GEO_MIN_LAT);
    *max_lat = GEO_MIN(lat + lat_delta, GEO_MAX_LAT);
    *min_lng = GEO_MAX(lng - lng_delta, GEO_MIN_LNG);
    *max_lng = GEO_MIN(lng + lng_delta, GEO_MAX_LNG);
}

// =========================================================
// Index Management
// =========================================================

GeoIndex* geo_index_create(size_t initial_capacity) {
    GeoIndex *index = (GeoIndex*)malloc(sizeof(GeoIndex));
    if (GEO_UNLIKELY(!index)) return NULL;
    
    if (initial_capacity == 0) initial_capacity = 1024;
    
    index->records = (GeoRecord*)aligned_alloc(16, initial_capacity * sizeof(GeoRecord));
    if (GEO_UNLIKELY(!index->records)) {
        free(index);
        return NULL;
    }
    
    index->count = 0;
    index->capacity = initial_capacity;
    index->sorted = false;
    
    return index;
}

void geo_index_destroy(GeoIndex *index) {
    if (index) {
        free(index->records);
        free(index);
    }
}

static bool geo_index_grow(GeoIndex *index) {
    size_t new_capacity = index->capacity * 2;
    GeoRecord *new_records = (GeoRecord*)aligned_alloc(16, new_capacity * sizeof(GeoRecord));
    if (GEO_UNLIKELY(!new_records)) return false;
    
    memcpy(new_records, index->records, index->count * sizeof(GeoRecord));
    free(index->records);
    index->records = new_records;
    index->capacity = new_capacity;
    
    return true;
}

bool geo_index_add(GeoIndex *index, uint64_t id, double lat, double lng) {
    if (GEO_UNLIKELY(!index)) return false;
    
    if (GEO_UNLIKELY(index->count >= index->capacity)) {
        if (!geo_index_grow(index)) return false;
    }
    
    index->records[index->count].id = id;
    index->records[index->count].z = geo_encode(lat, lng);
    index->count++;
    index->sorted = false;
    
    return true;
}

bool geo_index_add_batch(GeoIndex *index, const uint64_t *ids,
                         const double *lats, const double *lngs, size_t count) {
    if (GEO_UNLIKELY(!index || !ids || !lats || !lngs)) return false;
    
    while (index->count + count > index->capacity) {
        if (!geo_index_grow(index)) return false;
    }
    
    for (size_t i = 0; i < count; i++) {
        index->records[index->count + i].id = ids[i];
        index->records[index->count + i].z = geo_encode(lats[i], lngs[i]);
    }
    
    index->count += count;
    index->sorted = false;
    
    return true;
}

int geo_compare_records_by_z(const void *a, const void *b) {
    const GeoRecord *r1 = (const GeoRecord*)a;
    const GeoRecord *r2 = (const GeoRecord*)b;
    if (r1->z < r2->z) return -1;
    if (r1->z > r2->z) return 1;
    return 0;
}

int geo_compare_records_by_id(const void *a, const void *b) {
    const GeoRecord *r1 = (const GeoRecord*)a;
    const GeoRecord *r2 = (const GeoRecord*)b;
    if (r1->id < r2->id) return -1;
    if (r1->id > r2->id) return 1;
    return 0;
}

void geo_index_build(GeoIndex *index) {
    if (!index || index->sorted) return;
    qsort(index->records, index->count, sizeof(GeoRecord), geo_compare_records_by_z);
    index->sorted = true;
}

void geo_index_clear(GeoIndex *index) {
    if (index) {
        index->count = 0;
        index->sorted = false;
    }
}

// =========================================================
// Binary Search - Optimized with prefetch
// =========================================================

size_t geo_lower_bound(const GeoRecord *records, size_t n, uint64_t key) {
    if (n == 0) return 0;
    
    size_t lo = 0;
    size_t hi = n;
    
    while (lo < hi) {
        size_t mid = lo + ((hi - lo) >> 1);
        
        // Prefetch likely next access
        if (hi - lo > 8) {
            GEO_PREFETCH(&records[lo + ((mid - lo) >> 1)]);
            GEO_PREFETCH(&records[mid + ((hi - mid) >> 1)]);
        }
        
        if (records[mid].z < key) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    
    return lo;
}

size_t geo_upper_bound(const GeoRecord *records, size_t n, uint64_t key) {
    if (n == 0) return 0;
    
    size_t lo = 0;
    size_t hi = n;
    
    while (lo < hi) {
        size_t mid = lo + ((hi - lo) >> 1);
        
        // Prefetch likely next access
        if (hi - lo > 8) {
            GEO_PREFETCH(&records[lo + ((mid - lo) >> 1)]);
            GEO_PREFETCH(&records[mid + ((hi - mid) >> 1)]);
        }
        
        if (records[mid].z <= key) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    
    return lo;
}

// =========================================================
// Range Building
// =========================================================

int geo_build_ranges(double lat, double lng, double radius_km,
                     ZRange *out, int max_ranges) {
    if (!out || max_ranges < 1) return 0;
    
    double min_lat, max_lat, min_lng, max_lng;
    geo_bounding_box(lat, lng, radius_km, &min_lat, &max_lat, &min_lng, &max_lng);
    
    // Simple approach: single bounding box range
    uint64_t z_min = geo_encode(min_lat, min_lng);
    uint64_t z_max = geo_encode(max_lat, max_lng);
    
    // Ensure min <= max
    if (z_min > z_max) {
        uint64_t tmp = z_min;
        z_min = z_max;
        z_max = tmp;
    }
    
    out[0].min = z_min;
    out[0].max = z_max;
    
    if (max_ranges >= 4) {
        // Add corner ranges for better coverage
        uint64_t corners[4] = {
            geo_encode(min_lat, min_lng),
            geo_encode(min_lat, max_lng),
            geo_encode(max_lat, min_lng),
            geo_encode(max_lat, max_lng)
        };
        
        uint64_t center = geo_encode(lat, lng);
        
        int count = 0;
        for (int i = 0; i < 4 && count < max_ranges; i++) {
            uint64_t a = GEO_MIN(center, corners[i]);
            uint64_t b = GEO_MAX(center, corners[i]);
            
            // Check for overlap with existing ranges
            bool overlaps = false;
            for (int j = 0; j < count; j++) {
                if (a <= out[j].max && b >= out[j].min) {
                    // Merge ranges
                    out[j].min = GEO_MIN(out[j].min, a);
                    out[j].max = GEO_MAX(out[j].max, b);
                    overlaps = true;
                    break;
                }
            }
            
            if (!overlaps) {
                out[count].min = a;
                out[count].max = b;
                count++;
            }
        }
        
        return count > 0 ? count : 1;
    }
    
    return 1;
}

int geo_build_ranges_adaptive(double lat, double lng, double radius_km,
                              ZRange *out, int max_ranges, int precision_level) {
    if (!out || max_ranges < 1) return 0;
    
    double min_lat, max_lat, min_lng, max_lng;
    geo_bounding_box(lat, lng, radius_km, &min_lat, &max_lat, &min_lng, &max_lng);
    
    // Subdivide based on precision level
    int divisions = 1 << precision_level;  // 2^precision_level
    if (divisions * divisions > max_ranges) {
        divisions = (int)sqrt((double)max_ranges);
    }
    
    double lat_step = (max_lat - min_lat) / divisions;
    double lng_step = (max_lng - min_lng) / divisions;
    
    int count = 0;
    for (int i = 0; i < divisions && count < max_ranges; i++) {
        for (int j = 0; j < divisions && count < max_ranges; j++) {
            double cell_min_lat = min_lat + i * lat_step;
            double cell_max_lat = min_lat + (i + 1) * lat_step;
            double cell_min_lng = min_lng + j * lng_step;
            double cell_max_lng = min_lng + (j + 1) * lng_step;
            
            uint64_t z1 = geo_encode(cell_min_lat, cell_min_lng);
            uint64_t z2 = geo_encode(cell_max_lat, cell_max_lng);
            
            out[count].min = GEO_MIN(z1, z2);
            out[count].max = GEO_MAX(z1, z2);
            count++;
        }
    }
    
    return count;
}

// =========================================================
// Result Management
// =========================================================

GeoSearchResult* geo_result_create(size_t initial_capacity) {
    GeoSearchResult *result = (GeoSearchResult*)malloc(sizeof(GeoSearchResult));
    if (GEO_UNLIKELY(!result)) return NULL;
    
    if (initial_capacity == 0) initial_capacity = 64;
    
    result->results = (GeoRecord*)malloc(initial_capacity * sizeof(GeoRecord));
    if (GEO_UNLIKELY(!result->results)) {
        free(result);
        return NULL;
    }
    
    result->count = 0;
    result->capacity = initial_capacity;
    
    return result;
}

void geo_result_destroy(GeoSearchResult *result) {
    if (result) {
        free(result->results);
        free(result);
    }
}

bool geo_result_add(GeoSearchResult *result, const GeoRecord *record) {
    if (GEO_UNLIKELY(!result || !record)) return false;
    
    if (GEO_UNLIKELY(result->count >= result->capacity)) {
        size_t new_capacity = result->capacity * 2;
        GeoRecord *new_results = (GeoRecord*)realloc(result->results, 
                                                      new_capacity * sizeof(GeoRecord));
        if (GEO_UNLIKELY(!new_results)) return false;
        result->results = new_results;
        result->capacity = new_capacity;
    }
    
    result->results[result->count++] = *record;
    return true;
}

// Context for distance sorting
typedef struct {
    double lat;
    double lng;
} DistanceContext;

static DistanceContext g_dist_ctx;

static int compare_by_distance(const void *a, const void *b) {
    const GeoRecord *r1 = (const GeoRecord*)a;
    const GeoRecord *r2 = (const GeoRecord*)b;
    
    GeoPoint p1 = geo_decode(r1->z);
    GeoPoint p2 = geo_decode(r2->z);
    
    double d1 = geo_haversine_km(g_dist_ctx.lat, g_dist_ctx.lng, p1.lat, p1.lng);
    double d2 = geo_haversine_km(g_dist_ctx.lat, g_dist_ctx.lng, p2.lat, p2.lng);
    
    if (d1 < d2) return -1;
    if (d1 > d2) return 1;
    return 0;
}

void geo_result_sort_by_distance(GeoSearchResult *result, double lat, double lng) {
    if (!result || result->count == 0) return;
    
    g_dist_ctx.lat = lat;
    g_dist_ctx.lng = lng;
    
    qsort(result->results, result->count, sizeof(GeoRecord), compare_by_distance);
}

// =========================================================
// Search Functions
// =========================================================

GeoSearchResult* geo_search_radius(const GeoIndex *index,
                                   double lat, double lng, double radius_km,
                                   GeoSearchStats *stats) {
    if (GEO_UNLIKELY(!index || !index->sorted)) return NULL;
    
    double start_time = geo_get_time_ms();
    
    GeoSearchResult *result = geo_result_create(64);
    if (GEO_UNLIKELY(!result)) return NULL;
    
    // Get bounding box
    double min_lat, max_lat, min_lng, max_lng;
    geo_bounding_box(lat, lng, radius_km, &min_lat, &max_lat, &min_lng, &max_lng);
    
    // Build ranges
    ZRange ranges[16];
    int range_count = geo_build_ranges(lat, lng, radius_km, ranges, 16);
    
    uint64_t scanned = 0;
    uint64_t matched = 0;
    
    // Track seen IDs to avoid duplicates
    uint64_t last_id = UINT64_MAX;
    
    for (int r = 0; r < range_count; r++) {
        size_t start = geo_lower_bound(index->records, index->count, ranges[r].min);
        
        for (size_t i = start; i < index->count && index->records[i].z <= ranges[r].max; i++) {
            scanned++;
            
            // Skip duplicates
            if (index->records[i].id == last_id) continue;
            last_id = index->records[i].id;
            
            GeoPoint p = geo_decode(index->records[i].z);
            
            // Quick bounding box check
            if (p.lat < min_lat || p.lat > max_lat ||
                p.lng < min_lng || p.lng > max_lng) {
                continue;
            }
            
            // Precise distance check
            double dist = geo_haversine_km(lat, lng, p.lat, p.lng);
            if (dist <= radius_km) {
                geo_result_add(result, &index->records[i]);
                matched++;
            }
        }
    }
    
    if (stats) {
        stats->records_scanned = scanned;
        stats->records_matched = matched;
        stats->ranges_checked = range_count;
        stats->search_time_ms = geo_get_time_ms() - start_time;
    }
    
    return result;
}

GeoSearchResult* geo_search_bbox(const GeoIndex *index,
                                 double min_lat, double max_lat,
                                 double min_lng, double max_lng,
                                 GeoSearchStats *stats) {
    if (GEO_UNLIKELY(!index || !index->sorted)) return NULL;
    
    double start_time = geo_get_time_ms();
    
    GeoSearchResult *result = geo_result_create(64);
    if (GEO_UNLIKELY(!result)) return NULL;
    
    // Encode corners
    uint64_t z_min = geo_encode(min_lat, min_lng);
    uint64_t z_max = geo_encode(max_lat, max_lng);
    
    if (z_min > z_max) {
        uint64_t tmp = z_min;
        z_min = z_max;
        z_max = tmp;
    }
    
    uint64_t scanned = 0;
    uint64_t matched = 0;
    
    size_t start = geo_lower_bound(index->records, index->count, z_min);
    
    for (size_t i = start; i < index->count && index->records[i].z <= z_max; i++) {
        scanned++;
        
        GeoPoint p = geo_decode(index->records[i].z);
        
        if (p.lat >= min_lat && p.lat <= max_lat &&
            p.lng >= min_lng && p.lng <= max_lng) {
            geo_result_add(result, &index->records[i]);
            matched++;
        }
    }
    
    if (stats) {
        stats->records_scanned = scanned;
        stats->records_matched = matched;
        stats->ranges_checked = 1;
        stats->search_time_ms = geo_get_time_ms() - start_time;
    }
    
    return result;
}

GeoSearchResult* geo_search_knn(const GeoIndex *index,
                                double lat, double lng, size_t k,
                                double max_radius_km,
                                GeoSearchStats *stats) {
    if (GEO_UNLIKELY(!index || !index->sorted || k == 0)) return NULL;
    
    // Start with a small radius and expand if needed
    double radius = max_radius_km / 10.0;
    if (radius < 1.0) radius = 1.0;
    
    GeoSearchResult *result = NULL;
    
    while (radius <= max_radius_km) {
        if (result) geo_result_destroy(result);
        
        result = geo_search_radius(index, lat, lng, radius, stats);
        
        if (result && result->count >= k) {
            break;
        }
        
        radius *= 2.0;
    }
    
    if (result && result->count > 0) {
        geo_result_sort_by_distance(result, lat, lng);
        
        // Trim to k results
        if (result->count > k) {
            result->count = k;
        }
    }
    
    return result;
}

// =========================================================
// Utility Functions
// =========================================================

double geo_get_time_ms(void) {
#ifdef __APPLE__
    static mach_timebase_info_data_t timebase = {0};
    if (timebase.denom == 0) {
        mach_timebase_info(&timebase);
    }
    uint64_t time = mach_absolute_time();
    return (double)(time * timebase.numer / timebase.denom) / 1000000.0;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1000000.0;
#endif
}

bool geo_is_valid_lat(double lat) {
    return lat >= GEO_MIN_LAT && lat <= GEO_MAX_LAT && !isnan(lat) && !isinf(lat);
}

bool geo_is_valid_lng(double lng) {
    return lng >= GEO_MIN_LNG && lng <= GEO_MAX_LNG && !isnan(lng) && !isinf(lng);
}

bool geo_is_valid_point(double lat, double lng) {
    return geo_is_valid_lat(lat) && geo_is_valid_lng(lng);
}

double geo_clamp_lat(double lat) {
    return GEO_CLAMP(lat, GEO_MIN_LAT, GEO_MAX_LAT);
}

double geo_clamp_lng(double lng) {
    return GEO_CLAMP(lng, GEO_MIN_LNG, GEO_MAX_LNG);
}

double geo_wrap_lng(double lng) {
    while (lng > 180.0) lng -= 360.0;
    while (lng < -180.0) lng += 360.0;
    return lng;
}
