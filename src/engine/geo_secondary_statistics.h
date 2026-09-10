#ifndef GEO_SECONDARY_STATISTICS_H
#define GEO_SECONDARY_STATISTICS_H

#include "geobolt/geobolt.h"
#include "geo_rocks_bridge.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define GEO_SECONDARY_HISTOGRAM_MAX_BINS 64U

typedef struct {
    uint64_t counts[GEO_SECONDARY_HISTOGRAM_MAX_BINS];
    uint32_t distinct_values[GEO_SECONDARY_HISTOGRAM_MAX_BINS];
    uint32_t boundary_offsets[GEO_SECONDARY_HISTOGRAM_MAX_BINS + 1U];
    unsigned char *boundaries;
    uint32_t boundary_size;
    uint32_t bin_count;
} GeoSecondaryHistogram;

void geo_secondary_histogram_destroy(GeoSecondaryHistogram *histogram);

int geo_secondary_compare_encoded(const unsigned char *first,
                                  size_t first_size,
                                  const unsigned char *second,
                                  size_t second_size);
bool geo_secondary_encoded_key_valid(GeoDatabaseIndexType type,
                                     const unsigned char *value,
                                     size_t value_size);

GeoDatabaseStatus geo_secondary_histogram_build(GeoRocksDatabase *rocks,
                                                uint64_t index_id,
                                                GeoDatabaseIndexType type,
                                                uint64_t entry_count,
                                                GeoSecondaryHistogram *histogram);
GeoDatabaseStatus geo_secondary_histogram_load_metadata(GeoRocksDatabase *rocks,
                                                        uint64_t index_id,
                                                        GeoDatabaseIndexType type,
                                                        uint32_t bin_count,
                                                        uint32_t boundary_size,
                                                        GeoSecondaryHistogram *histogram);
GeoDatabaseStatus geo_secondary_histogram_load_counts(GeoRocksDatabase *rocks,
                                                      uint64_t index_id,
                                                      uint64_t entry_count,
                                                      GeoSecondaryHistogram *histogram);
GeoDatabaseStatus geo_secondary_histogram_put_metadata(GeoRocksBatch *batch,
                                                       uint64_t index_id,
                                                       GeoDatabaseIndexType type,
                                                       const GeoSecondaryHistogram *histogram);
GeoDatabaseStatus geo_secondary_histogram_put_all_counts(GeoRocksBatch *batch,
                                                         uint64_t index_id,
                                                         const GeoSecondaryHistogram *histogram);
GeoDatabaseStatus geo_secondary_histogram_put_count(GeoRocksBatch *batch,
                                                    uint64_t index_id,
                                                    uint32_t bin,
                                                    uint64_t count);
GeoDatabaseStatus geo_secondary_histogram_delete(GeoRocksBatch *batch, uint64_t index_id);

uint32_t geo_secondary_histogram_find_bin(const GeoSecondaryHistogram *histogram,
                                         const unsigned char *encoded_value,
                                         size_t encoded_size);
uint64_t geo_secondary_histogram_estimate(const GeoSecondaryHistogram *histogram,
                                         GeoDatabaseIndexOperator operation,
                                         const unsigned char *lower,
                                         size_t lower_size,
                                         const unsigned char *upper,
                                         size_t upper_size);

#endif
