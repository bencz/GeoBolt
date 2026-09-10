#ifndef GEO_SPATIAL_MEMTABLE_H
#define GEO_SPATIAL_MEMTABLE_H

#include "geobolt/geobolt.h"
#include "geo_index_private.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct GeoSpatialMemtable GeoSpatialMemtable;

GeoSpatialMemtable *geo_spatial_memtable_create(size_t maximum_operations);
void geo_spatial_memtable_destroy(GeoSpatialMemtable *memtable);

bool geo_spatial_memtable_can_append(const GeoSpatialMemtable *memtable, size_t operation_count);
bool geo_spatial_memtable_append(GeoSpatialMemtable *memtable,
                                 const GeoDatabaseMutation *operations,
                                 size_t operation_count,
                                 uint64_t first_sequence);
bool geo_spatial_memtable_should_rotate(const GeoSpatialMemtable *memtable);
bool geo_spatial_memtable_rotate(GeoSpatialMemtable *memtable);
bool geo_spatial_memtable_has_active(const GeoSpatialMemtable *memtable);
bool geo_spatial_memtable_has_frozen(const GeoSpatialMemtable *memtable);

bool geo_spatial_memtable_export_frozen(const GeoSpatialMemtable *memtable,
                                        GeoDatabaseMutation **operations,
                                        size_t *operation_count,
                                        uint64_t *first_sequence,
                                        uint64_t *next_sequence);
void geo_spatial_memtable_release_frozen(GeoSpatialMemtable *memtable);

bool geo_spatial_memtable_allows_persisted_record(const GeoRecord *record, void *context);
bool geo_spatial_memtable_search_radius(const GeoSpatialMemtable *memtable,
                                        double latitude,
                                        double longitude,
                                        double radius_km,
                                        GeoSearchResult *result,
                                        size_t *count,
                                        GeoRecordFilter filter,
                                        void *filter_context,
                                        GeoSearchStats *stats);

size_t geo_spatial_memtable_active_operations(const GeoSpatialMemtable *memtable);
size_t geo_spatial_memtable_frozen_operations(const GeoSpatialMemtable *memtable);

#endif
