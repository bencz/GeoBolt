#ifndef GEO_OBJECT_CACHE_H
#define GEO_OBJECT_CACHE_H

#include "geobolt/geodoc.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct GeoObjectCache GeoObjectCache;
typedef struct GeoObjectCacheEntry GeoObjectCacheEntry;

typedef struct {
    GeoDocView document;
    uint64_t sequence;
    uint64_t morton_code;
} GeoObjectCacheView;

/*
 * Cache lifetime and writer replacement are protected by the database visibility gate.
 * Concurrent readers may only populate empty slots. Consequently, a returned view remains
 * valid until the caller releases its visibility read lock.
 */
GeoObjectCache *geo_object_cache_create(size_t byte_limit);
void geo_object_cache_destroy(GeoObjectCache *cache);

GeoObjectCacheEntry *geo_object_cache_entry_create(uint64_t object_id,
                                                   uint64_t sequence,
                                                   uint64_t morton_code,
                                                   GeoDocView document);
void geo_object_cache_entry_destroy(GeoObjectCacheEntry *entry);

bool geo_object_cache_find(const GeoObjectCache *cache, uint64_t object_id, GeoObjectCacheView *view);

/* Takes ownership of entry, whether or not the entry can be retained. */
void geo_object_cache_fill(GeoObjectCache *cache, GeoObjectCacheEntry *entry);

/* Writer-gate-only operations. publish takes ownership of entry. */
void geo_object_cache_publish(GeoObjectCache *cache, GeoObjectCacheEntry *entry);
void geo_object_cache_remove(GeoObjectCache *cache, uint64_t object_id);

#endif
