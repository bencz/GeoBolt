#include "geo_object_cache.h"

#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

#define GEO_OBJECT_CACHE_MINIMUM_SLOTS 16U
#define GEO_OBJECT_CACHE_BYTES_PER_SLOT 256U
#define GEO_OBJECT_CACHE_WAYS 8U

_Static_assert((GEO_OBJECT_CACHE_WAYS & (GEO_OBJECT_CACHE_WAYS - 1U)) == 0U,
               "object-cache associativity must be a power of two");
_Static_assert(GEO_OBJECT_CACHE_MINIMUM_SLOTS % GEO_OBJECT_CACHE_WAYS == 0U,
               "the minimum object-cache capacity must contain complete sets");

struct GeoObjectCacheEntry {
    size_t allocation_size;
    uint64_t object_id;
    uint64_t sequence;
    uint64_t morton_code;
    uint32_t document_size;
    unsigned char document[];
};

struct GeoObjectCache {
    _Atomic(GeoObjectCacheEntry *) *slots;
    atomic_size_t used_bytes;
    size_t data_limit;
    size_t capacity;
    size_t set_mask;
};

static uint64_t object_cache_hash(uint64_t value)
{
    value ^= value >> 30U;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27U;
    value *= UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31U);
}

static size_t object_cache_capacity(size_t byte_limit)
{
    size_t requested = byte_limit / GEO_OBJECT_CACHE_BYTES_PER_SLOT;
    size_t capacity = GEO_OBJECT_CACHE_MINIMUM_SLOTS;

    while (capacity <= requested / 2U) {
        capacity *= 2U;
    }

    return capacity;
}

static bool object_cache_reserve(GeoObjectCache *cache, size_t bytes)
{
    size_t used = atomic_load_explicit(&cache->used_bytes, memory_order_relaxed);

    for (;;) {
        if (used > cache->data_limit || bytes > cache->data_limit - used) {
            return false;
        }

        if (atomic_compare_exchange_weak_explicit(&cache->used_bytes,
                                                  &used,
                                                  used + bytes,
                                                  memory_order_relaxed,
                                                  memory_order_relaxed)) {
            return true;
        }
    }
}

static void object_cache_release(GeoObjectCache *cache, size_t bytes)
{
    (void) atomic_fetch_sub_explicit(&cache->used_bytes, bytes, memory_order_relaxed);
}

GeoObjectCache *geo_object_cache_create(size_t byte_limit)
{
    size_t capacity = object_cache_capacity(byte_limit);

    if (capacity > SIZE_MAX / sizeof(_Atomic(GeoObjectCacheEntry *))) {
        return NULL;
    }

    size_t slot_bytes = capacity * sizeof(_Atomic(GeoObjectCacheEntry *));

    if (byte_limit <= sizeof(GeoObjectCache) || slot_bytes > byte_limit - sizeof(GeoObjectCache)) {
        return NULL;
    }

    GeoObjectCache *cache = malloc(sizeof(*cache));

    if (!cache) {
        return NULL;
    }

    cache->slots = malloc(slot_bytes);

    if (!cache->slots) {
        free(cache);
        return NULL;
    }

    for (size_t index = 0U; index < capacity; ++index) {
        atomic_init(cache->slots + index, NULL);
    }

    atomic_init(&cache->used_bytes, 0U);
    cache->data_limit = byte_limit - sizeof(*cache) - slot_bytes;
    cache->capacity = capacity;
    cache->set_mask = capacity / GEO_OBJECT_CACHE_WAYS - 1U;
    return cache;
}

void geo_object_cache_destroy(GeoObjectCache *cache)
{
    if (!cache) {
        return;
    }

    for (size_t index = 0U; index < cache->capacity; ++index) {
        GeoObjectCacheEntry *entry = atomic_load_explicit(cache->slots + index, memory_order_relaxed);
        free(entry);
    }

    free(cache->slots);
    free(cache);
}

GeoObjectCacheEntry *geo_object_cache_entry_create(uint64_t object_id,
                                                   uint64_t sequence,
                                                   uint64_t morton_code,
                                                   GeoDocView document)
{
    if (!object_id || (document.size && !document.data) || document.size > UINT32_MAX ||
        document.size > SIZE_MAX - sizeof(GeoObjectCacheEntry)) {
        return NULL;
    }

    size_t allocation_size = sizeof(GeoObjectCacheEntry) + document.size;
    GeoObjectCacheEntry *entry = malloc(allocation_size);

    if (!entry) {
        return NULL;
    }

    entry->allocation_size = allocation_size;
    entry->object_id = object_id;
    entry->sequence = sequence;
    entry->morton_code = morton_code;
    entry->document_size = (uint32_t) document.size;
    if (document.size) {
        memcpy(entry->document, document.data, document.size);
    }
    return entry;
}

void geo_object_cache_entry_destroy(GeoObjectCacheEntry *entry)
{
    free(entry);
}

bool geo_object_cache_find(const GeoObjectCache *cache, uint64_t object_id, GeoObjectCacheView *view)
{
    if (!cache || !object_id || !view) {
        return false;
    }

    size_t first_slot = ((size_t) object_cache_hash(object_id) & cache->set_mask) * GEO_OBJECT_CACHE_WAYS;

    for (size_t way = 0U; way < GEO_OBJECT_CACHE_WAYS; ++way) {
        GeoObjectCacheEntry *entry = atomic_load_explicit(cache->slots + first_slot + way, memory_order_acquire);

        if (entry && entry->object_id == object_id) {
            *view = (GeoObjectCacheView) {
                .document = { .data = entry->document, .size = entry->document_size },
                .sequence = entry->sequence,
                .morton_code = entry->morton_code,
            };
            return true;
        }
    }

    return false;
}

void geo_object_cache_fill(GeoObjectCache *cache, GeoObjectCacheEntry *entry)
{
    if (!entry) {
        return;
    }
    if (!cache || !object_cache_reserve(cache, entry->allocation_size)) {
        free(entry);
        return;
    }

    size_t first_slot = ((size_t) object_cache_hash(entry->object_id) & cache->set_mask) * GEO_OBJECT_CACHE_WAYS;

    for (size_t way = 0U; way < GEO_OBJECT_CACHE_WAYS; ++way) {
        GeoObjectCacheEntry *expected = atomic_load_explicit(cache->slots + first_slot + way, memory_order_acquire);

        if (expected && expected->object_id == entry->object_id) {
            object_cache_release(cache, entry->allocation_size);
            free(entry);
            return;
        }

        if (!expected) {
            if (atomic_compare_exchange_strong_explicit(cache->slots + first_slot + way,
                                                        &expected,
                                                        entry,
                                                        memory_order_release,
                                                        memory_order_acquire)) {
                return;
            }

            /* A failed CAS returns the entry published by the competing reader. Stop here
             * when both misses resolved the same object instead of consuming another way. */
            if (expected && expected->object_id == entry->object_id) {
                object_cache_release(cache, entry->allocation_size);
                free(entry);
                return;
            }
        }
    }

    object_cache_release(cache, entry->allocation_size);
    free(entry);
}

void geo_object_cache_publish(GeoObjectCache *cache, GeoObjectCacheEntry *entry)
{
    if (!entry) {
        return;
    }
    if (!cache) {
        free(entry);
        return;
    }

    uint64_t hash = object_cache_hash(entry->object_id);
    size_t first_slot = ((size_t) hash & cache->set_mask) * GEO_OBJECT_CACHE_WAYS;
    size_t target_way = GEO_OBJECT_CACHE_WAYS;
    size_t reclaimed = 0U;

    for (size_t way = 0U; way < GEO_OBJECT_CACHE_WAYS; ++way) {
        GeoObjectCacheEntry *old_entry = atomic_load_explicit(cache->slots + first_slot + way, memory_order_relaxed);

        if (old_entry && old_entry->object_id == entry->object_id) {
            atomic_store_explicit(cache->slots + first_slot + way, NULL, memory_order_relaxed);
            reclaimed += old_entry->allocation_size;
            free(old_entry);

            if (target_way == GEO_OBJECT_CACHE_WAYS) {
                target_way = way;
            }
        } else if (!old_entry && target_way == GEO_OBJECT_CACHE_WAYS) {
            target_way = way;
        }
    }

    if (target_way == GEO_OBJECT_CACHE_WAYS) {
        target_way = (size_t) (hash >> 32U) & (GEO_OBJECT_CACHE_WAYS - 1U);
        GeoObjectCacheEntry *victim = atomic_load_explicit(cache->slots + first_slot + target_way, memory_order_relaxed);

        atomic_store_explicit(cache->slots + first_slot + target_way, NULL, memory_order_relaxed);
        reclaimed += victim->allocation_size;
        free(victim);
    }

    size_t used = atomic_load_explicit(&cache->used_bytes, memory_order_relaxed);
    size_t available = used <= cache->data_limit && reclaimed <= used
                           ? cache->data_limit - used + reclaimed
                           : 0U;

    if (entry->allocation_size > available) {
        object_cache_release(cache, reclaimed);
        free(entry);
        return;
    }

    atomic_store_explicit(cache->slots + first_slot + target_way, entry, memory_order_release);
    atomic_store_explicit(&cache->used_bytes, used - reclaimed + entry->allocation_size, memory_order_relaxed);
}

void geo_object_cache_remove(GeoObjectCache *cache, uint64_t object_id)
{
    if (!cache || !object_id) {
        return;
    }

    size_t first_slot = ((size_t) object_cache_hash(object_id) & cache->set_mask) * GEO_OBJECT_CACHE_WAYS;

    for (size_t way = 0U; way < GEO_OBJECT_CACHE_WAYS; ++way) {
        GeoObjectCacheEntry *entry = atomic_load_explicit(cache->slots + first_slot + way, memory_order_relaxed);

        if (entry && entry->object_id == object_id) {
            atomic_store_explicit(cache->slots + first_slot + way, NULL, memory_order_release);
            object_cache_release(cache, entry->allocation_size);
            free(entry);
        }
    }
}
