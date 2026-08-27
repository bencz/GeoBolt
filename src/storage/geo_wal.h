#ifndef GEO_WAL_H
#define GEO_WAL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct GeoWal GeoWal;

typedef bool (*GeoWalReplayFunction)(void *context,
                                     uint64_t first_sequence,
                                     uint32_t entry_count,
                                     const void *payload,
                                     uint32_t payload_size);

GeoWal *geo_wal_open(const char *directory,
                     size_t segment_size,
                     GeoWalReplayFunction replay,
                     void *replay_context,
                     uint64_t *next_sequence);
void geo_wal_close(GeoWal *wal);

bool geo_wal_append(GeoWal *wal,
                    uint64_t first_sequence,
                    uint32_t entry_count,
                    const void *payload,
                    uint32_t payload_size);
bool geo_wal_sync(GeoWal *wal);
bool geo_wal_checkpoint(GeoWal *wal, uint64_t next_sequence);

#endif
