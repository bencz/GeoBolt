# GeoBolt Index

GeoBolt Index is a high-performance spatial indexing library focused on lightning-fast
geospatial lookups and analysis. Built for production workloads with **SIMD optimizations**
for ARM64 (NEON) and x86-64 (AVX2/SSE), it provides sub-centimeter precision using
64-bit Morton (Z-order) codes.

## Features

- **Sub-centimeter precision** – 32-bit per-axis normalization keeps errors under 1 cm
- **⚡ SIMD Optimized** – ARM64 NEON and x86-64 AVX2/SSE implementations with 10x+ speedups
- **Thread-safe reads** – Safe for concurrent search operations after index build
- **Cache-friendly** – 16-byte aligned records with prefetch hints for maximum throughput
- **Battle-tested** – 48 comprehensive tests covering correctness, stress, and performance
- **Proven scale** – Tested with 10+ million points, 40K+ searches/second

## Quick Start

### Requirements

- C11 compiler (clang recommended)
- make
- macOS or Linux
- pthread (for multi-threaded demo)

### Build & Run

```bash
# Build everything
make all

# Run comprehensive benchmark (10M points)
make demo

# Run test suite
make test

# Show all targets
make help
```

## API Reference

### Core Functions

```c
#include "geo_index.h"

// Create and populate index
GeoIndex *index = geo_index_create(1000000);
geo_index_add(index, id, lat, lng);
geo_index_build(index);  // Sort by Z-order - makes index read-only

// Search operations (thread-safe after build)
GeoSearchResult *r = geo_search_radius(index, lat, lng, radius_km, &stats);
GeoSearchResult *r = geo_search_knn(index, lat, lng, k, max_km, &stats);
GeoSearchResult *r = geo_search_bbox(index, min_lat, max_lat, min_lng, max_lng, &stats);

// Cleanup
geo_result_destroy(r);
geo_index_destroy(index);
```

### SIMD Batch Functions

```c
#include "geo_index.h"  // SIMD included automatically

// Batch encoding/decoding (up to 1.5x faster)
geo_simd_encode_batch(lats, lngs, z_codes, count);
geo_simd_decode_batch(z_codes, lats, lngs, count);

// Batch distance calculations (up to 5x faster)
geo_simd_haversine_batch(lat1, lng1, lats, lngs, distances, count);
geo_simd_fast_distance_batch(lat1, lng1, lats, lngs, distances, count);

// Batch filtering (up to 33x faster)
geo_simd_filter_radius(lats, lngs, count, center_lat, center_lng, radius_km, mask);
geo_simd_filter_bbox(lats, lngs, count, min_lat, max_lat, min_lng, max_lng, mask);

// SIMD info
bool available = geo_simd_available();
const char *name = geo_simd_get_name();  // "ARM64 NEON", "x86-64 AVX2", etc.
```

### Key Types

```c
typedef struct { double lat, lng; } GeoPoint;
typedef struct { uint64_t id, z; } GeoRecord;  // 16-byte aligned
typedef struct { GeoRecord *results; size_t count, capacity; } GeoSearchResult;
typedef struct { uint64_t records_scanned, records_matched; double search_time_ms; } GeoSearchStats;
```

## Performance Benchmarks

### SIMD vs Scalar (ARM64 Apple Silicon M1)

| Operation | Scalar | SIMD | Speedup |
|-----------|--------|------|---------|
| Encode (1M pts) | 2.08 ms | 1.43 ms | **1.45x** |
| Decode (1M pts) | 1.16 ms | 0.97 ms | **1.20x** |
| Haversine (1M) | 10.98 ms | 1.92 ms | **5.72x** |
| Filter Radius (1M) | 6.85 ms | 0.19 ms | **36.25x** |

**Average SIMD Speedup: 11x**

### Large Scale Performance

| Metric | Result |
|--------|--------|
| Index 10M points | 55 ms insert + 367 ms build |
| Radius search (50km, 10M pts) | 0.086 ms avg |
| Multi-thread speedup (8 threads) | **5.79x** |
| Searches per second | **41,780** |
| Encoding precision | < 0.01 meters |

### Thread Safety

```
✓ Reads are thread-safe after geo_index_build()
✗ Writes (geo_index_add) are NOT thread-safe
```

The library is designed for a **build-once, query-many** pattern:
1. Build the index in a single thread
2. Call `geo_index_build()` to finalize
3. Query from multiple threads concurrently

## Build Targets

```bash
make all            # Build everything
make lib            # Static library only
make test           # Run test suite (48 tests)
make test-debug     # Run with AddressSanitizer
make demo           # Comprehensive benchmark (10M points)
make main           # Basic API demo
make benchmark      # Performance benchmarks only
make simd-benchmark # SIMD-specific benchmarks
make clean          # Remove build artifacts
make help           # Show all options
```

## Precision Analysis

| Location | Encoding Error |
|----------|----------------|
| São Paulo, Brazil | 0.0021 m |
| New York, USA | 0.0058 m |
| Tokyo, Japan | 0.0049 m |
| North Pole | 0.0029 m |
| Prime Meridian | 0.0052 m |

**Sub-centimeter precision worldwide!**

## Roadmap

- [ ] Hilbert curve alternative encoding
- [ ] Streaming ingestion with chunked merges
- [ ] Memory-mapped persistent index
- [ ] Vincenty distance formula option
- [ ] Python/Rust/Go bindings

## Contributing

Contributions welcome! Areas of interest:
- Additional SIMD optimizations
- Alternative distance metrics
- Language bindings
- Documentation improvements ( I'm using AI to improve the documentation... )
