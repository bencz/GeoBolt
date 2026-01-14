# GeoBolt Index

GeoBolt Index is a high-performance spatial indexing library focused on lightning-fast
geospatial lookups and analysis. It provides a production-ready data structure built on
64-bit Morton (Z-order) codes with precision near one meter, plus a rigorous C test
suite that benchmarks millions of operations per second.

## Why GeoBolt?

- **Sub-meter precision** – 32-bit per-axis normalization keeps errors under 1 cm in
  practice.
- **Cache-friendly layout** – 16-byte-aligned records and prefetch hints ensure high
  throughput even on large datasets (tested up to 5 million points).
- **Complete toolkit** – encoding/decoding, radius queries, bounding boxes, KNN search,
  and validation utilities all bundled in a single static library.
- **Battle-tested** – 39 deterministic tests cover correctness, numerical stability,
  edge cases (poles, antimeridian), stress, and performance regressions.

## Repository Layout

```
.
├── geo_index.h        # Public API
├── geo_index.c        # Implementation
├── main.c             # Demo program
├── test_geo_index.c   # Comprehensive test suite
└── Makefile           # Build + test orchestration
```

## Getting Started

### Requirements

- clang (or another C11 compiler)
- make
- macOS or Linux (timing utilities adapt automatically)

### Build Everything

```bash
make all
```

This produces:
- `build/libgeoindex.a` – static library
- `bin/main` – demo showcasing the API
- `bin/test_geo_index` – performance + regression tests

### Run the Demo

```bash
make main
```

Sample output:
```
=== GeoIndex Demo ===
...
FOUND id=10 lat=-23.5505200 lng=-46.6333090 dist=0.000 km
```

### Run the Test Suite

```bash
make test
```

- Executes all 39 tests (correctness + stress + benchmarks)
- Automatically fails if the binary exits with an error or hangs longer than 120s

For AddressSanitizer + UBSan builds:

```bash
make test-debug
```

## Key APIs

Include the header and link against the static library:

```c
#include "geo_index.h"

GeoIndex *index = geo_index_create(1000);
geo_index_add(index, 42, -23.55, -46.63);
geo_index_build(index);

GeoSearchStats stats;
GeoSearchResult *hits = geo_search_radius(index, -23.55, -46.63, 5.0, &stats);
...
geo_result_destroy(hits);
geo_index_destroy(index);
```

Highlights:

| Function | Description |
|----------|-------------|
| `geo_encode/geo_decode` | Convert latitude/longitude to/from 64-bit Morton codes |
| `geo_search_radius` | Deduplicated radius filtering with bounding-box pruning |
| `geo_search_knn` | Adaptive KNN expanding search radius up to a max distance |
| `geo_build_ranges` | Builds Morton ranges for efficient filtering |
| `geo_fast_distance_km` | Low-cost equirectangular approximation |
| `geo_get_time_ms` | Cross-platform high-resolution timer |

## Performance Snapshot

| Benchmark | Result |
|-----------|--------|
| Encode | 318 M ops/sec |
| Decode | 1001 M ops/sec |
| Binary searches | 8.17 M ops/sec |
| Radius search (1M pts, 100 queries, 100 km) | 0.035 ms avg |
| 5M record build | 375 ms total |

*(Measured on Apple Silicon, `-O3 -march=native -flto -ffast-math`)*

## Roadmap Ideas

1. SIMD-powered Hilbert curve encoder.
2. Streaming ingestion with chunked merges.
3. Optional persistent index serialization (memory-mapped file).

## Contributing

PRs and issues are welcome once the GitHub repository is published. Suggested areas:
- Additional heuristics for adaptive range subdivision.
- Alternative distance metrics (Vincenty, ECEF).
- Bindings for Rust/Python/Go.

## License

Choose the license that best matches your distribution goals (MIT, Apache-2.0, etc.).
Add the corresponding `LICENSE` file before publishing on GitHub.
