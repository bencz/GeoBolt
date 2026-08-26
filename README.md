# GeoBolt

GeoBolt is a compact C11 geospatial index optimized for build-once, query-many workloads. Each coordinate pair is normalized to two
32-bit integers and interleaved into one 64-bit Morton code. Together with the 64-bit application ID, every indexed record occupies
exactly 16 bytes.

## Highlights

- 64-bit Morton encoding with sub-centimeter quantization error.
- Single-pass, budget-aware Morton-prefix coverage with cardinality-aware depth selection instead of fixed-depth retries or one large
  false-positive interval.
- Sparse six-level density hierarchy with persisted Morton bounds/counts and a cost-based cover planner for spatial hotspots.
- Exact spherical radius filtering and antimeridian-aware bounding boxes, including a fused Morton-to-bbox SIMD kernel.
- Interior-cell fast paths and allocation-free count-only radius/bounding-box queries.
- Six-pass 11-bit LSD radix build for large indexes, including a lock-free parallel scatter path and a persistent streaming-build pool.
- Best-first hierarchical kNN traversal with exact spherical cell lower bounds and a bounded result heap.
- Reusable result buffers for allocation-sensitive query loops.
- ID-only 8-byte results and allocation-free caller-owned fill operations.
- Persistent SoA batch executor with dynamic chunk scheduling, two-pass contiguous materialization, immutable segments, replicas, and
  shards.
- Native read-only persistence with direct `mmap` queries on Unix-like systems.
- External-memory streaming builds with bounded chunks, sorted runs, k-way merge, and atomic publication.
- Versioned immutable segment sets with checksummed manifests, aggregate queries, and direct Morton compaction.
- Linearizable batch deletion and ID reinsertion through a durable mutation WAL, with tombstone reclamation during compaction.
- Cache-line-sharded QSBR reader counters for segment publication and compaction; query readers do not acquire an rwlock.
- Selective runtime AVX-512/AVX2+FMA dispatch on x86-64, NEON on ARM64, and a portable scalar backend.
- SIMD batch validation and direct coordinate/ID encoding into the final 16-byte record layout.
- Packed 64-bit match masks and branch-efficient result compaction for sparse and dense query results.
- Sanitizer-friendly allocation and vector-tail handling; SIMD loads do not require aligned user buffers.

## Build

Requirements are a C11 compiler, `make`, and a Unix-like environment for memory mapping.

```bash
make all
make test
make test-debug
make test-scalar
make sample
make demo
make benchmark-tombstones
```

Compiler and linker settings remain overridable:

```bash
make clean
make CC=clang CFLAGS_OPT="-O3 -flto"
```

The production defaults deliberately do not enable `-ffast-math`. Query validation depends on defined IEEE behavior for invalid
coordinates and on controlled numerical error at poles, the antimeridian, and predicate boundaries.

Optimized executables automatically link jemalloc when it is available through `pkg-config` or a compiler link probe. The fallback is
the system allocator, so jemalloc is not a required dependency. Sanitizer builds always retain the sanitizer allocator/interceptors.

```bash
make allocator-info
make USE_JEMALLOC=0 all    # Force the system allocator for an A/B comparison.
make USE_JEMALLOC=1 all    # Force -ljemalloc when detection metadata is unavailable.
```

Changing `USE_JEMALLOC` forces optimized executables to relink, preventing a stale binary from silently retaining the previous
allocator selection.

Source organization is intentional: correctness and concurrency validation lives in `tests/`, capacity workloads live in
`benchmarks/`, and complete public-API examples live in `samples/`. Production modules remain at the repository root until a separate
source/include hierarchy provides a measurable build or maintenance benefit.

## Core API

```c
#include "geo_index.h"

GeoIndex *index = geo_index_create(1000000);

geo_index_add(index, id, latitude, longitude);
geo_index_build(index);

GeoSearchStats stats;
GeoSearchResult *nearby = geo_search_radius(index, latitude, longitude, radius_km, &stats);
GeoSearchResult *nearest = geo_search_knn(index, latitude, longitude, k, max_radius_km, NULL);

geo_result_destroy(nearby);
geo_result_destroy(nearest);
geo_index_destroy(index);
```

For large in-memory builds, `geo_index_build_parallel(index, 0)` selects up to eight workers automatically. Passing an explicit worker
count overrides the automatic choice. Histograms are private to each worker and each scatter writes to precomputed disjoint ranges, so
the radix passes require no locks or atomics.

For repeated queries, retain one result object and avoid repeated heap allocation:

```c
GeoSearchResult *result = geo_result_create(256);

for (size_t i = 0; i < query_count; ++i) {
    if (!geo_search_radius_reuse(index,
                                 query_latitudes[i],
                                 query_longitudes[i],
                                 radius_km,
                                 result,
                                 NULL)) {
        break;
    }

    consume(result->results, result->count);
}

geo_result_destroy(result);
```

When only cardinality is needed, use the count-only APIs. They execute the exact predicate without allocating or copying records.
Bounding-box cells that are completely inside the query are counted directly from their Morton ranges; only boundary cells are decoded
and filtered.

```c
size_t nearby_count;

if (geo_search_radius_count(index, latitude, longitude, radius_km, &nearby_count, NULL)) {
    consume_count(nearby_count);
}
```

When Morton codes are not needed downstream, ID-only results cut output storage from 16 to 8 bytes per match. The caller-owned variant
is the allocation-free fill primitive used by the batch executor.

```c
GeoIdResult *ids = geo_id_result_create(256);

geo_search_radius_ids_reuse(index, latitude, longitude, radius_km, ids, NULL);
consume_ids(ids->ids, ids->count);

geo_id_result_destroy(ids);
```

Large query batches use SoA inputs and a persistent worker pool. The executor counts first, computes a prefix sum, and fills one
contiguous allocation in the second pass, so worker threads never call `realloc` while materializing results.

```c
GeoQueryExecutorConfig config = {
    .thread_count = 8,
    .scheduling_chunk = 16,
    .pin_workers = true,
};
GeoQueryExecutor *executor = geo_query_executor_create(index, &config);
GeoBatchIdResult *batch = geo_batch_id_result_create(query_count, 0);

geo_query_executor_search_radius_ids(executor,
                                     query_latitudes,
                                     query_longitudes,
                                     query_radii_km,
                                     query_count,
                                     batch);

for (size_t query = 0; query < batch->query_count; ++query) {
    consume_ids(batch->ids + batch->offsets[query], batch->offsets[query + 1] - batch->offsets[query]);
}

geo_batch_id_result_destroy(batch);
geo_query_executor_destroy(executor);
```

`geo_query_executor_create_replicated()` assigns complete core-local index replicas to workers. For disjoint partitions,
`geo_query_executor_create_sharded()` maintains one cache-line-isolated work queue per shard, computes partial counts independently,
and concatenates each query's results in shard order. `worker_cpus` may provide an explicit topology-aware CPU mapping; otherwise pinned
workers use sequential logical CPU IDs.

`geo_query_executor_create_segment_set()` runs the same scheduler over an immutable segment collection. One QSBR snapshot remains active
across count, prefix-sum, and fill, so removals, reinsertions, or compaction cannot make the two materialization passes disagree. Worker
kernels consume that snapshot directly and never re-enter the reader gate.

## Persistent indexes

Only a built, sorted index can be saved. A mapped index is read-only and uses the same search APIs without copying its record array.

```c
geo_index_save(index, "places.geobolt");

GeoIndex *mapped = geo_index_open_mmap("places.geobolt");

if (mapped && geo_index_is_read_only(mapped)) {
    GeoSearchResult *result = geo_search_radius(mapped, latitude, longitude, 25.0, NULL);
    geo_result_destroy(result);
}

geo_index_destroy(mapped);
```

Format version 5 stores the records, fixed-width prefix offsets, and sparse density hierarchy as independently checksummed sections.
Its 64-bit polynomial checksum folds four words per iteration, exposing independent integer-multiply chains to both x86-64 and ARM64
while remaining identical between streamed writes and one-shot mmap validation.
Each hierarchy node contains an aligned Morton prefix and exact begin/end offsets. Only children of cells exceeding four times the
baseline density are stored, preventing the exponential memory growth of a dense quadtree. Opening a mapped index validates the exact
file size, section checksums, prefix bounds, density metadata, and Morton ordering before publishing the view. Older formats are rejected;
the project deliberately does not carry an unmeasured compatibility path.

The representation remains native-endian and records its format version and `GeoRecord` size. It is intended for fast reopening on
compatible systems, not yet as a cross-endian interchange format.

## External-memory streaming builds

The streaming builder keeps at most one mutable chunk in memory. Full chunks are SIMD-encoded, radix-sorted, and written as temporary
runs. A 32-way hierarchical merge compacts runs as they accumulate, so descriptor and final-merge heap usage grow logarithmically rather
than linearly with the number of chunks. The final merge produces a directly mappable index, which is published atomically only after its
data and metadata have been flushed.

Parallel builders create their radix workers, scratch array, and cache-line-aligned per-worker histograms once. Every full chunk reuses
those resources; the final partial chunk automatically reduces its active worker count or uses the serial sorter when parallel startup
would dominate.

When the complete input fits in the configured chunk, publication writes the sorted records directly from memory into the persisted
format. No temporary run is created or read back. Multi-run builds use 64 KiB input blocks per run, avoid redundant seek/flush cycles,
and compute section checksums only for the final persisted output rather than for disposable intermediate merges.

Merge cursors read blocks instead of individual records, and replacement of the heap minimum requires only one sift operation per
record. On Linux, temporary runs and final output are preallocated and marked for sequential access to reduce fragmentation and improve
kernel readahead on datasets larger than memory.

The streaming density writer defers child-hierarchy accounting while a Morton root remains below the persisted density threshold.
Uniform roots therefore pay only the root transition and a short Morton buffer instead of updating six child levels that would later be
discarded. Once a root crosses the threshold, its buffered Morton values are replayed exactly once and subsequent records use the normal
hot-root path. Deferred storage is capped at 131,072 records; larger thresholds retain the constant-memory eager implementation.

```c
GeoStreamBuilder *builder = geo_stream_builder_create("planet.geobolt", "/fast-temporary-storage", 1000000);

while (read_next_batch(ids, latitudes, longitudes, &count)) {
    if (!geo_stream_builder_add_batch(builder, ids, latitudes, longitudes, count)) {
        break;
    }
}

GeoStreamBuildStats build_stats;
bool published = geo_stream_builder_finish(builder, &build_stats);

geo_stream_builder_destroy(builder);
```

CPU-heavy builds can use `geo_stream_builder_create_parallel()` with an explicit worker count, or zero for automatic selection.

If an upstream stage already produces compact Morton keys, `geo_stream_builder_add_records()` accepts contiguous 16-byte
`GeoRecord { id, z }` values directly and bypasses coordinate validation and encoding.

The chunk capacity controls the dominant in-memory allocation. `GeoStreamBuildStats` reports initial runs, intermediate merge count,
peak open runs, and timings for chunk construction, hierarchical compaction, and final merge. Temporary runs should be placed on local
NVMe for production-scale builds; the final atomic rename occurs inside the destination directory. `runs_created` is zero when the
single in-memory chunk can be published directly.

## Continuous ingestion with immutable segments

A segment set publishes independently built immutable indexes through a small checksummed manifest. Adding a segment opens and validates
it first, then atomically publishes the next manifest generation. Readers reopening the manifest therefore see either the previous
complete generation or the next complete generation.

```c
GeoSegmentSet *segments = geo_segment_set_create("active.manifest");

geo_segment_set_enable_background_compaction(segments, "/var/lib/geobolt/segments", 8);

geo_segment_set_add_file(segments, "segment-0001.geobolt");
geo_segment_set_add_file(segments, "segment-0002.geobolt");

// A live update is an upsert: older physical copies of every ID in the new segment become invisible atomically.
geo_segment_set_upsert_file(segments, "vehicle-updates-0003.geobolt");

uint64_t deleted_ids[] = { 1042, 1097, 1201 };
geo_segment_set_remove_ids(segments, deleted_ids, 3);

size_t nearby_count;
geo_segment_set_search_radius_count(segments, latitude, longitude, radius_km, &nearby_count, NULL);

GeoSegmentCompactionStats compact_stats;
geo_segment_set_compact(segments, "segment-compacted.geobolt", &compact_stats);

// Explicit worker control is available for capacity planning; the ordinary API chooses automatically.
geo_segment_set_compact_with_workers(segments, "segment-compacted-manual.geobolt", 8, &compact_stats);

geo_segment_set_destroy(segments);
```

Deletion is logical on publication and physical on compaction. A mutation is first appended and synced to the manifest-adjacent
WAL; the new manifest generation then commits the exact WAL prefix. Two crash-safe WAL slots allow checkpoints to rewrite only active
ID states and publish the replacement through the manifest before reclaiming the old slot. Checkpointing is automatic when historical
operations exceed four times the active mutation cardinality, and remains explicitly callable for operational control. A crash can
therefore leave only an ignored WAL tail or obsolete slot, never a manifest that references an incomplete deletion. Queries consult a
sparse ID exception table only when removals or reinsertions exist. Untouched IDs retain the original 16-byte record and the
mutation-free query path has no hash lookup. Each in-memory exception keeps a naturally aligned 16-byte payload: occupancy, state, and
generation are encoded in one 64-bit metadata word, placing four payloads in a 64-byte cache line without packed or unaligned
structures. A separate one-byte-per-slot control array stores nonzero hash fingerprints. Lookups inspect controls eight at a time with
portable SWAR operations and touch the 16-byte payload only when a fingerprint matches, substantially reducing random cache traffic for
the common negative lookup.

`geo_segment_set_upsert_file()` applies last-operation-wins semantics to every ID in its segment, including an ID that is already live;
`geo_segment_set_add_file()` retains append-only semantics for disjoint IDs. The exception stores the segment generation, so an older
physical copy cannot reappear. Compaction writes only the visible generation, removes tombstoned records, carries LIVE generations
across partial merges, and collapses exceptions back into the common lookup-free state only after an integral generation-stable merge.

Radius, bounding-box, count-only, and exact kNN queries operate across every active segment. The `_reuse` result paths append directly
into the caller buffer and do not allocate when its capacity is sufficient. Large compactions sample the already sorted inputs to form
balanced contiguous Morton ranges, then one visibility pass computes exact per-range output offsets and the persisted high-16 prefix
histogram. Independent workers merge those ranges into disjoint `pwrite` offsets. Each partition computes a zero-seed polynomial
checksum that is combined in output order, so checksum generation requires no serial reread. Density metadata is constructed from one
sequential page-cache pass after the parallel write. Small merges retain the lower-overhead single-pass path. Neither path decodes,
re-encodes, nor radix-sorts mmap records. After the compacted file is durable, one manifest rename makes it active. Superseded segment
files remain available for an explicit retention policy.

Compaction specializes visibility work by mutation cardinality. With no mutations, the physical record count is already exact and the
precount pass is skipped. Up to eight active exceptions use the cache-resident small-entry path directly. Larger exception sets build a
temporary contiguous bitmap using one bit per physical record; each 64-bit word is assembled and written once, then reused by the merge
instead of repeating a hash-table lookup for every record. The bitmap is released before publication and never becomes persisted state.

Radius and bounding-box queries prepare their immutable Morton range plan once and reuse it across every segment. This avoids repeating
cover construction and scalar trigonometric setup for each immutable file. The cost planner evaluates aggregate candidate counts across
all active segments before selecting a shared depth. Read sections use only the QSBR gate and active-reader atomics.

The built-in scheduler starts one background worker when the active segment count exceeds its configured limit. Size pressure selects
up to four of the smallest segments from the same approximate size tier, merges only those inputs, and preserves larger segments
untouched. This avoids repeatedly rewriting the entire historical index as new level-0 segments arrive. If the active count is still
above the limit, the worker performs additional tiered cycles before it exits.

Mutation pressure uses a different policy because retaining a large exception table taxes every query. When at least 4,096 active
mutations cover one eighth of the physical records, the worker performs one integral rewrite, including the single-segment case, and
returns the set to its mutation-free query path. A second pass is allowed if a concurrent publisher introduced enough mutations during
the first rewrite; the bounded retry prevents sustained writes from causing unbounded background write amplification. Every output
receives a unique immutable name inside the configured directory.

`geo_segment_set_configure_background_compaction()` accepts a `GeoSegmentCompactionPolicy` when storage bandwidth or query SLAs require
a different crossover. Its ratio is evaluated with overflow-safe integer arithmetic; no floating-point calculation enters publication
or scheduling decisions. The simpler enable function installs the documented 4,096-entry, 1/8-ratio, two-pass defaults.

Use `geo_segment_set_wait_for_background_compaction()` during controlled shutdown or when aggregate worker statistics are required.
The long merge phase holds a QSBR read epoch while publishers are serialized separately. Readers are gated only for the short
manifest/snapshot handoff that replaces the active mmap array and reclaims the old one.

## SIMD API

Batch functions accept unaligned arrays and handle any element count, including scalar tails:

```c
geo_simd_encode_batch(latitudes, longitudes, morton_codes, count);
geo_simd_decode_batch(morton_codes, latitudes, longitudes, count);
geo_simd_haversine_batch(center_latitude, center_longitude, latitudes, longitudes, distances, count);
geo_simd_filter_radius(latitudes, longitudes, count, center_latitude, center_longitude, radius_km, mask);
```

On x86-64, the translation unit is safe to run on machines without AVX2: one initialization selects an immutable backend vtable.
AVX-512F+DQ+VL handles decode and integer bounding-box kernels, AVX2+FMA handles the trigonometric pipeline, and shared scalar kernels
provide the final fallback. Dispatch occurs once per batch, never inside an element loop, and capability detection is not repeated by
each public operation.
Keeping the radius pipeline at one vector width avoids the AVX-512-to-AVX2 frequency and transition penalty observed on the development
host. Radius scans extract and decode Morton members in narrow vector batches, then evaluate the prepared Haversine predicate over the
L1-resident coordinate arrays. Define `GEO_SIMD_FORCE_SCALAR` to compile the portable backend explicitly.

## Current validation

The current x86-64 validation run completed 69/69 optimized tests with Clang 22 and GCC, plus 69/69
AddressSanitizer/UndefinedBehaviorSanitizer tests. The forced scalar backend completed 58/58 applicable tests. ThreadSanitizer completed
the affected segment paths, including automatic tombstone reclamation, simultaneous insert/remove, writes during an active merge, and
a writer waiting behind a two-pass query snapshot. Clang's static analyzer completed the changed translation units
without diagnostics. The ARM64 source shares the same scalar tail semantics, but the changes in this revision have not yet been
validated with an AArch64 sysroot or ARM64 hardware.

The following x86-64 ranges came from four consecutive runs on an otherwise idle development host. Deployment capacity measurements
should additionally pin cores and control CPU frequency:

| Operation | Observed speedup over scalar |
|---|---:|
| Morton encode | 2.01–2.14x |
| Morton decode | 6.30–7.33x |
| Haversine batch | 2.32–2.41x |
| Radius predicate | 13.18–14.78x |

In a separate idle-host audit with ten million input points, warm-buffer direct AoS encoding completed in 24.23–24.80 ms
(approximately 403–413 million records/second). The complete public ingestion path, including validation and first-touch of newly
allocated record pages, completed in 55.18–55.93 ms. Parallel build completed in 118.84–121.53 ms.

Radius-query planning now chooses its Morton cover depth from the index cardinality, targeting approximately 64 globally averaged
records per cell. In the 100,000-record, 10,000-query stress test, replacing the former fixed depth of 24 reduced total exact-query
latency from approximately 77–84 ms to 8–14 ms across repeated runs. The fastest idle runs completed in 7.74–9.39 ms, exceeding one
million single-threaded queries per second. All result checksums remained identical across every tested depth.

The persisted sparse hierarchy measures actual local density at up to six additional quadtree levels. Uniform indexes record no child
nodes and take a predictable branch-only fast path. For skewed indexes, the cost planner evaluates the real candidate count and range
overhead of every available depth, caches the selected binary-search bounds, and skips the exact predicate entirely when a global query
can copy or count the complete index. In a synthetic 500,000-record hotspot, 2,000 exact 500 m materialized
queries improved from 3.075–3.154 seconds to 0.473–0.477 seconds, a 6.45–6.67x speedup, with identical checksums. The center query
scanned 96,385 candidates instead of all 500,000 records. On a separate globally uniform 100,000-record workload, the median time for
10,000 queries changed from 7.880 ms to 7.912 ms, a difference of approximately 0.4% and within run-to-run noise. The hotspot-specific
speedup should not be extrapolated to globally uniform data.

In a separate idle-host batch audit, 100,000 exact 50 km ID queries over one million records took 149.07–155.25 ms sequentially. The
persistent eight-worker executor completed both count and fill passes in 41.02–44.04 ms, a 3.47–3.66x speedup, and produced the same
4,640,223-ID checksum. This workload includes the cost of one contiguous result allocation but excludes worker creation because the pool
is intentionally retained across batches.

An alternating allocator comparison found no material single-threaded query benefit from jemalloc in the already-amortized result-buffer
workload. Median time for 10,000 allocating queries changed from 8.265 ms with the system allocator to 8.360 ms with jemalloc; the reused
buffer path changed from 7.790 ms to 7.680 ms. Automatic jemalloc selection remains useful for production concurrency and fragmentation,
while `USE_JEMALLOC=0` keeps controlled comparisons straightforward.

For eight immutable segments containing two million records in total, preparing one cardinality-aware query plan and sharing it across
all segments completed in approximately 10.25–10.84 ms, compared with 16.0–16.35 ms when rebuilding the plan per segment. This is a
1.48–1.60x remaining benefit from plan sharing; cardinality-aware planning also reduced the absolute shared-plan time from the former
fixed-depth result of approximately 32–34 ms. All result counts were identical.

The ten-million-record demo is intentionally more demanding: 10,000 exact 50 km queries materialized approximately 4.7 million result
records. The latest idle-host run completed in 77.87 ms with one query thread and 16.47 ms with eight threads, or approximately 128,000
and 607,000 queries per second respectively, with an identical 4,713,975-result checksum. Dense-query latency is increasingly dominated
by exact filtering and copying the returned records rather than Morton-plan construction.

For five million random Morton keys, the parallel radix build measured approximately 126 ms with one worker and 62–64 ms with 4–8
workers. At that point additional workers were limited by memory bandwidth. A five-million-record external build with four sorting
workers completed in approximately 158–170 ms from doubles, or approximately 126 ms when ingesting pre-encoded `GeoRecord` values, on
the local temporary filesystem. These storage figures are host-specific and are not substitutes for NVMe production benchmarks.

For repeated 600,000-record radix chunks, an idle-host alternating comparison used Clang 22, jemalloc, four workers, 11 measured chunks
per path, and identical `id`/Morton checksums. Reusing the worker pool, scratch array, and histograms completed in 49.00–49.68 ms across
three runs; reconstructing those resources for every chunk took 60.49–61.83 ms. The measured sorting-stage speedup was 1.22–1.26x.

For a three-million-record pre-encoded streaming build that fits in one chunk, an alternating ten-round comparison measured the former
eager density writer at a 113.60 ms median end-to-end and 59.13 ms median final-emission time. Deferred cold-root accounting reduced the
medians to 73.94 ms and 20.30 ms respectively, a 1.54x end-to-end and 2.91x emission-stage speedup. Every round used the same input
checksum and passed mmap checksum, ordering, density, and full-count validation. With 30 external runs of 100,000 records, median total
time fell from 179.54 ms to 138.78 ms and median final merge from 119.97 ms to 78.63 ms, a 1.29x and 1.53x speedup respectively. An
80%-hotspot distribution retained the complete density hierarchy and completed in 82.34–86.73 ms for direct publication or
153.18–157.07 ms with 30 runs.

A power-of-two winner tree was also evaluated for the 30-way merge. Although it halves the theoretical comparisons per level, its
index indirection increased the final-merge median from 82.89 ms to 85.61 ms across seven alternating pairs. The prototype was removed;
the compact binary heap and its single sift-after-replacement remain faster for the current 32-run fan-in.

Three further heap variants were rejected after alternating comparisons: a true loser tree raised typical final-merge time from
82–87 ms to 94–97 ms; separating a 16-byte hot Morton/run key from the cold ID raised the median from 83.019 ms to 84.460 ms; and a
Floyd bottom-up replacement raised typical time to 94–97 ms. The retained binary heap has fewer dependent loads and stops its sift as
soon as the replacement reaches its actual level.

Format version 5 replaced the serial per-word checksum dependency with four-word polynomial folding. On the same three-million-record,
30-run uniform workload, 11 alternating pairs with Clang 22.1.8, jemalloc, four sorting workers, identical input checksums, and mmap
validation reduced median final merge from 85.075 ms to 81.626 ms and median end-to-end time from 147.343 ms to 144.108 ms. Across seven
alternating 80%-hotspot pairs, final-merge median fell from 96.752 ms to 93.754 ms. These reported-idle-host improvements are 4.1%, 2.2%,
and 3.1% respectively; storage and CPU-frequency variance still require production-local capacity measurements.

In a local 1.6-million-record experiment, compacting 16 mmap segments took approximately 33 ms. A repeated 50 km count query improved
from approximately 395–433 microseconds across all 16 segments to approximately 25 microseconds after compaction, a 16–17x latency
reduction. This is why production deployments should bound the number of active level-0 segments.

The partitioned compaction benchmark alternated serial and eight-worker runs over the same eight one-million-record mmap segments,
using Clang 22.1.8, jemalloc, 32 Morton partitions, validated record cardinality, and a checksummed reopen after every output. Across five
idle-host uniform pairs, serial compaction averaged 174.574 ms and parallel compaction 97.748 ms, a 1.786x speedup. With 80% of records
inside one São Paulo hotspot, adaptive sampled boundaries retained a 1.572x speedup: 199.950 ms versus 127.202 ms. These measurements
include visibility counting, prefix and density metadata construction, durable output publication, and manifest replacement; CPU
frequency and affinity were not pinned.

For one million records with 100,000 active tombstones, ten warmed global count queries perform ten million visibility checks and
completed in 52.67–58.45 ms across three consecutive runs on the reported-idle x86-64 host with Clang 22.1.8, jemalloc, one query
thread, and identical nine-million result checksums. Before the fingerprint control array, consecutive runs of the same workload took
102.43–120.05 ms without an explicit in-process warmup; this earlier range is useful directional evidence, but was not an alternating,
frequency-pinned comparison. Rewriting the same single segment and physically reclaiming 100,000 records took 27.83–28.19 ms end to end
on the local temporary filesystem, compared with the previously measured 34.16–34.84 ms before visibility-bitmap reuse. A mutation-free
rewrite of the resulting 900,000 records took 17.55–18.23 ms. For this deliberately global query shape, reclamation amortizes after
roughly five subsequent queries; narrower workloads and production storage require their own policy crossover.

The bounding-box pipeline now compares the compacted 32-bit latitude and longitude fields directly against exact normalized inclusive
thresholds. On this host, its AVX2 kernel processed approximately 0.95–0.99 billion points per second, about 2.3x the former floating-point
decode/filter path. AVX-512 reached approximately 1.35–1.39 billion points per second with byte masks and 1.82–1.86 billion points per
second for count-only filtering.

Materialized queries use packed 64-bit match masks. Full words become contiguous copies; partial words are traversed with trailing-zero
count and `bits &= bits - 1`, avoiding one unpredictable branch per candidate. The measured result-compaction improvement was
approximately 1.43–1.59x in the earlier controlled experiment. A fully fused Morton-decode/Haversine kernel was also tested with
pre-touched buffers and alternating execution order. It reached only 0.858–0.985x the throughput of the two-stage executor because the
combined integer decode and trigonometric polynomials increase register pressure, so the faster two-stage pipeline remains active.

An on-query circle-aware cell classifier was also evaluated as a way to copy cells proven to be inside the radius without exact
per-record filtering. Its additional cell refinement and bound calculations increased the 10,000-query reuse time from 7.82 ms to
11.02 ms and slightly regressed the 1,000 km query. The prototype was removed. A future interior-cell path should use bounds and density
metadata persisted at build time, avoiding repeated geometric work in the query hot path.

The global SIMD distance comparison includes poles, the antimeridian, and an antipodal point. Its observed maximum absolute difference
from scalar Haversine was approximately `6e-9 km` on the current x86-64 host.

## Threading model

Concurrent reads are safe after `geo_index_build()` as long as no thread mutates or destroys that index. Segment-set readers increment
one of 128 cache-line-isolated counters selected from thread-local state. Publishers close the reader gate, wait for every reader shard
to quiesce, publish new arrays, reclaim the old snapshot, and reopen the gate. Queries may run concurrently with each other and with
manual or scheduled compaction; segment publishers are internally serialized. Segment-set destruction waits for its scheduled
compaction worker, but all external callers must still be quiescent. Each query returns independent storage, while the `_reuse`
variants require one result buffer per concurrently executing caller.

## Data layout policy

`GeoRecord` is naturally aligned, exactly 16 bytes, and has no internal padding: a 64-bit application ID followed by the packed 64-bit
Morton coordinate. Compile-time size, alignment, and field-offset assertions protect both the hot array layout and mmap ABI. Persisted
file and manifest headers have equivalent size/offset assertions. Density metadata is index-level and sparse: it adds no per-record hot
field and keeps `GeoRecord` at 16 bytes. Cache-line-aligned batch queues prevent atomic counters for different NUMA shards from sharing
a coherence line.

The code intentionally does not use `#pragma pack` or packed attributes for traversed records. Saving a few header bytes would not affect
throughput, while unaligned 64-bit/vector access can make ARM64 slower and complicate direct mmap access. Unions are likewise not used
for type punning; bit operations and intrinsics express the representation without violating strict aliasing. `volatile` is reserved for
benchmark sinks—thread synchronization uses C11 atomics and pthread primitives.

## Code style

The repository uses four-space indentation, braces for every control-flow body, logical blank-line separation, and a 140-column limit.
The checked-in `.clang-format` documents the intended layout.

## Remaining directions

- Portable, checksummed, cross-endian persistence format.
- Multi-level compaction with record-count/byte-size policies and retention-aware garbage collection.
- Optional Hilbert ordering for workloads where its improved locality offsets the higher encoding cost.
- ARM SVE2 and additional fused kernels where architecture-specific benchmarks justify them.
- Automatic Linux NUMA topology discovery and memory-policy binding on top of the explicit CPU/replica mapping API.
- Per-worker prepared-plan caches for repeated query centers and radius classes.
- Language bindings and a stable opaque-ABI layer.
