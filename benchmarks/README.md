# Benchmarks

Benchmarks are capacity and regression tools, not correctness tests. Record the compiler, CPU/backend, allocator, thread count, dataset,
checksums, warmup policy, and whether the host was idle when publishing results.

- `make demo` runs the large end-to-end ingestion, SIMD, query, and parallel-query workload.
- `make benchmark-tombstones` isolates visibility filtering before and after physical compaction. Use `BENCHMARK_ARGS` to override its
  record count, query iterations, and removal cardinality. Each state receives one unmeasured warmup query before the reported rounds;
  both internal compaction time and complete public-call wall time are emitted. A second mutation-free rewrite reports the underlying
  merge-and-persistence cost separately from mutation visibility processing.
- `make benchmark-stream` measures pre-encoded streaming ingestion, radix chunk construction, intermediate merges, and final
  publication. Its arguments are record count, chunk capacity, sort threads, rounds, and an optional `uniform` or `hotspot`
  distribution. The hotspot shape places 80% of records under one 16-bit Morton prefix to exercise persisted density construction.
- `make benchmark-compaction` performs alternating serial and parallel durable compactions over the same immutable source segments.
  Its arguments are segment count, records per segment, parallel workers, rounds, and an optional `uniform` or `hotspot` distribution.
  Every output is reopened through the checksummed mmap loader and its global visible cardinality is verified before timing is accepted.
- `make benchmark-server` measures persistent TCP round trips, exact radius-count commands, and concurrent durable single-object
  upserts independently. Upserts use disjoint IDs per client so duplicate-ID fallback does not distort group-commit throughput. Its
  arguments are client
  connections, operations per client, and indexed objects. The default is a short regression profile; pass larger values explicitly for
  capacity runs. Per-client random seeds are fixed independently of memory addresses. Every connection warms up with PING and a radius
  query before timing; insertion setup is outside timing and measured upserts use new IDs. The benchmark reports a deterministic query
  checksum, verifies the final global visible object count, and prints aggregate operations/second, wall-time per aggregate operation,
  server backpressure, and peak retained payload bytes. `aggregate_us_per_op` is inverse throughput, not individual request latency.
  The write phase also warms each connection with 32 independent upserts before timing. An optional fourth argument, `metadata`, enables
  seven typed secondary indexes and two alternating GeoDoc versions with a 1 KiB unindexed payload. Each object receives two separate
  requests, exercising insertion followed by replacement of every indexed value; operations per client must be even. Final verification
  checks all seven typed indexes and prints their deterministic ID checksum. This is a bounded regression workload, not the requested
  million-actor sustained simulation. Example: `bin/benchmark_server 4 16000 10000 metadata`.
  Set `GEOBOLT_BENCHMARK_DIRECTORY` to an existing directory on the intended storage filesystem; the benchmark creates and cleans up
  its own unique child directory there. It defaults to `/tmp`, which can be **tmpfs**: such runs measure the synchronized API path on
  volatile memory, not disk durability cost. Record filesystem/device information separately when comparing disk-backed workloads.
- `make benchmark-query-planner` measures the spatial-driven conjunction case that historically materialized a broad typed predicate
  before searching a small radius. Its arguments are object count, query count, and canonical object-cache MiB; zero disables that
  cache for a controlled RocksDB fallback comparison. The workload reports elapsed time, secondary-index entries scanned, post-cache
  canonical lookups, selected Morton ranges, estimated spatial candidates, result cardinality, and checksum. Its default 100,000-object
  grid deliberately gives the metadata predicate 100% selectivity and keeps the geographic result sparse.
- `make pgo-generate` builds instrumentation; train it with representative binaries before `make pgo-merge` and `make pgo-use`.

The dedicated benchmarks print the compiler version, selected runtime SIMD backend, and linked allocator before workload metadata and
checksums. Whether the host was idle, CPU affinity, frequency policy, and storage topology remain operator-controlled and must be noted
with any published result.
