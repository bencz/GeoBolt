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
- `make benchmark-server` measures persistent TCP round trips and exact radius-count commands independently. Its arguments are client
  connections, operations per client, and indexed objects. The default is a short regression profile; pass larger values explicitly for
  capacity runs. It reports aggregate operations/second, mean round-trip latency, server backpressure, and peak retained payload bytes.
- `make pgo-generate` builds instrumentation; train it with representative binaries before `make pgo-merge` and `make pgo-use`.

The dedicated benchmarks print the compiler version, selected runtime SIMD backend, and linked allocator before workload metadata and
checksums. Whether the host was idle, CPU affinity, frequency policy, and storage topology remain operator-controlled and must be noted
with any published result.
