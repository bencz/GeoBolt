# Tests

This directory contains correctness, persistence, concurrency, edge-case, and stress validation. Tests must protect a meaningful public
or internal invariant; performance measurements belong in `benchmarks/`.

Run the complete optimized suite with `make test`, the portable fallback with `make test-scalar`, or sanitizers with
`ASAN_OPTIONS=detect_leaks=0 make test-debug`. `TEST_FILTER=<substring>` selects a focused subset without changing the compiled code.

`test_geo_index.c` validates the low-level index, SIMD backends, immutable segments, and concurrent compaction. `test_database.c`
validates the single-server commit boundary. Its recovery case restores an older durable state watermark while retaining the newer WAL
and already-published segment, reproducing a real crash window and proving that replay remains idempotent.

`test_server.c` covers protocol corruption, authentication, persistent connections, concurrent clients, bounded work-queue rejection,
global payload-memory backpressure, durability across restart, and reactor shutdown. `test_daemon.c` launches the real `geoboltd`
executable, commits through the driver, delivers `SIGTERM`, and reopens the database to verify the process lifecycle and durable result.

## Uber-like actor soak

`make soak-uber` runs a five-minute concurrent workload by default. A fixed publisher pool schedules thousands of independent vehicle
actors by randomized deadlines; it does not create one operating-system thread per vehicle. Active vehicles move through random
micro-updates or deletes, inactive vehicles reinsert at new locations, publishers create immutable microbatch segments concurrently,
query workers continuously mix radius, ID, bounding-box, and kNN operations, and normal background compaction remains enabled.

An independent auditor periodically closes only the workload publication gate, then compares a global visible count and randomized
exact-radius ID results against the in-memory vehicle oracle. This protects delete/reinsert last-operation-wins semantics while ordinary
queries, publishers, and compaction remain concurrent between audits.

Arguments are positional:

```text
bin/soak_uber [seconds] [vehicles] [publishers] [query-threads] [batch]
              [min-delay-ms] [max-delay-ms] [max-segments] [report-seconds] [audit-seconds]
```

Use `make soak-uber-smoke` for a ten-second functional run and `make soak-uber-tsan` for a short ThreadSanitizer campaign. Override
`SOAK_ARGS`, `SOAK_SMOKE_ARGS`, or `SOAK_TSAN_ARGS` when exercising a different topology. The process reports event/query
throughput, latency histograms, CPU time, RSS and peak RSS, memory-growth slope, page faults, context switches, kernel-accounted I/O,
active segment pressure, and background-compaction work. Successful runs remove their isolated temporary directory; failed runs retain
it and print the path for inspection.
