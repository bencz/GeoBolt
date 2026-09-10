# Tests

This directory contains correctness, persistence, concurrency, edge-case, and stress validation. Tests must protect a meaningful public
or internal invariant; performance measurements belong in `benchmarks/`.

Run the complete optimized suite with `make test`, the portable fallback with `make test-scalar`, or sanitizers with
`make test-debug`. All three profiles cover index, database, storage, server, and daemon behavior;
`TEST_FILTER=<substring>` selects a focused low-level index subset without changing the compiled code.
Disable LeakSanitizer only after its documented ptrace incompatibility occurs; keep AddressSanitizer and UndefinedBehaviorSanitizer
enabled. Running outside a tracing sandbox can preserve leak checking as well.

`test_geo_index.c` validates the low-level index, SIMD backends, immutable segments, and concurrent compaction. `test_database.c`
validates the RocksDB-backed single-server commit boundary, both spatial-watermark crash windows, derived-index reconstruction,
conditional mutations, metadata, durable server-generated IDs, idempotency replay/conflict across reopen, active/frozen spatial
memtables, simultaneous inserts/deletes/readers during asynchronous flush, persistent group commit, and every typed secondary encoding.
The database publication scenario reopens with an empty process-local cache, then keeps three readers continuously querying while
twelve full 256-object delete/reinsert replacement rounds commit. It intentionally contains no scheduler yield or reduced-cardinality
escape hatch, so duplicate cold-cache fills, writer starvation, and partially visible batches remain observable regressions.
The secondary-index scenario protects build over existing objects, equality/range ordering, embedded-zero bytes, atomic key moves,
delete, strict type rejection, exact cardinality, persisted histogram ordering, progressive intersections, zero-cardinality reopen,
histogram-corruption rejection, direct GeoDoc filtering of spatial candidates, both cost-planned execution paths, batched canonical
reads, reopen, and drop.
`test_storage.c` separately verifies pinned `MultiGet` reuse, mixed found/not-found results, and snapshot isolation.

`test_server.c` covers protocol corruption, authentication, persistent connections, generated IDs, remote idempotency, typed-index
create/list/query/drop, conjunctive radius planning over protocol v4, concurrent clients, bounded work-queue rejection,
payload backpressure, durability across restart, and reactor shutdown. `test_daemon.c` launches the real `geoboltd` executable,
commits through the driver, delivers `SIGTERM`, and reopens the database to verify the process lifecycle and durable result.

The server disconnect scenario uses eight concurrent clients and 128 connections per client, mixing TCP resets with three pipelined
PING frames followed by write half-close. Payloads include zero, one, and 4,096 bytes. Every half-closed connection must return all
responses in order before EOF; all connection and request-payload accounting must return to zero before server destruction. Separate
cases truncate a mutation header/body, verify that only the complete mutation commits, and reject impossible declared operation counts
before allocating the decoded operation array. These scenarios run unchanged in optimized, scalar, ASan/UBSan, and TSan builds.
An 8 MiB PING response with a bounded receive window also remains unread while a second client completes a PING. The first client then
verifies the complete response and EOF, exercising partial-send resumption and reactor progress after write half-close.

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

## Failure-boundary regression coverage

On Linux, `make test`, `make test-debug`, and `make test-scalar` also run the corresponding `test_database_failures` executable.
It uses GNU/ELF linker interposition to inject an iterator allocation failure and to terminate a subprocess after a durable `BUILDING`
definition. The exact production index caller is compiled without LTO only in this executable so interposition cannot be optimized
away; no production fault hooks or substitute index algorithms are introduced. Tests verify failure status, absence of a false `READY`
index, recovery of a rejected/interrupted build, canonical object preservation, retry, and correct results after another reopen.
An injected rollback I/O failure additionally checks write quarantine and recovery of the remaining durable build definition.
The same executable observes the real coordinator's timed waits and asserts that a full one-operation group never enters its collection
timer. It checks scheduling behavior rather than a fragile wall-time threshold; the production coordinator caller is also compiled
without LTO in this executable, and all other implementations are unchanged.
Run this subset with `make test-failures`, `make test-failures-debug`, or `make test-failures-scalar`.

The shared-client server regression also exercises concurrent commands after server shutdown. Run it under ThreadSanitizer to cover
the descriptor-state transition as well as successful request serialization.
