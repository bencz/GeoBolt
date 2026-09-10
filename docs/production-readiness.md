# Single-server production readiness

This is an implementation backlog with explicit acceptance criteria, not a production certification. HA, replication, and distributed
transactions remain outside the current milestone. Storage-specific work is detailed in [storage-roadmap.md](storage-roadmap.md).

## Existing foundation

The implementation already includes the RocksDB canonical transaction, spatial active/frozen memtables, asynchronous flush and Morton
compaction, group commit, generic IDs and GeoDoc, typed secondary indexes, persisted cardinalities/histograms, idempotency, a versioned
binary protocol, persistent workers, authentication, client serialization, and connection/request admission bounds.

The reactor now preserves complete requests and ordered responses after write half-close, delays connection reclamation until the
event batch and worker ownership have ended, and validates declared object-operation counts before allocation. Its completion queue
coalesces notifications and attempts socket writes directly. See [reactor validation](reactor-validation.md) for regression evidence.

## Release blockers and next acceptance criteria

The latest audit's iterator-OOM false publication, rejected-index reopen failure, and shared-client disconnect data race have targeted
fixes and regression coverage. GeoDoc index maintenance now reuses validated views instead of revalidating every document per index.
Group collection now uses monotonic, capacity-aware waiting; full groups and shutdown do not wait for an unnecessary collection window.
These changes do not close the resource, operational, or sustained-load criteria below.
The typed-metadata [TCP comparison](metadata-write-validation.md) shows a repeatable tmpfs gain but an unresolved time-dependent
Btrfs/NVMe slowdown in both variants. Durable disk throughput improvement is not established.

1. **Time and memory bounds across the network path.** Add monotonic authentication, incomplete-frame, idle, and blocked-response
   deadlines. Bound response buffers and retained query results globally, in addition to request bodies. Verify slow senders/readers,
   connection churn, saturation, bounded RSS, and continued progress of ordinary clients.
2. **Bounded query execution.** Add result/iterator streaming with explicit continuation and consistency semantics for broad predicates.
   Carry limits through the planner and protocol before materializing a complete result. Validate cancellation, capacity exhaustion,
   exact results, and allocation behavior of reusable query paths, including cache misses.
3. **Operational durability.** Provide a supported consistent backup/restore procedure and executable verification. Exercise process
   termination at commit/publication boundaries, storage-full and I/O failures, restart, and acknowledged-write preservation. Existing
   watermark and corruption tests are useful coverage but do not replace the complete failure campaign.
4. **Representative sustained client/server load.** Implement the requested million-actor workload with independent per-actor requests,
   jitter, metadata/location updates, deletes/reinserts, and typed-plus-spatial searches. Measure offered and completed load, backlog,
   individual p50/p95/p99/p99.9 latency, CPU, RSS, and write amplification over five to ten minutes. The existing `soak_uber` exercises
   the index/segment layer directly and must not be presented as this network/database workload.
5. **Deployable operation.** Expose health/readiness and actionable runtime metrics, make database resource budgets configurable from the
   daemon, and document a supervised installation and restore drill. The current token-authenticated protocol requires trusted
   transport; TLS termination and its deployment configuration need explicit validation before untrusted-network use.

## Following performance and feature work

- Remove the writer pause during secondary-index creation using snapshot build plus ordered catch-up and crash-safe publication.
- Refresh histogram boundaries autonomously under distribution drift, with bounded maintenance work.
- Implement canonical GeoDoc partial patches and transactional secondary-key maintenance.
- Profile group-commit waiting, durable sync, response materialization, and query planning under the sustained mixed workload before
  selecting the next algorithmic or SIMD change.
- Redesign the active spatial overlay for urban hotspots: adaptive Morton subdivision, contiguous candidate blocks, and efficient
  rejection of obsolete versions. Preserve exact geographic predicates and last-operation-wins semantics.
- Extend group commit to conditional and idempotent requests with per-request outcomes and dependency ordering for repeated IDs;
  monotonic, capacity-aware collection is already implemented.
- Shorten the global visibility exclusion around durable sync using a consistent canonical/spatial generation, not independent
  snapshots that can expose mismatched metadata and location.
- Complete reusable query allocation behavior, including cache-miss admission, GeoDoc validation scratch, and cache replacement.
- Use live per-cell overlay estimates and hybrid predicate evaluation; avoid scanning every broad secondary range after a selective driver.
- Remove duplicate response materialization, integrate filtered kNN into the database/client/server path, and expose maintenance budgets.
- Finish exception containment at every fallible C++ bridge boundary and audit ignored synchronization errors.
- Run the changed portable paths on ARM64 hardware and validate native ISA selection; x86 and forced scalar results do not certify ARM.

Every completed item requires implementation, meaningful failure/regression coverage, updated contracts and docs, and appropriately
controlled measurements. A throughput result alone does not close any of these acceptance criteria.
