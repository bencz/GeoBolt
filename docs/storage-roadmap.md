# Storage and metadata roadmap

The current engine persists generic 64-bit object IDs and packed 64-bit Morton locations. Protocol version 1 supports durable upsert and
delete batches plus aggregate radius counts. It does not yet provide arbitrary metadata, conditional insert/update, point lookup, or
result-set materialization through the network API.

## Required next storage milestone

Before adding those commands, location and metadata must share one commit sequence and one crash-recovery decision. The storage layer
needs a mutable memtable, frozen memtables, group commit, asynchronous flush, and backpressure based on dirty bytes. Building one
immutable segment and synchronizing state for every unitary socket request is correct but not the target architecture for sustained
high-rate ingestion.

Metadata must remain outside the 16-byte hot spatial record. The canonical object value can contain location, user bytes, version, and
flags, while spatial segments retain only `{ object_id, morton_code }` and derive all secondary-index state asynchronously.

## RocksDB decision boundary

RocksDB is a plausible primary embedded storage engine because it already supplies WAL, batch commit, mutable/immutable memtables, SST
flush, checksums, snapshots, and concurrent compaction. A clean integration would store canonical objects, spatial deltas, and system
watermarks in one atomic RocksDB write batch; GeoBolt would build its specialized immutable spatial segments from that durable delta
stream.

RocksDB must not be added as an auxiliary dual-write database beside an independent GeoBolt WAL. That would require a transaction
coordinator to prevent location and metadata from committing independently, while also duplicating block caches, background I/O, and
compaction scheduling.

The choice between a native engine and RocksDB requires an equal-workload benchmark covering sustained writes, point reads, mixed
spatial queries, recovery time, write amplification, resident memory, and p50/p95/p99/p99.9 latency. Storage-engine selection must be
explicit at build/deployment time; unlike an allocator, it must never change silently based on host library availability.

DuckDB remains a possible analytical/export integration. Its columnar OLAP execution is useful for aggregations and ad-hoc SQL, but it
is not the transactional hot store for small concurrent object mutations.
