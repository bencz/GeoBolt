# Storage architecture and roadmap

## Implemented ownership boundary

RocksDB is GeoBolt's canonical single-server transactional engine. The integration is mandatory and explicit at build time; GeoBolt
does not silently choose a different storage engine. RocksDB owns its WAL, mutable and frozen memtables, SST creation, checksums,
snapshots, block cache, write-buffer accounting, flush, and compaction.

GeoBolt owns the domain-specific layers RocksDB does not provide:

- canonical object encoding with a generic nonzero 64-bit ID, commit sequence, packed 64-bit Morton location, and optional GeoDoc;
- exact geographic predicates, density-aware Morton covers, SIMD kernels, and result materialization;
- the derived immutable Morton segment set and its tombstone/last-write-wins visibility map;
- the commit coordinator that combines concurrent API calls into one canonical batch;
- the specialized active/frozen spatial memtable and its persistent asynchronous flush worker;
- typed, persisted GeoDoc equality/range indexes and their transactional cardinality catalog;
- wire protocol, admission control, authentication, driver, and server lifecycle.

There is no auxiliary native database WAL. One synchronized RocksDB `WriteBatch` contains object puts/deletes, append-only spatial
deltas, and catalog watermarks. This removes dual-write ambiguity: a batch is either absent or recoverable as one canonical decision.

## Column families

| Column family | Current role |
|---|---|
| `catalog` | format version, sequences, spatial watermark, ID allocators, typed index definitions, cardinalities, and histograms |
| `objects` | canonical object values keyed by big-endian `object_id` |
| `spatial_delta` | sequence-ordered durable changes used to catch the Morton index up after interruption |
| `idempotency` | durable request fingerprints, expiry, replay decisions, and compaction-filter reclamation |
| `secondary_index` | ordered typed values followed by generic object IDs; one shared architecture for every scalar index type |

All column families share one LRU block cache and one global write-buffer manager. `write_buffer_bytes` is therefore a database-wide
budget, not a value multiplied by the number of column families. Writes stall at the RocksDB boundary when that bounded budget is
exhausted instead of allowing unbounded RSS growth.

The current canonical database format is version 4. GeoBolt rejects another version instead of attempting an implicit in-place
migration; API and on-disk compatibility are not constraints in the current development milestone.

## Object and metadata model

The 16-byte hot spatial record remains `{ object_id, morton_code }`. Metadata never widens this record. A canonical object value uses an
explicit little-endian, versioned header followed by an optional canonical GeoDoc. GeoDoc supports null, boolean, signed/unsigned
64-bit integers, finite binary64, UTF-8 strings, arbitrary bytes, arrays, and nested objects. Object keys are canonicalized, duplicate
keys and invalid UTF-8 are rejected, and every document is checksummed and structurally validated before use.

Protocol v4 and both embedded/remote C APIs implement unconditional upsert, conditional insert, conditional full-object update, point
GET, delete, and atomic mixed batches. A full-object update replaces both location and metadata. IDs are generic and user supplied;
their meaning is never assumed to be a vehicle or another specific domain. Applications may alternatively request a server-generated
ID from a durable monotonic range; allocation and insertion share one canonical transaction.

Secondary definitions name a JSON Pointer and one strict physical type: boolean, signed/unsigned 64-bit integer, binary64, datetime,
string, or bytes. Missing/null fields are sparse; a present value of another type rejects the complete object transaction. Signed and
floating values use order-preserving bit transforms, while strings/bytes use a zero-safe terminator encoding, so RocksDB `Seek` serves
equality and inclusive/exclusive ranges without a metadata scan. The object update batch deletes the old key, inserts the new key, and
persists exact cardinality changes atomically with canonical object and spatial state. Numeric and common short string keys stay in
stack storage during preparation, avoiding per-index heap traffic on the usual write path.

Every nonempty built index also has an adaptive equi-depth histogram with at most 64 bins. The small index definition stores only the
bin count and boundary-byte count. Checksummed bin boundaries and build-time distinct counts live in one separate immutable catalog
record; current bin counts live in independent fixed-size records. A metadata update that remains in the same bin writes no histogram
record. Crossing a boundary changes only the affected counters in the canonical object batch, so selectivity stays transactionally
consistent without repeatedly serializing all boundaries. Static boundaries remain valid after deletes reduce cardinality to zero.

Index creation first commits a `BUILDING` definition and consumes a never-reused index ID. It scans canonical objects in bounded RocksDB
batches, then synchronously publishes the checksummed `READY` definition; restart detects and rebuilds any interrupted `BUILDING`
definition after clearing its private key prefix. Each acquisition is checked before another operation may overwrite its error status.
Failed construction rolls back the definition, secondary keys, and histogram records in one synchronized batch; the consumed ID is not
reused. Recovery also rolls back an interrupted, type-incompatible build instead of preventing access to valid canonical objects.
If durable rollback itself fails, the database rejects further writes in failed state until reopen reconciles the persisted definition.
Canonical corruption and storage failures still fail recovery. Final publication acquires its catalog lock and checks capacity before
persisting `READY`, so no subsequent fallible acquisition can leave the durable index unmaintained in memory.

Internal index maintenance and object serialization share validated immutable `GeoDocView` values. Input validation and persisted-object
decode retain structural and checksum verification; each individual secondary index no longer repeats the complete document validation.
Callers own mutation bytes and must leave them immutable until the synchronous write returns.

The current builder intentionally holds the single-server writer mutex for the scan,
so reads continue but writes pause. Replacing that pause with snapshot construction plus sequence-ordered delta catch-up remains the
next index-build optimization; the current behavior is explicit rather than exposing a partially maintained index.

## Recovery and derived state

`spatial_applied_sequence` may lag `next_sequence`, because the canonical transaction is made durable before the derived Morton
publication. Manifest format 5 stores a checksummed monotonic durable watermark and advances it atomically with segment publication.
Startup compares both watermarks before replay: a manifest-ahead state advances the catalog without duplicating the already-published
range; a catalog-ahead state rebuilds the derived index from canonical objects. Otherwise startup seeks directly to the common applied
sequence and replays contiguous deltas in bounded chunks.

The Morton manifest is a cache, not authority. If it is absent or structurally invalid, GeoBolt recreates it and scans canonical objects
in bounded chunks. Existing objects become new immutable Morton segments, then the canonical `next_sequence` is committed to the
manifest before the catalog acknowledgement. This ordering also prevents a second restart from rebuilding an already-complete cache.

## Spatial write path

Canonical commits append compact `{object_id, morton_code, operation}` entries to a preallocated spatial memtable instead of building
and synchronizing one `.gbi` file per API call. Each generation maintains an open-addressed last-write-wins ID table and 4,096 Morton
prefix buckets. Radius queries therefore suppress stale persisted copies with one hash lookup per candidate and inspect only overlay
buckets intersecting the geographic bounding box.

Group collection uses a monotonic deadline. It consumes compatible requests as they arrive and stops when the operation budget is full,
when a non-groupable request forms an ordering boundary, or when shutdown requires draining. A full group does not wait out the collection
window. Detached requests and the producer queue have separate tails while the collector releases its mutex to wait. Conditional and
idempotent requests still form individual transactions; extending their grouping requires per-request outcome/dependency handling.

At `spatial_memtable_max_operations`, active becomes immutable frozen and a replacement active generation is installed before the
persistent worker starts sorting or writing. The worker collapses repeated IDs, builds one Morton segment outside the visibility gate,
then publishes the segment and its final sequence watermark atomically. New commits and queries continue against active plus frozen
during the build. If both generations are full, writers apply bounded backpressure until frozen publication completes; memory does not
grow without limit. Checkpoint and orderly shutdown drain both generations and atomically remove already-applied spatial deltas from
RocksDB together with the catalog watermark.

## Remaining single-server work

The next storage features stay single-server; HA, replication, leader election, and distributed transactions remain explicitly out of
scope for this milestone.

- Nonblocking secondary-index builds using a RocksDB snapshot plus sequence-ordered delta catch-up; the persisted `BUILDING` recovery
  state, typed encoding, transactional maintenance, query path, and catalog are already implemented.
- Autonomous histogram reanalysis when accumulated distribution drift exceeds a bounded threshold. Persisted counts are already exact;
  refresh will replace static boundaries without making query correctness depend on an estimate.
- Bounded iterator/result streaming for very broad typed predicates. The planner already orders predicates by persisted selectivity,
  stops on an empty driver, filters every later scan against the progressively shrinking ID set, and issues canonical `MultiGet` reads
  in bounded reusable chunks when metadata is the cheaper driver.
- Partial GeoDoc `SET`/`UNSET` patches that preserve canonical encoding without forcing clients to resend an entire document.
- Equal-workload benchmarks for sustained ingestion, point reads, mixed spatial queries, recovery, write amplification, RSS, and
  p50/p95/p99/p99.9 latency. DuckDB remains a possible analytical/export integration, not the transactional hot store.
