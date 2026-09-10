# Single-server architecture

## Process boundary

`geoboltd` is one foreground process that owns exactly one database root. The current milestone is deliberately single-server: there is
no replication, leader election, distributed transaction, or cross-node failover protocol. A supervisor such as systemd or Kubernetes
owns process restart.

The process contains four long-lived execution domains:

1. One `epoll` reactor owns listening, accepting, connection parsing, response writes, and connection lifetime.
2. A bounded persistent job queue executes commands that can block or consume CPU. No worker is created per connection or request.
3. A persistent commit coordinator coalesces concurrent mutations, commits them through RocksDB, and publishes immutable spatial
   segments to concurrent readers.
4. The segment subsystem owns its background compaction worker and reusable parallel compaction pool.

A connection has at most one in-flight request. This keeps request ordering deterministic and bounds connection state without requiring a
per-connection response reorder buffer. Applications obtain parallelism through multiple persistent connections; a single C-driver
handle may be shared safely but intentionally serializes its complete request/response transactions.

The reactor consumes complete queued frames after a peer's write half-close and sends their ordered responses before closing on EOF.
An incomplete final header/body is discarded without dispatch. While a worker owns a connection, only error/hangup events remain
enabled; the reactor does not inspect worker-owned response fields. Completed jobs are published through the completion mutex, with an
eventfd notification on the empty-to-nonempty queue transition. Notifications retry interruption and tolerate an already readable,
saturated eventfd. The reactor attempts response writes directly on completion and enables `EPOLLOUT` only when the socket blocks.

Connections closed during an event batch are reclaimed after every pointer returned by that `epoll_wait` has been consumed. A running
worker additionally retains its connection until completion. The retirement list reuses the completion link after worker ownership
ends, preserving the 192-byte connection structure and preventing stale event pointers from reaching freed memory.

## Admission control

Resource growth is bounded at three independent boundaries:

- `max_connections` limits live accepted sockets.
- `work_queue_capacity` limits queued commands in addition to the configured active workers.
- `max_inflight_payload_bytes` limits request bodies retained across every connection.

The reactor reserves payload bytes before allocation. If a declared body would exceed the global budget, the server sends `BUSY` and
closes that connection after the response because the rejected body remains unread in the TCP stream. A full worker queue returns
`BUSY` without blocking the reactor and leaves the connection reusable after the response.

`max_frame_size` is both a protocol safety bound and a per-request memory bound. It must not exceed the compiled absolute limit of
64 MiB, and the global payload budget must be at least one maximum-sized frame.

## Authentication and transport

Every C-driver connection performs an explicit `AUTH` exchange before `geo_client_connect` succeeds. The server compares token contents
without a content-dependent early exit and erases its owned token buffer during destruction. Authentication is not encryption: the
current listener must remain on a trusted network or behind a trusted TLS transport.

Both endpoints enable `TCP_NODELAY`. The driver uses a persistent blocking socket with configured connect/read/write deadlines. The
server uses nonblocking sockets, `accept4`, and `epoll`; database work never runs on the reactor.

## Durability and visibility

RocksDB is the only canonical transactional store. One synchronized `WriteBatch` atomically commits object values, GeoDoc metadata,
spatial deltas, typed secondary-key deletes/inserts, exact index cardinalities, sequence counters, and operation counters. There is no
independent GeoBolt WAL and therefore no dual-write recovery window.

The Morton segment set is a derived read index. The engine holds the visibility gate across the RocksDB linearization point and spatial
memtable publication. The asynchronous frozen-generation worker later prepares an immutable segment outside that gate and commits the
derived application watermark in the same checksummed manifest replacement that publishes the segment. The catalog watermark is then
advanced as a recoverable acknowledgement. If a crash
leaves the manifest ahead, startup advances the catalog without replaying an already-published range. If the catalog is ahead of a
valid manifest, the derived index is rebuilt from canonical objects. A missing or invalid derived manifest follows the same bounded
rebuild path; canonical RocksDB state is never discarded to repair a cache.

The persistent commit coordinator groups concurrent unconditional upserts/deletes into one RocksDB batch. Conditional inserts and
updates remain individually ordered so their existence checks are linearizable. Collection uses a monotonic deadline and stops early
at a full group, an incompatible request, or shutdown; it never leaves the producer tail pointing into a detached group while waiting.
After the canonical sync, compact spatial operations
enter a bounded active memtable with an ID hash and 12-bit Morton-prefix buckets. A persistent worker freezes full generations,
collapses repeated IDs, and publishes one immutable segment asynchronously while the replacement active generation accepts writes.
Queries combine persisted, frozen, and active state under one visibility gate; overlay hashes suppress stale segment records without
allocating a count-only result. RocksDB owns its mutable/frozen memtables, WAL, SST flush, checksums, and compaction. GeoBolt owns the
domain-specific spatial memtable, Morton segment compaction, and exact spatial visibility.

Spatial segment contents and the spatial mutation log are reconstructible cache state, but publication still follows strict durability
ordering: segment and mutation bytes are synchronized before a manifest may reference them. If corruption is detected, GeoBolt rejects
the complete derived snapshot and rebuilds it from canonical RocksDB objects and spatial deltas. An explicit database checkpoint drains
active and frozen spatial generations, publishes the compact mutation snapshot, advances the catalog acknowledgement, removes its
applied spatial-delta prefix in the same RocksDB batch, and synchronously flushes RocksDB.

The coordinator lives directly in `src/engine/geo_commit_coordinator.c`; persistence/recovery remains in `src/engine`, canonical object
formats in `src/db/object`, metadata encoding in `src/db/metadata`, and the confined C++ RocksDB adapter in `src/storage/rocksdb`.
Server-generated IDs use the high half of the 64-bit space by convention. Their monotonic allocator is persisted in the catalog
transaction that inserts the object, so reopen cannot reuse a successfully committed ID. Caller-supplied IDs remain valid across the
full nonzero range; the
allocator probes past an occupied generated-range value instead of overwriting it.

Client-supplied idempotency keys are also canonical state. The first successful object batch atomically stores a semantic fingerprint,
sequence, and expiry beside the mutations. A matching retry returns the stored decision without another write or Morton publication;
different content under a live key is rejected. Expired entries stop deduplicating immediately and are reclaimed by the idempotency
column family's RocksDB compaction filter.

Typed secondary indexes share one catalog and one ordered RocksDB column family. Index IDs partition key ranges; the encoded scalar
value sorts next and the generic object ID is the final stable tie-breaker. Updates read the prior canonical object only when at least
one typed index exists, derive old/new keys, eliminate unchanged values, and commit every required key transition in the object's
canonical batch. This preserves last-write-wins semantics across update, delete, conditional failure, crash, and reopen without stale
query candidates. Index definition publication uses `BUILDING` and `READY` states, and startup completes an interrupted build before
accepting traffic.

Each built index persists at most 64 adaptive equi-depth bins. Immutable boundaries/distinct-count baselines are separated from mutable
per-bin counters, preventing definition-sized rewrites on the object write path. The same canonical RocksDB batch updates an object's
secondary keys, exact index cardinality, and only the histogram counters whose bins changed. Checksums, monotonic encoded boundaries,
type-specific key validation, exact counter totals, and zero-cardinality reopen are verified before a catalog becomes visible.

Conjunctive radius queries use a connection- or caller-owned reusable workspace. Persisted estimates are compared with the selected
multi-resolution Morton-cover estimate before any broad secondary iterator opens. When metadata is selective, typed scans run in
estimated-cardinality order and intersect through an epoch-tagged open-addressed ID set; each later scan materializes only IDs still in
the progressively shrinking conjunction, and an empty driver stops the remaining work. Bounded pinned RocksDB `MultiGet` then feeds
the exact spherical predicate.

When geography is selective, Morton segments and active/frozen generations execute the exact spherical predicate first. Only those
geographic matches proceed to canonical resolution. An eight-way set-associative, byte-bounded process-local cache stores the current
Morton code, sequence, and already validated GeoDoc for generic objects; it is not tied to a particular metadata type or workload.
Cache misses use bounded pinned RocksDB `MultiGet`, validate the canonical object once, and may populate an empty cache way. Both
secondary-driven and spatial-driven plans share this resolver.

Writers prepare replacement cache values outside the visibility write gate, then publish updates or remove deleted IDs only after the
synchronized RocksDB batch, secondary-index publication, and spatial-memtable append have succeeded. Readers may fill empty ways
concurrently but never reclaim entries. Writer replacement and reclamation run only while the exclusive visibility gate proves that no
query retains a view. A zero `object_cache_bytes` setting disables the layer and preserves the canonical RocksDB path. The default is
128 MiB, including the lookup table and requested entry bytes; allocator bookkeeping may add a small implementation-dependent amount.

The visibility gate adds writer admission ahead of the POSIX read/write lock. Once any writer announces that it is waiting, new readers
sleep on a condition instead of repeatedly barging ahead; only readers that had already entered the fast path may complete first. This
bounds writer acquisition under continuous query traffic while retaining one atomic-load fast-path check for the uncontended reader.

All typed predicates are evaluated directly against validated GeoDoc values, so a broad metadata filter never materializes its complete
secondary range. Both plans reuse the same query workspace and hold the database visibility read gate, so index maintenance, canonical
objects, cache entries, and spatial generations belong to one observable commit boundary. Every query releases pinned RocksDB values
before leaving that gate; reusable workspaces retain allocation capacity but never retain storage-owned memory across database close or
reopen.

For a single immutable segment, radius planning uses index-aware range bounds and reuses them during execution. Density refinements are
evaluated from finest to coarsest. A lower bound composed from the current candidate cardinality and the cheapest possible one-range
cover stops construction once no remaining coarser cover can beat the selected cost; this removes dominated recursive Morton-cover and
binary-search work without changing exact spherical filtering.

## Shutdown contract

`SIGINT` and `SIGTERM` are blocked in process threads and consumed by a dedicated signal-monitor thread. The signal requests reactor
termination through an `eventfd`, avoiding async-signal-unsafe handlers. After the reactor returns, destruction stops and joins the
persistent worker queue, closes connections, drains the commit coordinator, stops storage background work, and closes RocksDB. A
synchronized WriteBatch remains recoverable even if its client disconnects before receiving the response.

Library callers must request stop, join the thread executing `geo_server_run`, and only then call `geo_server_destroy`. Destroying a
server concurrently with its reactor is outside the API contract.
