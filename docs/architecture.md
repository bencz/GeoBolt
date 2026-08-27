# Single-server architecture

## Process boundary

`geoboltd` is one foreground process that owns exactly one database root. The current milestone is deliberately single-server: there is
no replication, leader election, distributed transaction, or cross-node failover protocol. A supervisor such as systemd or Kubernetes
owns process restart.

The process contains four long-lived execution domains:

1. One `epoll` reactor owns listening, accepting, connection parsing, response writes, and connection lifetime.
2. A bounded persistent job queue executes commands that can block or consume CPU. No worker is created per connection or request.
3. The database serializes durable mutations and publishes immutable spatial segments to concurrent readers.
4. The segment subsystem owns its background compaction worker and reusable parallel compaction pool.

A connection has at most one in-flight request. This keeps request ordering deterministic and bounds connection state without requiring a
per-connection response reorder buffer. Applications obtain parallelism through multiple persistent connections; a single C-driver
handle may be shared safely but intentionally serializes its complete request/response transactions.

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

A mutation batch receives a contiguous sequence range. Its immutable spatial segment is prepared, then the checksummed WAL frame is
appended and synchronized before the mutation becomes visible. The segment manifest and mutation visibility state are published before
the durable state watermark advances. Recovery can therefore replay a committed frame idempotently after a crash between publication
steps.

Automatic checkpoints bound the WAL and segment compaction bounds immutable-segment fanout. The current implementation still creates a
segment for each committed batch; the storage roadmap describes the pending memtable/group-commit milestone.

## Shutdown contract

`SIGINT` and `SIGTERM` are blocked in process threads and consumed by a dedicated signal-monitor thread. The signal requests reactor
termination through an `eventfd`, avoiding async-signal-unsafe handlers. After the reactor returns, destruction stops and joins the
persistent worker queue, closes connections, stops storage background work, and closes the WAL/database. A committed WAL frame remains
recoverable even if its client disconnects before receiving the response.

Library callers must request stop, join the thread executing `geo_server_run`, and only then call `geo_server_destroy`. Destroying a
server concurrently with its reactor is outside the API contract.
