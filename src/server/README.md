# Server boundary

`geoboltd` is the foreground single-server process intended for systemd, containers, and direct supervision. The Linux reactor owns
nonblocking accept/read/write state through `epoll`; a persistent bounded worker pool executes database commands outside that reactor.
Each connection permits one in-flight command, which preserves request order and bounds per-connection working state.

Write half-close drains complete queued requests and their ordered responses before EOF closes the connection; incomplete frames never
reach a worker. Connection memory is reclaimed after the current epoll event batch and any running worker have released ownership.
Worker completions coalesce eventfd notifications, and the reactor tries immediate writes before subscribing to socket writability.

Admission control has three independent bounds: connection count, queued jobs, and globally retained request-payload bytes. A saturated
job queue returns `BUSY`; a payload that would exceed the global memory budget is rejected before allocation and the connection is
closed after its `BUSY` response because its unread body cannot safely remain in the stream. Runtime counters expose accepted/rejected
connections, protocol/authentication errors, backpressure, current payload bytes, and the observed payload-byte high-water mark.

The binary protocol is versioned, endian-independent, length-delimited, and checksums both headers and payloads. Authentication tokens
are compared without content-dependent early exit, copied into server-owned memory, and erased at shutdown. The protocol does not yet
provide encryption; never expose a token-authenticated listener across an untrusted network without a trusted TLS tunnel or equivalent
network boundary.

`SIGINT` and `SIGTERM` stop admission, terminate the reactor, complete already committed database work while destroying the persistent
worker pool, close the database cleanly, and return a process status suitable for a supervisor. HA, replication, and geo-replication
remain deliberately outside this single-server milestone.
