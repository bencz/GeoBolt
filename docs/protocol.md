# Binary protocol version 4

All integers and IEEE-754 binary64 bit patterns use network byte order. Native C structure layout is never sent on the wire. Every
request and response contains a fixed 48-byte header followed by exactly `payload_size` bytes.

## Header

| Offset | Size | Field | Contract |
|---:|---:|---|---|
| 0 | 4 | magic | `0x47424c54` (`GBLT`) |
| 4 | 2 | version | `4` |
| 6 | 2 | opcode | command identifier |
| 8 | 4 | flags | `0` for requests, `1` for responses |
| 12 | 4 | status | `0` in requests; response status in responses |
| 16 | 8 | request ID | copied unchanged into the response |
| 24 | 4 | payload size | body length, excluding the header |
| 28 | 4 | reserved | must be zero |
| 32 | 8 | payload checksum | checksum of the exact body bytes |
| 40 | 8 | header checksum | checksum of all 48 bytes with this field zeroed |

The checksum detects corruption and framing mistakes; it is not a MAC and does not provide confidentiality. A malformed header,
invalid request flags/status, or payload checksum mismatch is a fatal stream error and closes the connection.

A TCP write half-close means the client will send no additional bytes. Complete buffered frames still execute in order and receive
their responses; EOF after a partial header or body discards that incomplete frame without dispatch. A reset or full transport failure
can interrupt response delivery, so an already dispatched mutation still has the documented ambiguous-outcome semantics.

## Commands

| Value | Command | Request payload | Success payload |
|---:|---|---|---|
| 1 | `AUTH` | raw token bytes | empty |
| 2 | `PING` | arbitrary bytes | exact request bytes |
| 3 | `WRITE` | compact batch below | empty |
| 4 | `RADIUS_COUNT` | latitude, longitude, radius km as three binary64 values | one `u64` count |
| 5 | `STATS` | empty | 64-byte statistics block |
| 6 | `CHECKPOINT` | empty | empty |
| 7 | `WRITE_OBJECTS` | variable object batch below | empty |
| 8 | `GET_OBJECT` | one nonzero `u64` object ID | object response below |
| 9 | `INSERT_GENERATED` | generated-insert payload below | generated nonzero `u64` object ID |
| 10 | `WRITE_OBJECTS_IDEMPOTENT` | idempotency prefix plus a `WRITE_OBJECTS` payload | one `u64`: zero committed, one replayed |
| 11 | `CREATE_INDEX` | typed definition below | empty |
| 12 | `DROP_INDEX` | index name below | empty |
| 13 | `LIST_INDEXES` | empty | typed definition list below |
| 14 | `QUERY_INDEX` | typed equality/range predicate below | statistics and generic object IDs below |
| 15 | `QUERY_RADIUS` | radius plus one or more typed predicates below | planner statistics and matching records below |

Except for `AUTH`, a command before successful authentication returns `UNAUTHORIZED`. The C driver authenticates while connecting.

### Compact WRITE

The payload begins with a `u32 operation_count`, followed by exactly that many 24-byte entries. It is the low-overhead metadata-free
ingestion path.

| Entry offset | Size | Field |
|---:|---:|---|
| 0 | 8 | nonzero generic object ID |
| 8 | 8 | packed Morton coordinate; zero for delete |
| 16 | 4 | operation: `1` upsert, `2` delete |
| 20 | 4 | reserved zero |

### WRITE_OBJECTS

The payload begins with a `u32 operation_count`. Each operation then contains a 32-byte header followed immediately by
`document_size` bytes of canonical GeoDoc.

| Entry offset | Size | Field |
|---:|---:|---|
| 0 | 8 | nonzero generic object ID |
| 8 | 8 | packed Morton coordinate; zero for delete |
| 16 | 4 | `1` upsert, `2` delete, `3` conditional insert, `4` conditional full update |
| 20 | 4 | reserved zero |
| 24 | 4 | GeoDoc byte count |
| 28 | 4 | reserved zero |
| 32 | variable | exact GeoDoc bytes |

IDs must be unique within one batch. Delete requires a zero Morton code and no document. Insert returns `ALREADY_EXISTS` if its ID is
present; update returns `NOT_FOUND` if absent. Upsert and full update replace the complete prior location and metadata. The batch consumes
one contiguous commit sequence per entry and is one atomic RocksDB transaction.

### WRITE_OBJECTS_IDEMPOTENT

The payload starts with a `u32 key_size` and `u32 retention_seconds`, followed by the opaque key bytes and then one complete
`WRITE_OBJECTS` payload beginning at its `operation_count`. Keys contain 1–128 bytes and retention is 1–604800 seconds. The first
successful write stores its 128-bit semantic fingerprint and expiry in the same RocksDB transaction. A matching retry returns success
with replay flag one without consuming sequences or publishing a segment. Reusing a live key with different mutations returns
`IDEMPOTENCY_CONFLICT`. RocksDB compaction removes expired records autonomously.

### GET_OBJECT response

| Offset | Size | Field |
|---:|---:|---|
| 0 | 8 | object ID |
| 8 | 8 | last committed object sequence |
| 16 | 8 | packed Morton coordinate |
| 24 | 4 | GeoDoc byte count |
| 28 | 4 | reserved zero |
| 32 | variable | exact canonical GeoDoc bytes |

### INSERT_GENERATED

This command inserts a new object while allocating its generic ID from the server-owned durable range. The allocation counter and the
new object are committed in the same RocksDB batch. This opcode does not accept an idempotency key; use caller-generated IDs with
`WRITE_OBJECTS_IDEMPOTENT` when retry deduplication is required.

| Offset | Size | Field |
|---:|---:|---|
| 0 | 8 | packed Morton coordinate |
| 8 | 4 | GeoDoc byte count |
| 12 | 4 | reserved zero |
| 16 | variable | exact canonical GeoDoc bytes |

### STATS response

The response contains eight `u64` fields: committed operations, recovered spatial operations, explicit checkpoints, maintenance
failures, next sequence, physical Morton records, active Morton segments, and one reserved zero.

### Typed secondary-index commands

`CREATE_INDEX` begins with four `u32` values: type, name byte count, JSON Pointer byte count, and reserved zero. Name bytes and
JSON Pointer bytes follow without terminators. Names contain 1–63 bytes and pointers at most 511 bytes. Types are `1` boolean, `2`
signed 64-bit, `3` unsigned 64-bit, `4` finite binary64, `5` datetime represented by a signed 64-bit GeoDoc value, `6` UTF-8 string,
and `7` arbitrary bytes. Creation returns only after the index is durable and ready; interrupted `BUILDING` definitions are recovered
before startup accepts traffic.

`DROP_INDEX` contains a `u32` name byte count, a reserved-zero `u32`, and the raw name. `LIST_INDEXES` starts with a `u32` definition
count and reserved-zero `u32`. Each variable definition contains `u64 index_id`, `u64 entry_count`, `u32 type`, `u32 name_size`,
`u32 pointer_size`, reserved-zero `u32`, then name and pointer bytes. Cardinality is committed transactionally with object changes.

`QUERY_INDEX` begins with six `u32` values: operator, type, name size, lower-value size, upper-value size, and reserved zero. Name, lower,
and upper bytes follow. Operators are equality, less, less-or-equal, greater, greater-or-equal, and inclusive between (`1` through `6`).
Only between carries an upper value. Booleans occupy one byte; integer/datetime values and binary64 occupy eight network-order bytes;
string/bytes values are raw and may be empty. The response contains `u64 scanned_entries`, `u64 matched_entries`, `u32 result_count`,
reserved-zero `u32`, and `result_count` network-order object IDs. The normal frame limit bounds materialization; result streaming is not
silently emulated.

### Geographic and typed conjunction

`QUERY_RADIUS` starts with latitude, longitude, and radius km as three binary64 values, followed by `u32 predicate_count` and a
reserved-zero `u32`. Between 1 and 16 complete `QUERY_INDEX` predicate encodings follow. Persisted typed histograms and Morton density
are costed before execution. The metadata-driven plan intersects typed ranges progressively and then performs canonical object lookups.
The spatial-driven plan executes exact geography first, reads only the matching canonical objects in bounded `MultiGet` chunks, and
evaluates every typed predicate directly against GeoDoc. Both plans execute under one visibility snapshot.

The 80-byte success prefix contains, in order, seven `u64` counters (secondary entries scanned, final metadata matches,
estimated spatial candidates, canonical object lookups, spatial records scanned, spatial records matched, and Morton ranges checked),
one binary64 total query time, `u32 plan` (`1` spatial-driven or `2` secondary-driven), `u32 predicate_count`, `u32 result_count`, and a
reserved-zero `u32`. Each result is then a network-order `{ u64 object_id; u64 morton_code; }`. The configured frame limit remains a
hard materialization bound.

## Response statuses

| Value | Meaning |
|---:|---|
| 0 | success |
| 1 | invalid request or argument |
| 2 | authentication required or failed |
| 3 | bounded queue or payload budget is saturated |
| 4 | server allocation/internal failure |
| 5 | database failure |
| 6 | declared frame exceeds the configured maximum |
| 7 | protocol error |
| 8 | object not found |
| 9 | conditional insert conflict |
| 10 | live idempotency key was reused with different mutation content |

Version 4 does not yet provide result-set streaming, pipelining, TLS, or partial GeoDoc patches. Those features require explicit future
protocol contracts and cannot reuse reserved fields.
