# Binary protocol version 1

All integers and IEEE-754 binary64 bit patterns are encoded in network byte order. Native C structure layout is never sent on the wire.
Every request and response consists of a fixed 48-byte header followed by exactly `payload_size` bytes.

## Header

| Offset | Size | Field | Contract |
|---:|---:|---|---|
| 0 | 4 | magic | `0x47424c54` (`GBLT`) |
| 4 | 2 | version | `1` |
| 6 | 2 | opcode | command identifier |
| 8 | 4 | flags | `0` for requests, `1` for responses |
| 12 | 4 | status | `0` in requests; response status in responses |
| 16 | 8 | request ID | copied unchanged into the response |
| 24 | 4 | payload size | body length, excluding the header |
| 28 | 4 | reserved | must be zero |
| 32 | 8 | payload checksum | checksum of the exact body bytes |
| 40 | 8 | header checksum | checksum of all 48 bytes with this field zeroed |

The checksum detects corruption and framing mistakes; it is not a MAC and provides no authentication or confidentiality. A malformed
header, invalid request flags/status, or payload checksum mismatch is a fatal stream error and closes the connection.

## Opcodes

| Value | Command | Request payload | Success payload |
|---:|---|---|---|
| 1 | `AUTH` | raw token bytes | empty |
| 2 | `PING` | arbitrary bytes | exact request bytes |
| 3 | `WRITE` | batch described below | empty |
| 4 | `RADIUS_COUNT` | three binary64 values: latitude, longitude, radius km | one `u64` count |
| 5 | `STATS` | empty | 64-byte database statistics block |
| 6 | `CHECKPOINT` | empty | empty |

Except for `AUTH`, a command received before successful authentication returns `UNAUTHORIZED`. The driver always authenticates during
connection establishment.

### WRITE payload

The body begins with a big-endian `u32 operation_count`, followed by exactly that many 24-byte operations:

| Operation offset | Size | Field |
|---:|---:|---|
| 0 | 8 | nonzero generic object ID |
| 8 | 8 | packed Morton coordinate; ignored for delete |
| 16 | 4 | operation: `1` upsert, `2` delete |
| 20 | 4 | reserved, must be zero |

IDs must be unique within one batch. A batch is one durable database operation for publication purposes, but each contained mutation
consumes one sequence number and increments `committed_operations`.

### STATS response

The response contains eight big-endian `u64` fields: committed operations, recovered operations, checkpoints, maintenance failures,
next sequence, physical records, active segments, and one reserved zero field.

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

Version 1 intentionally has no metadata records, point `GET`, result-set streaming, pipelining, TLS, or compression. These features must
receive new explicitly versioned payload contracts; they must not be smuggled into reserved fields.
