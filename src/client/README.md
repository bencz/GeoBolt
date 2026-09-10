# Client boundary

The blocking C driver is exposed through `<geobolt/client.h>` and compiled into `libgeobolt.a`. It depends only on the public database
types and the private wire codec; it never reaches into engine or storage internals. One `GeoClient` is safe for concurrent commands:
its mutex serializes complete request/response transactions over the persistent TCP connection.
The connection-state check uses that same mutex. After an I/O or framing failure closes the stream, queued and subsequent commands
return `GEO_CLIENT_IO_ERROR`; destruction still requires all callers to be quiescent.

Connection establishment uses nonblocking `connect` plus `poll`, so `timeout_ms` applies to the connect attempt as well as subsequent
socket reads and writes. Small requests use one `sendmsg` for header and payload, and both endpoints enable `TCP_NODELAY` to avoid the
Nagle/delayed-ACK latency cliff. If I/O, allocation, or validation fails after a request is sent, the driver closes the descriptor rather
than risk reusing a desynchronized stream.

Every successful `geo_client_connect` has completed an explicit `AUTH` exchange, including deployments configured with an empty token.
When the server requires a token, omitting it or supplying the wrong value fails during connect instead of returning a handle whose first
application command is guaranteed to be rejected.

The current commands are authentication, ping, durable mutation batches, exact radius counts, database statistics, and an explicit
checkpoint. The server performs normal checkpoints autonomously; the explicit command exists for administration and validation, not as
a routine application responsibility. Callers must ensure no other thread is using a client while `geo_client_close` runs.
