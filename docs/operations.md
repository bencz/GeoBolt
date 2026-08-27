# Operating `geoboltd`

## Startup

Create a token file readable only by the server account, then run the daemon under a process supervisor:

```bash
install -m 600 /dev/null geobolt.token
printf '%s\n' 'replace-with-a-long-random-token' > geobolt.token
bin/geoboltd --data /var/lib/geobolt/primary --token-file geobolt.token --workers 16
```

The daemon rejects token files with group/other permission bits, embedded NUL bytes, an empty value, or more than 4096 bytes. One final
CR/LF sequence is stripped so a conventional newline-terminated secret file works as expected.

The default bind address is `127.0.0.1`. Binding to a non-loopback interface does not add encryption; use a trusted private network or a
TLS proxy until transport security is implemented in the protocol.

## Capacity configuration

- `--workers` controls concurrently executing database commands. Start near the number of cores available to the process, then measure
  the complete read/write mix; durable writes currently serialize in the engine.
- `--queue-capacity` absorbs short bursts but also increases queued latency. It is not a substitute for client-side retry with jitter.
- `--max-connections` bounds socket and per-connection state.
- `--max-frame` bounds one request or response body.
- `--max-inflight` bounds retained request bodies globally and must be at least `--max-frame`.
- `--backlog` controls the kernel listen queue, not the application worker queue.

Clients receiving `BUSY` should retry only idempotent operations automatically. A timed-out mutation is ambiguous: it may already be
durable even when its response was lost. Exactly-once retry requires an application idempotency key, which protocol version 1 does not
yet provide.

## Monitoring

The embedded server API exposes accepted, rejected, and active connections; completed requests; protocol and authentication failures;
backpressure rejections; current retained payload bytes; and peak retained payload bytes. Database stats expose sequence progress,
recovery, checkpoints, maintenance failures, physical records, and active segments.

Repeated growth in `backpressure_rejections` means arrival rate exceeds an explicit resource boundary. Determine whether the saturated
resource is connection slots, queued work, or retained payload memory before raising a limit. Raising queue depth can worsen p99 latency.

## Shutdown and recovery

Send `SIGTERM` or `SIGINT` and wait for a zero process exit. Do not use `SIGKILL` for routine rotation. A forced termination cannot lose a
successfully synchronized WAL frame, but clients whose responses were interrupted must treat their writes as ambiguous and reconcile
after reconnecting.

On restart, GeoBolt validates the state, manifest, immutable segment files, mutation checkpoint, and WAL sequence. It truncates only an
incomplete tail in the newest WAL segment; committed corruption causes startup failure rather than silent data loss.
