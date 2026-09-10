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

- `--workers` controls concurrently executing database commands. Concurrent writers feed the persistent group-commit coordinator; tune
  workers against the complete read/write mix and storage latency rather than core count alone.
- `--queue-capacity` absorbs short bursts but also increases queued latency. It is not a substitute for client-side retry with jitter.
- `--max-connections` bounds socket and per-connection state.
- `--max-frame` bounds one request or response body.
- `--max-inflight` bounds retained request bodies globally and must be at least `--max-frame`.
- `--backlog` controls the kernel listen queue, not the application worker queue.

Clients receiving `BUSY` should retry only idempotent operations automatically. A timed-out ordinary mutation is ambiguous: it may
already be durable even when its response was lost. Use `WRITE_OBJECTS_IDEMPOTENT` with a stable application key for retry-safe writes.
The first successful batch and its fingerprint are one transaction; matching retries report replay without another mutation. Choose a
retention window longer than the maximum retry horizon, up to seven days.

## Monitoring

The embedded server API exposes accepted, rejected, and active connections; completed requests; protocol and authentication failures;
backpressure rejections; current retained payload bytes; and peak retained payload bytes. Database stats expose sequence progress,
recovery, checkpoints, maintenance failures, physical records, active segments, and operation counts in both spatial memtable
generations.

`spatial_memtable_max_operations` is a hard per-generation bound and must be at least `group_commit_max_operations`. A full active
generation is swapped with a preallocated replacement and flushed asynchronously. If frozen has not completed when the replacement
also fills, writers wait for the persistent flush worker instead of allocating a third unbounded generation. A checkpoint drains both
generations before returning.

Embedded deployments can tune `GeoDatabaseConfig.object_cache_bytes`. The 128 MiB default keeps current generic objects close to the
query engine after validation; use a working-set measurement rather than sizing it from object count alone because GeoDoc sizes vary.
The configured budget includes cache slots and requested entry bytes. Set it to zero for a controlled fallback or when RocksDB's block
cache is intentionally the only read cache. `GeoDatabaseQueryStats.object_lookups` counts canonical RocksDB keys requested after cache
hits, making hit effectiveness observable without conflating it with the number of geographic candidates.

Repeated growth in `backpressure_rejections` means arrival rate exceeds an explicit resource boundary. Determine whether the saturated
resource is connection slots, queued work, or retained payload memory before raising a limit. Raising queue depth can worsen p99 latency.

## Shutdown and recovery

Send `SIGTERM` or `SIGINT` and wait for a zero process exit. Do not use `SIGKILL` for routine rotation. A forced termination cannot lose a
successfully synchronized RocksDB WriteBatch, but clients whose responses were interrupted must treat their writes as ambiguous and
reconcile with `GET` after reconnecting.

On restart, RocksDB validates and replays its WAL/SST state. GeoBolt validates the Morton manifest and segment files, replays durable
spatial deltas beyond the common applied watermark, accepts a checksummed manifest watermark that is ahead of the catalog after an
interrupted acknowledgement, or rebuilds the derived spatial index from canonical objects when its manifest is missing, invalid, or
behind the catalog. Canonical corruption fails startup rather than silently dropping data.
