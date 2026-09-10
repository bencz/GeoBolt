# Samples

Samples demonstrate supported public APIs with complete error handling. They must not access private structures or use benchmark-only
shortcuts.

`make sample` builds and runs `basic_usage.c` as `bin/basic_usage`.

`make client-sample` builds `client_usage.c` as `bin/client_usage`. It demonstrates insert with a caller ID, insert with a durable
server-generated generic ID, retry-safe idempotent ingestion, get/select, update, upsert, stats, and delete.
Start `geoboltd`, then run it with the token outside the process argument list:

```bash
GEOBOLT_TOKEN='secret' bin/client_usage 127.0.0.1 7447
```

`make metadata-sample` builds `client_metadata.c`. It stores a generic object whose example domain happens to be a vehicle, including
driver ID/name, plate, online state, and location in one atomic object. It then reads the returned GeoDoc, creates a reusable typed
boolean index on `/online`, executes an indexed equality query, and combines that predicate with a 5 km radius through the cost-based
client/server planner:

```bash
GEOBOLT_TOKEN='secret' bin/client_metadata 127.0.0.1 7447
```
