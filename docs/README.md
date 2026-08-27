# GeoBolt documentation

The documentation is versioned with the implementation. A change to a public API, wire frame, durability rule, configuration default,
operational procedure, or architecture boundary is incomplete until the corresponding document and sample are updated.

- [Single-server architecture](architecture.md) describes ownership, concurrency, backpressure, durability, and shutdown.
- [Binary protocol](protocol.md) defines the version 1 wire format and command payloads byte by byte.
- [Operations](operations.md) covers daemon startup, authentication, capacity limits, monitoring, and safe shutdown.
- [Storage roadmap](storage-roadmap.md) records the native-engine versus RocksDB boundary and the metadata/CRUD work that is not yet part
  of protocol version 1.

Runnable public-API programs live in [`samples/`](../samples/). Validation programs live in [`tests/`](../tests/); representative
capacity measurements live in [`benchmarks/`](../benchmarks/).
