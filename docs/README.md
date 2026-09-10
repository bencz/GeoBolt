# GeoBolt documentation

The documentation is versioned with the implementation. A change to a public API, wire frame, durability rule, configuration default,
operational procedure, or architecture boundary is incomplete until the corresponding document and sample are updated.

- [Single-server architecture](architecture.md) describes ownership, concurrency, backpressure, durability, and shutdown.
- [Binary protocol](protocol.md) defines the version 4 wire format and command payloads byte by byte.
- [Operations](operations.md) covers daemon startup, authentication, capacity limits, monitoring, and safe shutdown.
- [Production readiness](production-readiness.md) tracks remaining single-server release requirements and their acceptance criteria.
- [Reactor validation](reactor-validation.md) records lifecycle regressions, verification, and the reproducible TCP comparison.
- [Metadata write validation](metadata-write-validation.md) records failure-boundary fixes, the typed-metadata TCP comparison, and the
  distinction between tmpfs and disk-backed measurements, including the unresolved disk slowdown.
- [Storage architecture and roadmap](storage-roadmap.md) records the RocksDB/GeoBolt ownership boundary and remaining secondary-index
  work.

Runnable public-API programs live in [`samples/`](../samples/). Validation programs live in [`tests/`](../tests/); representative
capacity measurements live in [`benchmarks/`](../benchmarks/).
