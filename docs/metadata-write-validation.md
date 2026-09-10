# Metadata write-path validation

## Scope

The implementation preserves input and persisted-document validation while sharing immutable validated GeoDoc views across typed index
maintenance and canonical serialization. For seven indexes, this removes repeated full-document checksum/tree walks and their validation
allocations; it does not remove typed-field checks, old-key deletion, new-key insertion, or synchronized canonical writes.

This change accompanies fixes for iterator-OOM false index publication, type-rejected builds preventing reopen, and a shared-client
descriptor race on disconnect. Failed rollback quarantines further writes until recovery. The coordinator additionally uses monotonic,
capacity-aware collection. The comparisons below isolate the validated-view change; the later collection change is validated separately.

## Workload and environment

- Command: `bin/benchmark_server 4 16000 10000 metadata`.
- Four synchronous TCP clients, four server workers, one reactor, and the persistent commit/maintenance workers.
- 10,000 seeded globally distributed spatial objects; 64,000 individually submitted measured upserts. Each new ID receives two
  operations, first with the old metadata values and then replacements. No client-side write batching in the measured phase.
- Seven indexes: boolean, signed integer, unsigned integer, double, datetime, string, and bytes. Both GeoDoc versions include a 1 KiB
  unindexed byte field. Every index is checked after measurement, not just the geographic cardinality.
- Each connection warms PING/radius paths; each writer also performs 32 unmeasured individual upserts before the start barrier.
- Final cardinality: **42,000**; geographic query checksum: **2654983261273313675**; typed-index ID checksum: **14383500815061827456**.
  Every reported run completed **192,182** requests with **zero** backpressure rejections.
- Intel Core i9-11900K, 8 cores/16 threads; Clang 22.1.8, C11 `-O3 -flto -Werror`, x86 AVX-512 runtime backend,
  jemalloc 5.3, RocksDB 10.2.1.
- Whole-process affinity restricted to CPUs 0–7. Workers were **not individually pinned**; governor remained `powersave`, turbo/frequency
  were not fixed, and the host was **not exclusively reserved**. Pre-run load was 0.40/0.59/0.62. These are workload observations, not
  maximum-capacity certification or small-delta claims.

## tmpfs comparison

`/tmp` was confirmed to be **tmpfs**. Synchronous RocksDB APIs were used, but this storage is volatile and does not measure NVMe fsync
cost or power-loss durability. Five alternating before/after pairs produced:

| Pair | Before upserts/s | After upserts/s | Before user CPU seconds | After user CPU seconds |
|---|---:|---:|---:|---:|
| 1 | 17,633 | 24,517 | 2.62 | 1.70 |
| 2 | 17,784 | 24,705 | 2.54 | 1.65 |
| 3 | 17,739 | 24,574 | 2.55 | 1.72 |
| 4 | 17,775 | 24,598 | 2.61 | 1.63 |
| 5 | 17,694 | 24,497 | 2.60 | 1.69 |

Median write throughput increased from 17,739 to 24,574 upserts/s (**38.5%**). User CPU seconds are for the complete process lifecycle,
including setup, other benchmark phases, verification, and teardown; they are not isolated commit CPU time. Peak RSS remained in the
same approximate 135–138 MiB range; there is no demonstrated RSS reduction.

### Final-source regression check

After the collector changes and final correctness matrix, another five alternating pairs compared the original baseline against the
combined final source. The workload, warmups, affinity, compiler, allocator, and checksums were unchanged. No test/build jobs overlapped
this series; observed host load was 0.43/0.55/0.79, without exclusive reservation or fixed frequency. The operator had confirmed no
other heavy disk process. This remains a tmpfs CPU/software-path comparison, not a durable-storage capacity claim.

| Pair | Baseline upserts/s | Final upserts/s | Baseline user CPU seconds | Final user CPU seconds |
|---|---:|---:|---:|---:|
| 1 | 17,660 | 24,788 | 2.57 | 1.58 |
| 2 | 17,697 | 24,547 | 2.55 | 1.68 |
| 3 | 17,732 | 25,011 | 2.62 | 1.52 |
| 4 | 17,720 | 24,715 | 2.52 | 1.59 |
| 5 | 17,650 | 24,612 | 2.62 | 1.61 |

The combined final-source median is 24,715 versus 17,697 upserts/s (**39.7% higher**). The complete-process median user CPU time falls
from 2.57 to 1.59 seconds. Peak RSS ranges overlap (baseline 137,960–140,596 KiB; final 137,764–141,080 KiB). Every run preserves the
42,000-object cardinality and both checksums specified above. This confirms the combined change preserves the earlier large tmpfs
gain; it does not isolate or establish an additional small throughput benefit from the collector alone.

## Btrfs/NVMe comparison — inconclusive

The same benchmark source was linked against the before/after libraries and run with
`GEOBOLT_BENCHMARK_DIRECTORY=/home/bencz/programming/GeoBolt`. Each run created and removed its own unique child directory. This is
Btrfs on `/dev/nvme1n1p3`, a Samsung SSD 980 500 GB. The request semantics, seeds, checksums, affinity, and client counts were unchanged.

| Pair | Before upserts/s | After upserts/s | Before total seconds | After total seconds |
|---|---:|---:|---:|---:|
| 1 | 1,785 | 2,095 | 36.65 | 31.37 |
| 2 | 1,183 | 1,752 | 54.93 | 37.68 |
| 3 | 711 | 437 | 91.30 | 147.86 |
| 4 | 451 | 448 | 143.12 | 143.81 |
| 5 | 411 | 452 | 156.78 | 142.81 |

Both variants slowed substantially as the sequence progressed. The raw medians are 711 before and 452 after; this apparent **36.4%
decrease is not hidden**, but the strong time trend and overlapping ranges do not establish an attributable code regression or gain.
Do not use the tmpfs improvement as a claim of improved durable disk throughput. The disk slowdown remains an investigation item,
including synchronization latency, filesystem/device state, background I/O, and group occupancy at only four synchronous clients.

The operator confirmed that no other heavy disk/NVMe process was running. External load is therefore not an established explanation
for the decline. Filesystem/device behavior and GeoBolt's own sustained write/synchronization pattern remain to be isolated.

### Follow-up synchronization diagnostic

The final coordinator was also exercised with `strace -q -f -w -c` over the same TCP metadata workload, reduced to 4,000 operations per
client for diagnosis. This is an instrumented run, not a throughput comparison. Affinity, compiler, backend, allocator, and disk were
unchanged. It completed 48,182 requests without backpressure, with 18,000 visible objects, geographic checksum 15491120351879726942,
and typed-index checksum 13759338793768420832.

Across setup, warmup, measured requests, verification, and teardown, 8,104 `fdatasync` calls accumulated 10.027 seconds of elapsed wait
(1.237 ms mean); another 22 `fsync` calls accumulated 0.032 seconds. The tracing used wall-clock accounting (`-w`), not the default
system-CPU accounting. The measured write phase lasted approximately 11.56 seconds, but the synchronization totals cover the whole
process lifetime, so their ratio is not a precise commit-time percentage. Tracing itself perturbs scheduling and grouping.

The count suggests small commit groups with four synchronous clients, but includes setup and maintenance and is not a direct group
occupancy measurement. The large aggregate `futex`/socket wait totals include overlapping waits in many threads and must not be added
to derive request latency. Synchronization is a significant limiting resource in this diagnostic; it does not yet explain the earlier
time-dependent degradation. No WAL synchronization, durability guarantee, or logical request was removed to improve the numbers.

## Acceptance and remaining work

Meaningful regression coverage includes allocation failure at the iterator boundary, incompatible-type rollback, subprocess interruption
after durable `BUILDING`, rollback I/O failure and write quarantine, successful retry/reopen, all typed key updates, and simultaneous
shared-client requests after disconnect. Full-group collection is checked by observing actual timed-wait calls, not a timing threshold.

The final source passed the full optimized Clang and GCC suites, the forced-scalar suite, and ASan/UBSan with leak detection enabled.
These include index, database, storage, server, daemon, and injected index-failure tests. The index suite reports 71 cases normally
and 60 with scalar forced. Separate ThreadSanitizer database/server runs passed, including simultaneous mutations, group commit,
shared-client disconnection, and reactor connection churn. Clang static analysis of the changed C modules/tests/benchmark completed
without diagnostics. All builds retained `-Werror`; sanitizer timings are not performance evidence. ARM64 and Darwin were not executed.

The million-actor network soak, broad-query memory bounds, online index building, and larger commit/visibility redesign remain tracked
in [production-readiness.md](production-readiness.md). This benchmark is not a substitute for those features or validation workloads.
