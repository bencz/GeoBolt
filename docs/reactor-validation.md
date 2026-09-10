# Reactor lifecycle validation — 2026-09-10

## Changes and regression evidence

The previous reactor closed on `EPOLLRDHUP` even when a complete request still needed a response. It also freed closed connections while
other pointers from the same `epoll_wait` batch could still be consumed. Reclamation is now deferred until the batch ends and any
worker completes. The reactor avoids worker-owned response fields while a job is in flight. Zero-payload parsing no longer performs
pointer arithmetic on a null payload pointer. Declared object-operation counts must fit the received frame before decoder allocation.

`tests/test_server.c` now covers truncated mutation headers/bodies, impossible operation counts, an acknowledged conditional insert
after FIN, and eight concurrent workers each opening 128 connections. Each sends three pipelined frames (0, 1, and 4,096 bytes), then
either resets the connection or half-closes and verifies all responses, checksums, order, and EOF. Connection and payload accounting
must return to zero before server destruction. The same test fails with the prior reactor and passes with the current implementation.
An additional 8 MiB response is held unread while an independent client connects and completes a PING, then drained and verified byte
for byte through EOF. This covers blocked-send resumption with a 128 KiB receive-buffer request, not just small immediate responses.

The optimized full suite passed (71 index tests, database, storage, server, and daemon). The server suite also passed with GCC, forced
scalar, ASan/UBSan with LeakSanitizer enabled, and ThreadSanitizer. Clang static analysis of the changed C translation units reported no
diagnostics. No ARM64 hardware run was performed for this change.

## TCP comparison

Both variants use the current engine/decoder and the same updated benchmark. Only the reactor implementation differs. Benchmark seeds
are deterministic per client, each connection performs unmeasured PING/radius warmups, and every run validates final visible count.
This removes the old address-dependent query sequence and corrects the old inverse-throughput label that incorrectly implied request
latency. No workload or assertion was reduced to pass a test.

Environment: Intel Core i9-11900K, 8 cores/16 hardware threads, Clang 22.1.8, C11 `-O3 -flto -Werror`, x86-64 runtime AVX-512 dispatch,
jemalloc 5.3.0, RocksDB 10.2.1. Four client threads, four database workers, one reactor, plus normal storage/background workers. The
100,000-object dataset uses deterministic globally distributed coordinates; each measured phase runs 10,000 operations per client.
Radius queries use 25 km. Upserts use disjoint new IDs and retain synchronous durability. Five alternating current/prior pairs ran
sequentially, after correctness builds/tests completed.

Host load before the campaign was 0.16/0.58/0.74; exclusive host idleness was not established. The governor was `powersave` with turbo
enabled, frequency was not fixed, and workers were not individually pinned. These are exploratory observations, not capacity limits or
proof of a small throughput improvement. Timings from strace or sanitizers are excluded.

| Phase | Prior runs (operations/s) | Current runs (operations/s) | Prior median | Current median |
|---|---|---|---:|---:|
| PING | 209073, 208515, 203450, 219042, 214715 | 214884, 232300, 227854, 222470, 225758 | 209073 | 225758 |
| Radius count | 197779, 205666, 207831, 206599, 199304 | 207298, 216537, 209968, 208537, 211260 | 205666 | 209968 |
| Durable upsert | 36077, 36498, 36184, 36090, 36184 | 36377, 36390, 36429, 36497, 36493 | 36184 | 36429 |

Every run returned query checksum `16875619389353842661`, 140,000 final visible objects, 120,040 completed requests, and zero
backpressure rejections. These are client/server PING, radius-count, and metadata-free write phases, not the conjunctive engine-only
benchmark and not the planned independent million-actor workload.

## Work removed

A separate `strace -f -c -e trace=epoll_ctl,epoll_wait,write` comparison used four clients, 1,000 operations per phase/client and the same
100,000-object setup. Both returned checksum `8706293000384914404`, 104,000 visible objects and 12,040 completed requests.

| System call | Prior count | Current count |
|---|---:|---:|
| `epoll_ctl` | 36150 | 24110 |
| `epoll_wait` | 12253 | 8146 |
| `write` (all process descriptors) | 14159 | 6240 |

The 12,040-call reduction in `epoll_ctl` is one avoided registration per completed response, or 33.3% in this trace. Coalescing the
completion queue also reduces notifications; the `write` totals include storage writes and must not be described as eventfd-only.
Immediate response sending and notification coalescing preserve partial-write retry and mutex-based ownership transfer.

The prior reactor source was retained for the comparison at `/tmp/geobolt-reactor-review.MQCBo1/geo_server_before.c` (temporary artifact).
Its SHA-256 is `9363dde09d0f8e260a70da38e75084c9c6547f9a507df18ff4d8371e34e13741`. The measured current reactor source SHA-256 is
`88453d258977b156a2a1aca1316f0e40efe8c06e99b558b4af35556512e8c301`.
