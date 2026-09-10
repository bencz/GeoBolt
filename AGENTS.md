# GeoBolt Engineering Rules

These instructions apply to the entire repository. They are hard project rules, not suggestions. A change is incomplete if it violates
them, even when it compiles or appears faster in one benchmark run.

## Priorities

1. Preserve correctness, numerical behavior, memory safety, linearizability, and persistence guarantees.
2. Improve end-to-end throughput, tail latency, memory traffic, or build cost with evidence from a representative workload.
3. Keep the implementation readable enough to audit low-level behavior and maintain safely.
4. Prefer algorithmic and data-layout improvements over isolated instruction-level tricks.

GeoBolt targets very large production workloads. Optimize the whole pipeline and its limiting resource, not a synthetic inner loop in
isolation. Do not trade correctness or maintainability for an unmeasured theoretical gain.

API and on-disk compatibility are not required unless a task explicitly requires them. When breaking an API or format produces a better
design, update every in-repository caller, test, version marker, and document in the same change.

## Required engineering depth

- Work as a senior performance architect with deep C11, geographic-indexing, concurrency, persistence, memory-hierarchy, and systems
  knowledge. Do not stop after surface-level cleanup or the first locally faster implementation.
- Actively inspect algorithms, cardinalities, bit representations, masks, field layout, cache lines, branch behavior, vectorization,
  synchronization, allocation, NUMA placement, TLB pressure, and I/O passes for every performance-sensitive change.
- Consider advanced and lesser-known techniques when their preconditions fit: hierarchical bitmaps, broadword/SWAR operations,
  branchless masks, multiplicative hashing, radix decomposition, prefix compression, SIMD compaction, cache-line sharding, software
  pipelining, batched prefetch, and architecture-specific instruction selection.
- Use bitmasking and bitsets when they remove branches, compact intermediate state, enable SIMD, or avoid repeated geometric work. Keep
  scalar semantics explicit and handle every partial word and vector tail safely.
- Treat x86 SSE/AVX2/AVX-512 and ARM64 NEON/SVE/SVE2 as first-class backends. Understand their register pressure, frequency behavior,
  gather/scatter costs, mask semantics, transition costs, alignment constraints, and runtime feature requirements.
- Apply deep geographic-algorithm knowledge: exact spherical predicates, conservative cell bounds, multi-resolution covers, Morton and
  Hilbert locality, H3-inspired hierarchy, antimeridian handling, poles, and error-bounded approximations.
- API compatibility is not a design constraint for this project. Prefer the cleanest and fastest contract, then update all callers and
  persisted format versions atomically.
- Never use expertise as justification for opaque code or an unsupported performance claim. Advanced code must remain organized,
  documented around invariants, verified against a reference path, and retained only when representative measurements justify it.

## Mandatory code style

- Treat every compiler warning in first-party GeoBolt code as an error. Every supported Clang and GCC build profile, including C11,
  the confined C++ RocksDB bridge, optimized, debug, scalar, sanitizer, tests, samples, and benchmarks, must compile with `-Werror`.
- Never use `-w`, remove `-Werror`, or disable a warning category globally to make a build pass. Fix the cause. A local diagnostic
  suppression is allowed only when the construct is unavoidable, the exact compiler warning is documented, the suppression covers the
  smallest possible region, and both Clang and GCC validate it.
- Include third-party headers as system headers when appropriate so dependency diagnostics are not misattributed to GeoBolt. This does
  not permit suppressing warnings originating in first-party wrappers or bridge code.

- Use four spaces for indentation and never tabs.
- The line limit is 140 columns. There is no 80-column limit.
- Follow the checked-in `.clang-format` file. Function braces go on the following line; control-flow braces stay on the same line.
- Always use braces for `if`, `else`, loops, and similar control flow.
- Never place a control-flow body, function, or multiple statements on one line.
- Do not chain assignments or compress unrelated declarations and operations into one statement.
- Prefer one declaration per line when values have different meanings or lifetimes.
- Separate validation, preparation, execution, publication, cleanup, and error handling with logical blank lines.
- Use descriptive names. Single-letter names are acceptable only for conventional, tiny mathematical or loop scopes.
- Keep comments focused on invariants, non-obvious reasoning, ownership, numerical constraints, and measured tradeoffs. Do not narrate
  obvious syntax.
- Organize large translation units with clear sections and cohesive helpers. Split a module when responsibilities no longer form one
  auditable unit.
- Do not perform repository-wide mechanical formatting as part of an unrelated change.

Readable low-level code is a correctness feature. Dense, clever code that is difficult to review must be rewritten even when it is
technically valid.

Code aesthetics are not a performance objective. Never replace a measurably faster hot-path representation or algorithm merely because
an abstraction, data structure, or implementation looks cleaner. Intrinsics, explicit masks, specialized paths, unusual layouts, and
low-level control flow are acceptable when they produce a representative end-to-end gain and preserve correctness. Keep their
invariants and dispatch boundaries documented and their formatting auditable; do not confuse necessary low-level complexity with
condensed or disorganized code.

Do not duplicate algorithms, invariants, validation, or business semantics across scalar and ISA-specific paths. Share common behavior
through zero-cost `static inline` helpers, parameterized internal kernels, immutable dispatch tables, or controlled code generation. ISA
specialization may duplicate only the irreducible instruction-level kernel when intrinsics differ. A proposed abstraction must still be
checked in generated assembly and benchmarks so eliminating source duplication does not introduce indirect calls inside inner loops,
block inlining, or regress the complete workload.

## C and data-layout rules

- Use C11 and keep the portable scalar implementation valid. Guard POSIX, Linux, x86, and ARM-specific behavior explicitly.
- Avoid undefined behavior, implementation-defined type punning, signed overflow, invalid shifts, out-of-bounds vector tails, and
  unchecked size arithmetic.
- Check every size addition and multiplication that can be influenced by input, file metadata, or dataset cardinality.
- Preserve `GeoRecord` as a naturally aligned 16-byte `{ uint64_t id; uint64_t z; }` hot record unless a measured architectural change
  justifies replacing the format.
- Order structure fields deliberately. Group hot fields, place naturally aligned/wider fields before narrow flags where appropriate,
  and isolate independently written concurrent fields on separate cache lines when measurements or ownership require it.
- Add `_Static_assert` checks for persisted layouts, SIMD-dependent layouts, cache-line isolation, important sizes, alignment, and field
  offsets.
- Do not use `#pragma pack` or packed attributes on hot or mmap-traversed structures. Persisted headers should use explicit fixed-width
  layouts, padding, and serialization when portability is required.
- Do not use unions for aliasing. Use `memcpy`, explicit bit operations, or architecture intrinsics.
- `volatile` is not a synchronization primitive. Reserve it for hardware interfaces or benchmark sinks where its limited semantics are
  intentional.
- Add `restrict`, alignment assumptions, branch hints, prefetches, non-temporal stores, or huge-page advice only when their preconditions
  are guaranteed and representative benchmarks demonstrate a benefit.

## Performance work

- Start with profiling and cardinality/memory-traffic analysis. Rank work by expected end-to-end impact.
- Prefer fewer candidates, fewer passes, fewer cache misses, less synchronization, and less output copying before adding SIMD or assembly.
- Hot query paths must not allocate when the caller supplies reusable storage. Persistent executors must reuse worker threads and
  reusable scratch memory.
- Do not create threads per query, per batch, or per radix pass in a production hot path. Use a persistent pool when parallel startup
  cost can be material.
- Avoid global atomic counters that become coherence bottlenecks. Use sharded, per-worker, or chunked ownership where appropriate.
- Treat memory bandwidth, NUMA locality, TLB pressure, cache capacity, and false sharing as first-class constraints.
- jemalloc is optional. The Makefile may use it when detected and must retain a working system-allocator fallback. Code must not depend on
  jemalloc-specific behavior unless isolated behind an optional feature.
- Never enable `-ffast-math` or equivalent unsafe aggregate math flags. Preserve IEEE handling of NaN, infinities, poles, the
  antimeridian, antipodal inputs, and exact predicate boundaries.
- Individual safe flags such as contraction or architecture targets require numerical validation and a measured A/B result.
- Do not keep a speculative optimization merely because it is obscure or theoretically faster. Remove it if controlled measurements
  show a regression or no reliable benefit.
- Do not sacrifice measured hot-path performance for cosmetic elegance, abstraction purity, reduced source size, or a superficially
  prettier implementation. Prefer the fastest validated design that remains correct and auditable.

Every performance claim must state the workload, dataset shape, result cardinality/checksum, compiler, relevant CPU/backend, thread
count, allocator, and whether the host was idle. Use warmups and multiple alternating runs. Pin workers and control frequency when making
capacity or small-delta claims. Never present sanitizer, debug, or CPU-contended timings as production throughput.

## Production features and representative workloads

- Never add a production API, index type, storage path, planner rule, or special case solely to satisfy a test, sample, or benchmark.
  Implement the domain capability as a reusable, complete abstraction first; tests and benchmarks must consume that production path.
- Secondary indexing must use one typed architecture with shared catalog, encoding, maintenance, recovery, statistics, and planning.
  Boolean, signed/unsigned numeric, floating-point, datetime, and string indexes may have type-specific physical kernels, but must not
  become unrelated one-off subsystems created for a particular sample field.
- A workload that models independent actors must preserve independent operation semantics. Each actor schedules and sends its own
  request with realistic jitter; client-side batching must not merge those logical operations unless the modeled producer explicitly
  batches them in production.
- Internal database group commit is allowed and encouraged: it may coalesce already-independent requests after they reach the server,
  while preserving per-request results, ordering guarantees, durability, and linearization semantics.
- Do not equate one actor with one operating-system thread. Use event loops, timing wheels, coroutines, or sharded schedulers to represent
  very high actor cardinality without changing the number or independence of logical requests.

## SIMD and architecture-specific code

- Maintain a scalar reference path with equivalent semantics.
- x86 runtime dispatch must verify every ISA feature used by a target function and its tail calls. Do not execute AVX2, FMA, AVX-512, or
  BMI instructions based only on compile-time host capabilities.
- ARM64 NEON is the portable AArch64 baseline. SVE/SVE2 code must use runtime capability detection and vector-length-agnostic loops.
- SIMD functions must accept the documented alignment, handle zero length and every tail length, and never read or write past buffers.
- Avoid AVX-512 automatically when downclocking, transition costs, or register pressure make AVX2 faster for the complete pipeline.
- Keep polynomial approximations and fused predicates within explicit error bounds. Validate global coordinates, boundary predicates,
  poles, the antimeridian, and antipodal cases against the scalar reference.
- Inline assembly is a last resort. Prefer intrinsics; isolate assembly by architecture; document clobbers and ABI assumptions; provide a
  scalar/intrinsic fallback; and retain it only after end-to-end measurement.
- Inspect generated assembly when a change depends on vectorization, instruction selection, branch removal, or memory-access shape.

## Concurrency and NUMA

- Use C11 atomics or pthread synchronization with documented ownership and publication rules. Use the weakest correct memory order, but
  prefer clarity over an unproven relaxation.
- QSBR/RCU publication must prevent new readers, wait for the active epoch, publish a complete immutable snapshot, and reclaim old state
  only after quiescence.
- Keep writer critical sections short. Perform expensive allocation, sorting, merging, and file I/O outside the reader gate when safe.
- Mutations must have a defined linearization point. Insert, remove, reinsertion, compaction, reopen, and crash recovery must agree on
  last-operation-wins ordering.
- Pinning a thread is not sufficient NUMA support. Replica/shard memory placement, first-touch ownership, CPU topology, and cross-node
  result aggregation must be considered together.
- Destruction requires external callers to be quiescent unless the API explicitly supplies lifetime management.
- Any change to shared state, worker scheduling, publication, or reclamation requires a meaningful concurrent test and a ThreadSanitizer
  run on the affected path.

## Persistence and crash consistency

- Version every persisted format and assert its layout. Prefer fixed-width integer fields; do not introduce new persisted `size_t`, raw
  pointers, compiler-dependent enums, or implicit padding.
- Validate magic, version, endian marker, sizes, counts, offsets, monotonic ordering, arithmetic overflow, and file boundaries before
  exposing an mmap view.
- Durable replacement follows: create a temporary file in the destination directory, write, flush, `fsync` the file, rename atomically,
  and `fsync` the parent directory.
- A manifest may commit only data and WAL bytes already durable. Recovery must ignore uncommitted tails and reject incomplete committed
  state.
- Long-running WALs need bounded checkpoint/rewrite behavior. Compaction must reclaim tombstones, obsolete versions, and stale metadata
  without resurrecting IDs.
- Do not silently accept corrupted sorted order or metadata that can produce incorrect query results. Checksums are for corruption
  detection, not a replacement for structural validation.
- mmap ownership and lifetime must remain explicit. Never retain pointers after unmap or reclaim a mapped snapshot while readers can
  still access it.

## Error handling and ownership

- Make ownership transfer explicit at allocation, publication, and cleanup boundaries.
- Check system-call, pthread, allocation, I/O, and serialization results. An intentionally ignored result needs a comment explaining why
  failure is harmless.
- Leave externally visible state unchanged when an operation fails before its linearization point.
- Avoid partial initialization hazards. Cleanup functions must accept partially initialized objects when constructors use them.
- Do not silently fall back when the fallback changes correctness, durability, concurrency, or documented complexity.
- Public functions must define behavior for null pointers, zero counts, invalid coordinates, overflow, duplicate IDs, and read-only state.

## Tests and verification

Add tests for real failure modes and invariants, not for trivial getters, compiler behavior, or line coverage. A good test would fail for
a plausible regression and validates externally meaningful results.

Depending on the change, cover:

- exact results against a scalar or brute-force oracle;
- sparse, dense, globally uniform, and highly skewed datasets;
- zero/one/tail/vector-width counts and large cardinalities;
- poles, antimeridian crossing, antipodal points, invalid floating-point values, and exact boundaries;
- mmap reopen, malformed/truncated metadata, committed and uncommitted WAL tails, and compaction recovery;
- simultaneous insert/remove/reinsert with active readers, reopen, and physical reclamation;
- x86 runtime backends, ARM64 NEON, and the forced scalar backend;
- allocation-free reuse paths and exact capacity failures.

Before handing off a material C change, run the relevant subset and, when practical, the full matrix:

```bash
make clean
make -j
make test
make clean && make CC=gcc -j && make test CC=gcc
make clean && make test-scalar -j
make clean && ASAN_OPTIONS=detect_leaks=0 make test-debug -j
```

Run ThreadSanitizer separately for concurrency changes. LeakSanitizer may be disabled only when the execution environment uses ptrace
and LSan reports its known fatal incompatibility; ASan and UBSan must remain enabled. Also run the Clang static analyzer on changed C
translation units, `git diff --check`, and a line-length check against 140 columns.

Cross-compile ARM64 changes with a valid AArch64 sysroot or run them on ARM64 hardware. A host-only compile that mixes x86 headers with
an AArch64 target is not valid evidence.

## Change discipline

- Inspect existing code and dirty worktree state before editing. Preserve unrelated user changes.
- Keep each change scoped and reviewable. Do not combine broad cleanup with an unrelated optimization.
- Documentation is a hard deliverable, not deferred cleanup. Keep `docs/`, public-header contracts, README examples, persisted/wire format
  versions, operational procedures, benchmark methodology, and sample commands synchronized with every behavior or architecture change.
- Update README examples, format versions, comments, benchmark numbers, and test counts when behavior changes.
- Report what was measured, what merely compiled, what was not available on the host, and any remaining risk.
- Do not claim “maximum performance” from microbenchmarks alone. State the next limiting resource and the workload for which the result
  applies.

## Explicitly prohibited

- Condensed, difficult-to-read code.
- Duplicated algorithms, validation rules, ownership logic, or business semantics. Centralize them with zero-cost helpers, immutable
  vtables, parameterized kernels, or controlled generation; keep only irreducible ISA-specific instruction sequences separate.
- Stub functions, empty methods, no-op implementations, hard-coded success values, or deliberately incomplete behavior added only to
  compile or make a test pass. New APIs require a complete production implementation, error handling, ownership semantics, and
  meaningful behavioral coverage in the same change.
- `-ffast-math` or equivalent unsafe math modes.
- Packed hot records or unaligned mmap traversal for cosmetic space savings.
- `volatile` as thread synchronization.
- Per-query worker creation or avoidable per-query heap allocation in production hot paths.
- Unchecked persisted offsets, counts, sizes, or arithmetic.
- ISA-specific execution without correct runtime/compile-time guards and a portable fallback.
- Performance claims from a single run, a busy host, unequal workloads, or mismatched result checksums.
- Tests that exist only to increase test count and do not protect a meaningful invariant.
