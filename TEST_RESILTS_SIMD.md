==========================================
Running optimized tests...
==========================================
========================================
GEO INDEX TEST SUITE
========================================

--- MORTON CODE TESTS ---

[TEST 1] test_spread_compact_bits_roundtrip
  PASSED (0.022 ms)

[TEST 2] test_spread_bits_pattern
  PASSED (0.000 ms)

--- ENCODE/DECODE TESTS ---

[TEST 3] test_encode_decode_roundtrip
  PASSED (0.001 ms)

[TEST 4] test_encode_decode_edge_cases
  PASSED (0.001 ms)

[TEST 5] test_encode_clamping
  PASSED (0.000 ms)

[TEST 6] test_encode_ordering
    Z-diff near: 770360913, Z-diff far: 609689249138401650
  PASSED (0.001 ms)

--- DISTANCE TESTS ---

[TEST 7] test_haversine_known_distances
    Sao Paulo to Rio de Janeiro: 360.7 km (expected ~357.0 km)
    New York to Paris: 5837.2 km (expected ~5837.0 km)
    London to Moscow: 2500.5 km (expected ~2500.0 km)
    Tokyo to Singapore: 5311.2 km (expected ~5312.0 km)
  PASSED (0.006 ms)

[TEST 8] test_haversine_zero_distance
  PASSED (0.000 ms)

[TEST 9] test_haversine_symmetry
  PASSED (0.003 ms)

[TEST 10] test_haversine_triangle_inequality
  PASSED (0.004 ms)

[TEST 11] test_fast_distance_accuracy
    Precise: 1.2579 km, Fast: 1.2593 km, Error: 0.11%
  PASSED (0.001 ms)

--- INDEX TESTS ---

[TEST 12] test_index_create_destroy
  PASSED (0.009 ms)

[TEST 13] test_index_add_single
  PASSED (0.004 ms)

[TEST 14] test_index_add_batch
  PASSED (0.001 ms)

[TEST 15] test_index_auto_grow
  PASSED (0.016 ms)

[TEST 16] test_index_build_sorts
  PASSED (0.002 ms)

--- BINARY SEARCH TESTS ---

[TEST 17] test_lower_bound_basic
  PASSED (0.002 ms)

--- SEARCH TESTS ---

[TEST 18] test_search_radius_basic
    Found 2 results, scanned 2 records in 0.001 ms
  PASSED (0.004 ms)

[TEST 19] test_search_radius_finds_nearby
  PASSED (0.001 ms)

[TEST 20] test_search_radius_empty_result
  PASSED (0.001 ms)

[TEST 21] test_search_bbox_basic
    Found 3 results in South America bbox
  PASSED (0.002 ms)

[TEST 22] test_search_knn_basic
    Found 3 nearest neighbors
  PASSED (0.003 ms)

--- PRECISION TESTS ---

[TEST 23] test_precision_at_equator
    Error at equator: 0.0052 meters
  PASSED (0.001 ms)

[TEST 24] test_precision_at_poles
    Error near north pole: 0.0042 meters
  PASSED (0.001 ms)

[TEST 25] test_precision_random_points
    Random points: avg error = 0.0041 m, max error = 0.0096 m
  PASSED (0.948 ms)

--- EDGE CASE TESTS ---

[TEST 26] test_edge_antimeridian
    Points near antimeridian: 2
  PASSED (0.001 ms)

[TEST 27] test_edge_poles
    Points near north pole: 4
  PASSED (0.001 ms)

[TEST 28] test_edge_empty_index
  PASSED (0.001 ms)

[TEST 29] test_edge_single_point
  PASSED (0.008 ms)

--- VALIDATION TESTS ---

[TEST 30] test_validation_coordinates
  PASSED (0.000 ms)

[TEST 31] test_validation_clamp
  PASSED (0.000 ms)

[TEST 32] test_validation_wrap_lng
  PASSED (0.000 ms)

--- PERFORMANCE TESTS ---

[TEST 33] test_perf_encode_decode
    Encode: 1000000 ops in 6.91 ms (144.71 M ops/sec)
    Decode: 1000000 ops in 2.78 ms (359.93 M ops/sec)
  PASSED (9.694 ms)

[TEST 34] test_perf_index_build
    Build 1000 records: 0.09 ms
    Build 10000 records: 1.12 ms
    Build 100000 records: 11.00 ms
    Build 1000000 records: 75.87 ms
  PASSED (107.693 ms)

[TEST 35] test_perf_search_radius
    Radius 1 km: avg 0.001 ms, avg results 0.0
    Radius 10 km: avg 0.010 ms, avg results 2.6
    Radius 100 km: avg 0.034 ms, avg results 161.8
    Radius 1000 km: avg 0.331 ms, avg results 9381.4
  PASSED (112.108 ms)

[TEST 36] test_perf_binary_search
    1000000 binary searches in 103.27 ms (9.68 M ops/sec)
  PASSED (987.488 ms)

--- STRESS TESTS ---

[TEST 37] test_stress_large_dataset
    Creating index with 5000000 records...
    ┌─────────────────────────────────────────────────────────┐
    │ LARGE DATASET TEST (5000000 records)                    │
    ├─────────────────────────────────────────────────────────┤
    │ Scalar Insert:    21.93 ms (228.01 M ops/sec)           │
    │ SIMD Encode:       0.00 ms (   inf M ops/sec)           │
    │ Encode Speedup:    infx                                 │
    │ Build time:      373.92 ms                              │
    ├─────────────────────────────────────────────────────────┤
    │ Avg search (50km): 0.086 ms                             │
    └─────────────────────────────────────────────────────────┘
  PASSED (537.868 ms)

[TEST 38] test_stress_many_searches
    ┌─────────────────────────────────────────────────────────┐
    │ MANY SEARCHES TEST (10000 searches)                     │
    ├─────────────────────────────────────────────────────────┤
    │ Scalar Search:    31.51 ms (  317355 searches/sec)      │
    │ SIMD Distance:     0.00 ms (     inf batches/sec)       │
    │ Distance Speedup:   infx                                │
    └─────────────────────────────────────────────────────────┘
  PASSED (38.244 ms)

[TEST 39] test_stress_dense_area
    ┌─────────────────────────────────────────────────────────┐
    │ DENSE AREA TEST (100000 points, 5.0km radius)           │
    ├─────────────────────────────────────────────────────────┤
    │ Scalar Search:    1.135 ms, found 69416 points          │
    │ SIMD Filter:      0.167 ms, found 0 points              │
    │ Filter Speedup:   6.82x                                 │
    │ SIMD Haversine:   0.217 ms, found 67082 points          │
    │ Haversine Speedup:  5.23x                               │
    └─────────────────────────────────────────────────────────┘
  PASSED (8.036 ms)

--- SIMD TESTS ---
SIMD Implementation: ARM64 NEON

[TEST 40] test_simd_available
    SIMD Available: YES
    SIMD Implementation: ARM64 NEON
    Optimal Batch Size: 256
  PASSED (0.001 ms)

[TEST 41] test_simd_encode_correctness
    Compared 1000 encodings, 0 mismatches
  PASSED (0.010 ms)

[TEST 42] test_simd_decode_correctness
    Max lat diff: 0.0000000000, max lng diff: 0.0000000000
  PASSED (0.006 ms)

[TEST 43] test_simd_haversine_accuracy
    Max absolute diff: 11.6177 km
    Max relative error: 1.81%
  PASSED (0.002 ms)

--- SIMD BENCHMARKS (Scalar vs SIMD) ---

[TEST 44] test_simd_benchmark_encode
    Running SIMD encode benchmark...
    ┌─────────────────────────────────────────────────────────┐
    │ ENCODE BENCHMARK (100K points x 10 iterations)          │
    ├─────────────────────────────────────────────────────────┤
    │ Scalar:      2.01 ms  (496.74 M ops/sec)                │
    │ SIMD:        1.39 ms  (719.60 M ops/sec)                │
    │ Speedup:     1.45x                                      │
    └─────────────────────────────────────────────────────────┘
  PASSED (3.473 ms)

[TEST 45] test_simd_benchmark_decode
    Running SIMD decode benchmark...
    ┌─────────────────────────────────────────────────────────┐
    │ DECODE BENCHMARK (100K codes x 10 iterations)           │
    ├─────────────────────────────────────────────────────────┤
    │ Scalar:      1.12 ms  (891.80 M ops/sec)                │
    │ SIMD:        0.87 ms  (1144.87 M ops/sec)               │
    │ Speedup:     1.28x                                      │
    └─────────────────────────────────────────────────────────┘
  PASSED (2.006 ms)

[TEST 46] test_simd_benchmark_haversine
    Running SIMD haversine benchmark...
    ┌─────────────────────────────────────────────────────────┐
    │ HAVERSINE BENCHMARK (100K distances x 10 iterations)    │
    ├─────────────────────────────────────────────────────────┤
    │ Scalar:     10.83 ms  ( 92.35 M ops/sec)                │
    │ SIMD:        2.04 ms  (489.16 M ops/sec)                │
    │ Speedup:     5.30x                                      │
    └─────────────────────────────────────────────────────────┘
  PASSED (12.945 ms)

[TEST 47] test_simd_benchmark_filter_radius
    Running SIMD filter_radius benchmark...
    ┌─────────────────────────────────────────────────────────┐
    │ FILTER RADIUS BENCHMARK (100K points x 5 iterations)    │
    ├─────────────────────────────────────────────────────────┤
    │ Scalar:      6.80 ms  ( 73.58 M ops/sec)                │
    │ SIMD:        0.16 ms  (3052.65 M ops/sec)               │
    │ Speedup:    41.49x                                      │
    └─────────────────────────────────────────────────────────┘
  PASSED (7.005 ms)

[TEST 48] test_simd_comprehensive_benchmark

    ╔═════════════════════════════════════════════════════════╗
    ║     COMPREHENSIVE SIMD vs SCALAR BENCHMARK SUMMARY      ║
    ╠═════════════════════════════════════════════════════════╣
    ║ Architecture: ARM64 NEON                                ║
    ╠═══════════════╦═══════════════╦═══════════════╦═════════╣
    ║   Operation   ║  Scalar (ms)  ║   SIMD (ms)   ║ Speedup ║
    ╠═══════════════╬═══════════════╬═══════════════╬═════════╣
    ║ Encode        ║          2.06 ║          1.51 ║   1.37x ║
    ║ Decode        ║          1.15 ║          0.88 ║   1.31x ║
    ║ Haversine     ║         10.83 ║          1.97 ║   5.51x ║
    ║ Filter Radius ║          6.52 ║          0.24 ║  27.33x ║
    ╚═══════════════╩═══════════════╩═══════════════╩═════════╝

    Average Speedup: 8.88x
  PASSED (25.247 ms)

========================================
TEST SUMMARY
========================================
Total:  48
Passed: 48
Failed: 0
========================================
