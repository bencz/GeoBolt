==========================================
Running optimized tests...
==========================================
========================================
GEO INDEX TEST SUITE
========================================

--- MORTON CODE TESTS ---

[TEST 1] test_spread_compact_bits_roundtrip
  PASSED

[TEST 2] test_spread_bits_pattern
  PASSED

--- ENCODE/DECODE TESTS ---

[TEST 3] test_encode_decode_roundtrip
  PASSED

[TEST 4] test_encode_decode_edge_cases
  PASSED

[TEST 5] test_encode_clamping
  PASSED

[TEST 6] test_encode_ordering
    Z-diff near: 770360913, Z-diff far: 609689249138401650
  PASSED

--- DISTANCE TESTS ---

[TEST 7] test_haversine_known_distances
    Sao Paulo to Rio de Janeiro: 360.7 km (expected ~357.0 km)
    New York to Paris: 5837.2 km (expected ~5837.0 km)
    London to Moscow: 2500.5 km (expected ~2500.0 km)
    Tokyo to Singapore: 5311.2 km (expected ~5312.0 km)
  PASSED

[TEST 8] test_haversine_zero_distance
  PASSED

[TEST 9] test_haversine_symmetry
  PASSED

[TEST 10] test_haversine_triangle_inequality
  PASSED

[TEST 11] test_fast_distance_accuracy
    Precise: 1.2579 km, Fast: 1.2593 km, Error: 0.11%
  PASSED

--- INDEX TESTS ---

[TEST 12] test_index_create_destroy
  PASSED

[TEST 13] test_index_add_single
  PASSED

[TEST 14] test_index_add_batch
  PASSED

[TEST 15] test_index_auto_grow
  PASSED

[TEST 16] test_index_build_sorts
  PASSED

--- BINARY SEARCH TESTS ---

[TEST 17] test_lower_bound_basic
  PASSED

--- SEARCH TESTS ---

[TEST 18] test_search_radius_basic
    Found 2 results, scanned 2 records in 0.001 ms
  PASSED

[TEST 19] test_search_radius_finds_nearby
  PASSED

[TEST 20] test_search_radius_empty_result
  PASSED

[TEST 21] test_search_bbox_basic
    Found 3 results in South America bbox
  PASSED

[TEST 22] test_search_knn_basic
    Found 3 nearest neighbors
  PASSED

--- PRECISION TESTS ---

[TEST 23] test_precision at equator
    Error at equator: 0.0052 meters
  PASSED

[TEST 24] test_precision_at_poles
    Error near north pole: 0.0042 meters
  PASSED

[TEST 25] test_precision_random_points
    Random points: avg error = 0.0041 m, max error = 0.0096 m
  PASSED

--- EDGE CASE TESTS ---

[TEST 26] test_edge_antimeridian
    Points near antimeridian: 2
  PASSED

[TEST 27] test_edge_poles
    Points near north pole: 4
  PASSED

[TEST 28] test_edge_empty_index
  PASSED

[TEST 29] test_edge_single_point
  PASSED

--- VALIDATION TESTS ---

[TEST 30] test_validation_coordinates
  PASSED

[TEST 31] test_validation_clamp
  PASSED

[TEST 32] test_validation_wrap_lng
  PASSED

--- PERFORMANCE TESTS ---

[TEST 33] test_perf_encode_decode
    Encode: 1000000 ops in 2.08 ms (480.03 M ops/sec)
    Decode: 1000000 ops in 0.65 ms (1543.81 M ops/sec)
  PASSED

[TEST 34] test_perf_index_build
    Build 1000 records: 0.05 ms
    Build 10000 records: 0.47 ms
    Build 100000 records: 6.37 ms
    Build 1000000 records: 67.81 ms
  PASSED

[TEST 35] test_perf_search_radius
    Radius 1 km: avg 0.001 ms, avg results 0.0
    Radius 10 km: avg 0.010 ms, avg results 2.6
    Radius 100 km: avg 0.037 ms, avg results 161.8
    Radius 1000 km: avg 0.331 ms, avg results 9381.4
  PASSED

[TEST 36] test_perf_binary_search
    1000000 binary searches in 126.99 ms (7.87 M ops/sec)
  PASSED

--- STRESS TESTS ---

[TEST 37] test_stress_large_dataset
    Creating index with 5000000 records...
    Insert time: 56.19 ms (88.98 M ops/sec)
    Build time: 375.29 ms
    Avg search time (50km radius): 0.086 ms
  PASSED

[TEST 38] test_stress_many_searches
    10000 searches in 31.91 ms (313357.65 searches/sec)
  PASSED

[TEST 39] test_stress_dense_area
    Dense area: found 69416 of 100000 in 1.193 ms
  PASSED

========================================
TEST SUMMARY
========================================
Total:  39
Passed: 39
Failed: 0
========================================
