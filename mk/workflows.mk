# Run demo benchmark
demo: dirs $(DEMO_BIN)
	@echo "=========================================="
	@echo "Running comprehensive benchmark demo..."
	@echo "=========================================="
	@$(DEMO_BIN)

# Run tests (optimized)
test: dirs $(TEST_BIN) $(DATABASE_TEST_BIN) $(SERVER_TEST_BIN) $(DAEMON_TEST_BIN)
	@echo "=========================================="
	@echo "Running optimized tests..."
	@echo "=========================================="
	@perl -e 'alarm 120; exec @ARGV' $(TEST_BIN) $(TEST_FILTER) || (echo "Test timed out or failed!" && exit 1)
	@$(DATABASE_TEST_BIN)
	@$(SERVER_TEST_BIN)
	@$(DAEMON_TEST_BIN)

# Run tests (debug with sanitizers)
test-debug: dirs $(TEST_BIN_DEBUG) $(DATABASE_TEST_BIN_DEBUG) $(SERVER_TEST_BIN_DEBUG) $(DAEMON_TEST_BIN_DEBUG)
	@echo "=========================================="
	@echo "Running debug tests with sanitizers..."
	@echo "=========================================="
	@$(TEST_BIN_DEBUG) $(TEST_FILTER)
	@$(DATABASE_TEST_BIN_DEBUG)
	@$(SERVER_TEST_BIN_DEBUG)
	@GEOBOLTD_TEST_BINARY=./$(DAEMON_BIN_DEBUG) $(DAEMON_TEST_BIN_DEBUG)

test-scalar: dirs $(TEST_BIN_SCALAR) $(DATABASE_TEST_BIN_SCALAR) $(SERVER_TEST_BIN_SCALAR)
	@echo "=========================================="
	@echo "Running portable scalar-backend tests..."
	@echo "=========================================="
	@$(TEST_BIN_SCALAR) $(TEST_FILTER)
	@$(DATABASE_TEST_BIN_SCALAR)
	@$(SERVER_TEST_BIN_SCALAR)

# Build and run the public-API sample
sample: dirs $(SAMPLE_BIN)
	@$(SAMPLE_BIN)

client-sample: dirs $(CLIENT_SAMPLE_BIN)
	@echo "Run with: GEOBOLT_TOKEN=secret $(CLIENT_SAMPLE_BIN) HOST PORT"

# Quick benchmark
benchmark: dirs $(TEST_BIN)
	@echo "=========================================="
	@echo "Running performance benchmarks..."
	@echo "=========================================="
	@$(TEST_BIN) 2>&1 | grep -A 100 "PERFORMANCE TESTS"

benchmark-tombstones: dirs $(TOMBSTONE_BENCH_BIN)
	@$(TOMBSTONE_BENCH_BIN) $(BENCHMARK_ARGS)

benchmark-stream: dirs $(STREAM_BENCH_BIN)
	@$(STREAM_BENCH_BIN) $(BENCHMARK_ARGS)

benchmark-compaction: dirs $(COMPACTION_BENCH_BIN)
	@$(COMPACTION_BENCH_BIN) $(BENCHMARK_ARGS)

benchmark-server: dirs $(SERVER_BENCH_BIN)
	@$(SERVER_BENCH_BIN) $(BENCHMARK_ARGS)

soak-uber: dirs $(SOAK_BIN)
	@$(SOAK_BIN) $(SOAK_ARGS)

soak-uber-smoke: dirs $(SOAK_BIN)
	@$(SOAK_BIN) $(SOAK_SMOKE_ARGS)

soak-uber-tsan:
	@$(MAKE) clean
	@$(MAKE) USE_JEMALLOC=0 CFLAGS_OPT="-g -O1 -fno-omit-frame-pointer -fsanitize=thread" \
		LDFLAGS="-fsanitize=thread" $(TEST_BIN) $(DATABASE_TEST_BIN) $(SERVER_TEST_BIN) $(DAEMON_TEST_BIN) $(SOAK_BIN)
	@TSAN_OPTIONS=halt_on_error=1 $(TEST_BIN) writes_during_compaction
	@TSAN_OPTIONS=halt_on_error=1 $(DATABASE_TEST_BIN)
	@TSAN_OPTIONS=halt_on_error=1 $(SERVER_TEST_BIN)
	@TSAN_OPTIONS=halt_on_error=1 $(DAEMON_TEST_BIN)
	@TSAN_OPTIONS=halt_on_error=1 $(SOAK_BIN) $(SOAK_TSAN_ARGS)

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

-include $(BUILD_DIR)/*.d

# SIMD benchmark only
simd-benchmark: dirs $(TEST_BIN)
	@echo "=========================================="
	@echo "Running SIMD benchmarks..."
	@echo "Architecture: $(UNAME_M)"
	@echo "SIMD Flags: $(SIMD_FLAGS)"
	@echo "=========================================="
	@$(TEST_BIN) 2>&1 | grep -A 200 "SIMD"

allocator-info:
	@echo "Allocator:          $(ALLOCATOR_NAME)"
	@echo "USE_JEMALLOC:       $(USE_JEMALLOC)"
	@echo "jemalloc available: $(JEMALLOC_AVAILABLE)"
	@echo "jemalloc libraries: $(JEMALLOC_LIBS)"

pgo-generate:
	@mkdir -p $(PGO_DIR)
	@$(MAKE) clean
	@$(MAKE) CFLAGS_OPT="-O3 -flto -fprofile-instr-generate" LDFLAGS="$(LDFLAGS) -fprofile-instr-generate" all
	@echo "Run representative workloads with LLVM_PROFILE_FILE='$(PGO_RAW)', then execute 'make pgo-merge'."

pgo-merge:
	@command -v llvm-profdata >/dev/null 2>&1 || (echo "llvm-profdata is required" && exit 1)
	@llvm-profdata merge -output=$(PGO_DATA) $(PGO_DIR)/*.profraw

pgo-use:
	@test -f $(PGO_DATA) || (echo "Missing $(PGO_DATA); run pgo-generate, train, and pgo-merge first" && exit 1)
	@$(MAKE) clean
	@$(MAKE) CFLAGS_OPT="-O3 -flto -fprofile-instr-use=$(PGO_DATA)" all

# Help
help:
	@echo "GeoBolt Build System"
	@echo "====================="
	@echo ""
	@echo "Architecture: $(UNAME_M)"
	@echo "SIMD Source:  $(SIMD_SRC)"
	@echo "SIMD Flags:   $(SIMD_FLAGS)"
	@echo "Allocator:    $(ALLOCATOR_NAME)"
	@echo ""
	@echo "Targets:"
	@echo "  all            - Build library, tests, sample, and demo (default)"
	@echo "  lib            - Build static library only"
	@echo "  bin/geoboltd   - Build the single-server database daemon"
	@echo "  test           - Build and run optimized tests"
	@echo "  test-debug     - Build and run tests with sanitizers"
	@echo "  test-scalar    - Build and run the portable scalar backend"
	@echo "  sample         - Build and run the basic public-API example"
	@echo "  client-sample  - Build the remote C-driver example"
	@echo "  demo           - Run comprehensive benchmark demo (10M points)"
	@echo "  benchmark      - Run performance benchmarks"
	@echo "  benchmark-tombstones - Measure mutation-filter and post-compaction fast paths"
	@echo "  benchmark-stream - Measure radix chunks, hierarchical merges, and final publication"
	@echo "  soak-uber       - Run the configurable Uber-like concurrent actor soak (default: 5 minutes)"
	@echo "  soak-uber-smoke - Run a 10-second functional soak"
	@echo "  soak-uber-tsan  - Rebuild with ThreadSanitizer and run the short concurrent soak"
	@echo "  benchmark-compaction - Compare serial and partitioned Morton compaction"
	@echo "  benchmark-server - Measure persistent-connection protocol and query throughput"
	@echo "  simd-benchmark - Run SIMD-specific benchmarks"
	@echo "  allocator-info - Show allocator auto-detection details"
	@echo "  pgo-generate    - Build instrumented binaries for representative training"
	@echo "  pgo-merge       - Merge collected Clang profiles"
	@echo "  pgo-use         - Build optimized binaries using the merged profile"
	@echo "  clean          - Remove build artifacts"
	@echo "  help           - Show this help message"
	@echo ""
	@echo "Build Flags:"
	@echo "  Optimized:  $(CFLAGS_OPT)"
	@echo "  Debug:      $(CFLAGS_DEBUG)"
	@echo "  SIMD:       $(SIMD_FLAGS)"
	@echo "  Allocator:  $(ALLOCATOR_NAME)"
	@echo "  Test filter: $(TEST_FILTER)"
