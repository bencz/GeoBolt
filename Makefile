# GeoIndex Makefile
# Production-grade build system

CC = clang
AR ?= ar
CFLAGS ?= -Wall -Wextra -Werror -std=c11 -pedantic
CPPFLAGS ?= -D_POSIX_C_SOURCE=200809L -I.
CFLAGS_OPT ?= -O3 -flto
CFLAGS_DEBUG ?= -g -O1 -DDEBUG -fno-omit-frame-pointer -fsanitize=address,undefined
LDFLAGS ?=
LDLIBS ?= -lm
THREAD_LIBS ?= -pthread
USE_JEMALLOC ?= auto
TEST_FILTER ?=
PGO_DIR ?= profiles
PGO_RAW ?= $(PGO_DIR)/geobolt-%p.profraw
PGO_DATA ?= $(PGO_DIR)/geobolt.profdata
BENCHMARK_CPPFLAGS = -DGEO_BENCHMARK_ALLOCATOR=\"$(ALLOCATOR_NAME)\"
SOAK_ARGS ?=
SOAK_SMOKE_ARGS ?= 10 10000 2 4 512 10 250 16 2 2
SOAK_TSAN_ARGS ?= 15 10000 2 4 512 10 250 16 5 2

# jemalloc is enabled for optimized executables when it can be linked. Sanitizer
# builds deliberately keep their own allocator/interceptors. USE_JEMALLOC=0
# disables it for controlled comparisons; USE_JEMALLOC=1 forces -ljemalloc.
JEMALLOC_CANDIDATE_LIBS := $(shell \
	if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists jemalloc; then \
		pkg-config --libs jemalloc; \
	else \
		echo -ljemalloc; \
	fi)

JEMALLOC_AVAILABLE := $(shell \
	if printf 'int main(void) { return 0; }\n' | \
		CCACHE_DISABLE=1 $(CC) -x c - $(JEMALLOC_CANDIDATE_LIBS) -o /dev/null >/dev/null 2>&1; then \
		echo 1; \
	else \
		echo 0; \
	fi)

ifeq ($(USE_JEMALLOC),auto)
    JEMALLOC_ENABLED := $(JEMALLOC_AVAILABLE)
else ifeq ($(USE_JEMALLOC),1)
    JEMALLOC_ENABLED := 1
else ifeq ($(USE_JEMALLOC),0)
    JEMALLOC_ENABLED := 0
else
    $(error USE_JEMALLOC must be auto, 0, or 1)
endif

ifeq ($(JEMALLOC_ENABLED),1)
    JEMALLOC_LIBS := $(JEMALLOC_CANDIDATE_LIBS)
    OPT_LDLIBS := $(LDLIBS) $(JEMALLOC_LIBS)
    ALLOCATOR_NAME := jemalloc
else
    JEMALLOC_LIBS :=
    OPT_LDLIBS := $(LDLIBS)
    ALLOCATOR_NAME := system
endif

# Directories
SRC_DIR = .
BUILD_DIR = build
BIN_DIR = bin
TEST_DIR = tests
BENCHMARK_DIR = benchmarks
SAMPLE_DIR = samples

# Source files
LIB_SRC = geo_index.c
DENSITY_SRC = geo_index_density.c
STREAM_SRC = geo_index_stream.c
PARALLEL_SRC = geo_index_parallel.c
BATCH_SRC = geo_index_batch.c
SEGMENTS_SRC = geo_index_segments.c
IO_SRC = geo_index_io.c
SIMD_SRC_ARM64 = geo_index_simd_arm64.c
SIMD_SRC_X86 = geo_index_simd_x86.c
SIMD_SRC_SCALAR = geo_index_simd_scalar.c

# Object files
LIB_OBJ = $(BUILD_DIR)/geo_index.o
DENSITY_OBJ = $(BUILD_DIR)/geo_index_density.o
STREAM_OBJ = $(BUILD_DIR)/geo_index_stream.o
PARALLEL_OBJ = $(BUILD_DIR)/geo_index_parallel.o
BATCH_OBJ = $(BUILD_DIR)/geo_index_batch.o
SEGMENTS_OBJ = $(BUILD_DIR)/geo_index_segments.o
IO_OBJ = $(BUILD_DIR)/geo_index_io.o
LIB_OBJ_DEBUG = $(BUILD_DIR)/geo_index_debug.o
DENSITY_OBJ_DEBUG = $(BUILD_DIR)/geo_index_density_debug.o
STREAM_OBJ_DEBUG = $(BUILD_DIR)/geo_index_stream_debug.o
PARALLEL_OBJ_DEBUG = $(BUILD_DIR)/geo_index_parallel_debug.o
BATCH_OBJ_DEBUG = $(BUILD_DIR)/geo_index_batch_debug.o
SEGMENTS_OBJ_DEBUG = $(BUILD_DIR)/geo_index_segments_debug.o
IO_OBJ_DEBUG = $(BUILD_DIR)/geo_index_io_debug.o
LIB_OBJ_SCALAR = $(BUILD_DIR)/geo_index_scalar.o
DENSITY_OBJ_SCALAR = $(BUILD_DIR)/geo_index_density_scalar.o
STREAM_OBJ_SCALAR = $(BUILD_DIR)/geo_index_stream_scalar.o
PARALLEL_OBJ_SCALAR = $(BUILD_DIR)/geo_index_parallel_scalar.o
BATCH_OBJ_SCALAR = $(BUILD_DIR)/geo_index_batch_scalar.o
SEGMENTS_OBJ_SCALAR = $(BUILD_DIR)/geo_index_segments_scalar.o
IO_OBJ_SCALAR = $(BUILD_DIR)/geo_index_io_scalar.o
SCALAR_BACKEND_OBJ = $(BUILD_DIR)/geo_index_simd_forced_scalar.o

# SIMD object files (architecture-specific)
SIMD_OBJ_ARM64 = $(BUILD_DIR)/geo_index_simd_arm64.o
SIMD_OBJ_X86 = $(BUILD_DIR)/geo_index_simd_x86.o
SIMD_OBJ_SCALAR = $(BUILD_DIR)/geo_index_simd_scalar.o

# Detect architecture
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_M),arm64)
    SIMD_OBJ = $(SIMD_OBJ_ARM64)
    SIMD_SRC = $(SIMD_SRC_ARM64)
    SIMD_FLAGS =
else ifeq ($(UNAME_M),aarch64)
    SIMD_OBJ = $(SIMD_OBJ_ARM64)
    SIMD_SRC = $(SIMD_SRC_ARM64)
    SIMD_FLAGS =
else ifeq ($(UNAME_M),x86_64)
    SIMD_OBJ = $(SIMD_OBJ_X86)
    SIMD_SRC = $(SIMD_SRC_X86)
    # x86 runtime dispatch is implemented inside the SIMD translation unit.
    SIMD_FLAGS =
else
    SIMD_OBJ = $(SIMD_OBJ_SCALAR)
    SIMD_SRC = $(SIMD_SRC_SCALAR)
    SIMD_FLAGS =
endif

TEST_SRC = $(TEST_DIR)/test_geo_index.c
SAMPLE_SRC = $(SAMPLE_DIR)/basic_usage.c
DEMO_SRC = $(BENCHMARK_DIR)/demo_benchmark.c
TOMBSTONE_BENCH_SRC = $(BENCHMARK_DIR)/tombstone_visibility.c
STREAM_BENCH_SRC = $(BENCHMARK_DIR)/stream_build.c
COMPACTION_BENCH_SRC = $(BENCHMARK_DIR)/compaction_merge.c
SOAK_SRC = $(TEST_DIR)/soak_uber.c

# Targets
LIB_STATIC = $(BUILD_DIR)/libgeoindex.a
TEST_BIN = $(BIN_DIR)/test_geo_index
TEST_BIN_DEBUG = $(BIN_DIR)/test_geo_index_debug
TEST_BIN_SCALAR = $(BIN_DIR)/test_geo_index_scalar
SAMPLE_BIN = $(BIN_DIR)/basic_usage
DEMO_BIN = $(BIN_DIR)/demo_benchmark
BENCH_BIN = $(BIN_DIR)/benchmark
TOMBSTONE_BENCH_BIN = $(BIN_DIR)/benchmark_tombstones
STREAM_BENCH_BIN = $(BIN_DIR)/benchmark_stream
COMPACTION_BENCH_BIN = $(BIN_DIR)/benchmark_compaction
SOAK_BIN = $(BIN_DIR)/soak_uber

DIRECT_BUILD_TARGETS = $(LIB_OBJ) $(DENSITY_OBJ) $(STREAM_OBJ) $(PARALLEL_OBJ) $(BATCH_OBJ) $(SEGMENTS_OBJ) $(IO_OBJ) $(SIMD_OBJ) \
	$(LIB_OBJ_DEBUG) $(DENSITY_OBJ_DEBUG) $(STREAM_OBJ_DEBUG) $(PARALLEL_OBJ_DEBUG) $(BATCH_OBJ_DEBUG) $(SEGMENTS_OBJ_DEBUG) \
	$(IO_OBJ_DEBUG) $(SIMD_OBJ_DEBUG) $(LIB_OBJ_SCALAR) $(DENSITY_OBJ_SCALAR) $(STREAM_OBJ_SCALAR) $(PARALLEL_OBJ_SCALAR) \
	$(BATCH_OBJ_SCALAR) $(SEGMENTS_OBJ_SCALAR) $(IO_OBJ_SCALAR) $(SCALAR_BACKEND_OBJ) $(LIB_STATIC) $(TEST_BIN) \
	$(TEST_BIN_DEBUG) $(TEST_BIN_SCALAR) $(SAMPLE_BIN) $(DEMO_BIN) $(TOMBSTONE_BENCH_BIN) $(STREAM_BENCH_BIN) \
	$(COMPACTION_BENCH_BIN) $(SOAK_BIN)

.PHONY: all clean test test-debug test-scalar lib sample demo benchmark benchmark-tombstones benchmark-stream benchmark-compaction \
	simd-benchmark \
	soak-uber soak-uber-smoke soak-uber-tsan allocator-info pgo-generate pgo-merge pgo-use help dirs FORCE_OPT_LINK

# Default target
all: dirs lib $(TEST_BIN) $(SAMPLE_BIN) $(DEMO_BIN)

# Create directories
dirs:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

$(DIRECT_BUILD_TARGETS): | dirs

# Build static library (optimized)
$(LIB_OBJ): $(SRC_DIR)/geo_index.c $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_private.h $(SRC_DIR)/geo_index_simd.h \
		$(SRC_DIR)/geo_index_persistence.h $(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(DENSITY_OBJ): $(SRC_DIR)/$(DENSITY_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_persistence.h $(SRC_DIR)/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(STREAM_OBJ): $(SRC_DIR)/$(STREAM_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_persistence.h $(SRC_DIR)/geo_index_private.h \
		$(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(PARALLEL_OBJ): $(SRC_DIR)/$(PARALLEL_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(BATCH_OBJ): $(SRC_DIR)/$(BATCH_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SEGMENTS_OBJ): $(SRC_DIR)/$(SEGMENTS_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_persistence.h $(SRC_DIR)/geo_index_private.h \
		$(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(IO_OBJ): $(SRC_DIR)/$(IO_SRC) $(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

# Build SIMD object file
$(SIMD_OBJ): $(SRC_DIR)/$(SIMD_SRC) $(SRC_DIR)/geo_index_simd.h $(SRC_DIR)/geo_index_simd_scalar_kernels.h $(SRC_DIR)/geo_index.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) $(SIMD_FLAGS) -MMD -MP -c $< -o $@

$(LIB_STATIC): $(LIB_OBJ) $(DENSITY_OBJ) $(STREAM_OBJ) $(PARALLEL_OBJ) $(BATCH_OBJ) $(SEGMENTS_OBJ) $(IO_OBJ) $(SIMD_OBJ)
	$(AR) rcs $@ $^

lib: dirs $(LIB_STATIC)

# Build static library (debug)
$(LIB_OBJ_DEBUG): $(SRC_DIR)/geo_index.c $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_private.h \
		$(SRC_DIR)/geo_index_persistence.h $(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(DENSITY_OBJ_DEBUG): $(SRC_DIR)/$(DENSITY_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_persistence.h \
		$(SRC_DIR)/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(STREAM_OBJ_DEBUG): $(SRC_DIR)/$(STREAM_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_persistence.h \
		$(SRC_DIR)/geo_index_private.h $(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(PARALLEL_OBJ_DEBUG): $(SRC_DIR)/$(PARALLEL_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(BATCH_OBJ_DEBUG): $(SRC_DIR)/$(BATCH_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(SEGMENTS_OBJ_DEBUG): $(SRC_DIR)/$(SEGMENTS_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_persistence.h \
		$(SRC_DIR)/geo_index_private.h $(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(IO_OBJ_DEBUG): $(SRC_DIR)/$(IO_SRC) $(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

SIMD_OBJ_DEBUG = $(BUILD_DIR)/geo_index_simd_debug.o

$(SIMD_OBJ_DEBUG): $(SRC_DIR)/$(SIMD_SRC) $(SRC_DIR)/geo_index_simd.h $(SRC_DIR)/geo_index_simd_scalar_kernels.h \
		$(SRC_DIR)/geo_index_internal.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) $(SIMD_FLAGS) -MMD -MP -c $< -o $@

$(LIB_OBJ_SCALAR): $(SRC_DIR)/geo_index.c $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_private.h $(SRC_DIR)/geo_index_simd.h \
		$(SRC_DIR)/geo_index_persistence.h $(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(DENSITY_OBJ_SCALAR): $(SRC_DIR)/$(DENSITY_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_persistence.h \
		$(SRC_DIR)/geo_index_private.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(STREAM_OBJ_SCALAR): $(SRC_DIR)/$(STREAM_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_persistence.h \
		$(SRC_DIR)/geo_index_private.h $(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(PARALLEL_OBJ_SCALAR): $(SRC_DIR)/$(PARALLEL_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_private.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(BATCH_OBJ_SCALAR): $(SRC_DIR)/$(BATCH_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_private.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SEGMENTS_OBJ_SCALAR): $(SRC_DIR)/$(SEGMENTS_SRC) $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_persistence.h \
		$(SRC_DIR)/geo_index_private.h $(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(IO_OBJ_SCALAR): $(SRC_DIR)/$(IO_SRC) $(SRC_DIR)/geo_index_io.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SCALAR_BACKEND_OBJ): $(SRC_DIR)/$(SIMD_SRC_SCALAR) $(SRC_DIR)/geo_index_simd.h $(SRC_DIR)/geo_index_simd_scalar_kernels.h \
		$(SRC_DIR)/geo_index_internal.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

# Build test binary (optimized)
$(TEST_BIN): $(TEST_SRC) $(LIB_STATIC) FORCE_OPT_LINK
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) $< -L$(BUILD_DIR) -lgeoindex $(LDFLAGS) $(OPT_LDLIBS) $(THREAD_LIBS) -o $@

# Build test binary (debug with sanitizers)
$(TEST_BIN_DEBUG): $(TEST_SRC) $(LIB_OBJ_DEBUG) $(DENSITY_OBJ_DEBUG) $(STREAM_OBJ_DEBUG) $(PARALLEL_OBJ_DEBUG) \
		$(BATCH_OBJ_DEBUG) \
		$(SEGMENTS_OBJ_DEBUG) $(IO_OBJ_DEBUG) $(SIMD_OBJ_DEBUG)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) $< \
		$(LIB_OBJ_DEBUG) $(DENSITY_OBJ_DEBUG) $(STREAM_OBJ_DEBUG) $(PARALLEL_OBJ_DEBUG) $(BATCH_OBJ_DEBUG) \
		$(SEGMENTS_OBJ_DEBUG) $(SIMD_OBJ_DEBUG) \
		$(IO_OBJ_DEBUG) \
		$(LDFLAGS) $(LDLIBS) $(THREAD_LIBS) -o $@

$(TEST_BIN_SCALAR): $(TEST_SRC) $(LIB_OBJ_SCALAR) $(DENSITY_OBJ_SCALAR) $(STREAM_OBJ_SCALAR) \
		$(PARALLEL_OBJ_SCALAR) $(BATCH_OBJ_SCALAR) \
		$(SEGMENTS_OBJ_SCALAR) $(IO_OBJ_SCALAR) $(SCALAR_BACKEND_OBJ) FORCE_OPT_LINK
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) $< \
		$(LIB_OBJ_SCALAR) $(DENSITY_OBJ_SCALAR) $(STREAM_OBJ_SCALAR) $(PARALLEL_OBJ_SCALAR) $(BATCH_OBJ_SCALAR) \
		$(SEGMENTS_OBJ_SCALAR) \
		$(IO_OBJ_SCALAR) \
		$(SCALAR_BACKEND_OBJ) \
		$(LDFLAGS) $(OPT_LDLIBS) $(THREAD_LIBS) -o $@

# Build the public-API sample
$(SAMPLE_BIN): $(SAMPLE_SRC) $(LIB_STATIC) FORCE_OPT_LINK
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) $< -L$(BUILD_DIR) -lgeoindex $(LDFLAGS) $(OPT_LDLIBS) $(THREAD_LIBS) -o $@

# Build demo benchmark
$(DEMO_BIN): $(DEMO_SRC) $(LIB_STATIC) FORCE_OPT_LINK
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) $< -L$(BUILD_DIR) -lgeoindex $(LDFLAGS) $(OPT_LDLIBS) $(THREAD_LIBS) -o $@

$(TOMBSTONE_BENCH_BIN): $(TOMBSTONE_BENCH_SRC) $(BENCHMARK_DIR)/benchmark_common.h $(LIB_STATIC) FORCE_OPT_LINK
	$(CC) $(CPPFLAGS) $(BENCHMARK_CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) $< -L$(BUILD_DIR) -lgeoindex \
		$(LDFLAGS) $(OPT_LDLIBS) $(THREAD_LIBS) -o $@

$(STREAM_BENCH_BIN): $(STREAM_BENCH_SRC) $(BENCHMARK_DIR)/benchmark_common.h $(LIB_STATIC) FORCE_OPT_LINK
	$(CC) $(CPPFLAGS) $(BENCHMARK_CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) $< -L$(BUILD_DIR) -lgeoindex \
		$(LDFLAGS) $(OPT_LDLIBS) $(THREAD_LIBS) -o $@

$(COMPACTION_BENCH_BIN): $(COMPACTION_BENCH_SRC) $(BENCHMARK_DIR)/benchmark_common.h $(LIB_STATIC) FORCE_OPT_LINK
	$(CC) $(CPPFLAGS) $(BENCHMARK_CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) $< -L$(BUILD_DIR) -lgeoindex \
		$(LDFLAGS) $(OPT_LDLIBS) $(THREAD_LIBS) -o $@

$(SOAK_BIN): $(SOAK_SRC) $(BENCHMARK_DIR)/benchmark_common.h $(LIB_STATIC) FORCE_OPT_LINK
	$(CC) $(CPPFLAGS) $(BENCHMARK_CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) $< -L$(BUILD_DIR) -lgeoindex \
		$(LDFLAGS) $(OPT_LDLIBS) $(THREAD_LIBS) -o $@

# Run demo benchmark
demo: dirs $(DEMO_BIN)
	@echo "=========================================="
	@echo "Running comprehensive benchmark demo..."
	@echo "=========================================="
	@$(DEMO_BIN)

# Run tests (optimized)
test: dirs $(TEST_BIN)
	@echo "=========================================="
	@echo "Running optimized tests..."
	@echo "=========================================="
	@perl -e 'alarm 120; exec @ARGV' $(TEST_BIN) $(TEST_FILTER) || (echo "Test timed out or failed!" && exit 1)

# Run tests (debug with sanitizers)
test-debug: dirs $(TEST_BIN_DEBUG)
	@echo "=========================================="
	@echo "Running debug tests with sanitizers..."
	@echo "=========================================="
	@$(TEST_BIN_DEBUG) $(TEST_FILTER)

test-scalar: dirs $(TEST_BIN_SCALAR)
	@echo "=========================================="
	@echo "Running portable scalar-backend tests..."
	@echo "=========================================="
	@$(TEST_BIN_SCALAR) $(TEST_FILTER)

# Build and run the public-API sample
sample: dirs $(SAMPLE_BIN)
	@$(SAMPLE_BIN)

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

soak-uber: dirs $(SOAK_BIN)
	@$(SOAK_BIN) $(SOAK_ARGS)

soak-uber-smoke: dirs $(SOAK_BIN)
	@$(SOAK_BIN) $(SOAK_SMOKE_ARGS)

soak-uber-tsan:
	@$(MAKE) clean
	@$(MAKE) USE_JEMALLOC=0 CFLAGS_OPT="-g -O1 -fno-omit-frame-pointer -fsanitize=thread" \
		LDFLAGS="-fsanitize=thread" $(TEST_BIN) $(SOAK_BIN)
	@TSAN_OPTIONS=halt_on_error=1 $(TEST_BIN) writes_during_compaction
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
	@echo "GeoIndex Build System"
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
	@echo "  test           - Build and run optimized tests"
	@echo "  test-debug     - Build and run tests with sanitizers"
	@echo "  test-scalar    - Build and run the portable scalar backend"
	@echo "  sample         - Build and run the basic public-API example"
	@echo "  demo           - Run comprehensive benchmark demo (10M points)"
	@echo "  benchmark      - Run performance benchmarks"
	@echo "  benchmark-tombstones - Measure mutation-filter and post-compaction fast paths"
	@echo "  benchmark-stream - Measure radix chunks, hierarchical merges, and final publication"
	@echo "  soak-uber       - Run the configurable Uber-like concurrent actor soak (default: 5 minutes)"
	@echo "  soak-uber-smoke - Run a 10-second functional soak"
	@echo "  soak-uber-tsan  - Rebuild with ThreadSanitizer and run the short concurrent soak"
	@echo "  benchmark-compaction - Compare serial and partitioned Morton compaction"
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
