ifeq ($(origin CC), default)
    CC = clang
endif
ifeq ($(origin CXX), default)
    ifeq ($(notdir $(CC)),gcc)
        CXX = g++
    else
        CXX = clang++
    endif
endif
AR ?= ar
C_WARNINGS = -Wall -Wextra -Wpedantic -Werror
CXX_WARNINGS = -Wall -Wextra -Wpedantic -Werror
CFLAGS ?= $(C_WARNINGS) -std=c11
CXXFLAGS ?= $(CXX_WARNINGS) -std=c++17
CPPFLAGS ?= -D_POSIX_C_SOURCE=200809L -Iinclude -Isrc/client -Isrc/core -Isrc/db -Isrc/db/object -Isrc/engine -Isrc/protocol -Isrc/query \
    -Isrc/server -Isrc/storage -Isrc/storage/rocksdb -Isrc/runtime -Isrc/simd
ifeq ($(notdir $(CC)),gcc)
    DEFAULT_LTO_FLAG = -flto=auto
else
    DEFAULT_LTO_FLAG = -flto
endif
CFLAGS_OPT ?= -O3 $(DEFAULT_LTO_FLAG)
CXXFLAGS_OPT ?= -O3 $(DEFAULT_LTO_FLAG)
CFLAGS_DEBUG ?= -g -O1 -DDEBUG -fno-omit-frame-pointer -fsanitize=address,undefined
CXXFLAGS_DEBUG ?= -g -O1 -DDEBUG -fno-omit-frame-pointer -fsanitize=address,undefined
LDFLAGS ?=
LDLIBS ?= -lm $(ROCKSDB_LIBS) $(CXX_RUNTIME_LIBS)
THREAD_LIBS ?= -pthread
ROCKSDB_LIBS ?= $(shell pkg-config --libs rocksdb 2>/dev/null || echo -lrocksdb)
CXX_RUNTIME_LIBS ?= -lstdc++
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
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin
TEST_DIR = tests
BENCHMARK_DIR = benchmarks
SAMPLE_DIR = samples

# Source files
LIB_SRC = $(SRC_DIR)/core/geo_index.c
DENSITY_SRC = $(SRC_DIR)/query/geo_index_density.c
STREAM_SRC = $(SRC_DIR)/storage/geo_index_stream.c
PARALLEL_SRC = $(SRC_DIR)/runtime/geo_index_parallel.c
BATCH_SRC = $(SRC_DIR)/query/geo_index_batch.c
SEGMENTS_SRC = $(SRC_DIR)/storage/geo_index_segments.c
IO_SRC = $(SRC_DIR)/storage/geo_index_io.c
THREAD_POOL_SRC = $(SRC_DIR)/runtime/geo_thread_pool.c
ROCKS_BRIDGE_SRC = $(SRC_DIR)/storage/rocksdb/geo_rocks_bridge.cc
GEODOC_SRC = $(SRC_DIR)/db/metadata/geo_doc.c
DB_FORMAT_SRC = $(SRC_DIR)/db/object/geo_db_format.c
DATABASE_SRC = $(SRC_DIR)/engine/geo_database.c
COMMIT_COORDINATOR_SRC = $(SRC_DIR)/engine/geo_commit_coordinator.c
SPATIAL_MEMTABLE_SRC = $(SRC_DIR)/engine/geo_spatial_memtable.c
OBJECT_CACHE_SRC = $(SRC_DIR)/engine/geo_object_cache.c
VISIBILITY_GATE_SRC = $(SRC_DIR)/engine/geo_visibility_gate.c
SECONDARY_INDEX_SRC = $(SRC_DIR)/engine/geo_secondary_index.c
SECONDARY_STATISTICS_SRC = $(SRC_DIR)/engine/geo_secondary_statistics.c
QUERY_PLANNER_SRC = $(SRC_DIR)/engine/geo_query_planner.c
PROTOCOL_SRC = $(SRC_DIR)/protocol/geo_protocol.c
JOB_QUEUE_SRC = $(SRC_DIR)/runtime/geo_job_queue.c
CLIENT_SRC = $(SRC_DIR)/client/geo_client.c
SERVER_SRC = $(SRC_DIR)/server/geo_server.c
SERVER_COMMANDS_SRC = $(SRC_DIR)/server/geo_server_commands.c
DAEMON_SRC = $(SRC_DIR)/server/geoboltd.c
SIMD_SRC_ARM64 = $(SRC_DIR)/simd/geo_index_simd_arm64.c
SIMD_SRC_X86 = $(SRC_DIR)/simd/geo_index_simd_x86.c
SIMD_SRC_SCALAR = $(SRC_DIR)/simd/geo_index_simd_scalar.c

# Object files
LIB_OBJ = $(BUILD_DIR)/geo_index.o
DENSITY_OBJ = $(BUILD_DIR)/geo_index_density.o
STREAM_OBJ = $(BUILD_DIR)/geo_index_stream.o
PARALLEL_OBJ = $(BUILD_DIR)/geo_index_parallel.o
BATCH_OBJ = $(BUILD_DIR)/geo_index_batch.o
SEGMENTS_OBJ = $(BUILD_DIR)/geo_index_segments.o
IO_OBJ = $(BUILD_DIR)/geo_index_io.o
THREAD_POOL_OBJ = $(BUILD_DIR)/geo_thread_pool.o
ROCKS_BRIDGE_OBJ = $(BUILD_DIR)/geo_rocks_bridge.o
GEODOC_OBJ = $(BUILD_DIR)/geo_doc.o
DB_FORMAT_OBJ = $(BUILD_DIR)/geo_db_format.o
DATABASE_OBJ = $(BUILD_DIR)/geo_database.o
COMMIT_COORDINATOR_OBJ = $(BUILD_DIR)/geo_commit_coordinator.o
SPATIAL_MEMTABLE_OBJ = $(BUILD_DIR)/geo_spatial_memtable.o
OBJECT_CACHE_OBJ = $(BUILD_DIR)/geo_object_cache.o
VISIBILITY_GATE_OBJ = $(BUILD_DIR)/geo_visibility_gate.o
SECONDARY_INDEX_OBJ = $(BUILD_DIR)/geo_secondary_index.o
SECONDARY_STATISTICS_OBJ = $(BUILD_DIR)/geo_secondary_statistics.o
QUERY_PLANNER_OBJ = $(BUILD_DIR)/geo_query_planner.o
PROTOCOL_OBJ = $(BUILD_DIR)/geo_protocol.o
JOB_QUEUE_OBJ = $(BUILD_DIR)/geo_job_queue.o
CLIENT_OBJ = $(BUILD_DIR)/geo_client.o
SERVER_OBJ = $(BUILD_DIR)/geo_server.o
SERVER_COMMANDS_OBJ = $(BUILD_DIR)/geo_server_commands.o
LIB_OBJ_DEBUG = $(BUILD_DIR)/geo_index_debug.o
DENSITY_OBJ_DEBUG = $(BUILD_DIR)/geo_index_density_debug.o
STREAM_OBJ_DEBUG = $(BUILD_DIR)/geo_index_stream_debug.o
PARALLEL_OBJ_DEBUG = $(BUILD_DIR)/geo_index_parallel_debug.o
BATCH_OBJ_DEBUG = $(BUILD_DIR)/geo_index_batch_debug.o
SEGMENTS_OBJ_DEBUG = $(BUILD_DIR)/geo_index_segments_debug.o
IO_OBJ_DEBUG = $(BUILD_DIR)/geo_index_io_debug.o
THREAD_POOL_OBJ_DEBUG = $(BUILD_DIR)/geo_thread_pool_debug.o
ROCKS_BRIDGE_OBJ_DEBUG = $(BUILD_DIR)/geo_rocks_bridge_debug.o
GEODOC_OBJ_DEBUG = $(BUILD_DIR)/geo_doc_debug.o
DB_FORMAT_OBJ_DEBUG = $(BUILD_DIR)/geo_db_format_debug.o
DATABASE_OBJ_DEBUG = $(BUILD_DIR)/geo_database_debug.o
COMMIT_COORDINATOR_OBJ_DEBUG = $(BUILD_DIR)/geo_commit_coordinator_debug.o
SPATIAL_MEMTABLE_OBJ_DEBUG = $(BUILD_DIR)/geo_spatial_memtable_debug.o
OBJECT_CACHE_OBJ_DEBUG = $(BUILD_DIR)/geo_object_cache_debug.o
VISIBILITY_GATE_OBJ_DEBUG = $(BUILD_DIR)/geo_visibility_gate_debug.o
SECONDARY_INDEX_OBJ_DEBUG = $(BUILD_DIR)/geo_secondary_index_debug.o
SECONDARY_STATISTICS_OBJ_DEBUG = $(BUILD_DIR)/geo_secondary_statistics_debug.o
QUERY_PLANNER_OBJ_DEBUG = $(BUILD_DIR)/geo_query_planner_debug.o
PROTOCOL_OBJ_DEBUG = $(BUILD_DIR)/geo_protocol_debug.o
JOB_QUEUE_OBJ_DEBUG = $(BUILD_DIR)/geo_job_queue_debug.o
CLIENT_OBJ_DEBUG = $(BUILD_DIR)/geo_client_debug.o
SERVER_OBJ_DEBUG = $(BUILD_DIR)/geo_server_debug.o
SERVER_COMMANDS_OBJ_DEBUG = $(BUILD_DIR)/geo_server_commands_debug.o
LIB_OBJ_SCALAR = $(BUILD_DIR)/geo_index_scalar.o
DENSITY_OBJ_SCALAR = $(BUILD_DIR)/geo_index_density_scalar.o
STREAM_OBJ_SCALAR = $(BUILD_DIR)/geo_index_stream_scalar.o
PARALLEL_OBJ_SCALAR = $(BUILD_DIR)/geo_index_parallel_scalar.o
BATCH_OBJ_SCALAR = $(BUILD_DIR)/geo_index_batch_scalar.o
SEGMENTS_OBJ_SCALAR = $(BUILD_DIR)/geo_index_segments_scalar.o
IO_OBJ_SCALAR = $(BUILD_DIR)/geo_index_io_scalar.o
THREAD_POOL_OBJ_SCALAR = $(BUILD_DIR)/geo_thread_pool_scalar.o
ROCKS_BRIDGE_OBJ_SCALAR = $(BUILD_DIR)/geo_rocks_bridge_scalar.o
GEODOC_OBJ_SCALAR = $(BUILD_DIR)/geo_doc_scalar.o
DB_FORMAT_OBJ_SCALAR = $(BUILD_DIR)/geo_db_format_scalar.o
DATABASE_OBJ_SCALAR = $(BUILD_DIR)/geo_database_scalar.o
COMMIT_COORDINATOR_OBJ_SCALAR = $(BUILD_DIR)/geo_commit_coordinator_scalar.o
SPATIAL_MEMTABLE_OBJ_SCALAR = $(BUILD_DIR)/geo_spatial_memtable_scalar.o
OBJECT_CACHE_OBJ_SCALAR = $(BUILD_DIR)/geo_object_cache_scalar.o
VISIBILITY_GATE_OBJ_SCALAR = $(BUILD_DIR)/geo_visibility_gate_scalar.o
SECONDARY_INDEX_OBJ_SCALAR = $(BUILD_DIR)/geo_secondary_index_scalar.o
SECONDARY_STATISTICS_OBJ_SCALAR = $(BUILD_DIR)/geo_secondary_statistics_scalar.o
QUERY_PLANNER_OBJ_SCALAR = $(BUILD_DIR)/geo_query_planner_scalar.o
PROTOCOL_OBJ_SCALAR = $(BUILD_DIR)/geo_protocol_scalar.o
JOB_QUEUE_OBJ_SCALAR = $(BUILD_DIR)/geo_job_queue_scalar.o
CLIENT_OBJ_SCALAR = $(BUILD_DIR)/geo_client_scalar.o
SERVER_OBJ_SCALAR = $(BUILD_DIR)/geo_server_scalar.o
SERVER_COMMANDS_OBJ_SCALAR = $(BUILD_DIR)/geo_server_commands_scalar.o
SCALAR_BACKEND_OBJ = $(BUILD_DIR)/geo_index_simd_forced_scalar.o

DEBUG_LIBRARY_OBJECTS = $(LIB_OBJ_DEBUG) $(DENSITY_OBJ_DEBUG) $(STREAM_OBJ_DEBUG) $(PARALLEL_OBJ_DEBUG) \
	$(BATCH_OBJ_DEBUG) $(SEGMENTS_OBJ_DEBUG) $(IO_OBJ_DEBUG) $(THREAD_POOL_OBJ_DEBUG) \
	$(ROCKS_BRIDGE_OBJ_DEBUG) $(GEODOC_OBJ_DEBUG) $(DB_FORMAT_OBJ_DEBUG) $(DATABASE_OBJ_DEBUG) $(COMMIT_COORDINATOR_OBJ_DEBUG) \
	$(SPATIAL_MEMTABLE_OBJ_DEBUG) $(OBJECT_CACHE_OBJ_DEBUG) $(VISIBILITY_GATE_OBJ_DEBUG) $(SECONDARY_INDEX_OBJ_DEBUG) \
	$(SECONDARY_STATISTICS_OBJ_DEBUG) \
	$(QUERY_PLANNER_OBJ_DEBUG) \
	$(PROTOCOL_OBJ_DEBUG) \
	$(JOB_QUEUE_OBJ_DEBUG) \
	$(CLIENT_OBJ_DEBUG) $(SERVER_OBJ_DEBUG) \
	$(SERVER_COMMANDS_OBJ_DEBUG) $(SIMD_OBJ_DEBUG)

SCALAR_LIBRARY_OBJECTS = $(LIB_OBJ_SCALAR) $(DENSITY_OBJ_SCALAR) $(STREAM_OBJ_SCALAR) $(PARALLEL_OBJ_SCALAR) \
	$(BATCH_OBJ_SCALAR) $(SEGMENTS_OBJ_SCALAR) $(IO_OBJ_SCALAR) $(THREAD_POOL_OBJ_SCALAR) \
	$(ROCKS_BRIDGE_OBJ_SCALAR) $(GEODOC_OBJ_SCALAR) $(DB_FORMAT_OBJ_SCALAR) $(DATABASE_OBJ_SCALAR) \
	$(COMMIT_COORDINATOR_OBJ_SCALAR) $(SPATIAL_MEMTABLE_OBJ_SCALAR) $(OBJECT_CACHE_OBJ_SCALAR) $(VISIBILITY_GATE_OBJ_SCALAR) \
	$(SECONDARY_INDEX_OBJ_SCALAR) \
	$(SECONDARY_STATISTICS_OBJ_SCALAR) $(QUERY_PLANNER_OBJ_SCALAR) \
	$(PROTOCOL_OBJ_SCALAR) \
	$(JOB_QUEUE_OBJ_SCALAR) \
	$(CLIENT_OBJ_SCALAR) $(SERVER_OBJ_SCALAR) \
	$(SERVER_COMMANDS_OBJ_SCALAR) \
	$(SCALAR_BACKEND_OBJ)

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
DATABASE_TEST_SRC = $(TEST_DIR)/test_database.c
STORAGE_TEST_SRC = $(TEST_DIR)/test_storage.c
SERVER_TEST_SRC = $(TEST_DIR)/test_server.c
DAEMON_TEST_SRC = $(TEST_DIR)/test_daemon.c
TEST_SUPPORT_SRC = $(TEST_DIR)/test_support.c
SAMPLE_SRC = $(SAMPLE_DIR)/basic_usage.c
CLIENT_SAMPLE_SRC = $(SAMPLE_DIR)/client_usage.c
METADATA_SAMPLE_SRC = $(SAMPLE_DIR)/client_metadata.c
DEMO_SRC = $(BENCHMARK_DIR)/demo_benchmark.c
TOMBSTONE_BENCH_SRC = $(BENCHMARK_DIR)/tombstone_visibility.c
STREAM_BENCH_SRC = $(BENCHMARK_DIR)/stream_build.c
COMPACTION_BENCH_SRC = $(BENCHMARK_DIR)/compaction_merge.c
SERVER_BENCH_SRC = $(BENCHMARK_DIR)/server_roundtrip.c
QUERY_PLANNER_BENCH_SRC = $(BENCHMARK_DIR)/query_planner.c
SOAK_SRC = $(TEST_DIR)/soak_uber.c

# Targets
LIB_STATIC = $(BUILD_DIR)/libgeobolt.a
TEST_BIN = $(BIN_DIR)/test_geo_index
DATABASE_TEST_BIN = $(BIN_DIR)/test_database
STORAGE_TEST_BIN = $(BIN_DIR)/test_storage
SERVER_TEST_BIN = $(BIN_DIR)/test_server
DAEMON_TEST_BIN = $(BIN_DIR)/test_daemon
STORAGE_TEST_BIN_DEBUG = $(BIN_DIR)/test_storage_debug
SERVER_TEST_BIN_DEBUG = $(BIN_DIR)/test_server_debug
STORAGE_TEST_BIN_SCALAR = $(BIN_DIR)/test_storage_scalar
SERVER_TEST_BIN_SCALAR = $(BIN_DIR)/test_server_scalar
DAEMON_BIN = $(BIN_DIR)/geoboltd
DAEMON_BIN_DEBUG = $(BIN_DIR)/geoboltd_debug
DAEMON_BIN_SCALAR = $(BIN_DIR)/geoboltd_scalar
DAEMON_TEST_BIN_DEBUG = $(BIN_DIR)/test_daemon_debug
DAEMON_TEST_BIN_SCALAR = $(BIN_DIR)/test_daemon_scalar
TEST_BIN_DEBUG = $(BIN_DIR)/test_geo_index_debug
DATABASE_TEST_BIN_DEBUG = $(BIN_DIR)/test_database_debug
TEST_BIN_SCALAR = $(BIN_DIR)/test_geo_index_scalar
DATABASE_TEST_BIN_SCALAR = $(BIN_DIR)/test_database_scalar
SAMPLE_BIN = $(BIN_DIR)/basic_usage
CLIENT_SAMPLE_BIN = $(BIN_DIR)/client_usage
METADATA_SAMPLE_BIN = $(BIN_DIR)/client_metadata
DEMO_BIN = $(BIN_DIR)/demo_benchmark
BENCH_BIN = $(BIN_DIR)/benchmark
TOMBSTONE_BENCH_BIN = $(BIN_DIR)/benchmark_tombstones
STREAM_BENCH_BIN = $(BIN_DIR)/benchmark_stream
COMPACTION_BENCH_BIN = $(BIN_DIR)/benchmark_compaction
SERVER_BENCH_BIN = $(BIN_DIR)/benchmark_server
QUERY_PLANNER_BENCH_BIN = $(BIN_DIR)/benchmark_query_planner
SOAK_BIN = $(BIN_DIR)/soak_uber

DIRECT_BUILD_TARGETS = $(LIB_OBJ) $(DENSITY_OBJ) $(STREAM_OBJ) $(PARALLEL_OBJ) $(BATCH_OBJ) $(SEGMENTS_OBJ) $(IO_OBJ) \
	$(THREAD_POOL_OBJ) $(ROCKS_BRIDGE_OBJ) $(GEODOC_OBJ) $(DB_FORMAT_OBJ) $(DATABASE_OBJ) $(COMMIT_COORDINATOR_OBJ) $(PROTOCOL_OBJ) \
	$(SPATIAL_MEMTABLE_OBJ) $(OBJECT_CACHE_OBJ) $(VISIBILITY_GATE_OBJ) $(SECONDARY_INDEX_OBJ) $(SECONDARY_STATISTICS_OBJ) \
	$(QUERY_PLANNER_OBJ) \
	$(JOB_QUEUE_OBJ) $(CLIENT_OBJ) \
	$(SERVER_OBJ) \
	$(SERVER_COMMANDS_OBJ) \
	$(SIMD_OBJ) \
	$(LIB_OBJ_DEBUG) $(DENSITY_OBJ_DEBUG) $(STREAM_OBJ_DEBUG) $(PARALLEL_OBJ_DEBUG) $(BATCH_OBJ_DEBUG) $(SEGMENTS_OBJ_DEBUG) \
	$(IO_OBJ_DEBUG) $(THREAD_POOL_OBJ_DEBUG) $(ROCKS_BRIDGE_OBJ_DEBUG) $(GEODOC_OBJ_DEBUG) $(DB_FORMAT_OBJ_DEBUG) \
	$(DATABASE_OBJ_DEBUG) $(COMMIT_COORDINATOR_OBJ_DEBUG) \
	$(SPATIAL_MEMTABLE_OBJ_DEBUG) $(OBJECT_CACHE_OBJ_DEBUG) $(VISIBILITY_GATE_OBJ_DEBUG) $(SECONDARY_INDEX_OBJ_DEBUG) \
	$(SECONDARY_STATISTICS_OBJ_DEBUG) \
	$(QUERY_PLANNER_OBJ_DEBUG) \
	$(PROTOCOL_OBJ_DEBUG) \
	$(JOB_QUEUE_OBJ_DEBUG) \
	$(CLIENT_OBJ_DEBUG) $(SERVER_OBJ_DEBUG) $(SERVER_COMMANDS_OBJ_DEBUG) $(SIMD_OBJ_DEBUG) $(LIB_OBJ_SCALAR) \
	$(DENSITY_OBJ_SCALAR) \
	$(STREAM_OBJ_SCALAR) \
	$(PARALLEL_OBJ_SCALAR) $(BATCH_OBJ_SCALAR) $(SEGMENTS_OBJ_SCALAR) $(IO_OBJ_SCALAR) $(THREAD_POOL_OBJ_SCALAR) \
	$(ROCKS_BRIDGE_OBJ_SCALAR) $(GEODOC_OBJ_SCALAR) $(DB_FORMAT_OBJ_SCALAR) \
	$(DATABASE_OBJ_SCALAR) $(COMMIT_COORDINATOR_OBJ_SCALAR) $(SPATIAL_MEMTABLE_OBJ_SCALAR) $(OBJECT_CACHE_OBJ_SCALAR) \
	$(VISIBILITY_GATE_OBJ_SCALAR) \
	$(SECONDARY_INDEX_OBJ_SCALAR) \
	$(SECONDARY_STATISTICS_OBJ_SCALAR) \
	$(QUERY_PLANNER_OBJ_SCALAR) \
	$(PROTOCOL_OBJ_SCALAR) \
	$(CLIENT_OBJ_SCALAR) $(SERVER_OBJ_SCALAR) \
	$(SERVER_COMMANDS_OBJ_SCALAR) \
	$(SCALAR_BACKEND_OBJ) $(LIB_STATIC) $(TEST_BIN) $(DATABASE_TEST_BIN) $(STORAGE_TEST_BIN) $(SERVER_TEST_BIN) $(DAEMON_TEST_BIN) \
	$(TEST_BIN_DEBUG) $(DATABASE_TEST_BIN_DEBUG) $(STORAGE_TEST_BIN_DEBUG) $(SERVER_TEST_BIN_DEBUG) $(TEST_BIN_SCALAR) \
	$(DATABASE_TEST_BIN_SCALAR) $(STORAGE_TEST_BIN_SCALAR) $(SERVER_TEST_BIN_SCALAR) $(DAEMON_BIN) $(DAEMON_BIN_DEBUG) \
	$(DAEMON_BIN_SCALAR) $(DAEMON_TEST_BIN_DEBUG) $(DAEMON_TEST_BIN_SCALAR) $(SAMPLE_BIN) $(CLIENT_SAMPLE_BIN) \
	$(METADATA_SAMPLE_BIN) $(DEMO_BIN) \
	$(TOMBSTONE_BENCH_BIN) $(STREAM_BENCH_BIN) $(SERVER_BENCH_BIN) $(QUERY_PLANNER_BENCH_BIN) \
	$(COMPACTION_BENCH_BIN) $(SOAK_BIN)

.PHONY: all clean test test-debug test-scalar lib sample client-sample metadata-sample demo benchmark benchmark-tombstones \
	benchmark-stream \
	benchmark-compaction \
	benchmark-server \
	benchmark-query-planner \
	simd-benchmark \
	soak-uber soak-uber-smoke soak-uber-tsan allocator-info pgo-generate pgo-merge pgo-use help dirs FORCE_OPT_LINK
