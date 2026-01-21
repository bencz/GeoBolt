# GeoIndex Makefile
# Production-grade build system

CC = clang
CFLAGS = -Wall -Wextra -Werror -std=c11 -pedantic
CFLAGS_OPT = -O3 -march=native -flto -ffast-math
CFLAGS_DEBUG = -g -O0 -DDEBUG -fsanitize=address,undefined
LDFLAGS = -lm

# Directories
SRC_DIR = .
BUILD_DIR = build
BIN_DIR = bin

# Source files
LIB_SRC = geo_index.c
SIMD_SRC_ARM64 = geo_index_simd_arm64.c
SIMD_SRC_X86 = geo_index_simd_x86.c
SIMD_SRC_SCALAR = geo_index_simd_scalar.c

# Object files
LIB_OBJ = $(BUILD_DIR)/geo_index.o
LIB_OBJ_DEBUG = $(BUILD_DIR)/geo_index_debug.o

# SIMD object files (architecture-specific)
SIMD_OBJ_ARM64 = $(BUILD_DIR)/geo_index_simd_arm64.o
SIMD_OBJ_X86 = $(BUILD_DIR)/geo_index_simd_x86.o
SIMD_OBJ_SCALAR = $(BUILD_DIR)/geo_index_simd_scalar.o

# Detect architecture
UNAME_M := $(shell uname -m)
ifeq ($(UNAME_M),arm64)
    SIMD_OBJ = $(SIMD_OBJ_ARM64)
    SIMD_SRC = $(SIMD_SRC_ARM64)
    SIMD_FLAGS = -march=armv8-a+simd
else ifeq ($(UNAME_M),aarch64)
    SIMD_OBJ = $(SIMD_OBJ_ARM64)
    SIMD_SRC = $(SIMD_SRC_ARM64)
    SIMD_FLAGS = -march=armv8-a+simd
else ifeq ($(UNAME_M),x86_64)
    SIMD_OBJ = $(SIMD_OBJ_X86)
    SIMD_SRC = $(SIMD_SRC_X86)
    SIMD_FLAGS = -mavx2 -mfma
else
    SIMD_OBJ = $(SIMD_OBJ_SCALAR)
    SIMD_SRC = $(SIMD_SRC_SCALAR)
    SIMD_FLAGS =
endif

TEST_SRC = test_geo_index.c
MAIN_SRC = main.c
DEMO_SRC = demo_benchmark.c

# Targets
LIB_STATIC = $(BUILD_DIR)/libgeoindex.a
TEST_BIN = $(BIN_DIR)/test_geo_index
TEST_BIN_DEBUG = $(BIN_DIR)/test_geo_index_debug
MAIN_BIN = $(BIN_DIR)/main
DEMO_BIN = $(BIN_DIR)/demo_benchmark
BENCH_BIN = $(BIN_DIR)/benchmark

.PHONY: all clean test test-debug lib main benchmark help dirs

# Default target
all: dirs lib $(TEST_BIN) $(MAIN_BIN)

# Create directories
dirs:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

# Build static library (optimized)
$(LIB_OBJ): $(SRC_DIR)/geo_index.c $(SRC_DIR)/geo_index.h $(SRC_DIR)/geo_index_simd.h
	$(CC) $(CFLAGS) $(CFLAGS_OPT) -c $< -o $@

# Build SIMD object file
$(SIMD_OBJ): $(SRC_DIR)/$(SIMD_SRC) $(SRC_DIR)/geo_index_simd.h $(SRC_DIR)/geo_index.h
	$(CC) $(CFLAGS) $(CFLAGS_OPT) $(SIMD_FLAGS) -c $< -o $@

$(LIB_STATIC): $(LIB_OBJ) $(SIMD_OBJ)
	ar rcs $@ $^

lib: dirs $(LIB_STATIC)

# Build static library (debug)
$(LIB_OBJ_DEBUG): $(SRC_DIR)/geo_index.c $(SRC_DIR)/geo_index.h
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) -c $< -o $@

# Build test binary (optimized)
$(TEST_BIN): $(SRC_DIR)/test_geo_index.c $(LIB_STATIC)
	$(CC) $(CFLAGS) $(CFLAGS_OPT) $(SIMD_FLAGS) $< -L$(BUILD_DIR) -lgeoindex $(LDFLAGS) -o $@

# Build test binary (debug with sanitizers)
$(TEST_BIN_DEBUG): $(SRC_DIR)/test_geo_index.c $(LIB_OBJ_DEBUG)
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $< $(LIB_OBJ_DEBUG) $(LDFLAGS) -o $@

# Build main example
$(MAIN_BIN): $(SRC_DIR)/main.c $(LIB_STATIC)
	$(CC) $(CFLAGS) $(CFLAGS_OPT) $(SIMD_FLAGS) $< -L$(BUILD_DIR) -lgeoindex $(LDFLAGS) -o $@

# Build demo benchmark
$(DEMO_BIN): $(SRC_DIR)/demo_benchmark.c $(LIB_STATIC)
	$(CC) $(CFLAGS) $(CFLAGS_OPT) $(SIMD_FLAGS) $< -L$(BUILD_DIR) -lgeoindex $(LDFLAGS) -lpthread -o $@

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
	@perl -e 'alarm 120; exec @ARGV' $(TEST_BIN) || (echo "Test timed out or failed!" && exit 1)

# Run tests (debug with sanitizers)
test-debug: dirs $(TEST_BIN_DEBUG)
	@echo "=========================================="
	@echo "Running debug tests with sanitizers..."
	@echo "=========================================="
	@$(TEST_BIN_DEBUG)

# Build and run main example
main: dirs $(MAIN_BIN)
	@$(MAIN_BIN)

# Quick benchmark
benchmark: dirs $(TEST_BIN)
	@echo "=========================================="
	@echo "Running performance benchmarks..."
	@echo "=========================================="
	@$(TEST_BIN) 2>&1 | grep -A 100 "PERFORMANCE TESTS"

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# SIMD benchmark only
simd-benchmark: dirs $(TEST_BIN)
	@echo "=========================================="
	@echo "Running SIMD benchmarks..."
	@echo "Architecture: $(UNAME_M)"
	@echo "SIMD Flags: $(SIMD_FLAGS)"
	@echo "=========================================="
	@$(TEST_BIN) 2>&1 | grep -A 200 "SIMD"

# Help
help:
	@echo "GeoIndex Build System"
	@echo "====================="
	@echo ""
	@echo "Architecture: $(UNAME_M)"
	@echo "SIMD Source:  $(SIMD_SRC)"
	@echo "SIMD Flags:   $(SIMD_FLAGS)"
	@echo ""
	@echo "Targets:"
	@echo "  all            - Build library, tests, main, and demo (default)"
	@echo "  lib            - Build static library only"
	@echo "  test           - Build and run optimized tests"
	@echo "  test-debug     - Build and run tests with sanitizers"
	@echo "  main           - Build and run main example"
	@echo "  demo           - Run comprehensive benchmark demo (10M points)"
	@echo "  benchmark      - Run performance benchmarks"
	@echo "  simd-benchmark - Run SIMD-specific benchmarks"
	@echo "  clean          - Remove build artifacts"
	@echo "  help           - Show this help message"
	@echo ""
	@echo "Build Flags:"
	@echo "  Optimized:  $(CFLAGS_OPT)"
	@echo "  Debug:      $(CFLAGS_DEBUG)"
	@echo "  SIMD:       $(SIMD_FLAGS)"
