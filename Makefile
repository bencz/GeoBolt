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
LIB_OBJ = $(BUILD_DIR)/geo_index.o
LIB_OBJ_DEBUG = $(BUILD_DIR)/geo_index_debug.o

TEST_SRC = test_geo_index.c
MAIN_SRC = main.c

# Targets
LIB_STATIC = $(BUILD_DIR)/libgeoindex.a
TEST_BIN = $(BIN_DIR)/test_geo_index
TEST_BIN_DEBUG = $(BIN_DIR)/test_geo_index_debug
MAIN_BIN = $(BIN_DIR)/main
BENCH_BIN = $(BIN_DIR)/benchmark

.PHONY: all clean test test-debug lib main benchmark help dirs

# Default target
all: dirs lib $(TEST_BIN) $(MAIN_BIN)

# Create directories
dirs:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

# Build static library (optimized)
$(LIB_OBJ): $(SRC_DIR)/geo_index.c $(SRC_DIR)/geo_index.h
	$(CC) $(CFLAGS) $(CFLAGS_OPT) -c $< -o $@

$(LIB_STATIC): $(LIB_OBJ)
	ar rcs $@ $^

lib: dirs $(LIB_STATIC)

# Build static library (debug)
$(LIB_OBJ_DEBUG): $(SRC_DIR)/geo_index.c $(SRC_DIR)/geo_index.h
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) -c $< -o $@

# Build test binary (optimized)
$(TEST_BIN): $(SRC_DIR)/test_geo_index.c $(LIB_STATIC)
	$(CC) $(CFLAGS) $(CFLAGS_OPT) $< -L$(BUILD_DIR) -lgeoindex $(LDFLAGS) -o $@

# Build test binary (debug with sanitizers)
$(TEST_BIN_DEBUG): $(SRC_DIR)/test_geo_index.c $(LIB_OBJ_DEBUG)
	$(CC) $(CFLAGS) $(CFLAGS_DEBUG) $< $(LIB_OBJ_DEBUG) $(LDFLAGS) -o $@

# Build main example
$(MAIN_BIN): $(SRC_DIR)/main.c $(LIB_STATIC)
	$(CC) $(CFLAGS) $(CFLAGS_OPT) $< -L$(BUILD_DIR) -lgeoindex $(LDFLAGS) -o $@

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

# Help
help:
	@echo "GeoIndex Build System"
	@echo "====================="
	@echo ""
	@echo "Targets:"
	@echo "  all         - Build library, tests, and main (default)"
	@echo "  lib         - Build static library only"
	@echo "  test        - Build and run optimized tests"
	@echo "  test-debug  - Build and run tests with sanitizers"
	@echo "  main        - Build and run main example"
	@echo "  benchmark   - Run performance benchmarks"
	@echo "  clean       - Remove build artifacts"
	@echo "  help        - Show this help message"
	@echo ""
	@echo "Build Flags:"
	@echo "  Optimized:  $(CFLAGS_OPT)"
	@echo "  Debug:      $(CFLAGS_DEBUG)"
