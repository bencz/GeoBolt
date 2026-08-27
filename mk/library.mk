# Default target
all: dirs lib $(DAEMON_BIN) $(TEST_BIN) $(DATABASE_TEST_BIN) $(SERVER_TEST_BIN) $(DAEMON_TEST_BIN) $(SAMPLE_BIN) \
	$(CLIENT_SAMPLE_BIN) $(DEMO_BIN)

# Create directories
dirs:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

$(DIRECT_BUILD_TARGETS): | dirs

# Build static library (optimized)
$(LIB_OBJ): $(LIB_SRC) include/geobolt/geo_index.h src/core/geo_index_private.h include/geobolt/geo_index_simd.h \
		src/storage/geo_index_persistence.h src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(DENSITY_OBJ): $(DENSITY_SRC) include/geobolt/geo_index.h src/storage/geo_index_persistence.h src/core/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(STREAM_OBJ): $(STREAM_SRC) include/geobolt/geo_index.h src/storage/geo_index_persistence.h src/core/geo_index_private.h \
		src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(PARALLEL_OBJ): $(PARALLEL_SRC) include/geobolt/geo_index.h src/core/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(BATCH_OBJ): $(BATCH_SRC) include/geobolt/geo_index.h src/core/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SEGMENTS_OBJ): $(SEGMENTS_SRC) include/geobolt/geo_index.h src/storage/geo_index_persistence.h src/core/geo_index_private.h \
		src/storage/geo_index_io.h src/runtime/geo_thread_pool.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(IO_OBJ): $(IO_SRC) src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(THREAD_POOL_OBJ): $(THREAD_POOL_SRC) src/runtime/geo_thread_pool.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(WAL_OBJ): $(WAL_SRC) src/storage/geo_wal.h src/storage/geo_index_io.h src/storage/geo_index_persistence.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(DATABASE_OBJ): $(DATABASE_SRC) include/geobolt/geobolt.h include/geobolt/geo_index.h src/storage/geo_wal.h \
		src/storage/geo_index_io.h src/storage/geo_index_persistence.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(PROTOCOL_OBJ): $(PROTOCOL_SRC) src/protocol/geo_protocol.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(JOB_QUEUE_OBJ): $(JOB_QUEUE_SRC) src/runtime/geo_job_queue.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(CLIENT_OBJ): $(CLIENT_SRC) include/geobolt/client.h src/protocol/geo_protocol.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SERVER_OBJ): $(SERVER_SRC) include/geobolt/server.h include/geobolt/geobolt.h src/protocol/geo_protocol.h \
		src/runtime/geo_job_queue.h src/server/geo_server_internal.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SERVER_COMMANDS_OBJ): $(SERVER_COMMANDS_SRC) src/server/geo_server_internal.h include/geobolt/geobolt.h \
		src/protocol/geo_protocol.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

# Build SIMD object file
$(SIMD_OBJ): $(SIMD_SRC) include/geobolt/geo_index_simd.h src/simd/geo_index_simd_scalar_kernels.h include/geobolt/geo_index.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) $(SIMD_FLAGS) -MMD -MP -c $< -o $@

$(LIB_STATIC): $(LIB_OBJ) $(DENSITY_OBJ) $(STREAM_OBJ) $(PARALLEL_OBJ) $(BATCH_OBJ) $(SEGMENTS_OBJ) $(IO_OBJ) \
		$(THREAD_POOL_OBJ) $(WAL_OBJ) $(DATABASE_OBJ) $(PROTOCOL_OBJ) $(JOB_QUEUE_OBJ) $(CLIENT_OBJ) $(SERVER_OBJ) \
		$(SERVER_COMMANDS_OBJ) $(SIMD_OBJ)
	$(AR) rcs $@ $^

lib: dirs $(LIB_STATIC)

# Build static library (debug)
$(LIB_OBJ_DEBUG): $(LIB_SRC) include/geobolt/geo_index.h src/core/geo_index_private.h \
		src/storage/geo_index_persistence.h src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(DENSITY_OBJ_DEBUG): $(DENSITY_SRC) include/geobolt/geo_index.h src/storage/geo_index_persistence.h \
		src/core/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(STREAM_OBJ_DEBUG): $(STREAM_SRC) include/geobolt/geo_index.h src/storage/geo_index_persistence.h \
		src/core/geo_index_private.h src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(PARALLEL_OBJ_DEBUG): $(PARALLEL_SRC) include/geobolt/geo_index.h src/core/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(BATCH_OBJ_DEBUG): $(BATCH_SRC) include/geobolt/geo_index.h src/core/geo_index_private.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(SEGMENTS_OBJ_DEBUG): $(SEGMENTS_SRC) include/geobolt/geo_index.h src/storage/geo_index_persistence.h \
		src/core/geo_index_private.h src/storage/geo_index_io.h src/runtime/geo_thread_pool.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(IO_OBJ_DEBUG): $(IO_SRC) src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(THREAD_POOL_OBJ_DEBUG): $(THREAD_POOL_SRC) src/runtime/geo_thread_pool.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(WAL_OBJ_DEBUG): $(WAL_SRC) src/storage/geo_wal.h src/storage/geo_index_io.h src/storage/geo_index_persistence.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(DATABASE_OBJ_DEBUG): $(DATABASE_SRC) include/geobolt/geobolt.h include/geobolt/geo_index.h src/storage/geo_wal.h \
		src/storage/geo_index_io.h src/storage/geo_index_persistence.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(PROTOCOL_OBJ_DEBUG): $(PROTOCOL_SRC) src/protocol/geo_protocol.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(JOB_QUEUE_OBJ_DEBUG): $(JOB_QUEUE_SRC) src/runtime/geo_job_queue.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(CLIENT_OBJ_DEBUG): $(CLIENT_SRC) include/geobolt/client.h src/protocol/geo_protocol.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(SERVER_OBJ_DEBUG): $(SERVER_SRC) include/geobolt/server.h include/geobolt/geobolt.h src/protocol/geo_protocol.h \
		src/runtime/geo_job_queue.h src/server/geo_server_internal.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(SERVER_COMMANDS_OBJ_DEBUG): $(SERVER_COMMANDS_SRC) src/server/geo_server_internal.h include/geobolt/geobolt.h \
		src/protocol/geo_protocol.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

SIMD_OBJ_DEBUG = $(BUILD_DIR)/geo_index_simd_debug.o

$(SIMD_OBJ_DEBUG): $(SIMD_SRC) include/geobolt/geo_index_simd.h src/simd/geo_index_simd_scalar_kernels.h \
		src/core/geo_index_internal.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) $(SIMD_FLAGS) -MMD -MP -c $< -o $@

$(LIB_OBJ_SCALAR): $(LIB_SRC) include/geobolt/geo_index.h src/core/geo_index_private.h include/geobolt/geo_index_simd.h \
		src/storage/geo_index_persistence.h src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(DENSITY_OBJ_SCALAR): $(DENSITY_SRC) include/geobolt/geo_index.h src/storage/geo_index_persistence.h \
		src/core/geo_index_private.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(STREAM_OBJ_SCALAR): $(STREAM_SRC) include/geobolt/geo_index.h src/storage/geo_index_persistence.h \
		src/core/geo_index_private.h src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(PARALLEL_OBJ_SCALAR): $(PARALLEL_SRC) include/geobolt/geo_index.h src/core/geo_index_private.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(BATCH_OBJ_SCALAR): $(BATCH_SRC) include/geobolt/geo_index.h src/core/geo_index_private.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SEGMENTS_OBJ_SCALAR): $(SEGMENTS_SRC) include/geobolt/geo_index.h src/storage/geo_index_persistence.h \
		src/core/geo_index_private.h src/storage/geo_index_io.h src/runtime/geo_thread_pool.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(IO_OBJ_SCALAR): $(IO_SRC) src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(THREAD_POOL_OBJ_SCALAR): $(THREAD_POOL_SRC) src/runtime/geo_thread_pool.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(WAL_OBJ_SCALAR): $(WAL_SRC) src/storage/geo_wal.h src/storage/geo_index_io.h src/storage/geo_index_persistence.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(DATABASE_OBJ_SCALAR): $(DATABASE_SRC) include/geobolt/geobolt.h include/geobolt/geo_index.h src/storage/geo_wal.h \
		src/storage/geo_index_io.h src/storage/geo_index_persistence.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(PROTOCOL_OBJ_SCALAR): $(PROTOCOL_SRC) src/protocol/geo_protocol.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(JOB_QUEUE_OBJ_SCALAR): $(JOB_QUEUE_SRC) src/runtime/geo_job_queue.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(CLIENT_OBJ_SCALAR): $(CLIENT_SRC) include/geobolt/client.h src/protocol/geo_protocol.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SERVER_OBJ_SCALAR): $(SERVER_SRC) include/geobolt/server.h include/geobolt/geobolt.h src/protocol/geo_protocol.h \
		src/runtime/geo_job_queue.h src/server/geo_server_internal.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SERVER_COMMANDS_OBJ_SCALAR): $(SERVER_COMMANDS_SRC) src/server/geo_server_internal.h include/geobolt/geobolt.h \
		src/protocol/geo_protocol.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SCALAR_BACKEND_OBJ): $(SIMD_SRC_SCALAR) include/geobolt/geo_index_simd.h src/simd/geo_index_simd_scalar_kernels.h \
		src/core/geo_index_internal.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@
