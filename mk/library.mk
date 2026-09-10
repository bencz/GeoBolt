# Default target
all: dirs lib $(DAEMON_BIN) $(TEST_BIN) $(DATABASE_TEST_BIN) $(STORAGE_TEST_BIN) $(SERVER_TEST_BIN) $(DAEMON_TEST_BIN) $(SAMPLE_BIN) \
	$(CLIENT_SAMPLE_BIN) $(METADATA_SAMPLE_BIN) $(DEMO_BIN)

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

$(ROCKS_BRIDGE_OBJ): $(ROCKS_BRIDGE_SRC) src/storage/rocksdb/geo_rocks_bridge.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXXFLAGS_OPT) -MMD -MP -c $< -o $@

$(GEODOC_OBJ): $(GEODOC_SRC) include/geobolt/geodoc.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(DB_FORMAT_OBJ): $(DB_FORMAT_SRC) src/db/object/geo_db_format.h include/geobolt/geodoc.h \
		src/storage/rocksdb/geo_rocks_bridge.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(DATABASE_OBJ): $(DATABASE_SRC) include/geobolt/geobolt.h include/geobolt/geodoc.h include/geobolt/geo_index.h \
		src/engine/geo_database_internal.h src/db/object/geo_db_format.h src/storage/rocksdb/geo_rocks_bridge.h \
		src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(COMMIT_COORDINATOR_OBJ): $(COMMIT_COORDINATOR_SRC) src/engine/geo_database_internal.h include/geobolt/geobolt.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SPATIAL_MEMTABLE_OBJ): $(SPATIAL_MEMTABLE_SRC) src/engine/geo_spatial_memtable.h include/geobolt/geobolt.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(OBJECT_CACHE_OBJ): $(OBJECT_CACHE_SRC) src/engine/geo_object_cache.h include/geobolt/geodoc.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(VISIBILITY_GATE_OBJ): $(VISIBILITY_GATE_SRC) src/engine/geo_visibility_gate.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SECONDARY_INDEX_OBJ): $(SECONDARY_INDEX_SRC) src/engine/geo_secondary_index.h include/geobolt/geobolt.h \
		src/engine/geo_secondary_statistics.h src/db/object/geo_db_format.h src/storage/rocksdb/geo_rocks_bridge.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SECONDARY_STATISTICS_OBJ): $(SECONDARY_STATISTICS_SRC) src/engine/geo_secondary_statistics.h \
		src/db/object/geo_db_format.h src/storage/rocksdb/geo_rocks_bridge.h include/geobolt/geobolt.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(QUERY_PLANNER_OBJ): $(QUERY_PLANNER_SRC) src/engine/geo_query_planner.h src/engine/geo_database_internal.h \
		include/geobolt/geobolt.h src/core/geo_index_private.h
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
		$(THREAD_POOL_OBJ) $(ROCKS_BRIDGE_OBJ) $(GEODOC_OBJ) $(DB_FORMAT_OBJ) $(DATABASE_OBJ) $(COMMIT_COORDINATOR_OBJ) \
		$(SPATIAL_MEMTABLE_OBJ) $(OBJECT_CACHE_OBJ) $(VISIBILITY_GATE_OBJ) $(SECONDARY_INDEX_OBJ) $(SECONDARY_STATISTICS_OBJ) \
		$(QUERY_PLANNER_OBJ) \
		$(PROTOCOL_OBJ) \
		$(JOB_QUEUE_OBJ) $(CLIENT_OBJ) \
		$(SERVER_OBJ) \
		$(SERVER_COMMANDS_OBJ) $(SIMD_OBJ)
	$(RM) $@
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

$(ROCKS_BRIDGE_OBJ_DEBUG): $(ROCKS_BRIDGE_SRC) src/storage/rocksdb/geo_rocks_bridge.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXXFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(GEODOC_OBJ_DEBUG): $(GEODOC_SRC) include/geobolt/geodoc.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(DB_FORMAT_OBJ_DEBUG): $(DB_FORMAT_SRC) src/db/object/geo_db_format.h include/geobolt/geodoc.h \
		src/storage/rocksdb/geo_rocks_bridge.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(DATABASE_OBJ_DEBUG): $(DATABASE_SRC) include/geobolt/geobolt.h include/geobolt/geodoc.h include/geobolt/geo_index.h \
		src/engine/geo_database_internal.h src/db/object/geo_db_format.h src/storage/rocksdb/geo_rocks_bridge.h \
		src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(COMMIT_COORDINATOR_OBJ_DEBUG): $(COMMIT_COORDINATOR_SRC) src/engine/geo_database_internal.h include/geobolt/geobolt.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(SPATIAL_MEMTABLE_OBJ_DEBUG): $(SPATIAL_MEMTABLE_SRC) src/engine/geo_spatial_memtable.h include/geobolt/geobolt.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(OBJECT_CACHE_OBJ_DEBUG): $(OBJECT_CACHE_SRC) src/engine/geo_object_cache.h include/geobolt/geodoc.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(VISIBILITY_GATE_OBJ_DEBUG): $(VISIBILITY_GATE_SRC) src/engine/geo_visibility_gate.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(SECONDARY_INDEX_OBJ_DEBUG): $(SECONDARY_INDEX_SRC) src/engine/geo_secondary_index.h include/geobolt/geobolt.h \
		src/engine/geo_secondary_statistics.h src/db/object/geo_db_format.h src/storage/rocksdb/geo_rocks_bridge.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(SECONDARY_STATISTICS_OBJ_DEBUG): $(SECONDARY_STATISTICS_SRC) src/engine/geo_secondary_statistics.h \
		src/db/object/geo_db_format.h src/storage/rocksdb/geo_rocks_bridge.h include/geobolt/geobolt.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CFLAGS_DEBUG) -MMD -MP -c $< -o $@

$(QUERY_PLANNER_OBJ_DEBUG): $(QUERY_PLANNER_SRC) src/engine/geo_query_planner.h src/engine/geo_database_internal.h \
		include/geobolt/geobolt.h src/core/geo_index_private.h
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

$(ROCKS_BRIDGE_OBJ_SCALAR): $(ROCKS_BRIDGE_SRC) src/storage/rocksdb/geo_rocks_bridge.h
	$(CXX) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CXXFLAGS) $(CXXFLAGS_OPT) -MMD -MP -c $< -o $@

$(GEODOC_OBJ_SCALAR): $(GEODOC_SRC) include/geobolt/geodoc.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(DB_FORMAT_OBJ_SCALAR): $(DB_FORMAT_SRC) src/db/object/geo_db_format.h include/geobolt/geodoc.h \
		src/storage/rocksdb/geo_rocks_bridge.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(DATABASE_OBJ_SCALAR): $(DATABASE_SRC) include/geobolt/geobolt.h include/geobolt/geodoc.h include/geobolt/geo_index.h \
		src/engine/geo_database_internal.h src/db/object/geo_db_format.h src/storage/rocksdb/geo_rocks_bridge.h \
		src/storage/geo_index_io.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(COMMIT_COORDINATOR_OBJ_SCALAR): $(COMMIT_COORDINATOR_SRC) src/engine/geo_database_internal.h include/geobolt/geobolt.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SPATIAL_MEMTABLE_OBJ_SCALAR): $(SPATIAL_MEMTABLE_SRC) src/engine/geo_spatial_memtable.h include/geobolt/geobolt.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(OBJECT_CACHE_OBJ_SCALAR): $(OBJECT_CACHE_SRC) src/engine/geo_object_cache.h include/geobolt/geodoc.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(VISIBILITY_GATE_OBJ_SCALAR): $(VISIBILITY_GATE_SRC) src/engine/geo_visibility_gate.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SECONDARY_INDEX_OBJ_SCALAR): $(SECONDARY_INDEX_SRC) src/engine/geo_secondary_index.h include/geobolt/geobolt.h \
		src/engine/geo_secondary_statistics.h src/db/object/geo_db_format.h src/storage/rocksdb/geo_rocks_bridge.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(SECONDARY_STATISTICS_OBJ_SCALAR): $(SECONDARY_STATISTICS_SRC) src/engine/geo_secondary_statistics.h \
		src/db/object/geo_db_format.h src/storage/rocksdb/geo_rocks_bridge.h include/geobolt/geobolt.h
	$(CC) $(CPPFLAGS) -DGEO_SIMD_FORCE_SCALAR $(CFLAGS) $(CFLAGS_OPT) -MMD -MP -c $< -o $@

$(QUERY_PLANNER_OBJ_SCALAR): $(QUERY_PLANNER_SRC) src/engine/geo_query_planner.h src/engine/geo_database_internal.h \
		include/geobolt/geobolt.h src/core/geo_index_private.h
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
