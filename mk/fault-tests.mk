# GNU/ELF --wrap must observe the real bridge boundary. Compile only the selected
# production caller without LTO; the test does not replace an indexing algorithm.
ifeq ($(shell uname -s),Linux)

define INDEX_FAILURE_PROFILE
$(BUILD_DIR)/secondary_fault_$(1).o: $(SECONDARY_INDEX_SRC) $(wildcard src/engine/*.h) | dirs
	$$(CC) $$(CPPFLAGS) $$(CFLAGS) $(2) -fno-lto -MMD -MP -c $$< -o $$@

$(BUILD_DIR)/commit_fault_$(1).o: $(COMMIT_COORDINATOR_SRC) $(wildcard src/engine/*.h) | dirs
	$$(CC) $$(CPPFLAGS) $$(CFLAGS) $(2) -fno-lto -MMD -MP -c $$< -o $$@

$(BIN_DIR)/test_database_failures$(3): tests/test_database_failures.c $(TEST_SUPPORT_SRC) \
        $(BUILD_DIR)/secondary_fault_$(1).o $(BUILD_DIR)/commit_fault_$(1).o $(4) | dirs
	$$(CC) $$(CPPFLAGS) -Itests $$(CFLAGS) $(2) tests/test_database_failures.c $$(TEST_SUPPORT_SRC) \
		$(BUILD_DIR)/secondary_fault_$(1).o $(BUILD_DIR)/commit_fault_$(1).o $(4) $$(LDFLAGS) \
		-Wl,--wrap=geo_rocks_iterator_create,--wrap=geo_rocks_write,--wrap=pthread_cond_timedwait \
		$(5) $$(THREAD_LIBS) -o $$@

$(6): $(BIN_DIR)/test_database_failures$(3)
	@$(BIN_DIR)/test_database_failures$(3)

-include $(BUILD_DIR)/secondary_fault_$(1).d
-include $(BUILD_DIR)/commit_fault_$(1).d
endef

$(eval $(call INDEX_FAILURE_PROFILE,opt,$(CFLAGS_OPT),,$(LIB_STATIC),$(OPT_LDLIBS),test-failures))
$(eval $(call INDEX_FAILURE_PROFILE,debug,$(CFLAGS_DEBUG),_debug,\
    $(filter-out $(SECONDARY_INDEX_OBJ_DEBUG) $(COMMIT_COORDINATOR_OBJ_DEBUG),$(DEBUG_LIBRARY_OBJECTS)),\
    $(LDLIBS),test-failures-debug))
$(eval $(call INDEX_FAILURE_PROFILE,scalar,$(CFLAGS_OPT) -DGEO_SIMD_FORCE_SCALAR,_scalar,\
    $(filter-out $(SECONDARY_INDEX_OBJ_SCALAR) $(COMMIT_COORDINATOR_OBJ_SCALAR),$(SCALAR_LIBRARY_OBJECTS)),\
    $(OPT_LDLIBS),test-failures-scalar))

test: test-failures
test-debug: test-failures-debug
test-scalar: test-failures-scalar

.PHONY: test-failures test-failures-debug test-failures-scalar
endif
