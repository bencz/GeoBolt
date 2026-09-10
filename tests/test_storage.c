#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "geobolt/geodoc.h"
#include "geo_rocks_bridge.h"
#include "test_support.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_ASSERT(condition, message)                                                                                   \
    do {                                                                                                                  \
        if (!(condition)) {                                                                                               \
            fprintf(stderr, "FAIL: %s (%s:%d)\n", message, __FILE__, __LINE__);                                         \
            succeeded = false;                                                                                            \
            goto cleanup;                                                                                                 \
        }                                                                                                                 \
    } while (0)

static bool test_geodoc_nested_roundtrip_and_corruption(void)
{
    bool succeeded = true;
    GeoDocBuilder *builder = geo_doc_builder_create();
    GeoDocBuffer document = { 0 };
    GeoDocBuffer second_document = { 0 };
    uint32_t vehicle;
    uint32_t tags;

    TEST_ASSERT(builder != NULL, "GeoDoc builder must be created");
    TEST_ASSERT(geo_doc_builder_add_string(builder, 0U, "plate", "ABC1D23", 7U, NULL) == GEO_DOC_OK,
                "string field must be added");
    TEST_ASSERT(geo_doc_builder_add_object(builder, 0U, "vehicle", &vehicle) == GEO_DOC_OK,
                "nested object must be added");
    TEST_ASSERT(geo_doc_builder_add_uint64(builder, vehicle, "driver_id", 998877U, NULL) == GEO_DOC_OK,
                "nested integer must be added");
    TEST_ASSERT(geo_doc_builder_add_bool(builder, vehicle, "online", true, NULL) == GEO_DOC_OK,
                "nested boolean must be added");
    TEST_ASSERT(geo_doc_builder_add_array(builder, vehicle, "tags", &tags) == GEO_DOC_OK,
                "nested array must be added");
    TEST_ASSERT(geo_doc_builder_add_string(builder, tags, NULL, "priority", 8U, NULL) == GEO_DOC_OK,
                "array string must be added");
    TEST_ASSERT(geo_doc_builder_add_int64(builder, tags, NULL, -42, NULL) == GEO_DOC_OK,
                "array integer must be added");
    static const char invalid_utf8[] = { (char) 0xc0, (char) 0x80, '\0' };
    TEST_ASSERT(geo_doc_builder_add_string(builder, 0U, "invalid", invalid_utf8, 2U, NULL) == GEO_DOC_INVALID_ARGUMENT,
                "overlong UTF-8 must be rejected before serialization");
    TEST_ASSERT(geo_doc_builder_add_uint64(builder, vehicle, "driver_id", 1U, NULL) == GEO_DOC_DUPLICATE_KEY,
                "duplicate object keys must be rejected");
    TEST_ASSERT(geo_doc_builder_finish(builder, &document) == GEO_DOC_OK, "GeoDoc must serialize");
    TEST_ASSERT(geo_doc_builder_finish(builder, &second_document) == GEO_DOC_OK, "GeoDoc serialization must be repeatable");
    TEST_ASSERT(document.size == second_document.size && memcmp(document.data, second_document.data, document.size) == 0,
                "canonical serialization must be deterministic");

    GeoDocView view;
    GeoDocValue value;
    uint64_t driver_id;
    bool online;
    const void *tag;
    size_t tag_size;

    TEST_ASSERT(geo_doc_open(document.data, document.size, &view) == GEO_DOC_OK, "serialized GeoDoc must validate");
    TEST_ASSERT(geo_doc_find_pointer(view, "/vehicle/driver_id", &value) == GEO_DOC_OK,
                "nested JSON Pointer must resolve");
    TEST_ASSERT(geo_doc_value_uint64(value, &driver_id) == GEO_DOC_OK && driver_id == 998877U,
                "nested uint64 must round-trip exactly");
    TEST_ASSERT(geo_doc_find_pointer(view, "/vehicle/online", &value) == GEO_DOC_OK &&
                    geo_doc_value_bool(value, &online) == GEO_DOC_OK && online,
                "nested bool must round-trip exactly");
    TEST_ASSERT(geo_doc_find_pointer(view, "/vehicle/tags/0", &value) == GEO_DOC_OK &&
                    geo_doc_value_data(value, &tag, &tag_size) == GEO_DOC_OK && tag_size == 8U &&
                    memcmp(tag, "priority", tag_size) == 0,
                "array pointer must resolve exact string bytes");
    TEST_ASSERT(geo_doc_find_pointer(view, "/vehicle/missing", &value) == GEO_DOC_NOT_FOUND,
                "missing path must return NOT_FOUND");
    TEST_ASSERT(geo_doc_find_pointer(view, "/vehicle/~", &value) == GEO_DOC_INVALID_ARGUMENT,
                "truncated JSON Pointer escape must be rejected without releasing stack storage");

    ((unsigned char *) document.data)[document.size - 1U] ^= 1U;
    TEST_ASSERT(geo_doc_open(document.data, document.size, &view) == GEO_DOC_INVALID_FORMAT,
                "checksum mismatch must reject a corrupted GeoDoc");

cleanup:
    geo_doc_buffer_release(&second_document);
    geo_doc_buffer_release(&document);
    geo_doc_builder_destroy(builder);
    return succeeded;
}

static bool test_rocksdb_atomic_batch_snapshot_and_reopen(void)
{
    bool succeeded = true;
    char directory_template[] = "/tmp/geobolt-rocks-test-XXXXXX";
    char *directory = mkdtemp(directory_template);
    GeoRocksDatabase *database = NULL;
    GeoRocksBatch *batch = NULL;
    GeoRocksSnapshot *snapshot = NULL;
    GeoRocksIterator *iterator = NULL;
    GeoRocksMultiGet *multi_get = NULL;
    GeoRocksBuffer value = { 0 };
    GeoRocksStatus status;
    const unsigned char key1[] = { 0U, 0U, 0U, 0U, 0U, 0U, 0U, 1U };
    const unsigned char key2[] = { 0U, 0U, 0U, 0U, 0U, 0U, 0U, 2U };

    TEST_ASSERT(directory != NULL, "temporary RocksDB directory must be created");

    GeoRocksConfig config = geo_rocks_default_config();
    config.block_cache_bytes = 8U * 1024U * 1024U;
    config.write_buffer_bytes = 4U * 1024U * 1024U;
    config.background_jobs = 2;
    database = geo_rocks_open(directory, &config, &status);
    TEST_ASSERT(database != NULL && status.code == GEO_ROCKS_OK, "RocksDB must open with all GeoBolt column families");

    batch = geo_rocks_batch_create(database, 128U, &status);
    TEST_ASSERT(batch != NULL, "RocksDB batch must be created");
    TEST_ASSERT(geo_rocks_batch_put(batch, GEO_ROCKS_CF_OBJECTS, key1, sizeof(key1), "old", 3U, &status),
                "first object put must enter batch");
    TEST_ASSERT(geo_rocks_batch_put(batch, GEO_ROCKS_CF_CATALOG, "sequence", 8U, "1", 1U, &status),
                "catalog and object must share one batch");
    TEST_ASSERT(geo_rocks_write(database, batch, true, &status), "atomic synchronized write must commit");
    geo_rocks_batch_destroy(batch);
    batch = NULL;

    snapshot = geo_rocks_snapshot_create(database, &status);
    TEST_ASSERT(snapshot != NULL && geo_rocks_snapshot_sequence(snapshot) > 0U, "snapshot must expose a RocksDB sequence");

    batch = geo_rocks_batch_create(database, 128U, &status);
    TEST_ASSERT(batch != NULL, "replacement batch must be created");
    TEST_ASSERT(geo_rocks_batch_put(batch, GEO_ROCKS_CF_OBJECTS, key1, sizeof(key1), "new", 3U, &status),
                "replacement object must enter batch");
    TEST_ASSERT(geo_rocks_batch_put(batch, GEO_ROCKS_CF_OBJECTS, key2, sizeof(key2), "two", 3U, &status),
                "second object must enter batch");
    TEST_ASSERT(geo_rocks_write(database, batch, true, &status), "replacement batch must commit");
    geo_rocks_batch_destroy(batch);
    batch = NULL;

    TEST_ASSERT(geo_rocks_get(database, GEO_ROCKS_CF_OBJECTS, snapshot, key1, sizeof(key1), &value, &status),
                "snapshot read must find old value");
    TEST_ASSERT(value.size == 3U && memcmp(value.data, "old", 3U) == 0, "snapshot must preserve old value");
    geo_rocks_buffer_release(&value);
    TEST_ASSERT(geo_rocks_get(database, GEO_ROCKS_CF_OBJECTS, NULL, key1, sizeof(key1), &value, &status),
                "latest read must find replacement");
    TEST_ASSERT(value.size == 3U && memcmp(value.data, "new", 3U) == 0, "latest read must expose replacement");
    geo_rocks_buffer_release(&value);

    const uint64_t multi_keys[] = { 1U, 2U, 3U };
    const void *multi_value;
    size_t multi_value_size;

    multi_get = geo_rocks_multi_get_create(3U, &status);
    TEST_ASSERT(multi_get != NULL, "reusable RocksDB MultiGet workspace must be created");
    TEST_ASSERT(geo_rocks_multi_get_u64_be(database,
                                           GEO_ROCKS_CF_OBJECTS,
                                           snapshot,
                                           multi_keys,
                                           3U,
                                           multi_get,
                                           &status),
                "snapshot MultiGet must execute");
    TEST_ASSERT(geo_rocks_multi_get_result(multi_get, 0U, &multi_value, &multi_value_size) == GEO_ROCKS_OK &&
                    multi_value_size == 3U && memcmp(multi_value, "old", 3U) == 0,
                "snapshot MultiGet must retain the old first value without copying it");
    TEST_ASSERT(geo_rocks_multi_get_result(multi_get, 1U, &multi_value, &multi_value_size) == GEO_ROCKS_NOT_FOUND &&
                    multi_value == NULL && multi_value_size == 0U,
                "snapshot MultiGet must report a key created after the snapshot as missing");
    TEST_ASSERT(geo_rocks_multi_get_u64_be(database,
                                           GEO_ROCKS_CF_OBJECTS,
                                           NULL,
                                           multi_keys,
                                           3U,
                                           multi_get,
                                           &status),
                "latest MultiGet must reuse the pinned-value workspace");
    TEST_ASSERT(geo_rocks_multi_get_result(multi_get, 0U, &multi_value, &multi_value_size) == GEO_ROCKS_OK &&
                    multi_value_size == 3U && memcmp(multi_value, "new", 3U) == 0,
                "reused MultiGet must expose the replacement value");
    TEST_ASSERT(geo_rocks_multi_get_result(multi_get, 1U, &multi_value, &multi_value_size) == GEO_ROCKS_OK &&
                    multi_value_size == 3U && memcmp(multi_value, "two", 3U) == 0,
                "reused MultiGet must expose the second value");
    TEST_ASSERT(geo_rocks_multi_get_result(multi_get, 2U, &multi_value, &multi_value_size) == GEO_ROCKS_NOT_FOUND,
                "latest MultiGet must distinguish an actually missing key");
    geo_rocks_multi_get_release(multi_get);
    TEST_ASSERT(geo_rocks_multi_get_result(multi_get, 0U, &multi_value, &multi_value_size) == GEO_ROCKS_INVALID_ARGUMENT,
                "released MultiGet pins must no longer expose database-owned values");

    iterator = geo_rocks_iterator_create(database, GEO_ROCKS_CF_OBJECTS, NULL, &status);
    TEST_ASSERT(iterator != NULL, "object iterator must be created");
    geo_rocks_iterator_seek_first(iterator);
    size_t visited = 0U;

    while (geo_rocks_iterator_valid(iterator)) {
        size_t key_size;
        size_t value_size;

        TEST_ASSERT(geo_rocks_iterator_key(iterator, &key_size) != NULL && key_size == sizeof(key1),
                    "iterator key must remain valid until movement");
        TEST_ASSERT(geo_rocks_iterator_value(iterator, &value_size) != NULL && value_size == 3U,
                    "iterator value must remain valid until movement");
        visited++;
        geo_rocks_iterator_next(iterator);
    }

    TEST_ASSERT(geo_rocks_iterator_status(iterator, &status) && visited == 2U,
                "iterator must visit exactly both committed objects");
    geo_rocks_iterator_destroy(iterator);
    iterator = NULL;
    geo_rocks_snapshot_destroy(snapshot);
    snapshot = NULL;
    geo_rocks_close(database);
    database = NULL;

    config.create_if_missing = false;
    database = geo_rocks_open(directory, &config, &status);
    TEST_ASSERT(database != NULL, "existing RocksDB must reopen without creating state");
    TEST_ASSERT(geo_rocks_get(database, GEO_ROCKS_CF_OBJECTS, NULL, key2, sizeof(key2), &value, &status),
                "reopened database must retain synchronized batch");
    TEST_ASSERT(value.size == 3U && memcmp(value.data, "two", 3U) == 0, "reopened value must remain exact");

cleanup:
    geo_rocks_buffer_release(&value);
    geo_rocks_multi_get_destroy(multi_get);
    geo_rocks_iterator_destroy(iterator);
    geo_rocks_snapshot_destroy(snapshot);
    geo_rocks_batch_destroy(batch);
    geo_rocks_close(database);

    if (directory && !geobolt_test_remove_tree(directory)) {
        fprintf(stderr, "FAIL: RocksDB temporary tree could not be removed\n");
        succeeded = false;
    }

    return succeeded;
}

int main(void)
{
    printf("========================================\n");
    printf("GEOBOLT STORAGE TEST SUITE\n");
    printf("========================================\n\n");

    if (!test_geodoc_nested_roundtrip_and_corruption()) {
        return EXIT_FAILURE;
    }
    printf("[PASS] GeoDoc nested canonical round-trip and corruption rejection\n");

    if (!test_rocksdb_atomic_batch_snapshot_and_reopen()) {
        return EXIT_FAILURE;
    }
    printf("[PASS] RocksDB atomic batch, snapshot isolation, pinned MultiGet reuse, iteration, and reopen\n");
    return EXIT_SUCCESS;
}
