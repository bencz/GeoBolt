#include "geo_rocks_bridge.h"

#include <rocksdb/cache.h>
#include <rocksdb/compaction_filter.h>
#include <rocksdb/db.h>
#include <rocksdb/filter_policy.h>
#include <rocksdb/options.h>
#include <rocksdb/slice.h>
#include <rocksdb/statistics.h>
#include <rocksdb/table.h>
#include <rocksdb/write_batch.h>
#include <rocksdb/write_buffer_manager.h>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <memory>
#include <new>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

constexpr std::array<const char *, GEO_ROCKS_CF_COUNT> kColumnFamilyNames = {
    "catalog",
    "objects",
    "spatial_delta",
    "idempotency",
    "secondary_index",
};

constexpr size_t kIdempotencyValueSize = 48U;
constexpr size_t kIdempotencyExpirationOffset = 32U;

uint64_t load_u64_le(const char *bytes)
{
    uint64_t value = 0U;

    for (size_t index = 0; index < sizeof(value); ++index) {
        value |= static_cast<uint64_t>(static_cast<unsigned char>(bytes[index])) << (index * 8U);
    }

    return value;
}

class IdempotencyExpirationFilter final : public rocksdb::CompactionFilter {
public:
    explicit IdempotencyExpirationFilter(std::time_t now)
        : now_(now)
    {
    }

    bool Filter(int,
                const rocksdb::Slice &,
                const rocksdb::Slice &existing_value,
                std::string *,
                bool *) const override
    {
        if (existing_value.size() != kIdempotencyValueSize) {
            return false;
        }

        if (now_ < 0) {
            return false;
        }

        uint64_t expires_at = load_u64_le(existing_value.data() + kIdempotencyExpirationOffset);
        return expires_at != 0U && expires_at <= static_cast<uint64_t>(now_);
    }

    const char *Name() const override
    {
        return "GeoBoltIdempotencyExpirationFilter";
    }

private:
    std::time_t now_;
};

class IdempotencyExpirationFilterFactory final : public rocksdb::CompactionFilterFactory {
public:
    std::unique_ptr<rocksdb::CompactionFilter> CreateCompactionFilter(
        const rocksdb::CompactionFilter::Context &) override
    {
        return std::make_unique<IdempotencyExpirationFilter>(std::time(nullptr));
    }

    const char *Name() const override
    {
        return "GeoBoltIdempotencyExpirationFilterFactory";
    }
};

void clear_status(GeoRocksStatus *status)
{
    if (status == nullptr) {
        return;
    }

    status->code = GEO_ROCKS_OK;
    status->message[0] = '\0';
}

void set_status(GeoRocksStatus *status, GeoRocksStatusCode code, const char *message)
{
    if (status == nullptr) {
        return;
    }

    status->code = code;

    if (message == nullptr) {
        status->message[0] = '\0';
        return;
    }

    std::strncpy(status->message, message, sizeof(status->message) - 1U);
    status->message[sizeof(status->message) - 1U] = '\0';
}

GeoRocksStatusCode translate_status(const rocksdb::Status &status)
{
    if (status.ok()) {
        return GEO_ROCKS_OK;
    }
    if (status.IsNotFound()) {
        return GEO_ROCKS_NOT_FOUND;
    }
    if (status.IsInvalidArgument()) {
        return GEO_ROCKS_INVALID_ARGUMENT;
    }
    if (status.IsCorruption()) {
        return GEO_ROCKS_CORRUPTION;
    }
    if (status.IsBusy() || status.IsIncomplete() || status.IsTryAgain()) {
        return GEO_ROCKS_BUSY;
    }
    if (status.IsIOError()) {
        return GEO_ROCKS_IO_ERROR;
    }

    return GEO_ROCKS_INTERNAL_ERROR;
}

bool export_status(const rocksdb::Status &source, GeoRocksStatus *destination)
{
    if (source.ok()) {
        clear_status(destination);
        return true;
    }

    std::string message = source.ToString();
    set_status(destination, translate_status(source), message.c_str());
    return false;
}

rocksdb::Slice make_slice(const void *data, size_t size)
{
    return rocksdb::Slice(static_cast<const char *>(data), size);
}

bool valid_column_family(GeoRocksColumnFamily column_family)
{
    return column_family >= GEO_ROCKS_CF_CATALOG && column_family < GEO_ROCKS_CF_COUNT;
}

rocksdb::ColumnFamilyOptions make_column_family_options(const GeoRocksConfig &config,
                                                        const std::shared_ptr<rocksdb::Cache> &cache,
                                                        bool point_lookup,
                                                        bool expire_idempotency = false)
{
    rocksdb::ColumnFamilyOptions options;
    rocksdb::BlockBasedTableOptions table_options;

    table_options.block_cache = cache;
    table_options.cache_index_and_filter_blocks = true;
    table_options.cache_index_and_filter_blocks_with_high_priority = true;
    table_options.pin_l0_filter_and_index_blocks_in_cache = true;
    table_options.filter_policy.reset(rocksdb::NewBloomFilterPolicy(point_lookup ? 10.0 : 8.0, false));
    options.table_factory.reset(rocksdb::NewBlockBasedTableFactory(table_options));
    options.write_buffer_size = std::max<size_t>(config.write_buffer_bytes / 4U, 4U * 1024U * 1024U);
    options.max_write_buffer_number = 4;
    options.min_write_buffer_number_to_merge = 2;
    options.level_compaction_dynamic_level_bytes = true;

    if (expire_idempotency) {
        options.compaction_filter_factory = std::make_shared<IdempotencyExpirationFilterFactory>();
    }

    return options;
}

}  // namespace

struct GeoRocksDatabase {
    rocksdb::DB *database = nullptr;
    std::array<rocksdb::ColumnFamilyHandle *, GEO_ROCKS_CF_COUNT> column_families = {};
    rocksdb::ColumnFamilyHandle *default_column_family = nullptr;
    std::shared_ptr<rocksdb::Cache> block_cache;
    std::shared_ptr<rocksdb::WriteBufferManager> write_buffer_manager;
    std::shared_ptr<rocksdb::Statistics> statistics;
};

struct GeoRocksBatch {
    GeoRocksBatch(GeoRocksDatabase *database_value, size_t reserved_bytes)
        : database(database_value), batch(reserved_bytes)
    {
    }

    GeoRocksDatabase *database;
    rocksdb::WriteBatch batch;
};

struct GeoRocksSnapshot {
    GeoRocksDatabase *owner = nullptr;
    const rocksdb::Snapshot *snapshot = nullptr;
};

struct GeoRocksIterator {
    std::unique_ptr<rocksdb::Iterator> iterator;
};

struct GeoRocksMultiGet {
    explicit GeoRocksMultiGet(size_t capacity)
        : encoded_keys(capacity), keys(capacity), values(capacity), statuses(capacity)
    {
    }

    std::vector<std::array<char, sizeof(uint64_t)>> encoded_keys;
    std::vector<rocksdb::Slice> keys;
    std::vector<rocksdb::PinnableSlice> values;
    std::vector<rocksdb::Status> statuses;
    size_t result_count = 0U;
};

static void multi_get_reset_values(GeoRocksMultiGet *multi_get, size_t count)
{
    for (size_t index = 0U; index < count; ++index) {
        multi_get->values[index].Reset();
    }
}

extern "C" GeoRocksConfig geo_rocks_default_config(void)
{
    const unsigned hardware_threads = std::thread::hardware_concurrency();
    const unsigned background_jobs = hardware_threads == 0U ? 4U : std::min(hardware_threads, 16U);

    return GeoRocksConfig {
        256U * 1024U * 1024U,
        64U * 1024U * 1024U,
        static_cast<int>(background_jobs),
        true,
    };
}

extern "C" GeoRocksDatabase *geo_rocks_open(const char *directory,
                                              const GeoRocksConfig *config,
                                              GeoRocksStatus *status)
{
    clear_status(status);

    if (directory == nullptr || directory[0] == '\0') {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "RocksDB directory is empty");
        return nullptr;
    }

    try {
        const GeoRocksConfig selected = config == nullptr ? geo_rocks_default_config() : *config;

        if (selected.block_cache_bytes == 0U || selected.write_buffer_bytes == 0U || selected.background_jobs <= 0) {
            set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "RocksDB resource limits must be positive");
            return nullptr;
        }

        std::unique_ptr<GeoRocksDatabase> database(new GeoRocksDatabase());
        database->block_cache = rocksdb::NewLRUCache(selected.block_cache_bytes);
        database->write_buffer_manager = std::make_shared<rocksdb::WriteBufferManager>(selected.write_buffer_bytes,
                                                                                       nullptr,
                                                                                       true);
        database->statistics = rocksdb::CreateDBStatistics();

        rocksdb::DBOptions database_options;
        database_options.create_if_missing = selected.create_if_missing;
        database_options.create_missing_column_families = selected.create_if_missing;
        database_options.atomic_flush = true;
        database_options.max_background_jobs = selected.background_jobs;
        database_options.statistics = database->statistics;
        database_options.write_buffer_manager = database->write_buffer_manager;

        std::vector<rocksdb::ColumnFamilyDescriptor> descriptors;
        descriptors.reserve(GEO_ROCKS_CF_COUNT + 1U);
        descriptors.emplace_back(rocksdb::kDefaultColumnFamilyName,
                                 make_column_family_options(selected, database->block_cache, true));

        for (size_t index = 0; index < GEO_ROCKS_CF_COUNT; ++index) {
            const bool point_lookup = index == GEO_ROCKS_CF_CATALOG || index == GEO_ROCKS_CF_OBJECTS ||
                                      index == GEO_ROCKS_CF_IDEMPOTENCY;
            descriptors.emplace_back(kColumnFamilyNames[index],
                                     make_column_family_options(selected,
                                                                database->block_cache,
                                                                point_lookup,
                                                                index == GEO_ROCKS_CF_IDEMPOTENCY));
        }

        std::vector<rocksdb::ColumnFamilyHandle *> handles;
        rocksdb::Status open_status = rocksdb::DB::Open(database_options, directory, descriptors, &handles, &database->database);

        if (!export_status(open_status, status)) {
            return nullptr;
        }

        if (handles.size() != descriptors.size()) {
            set_status(status, GEO_ROCKS_INTERNAL_ERROR, "RocksDB returned an unexpected column family count");

            for (rocksdb::ColumnFamilyHandle *handle : handles) {
                delete handle;
            }

            delete database->database;
            database->database = nullptr;
            return nullptr;
        }

        database->default_column_family = handles[0];

        for (size_t index = 0; index < GEO_ROCKS_CF_COUNT; ++index) {
            database->column_families[index] = handles[index + 1U];
        }

        return database.release();
    } catch (const std::bad_alloc &) {
        set_status(status, GEO_ROCKS_OUT_OF_MEMORY, "RocksDB bridge allocation failed");
    } catch (const std::exception &exception) {
        set_status(status, GEO_ROCKS_INTERNAL_ERROR, exception.what());
    } catch (...) {
        set_status(status, GEO_ROCKS_INTERNAL_ERROR, "Unknown exception in RocksDB bridge");
    }

    return nullptr;
}

extern "C" void geo_rocks_close(GeoRocksDatabase *database)
{
    if (database == nullptr) {
        return;
    }

    for (rocksdb::ColumnFamilyHandle *handle : database->column_families) {
        delete handle;
    }

    delete database->default_column_family;
    delete database->database;
    delete database;
}

extern "C" GeoRocksBatch *geo_rocks_batch_create(GeoRocksDatabase *database,
                                                  size_t reserved_bytes,
                                                  GeoRocksStatus *status)
{
    clear_status(status);

    if (database == nullptr) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB batch database");
        return nullptr;
    }

    try {
        return new GeoRocksBatch(database, reserved_bytes);
    } catch (const std::bad_alloc &) {
        set_status(status, GEO_ROCKS_OUT_OF_MEMORY, "RocksDB batch allocation failed");
    } catch (const std::exception &exception) {
        set_status(status, GEO_ROCKS_INTERNAL_ERROR, exception.what());
    } catch (...) {
        set_status(status, GEO_ROCKS_INTERNAL_ERROR, "Unknown exception while creating RocksDB batch");
    }

    return nullptr;
}

extern "C" void geo_rocks_batch_destroy(GeoRocksBatch *batch)
{
    delete batch;
}

extern "C" bool geo_rocks_batch_put(GeoRocksBatch *batch,
                                     GeoRocksColumnFamily column_family,
                                     const void *key,
                                     size_t key_size,
                                     const void *value,
                                     size_t value_size,
                                     GeoRocksStatus *status)
{
    clear_status(status);

    if (batch == nullptr || !valid_column_family(column_family) || (key == nullptr && key_size != 0U) ||
        (value == nullptr && value_size != 0U)) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB batch put argument");
        return false;
    }

    return export_status(batch->batch.Put(batch->database->column_families[static_cast<size_t>(column_family)],
                                          make_slice(key, key_size),
                                          make_slice(value, value_size)),
                         status);
}

extern "C" bool geo_rocks_batch_delete(GeoRocksBatch *batch,
                                        GeoRocksColumnFamily column_family,
                                        const void *key,
                                        size_t key_size,
                                        GeoRocksStatus *status)
{
    clear_status(status);

    if (batch == nullptr || !valid_column_family(column_family) || (key == nullptr && key_size != 0U)) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB batch delete argument");
        return false;
    }

    return export_status(batch->batch.Delete(batch->database->column_families[static_cast<size_t>(column_family)],
                                             make_slice(key, key_size)),
                         status);
}

extern "C" bool geo_rocks_batch_delete_range(GeoRocksBatch *batch,
                                              GeoRocksColumnFamily column_family,
                                              const void *begin_key,
                                              size_t begin_key_size,
                                              const void *end_key,
                                              size_t end_key_size,
                                              GeoRocksStatus *status)
{
    clear_status(status);

    if (batch == nullptr || !valid_column_family(column_family) || (begin_key == nullptr && begin_key_size != 0U) ||
        (end_key == nullptr && end_key_size != 0U)) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB delete range argument");
        return false;
    }

    return export_status(batch->batch.DeleteRange(batch->database->column_families[static_cast<size_t>(column_family)],
                                                  make_slice(begin_key, begin_key_size),
                                                  make_slice(end_key, end_key_size)),
                         status);
}

extern "C" bool geo_rocks_write(GeoRocksDatabase *database,
                                 GeoRocksBatch *batch,
                                 bool synchronize,
                                 GeoRocksStatus *status)
{
    clear_status(status);

    if (database == nullptr || batch == nullptr || batch->database != database) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB write argument");
        return false;
    }

    rocksdb::WriteOptions options;
    options.sync = synchronize;
    options.disableWAL = false;

    return export_status(database->database->Write(options, &batch->batch), status);
}

extern "C" bool geo_rocks_get(GeoRocksDatabase *database,
                               GeoRocksColumnFamily column_family,
                               const GeoRocksSnapshot *snapshot,
                               const void *key,
                               size_t key_size,
                               GeoRocksBuffer *value,
                               GeoRocksStatus *status)
{
    clear_status(status);

    if (database == nullptr || !valid_column_family(column_family) || (key == nullptr && key_size != 0U) || value == nullptr ||
        (snapshot != nullptr && snapshot->owner != database)) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB get argument");
        return false;
    }

    value->data = nullptr;
    value->size = 0U;
    rocksdb::ReadOptions options;
    options.snapshot = snapshot == nullptr ? nullptr : snapshot->snapshot;
    std::string result;
    rocksdb::Status get_status = database->database->Get(options,
                                                         database->column_families[static_cast<size_t>(column_family)],
                                                         make_slice(key, key_size),
                                                         &result);

    if (!export_status(get_status, status)) {
        return false;
    }

    if (result.empty()) {
        return true;
    }

    void *storage = std::malloc(result.size());

    if (storage == nullptr) {
        set_status(status, GEO_ROCKS_OUT_OF_MEMORY, "RocksDB result allocation failed");
        return false;
    }

    std::memcpy(storage, result.data(), result.size());

    value->data = storage;
    value->size = result.size();
    return true;
}

extern "C" void geo_rocks_buffer_release(GeoRocksBuffer *buffer)
{
    if (buffer == nullptr) {
        return;
    }

    std::free(buffer->data);
    buffer->data = nullptr;
    buffer->size = 0U;
}

extern "C" GeoRocksMultiGet *geo_rocks_multi_get_create(size_t capacity, GeoRocksStatus *status)
{
    clear_status(status);

    if (capacity == 0U) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "RocksDB MultiGet capacity must be positive");
        return nullptr;
    }

    try {
        return new GeoRocksMultiGet(capacity);
    } catch (const std::bad_alloc &) {
        set_status(status, GEO_ROCKS_OUT_OF_MEMORY, "RocksDB MultiGet workspace allocation failed");
    } catch (const std::exception &exception) {
        set_status(status, GEO_ROCKS_INTERNAL_ERROR, exception.what());
    } catch (...) {
        set_status(status, GEO_ROCKS_INTERNAL_ERROR, "Unknown exception while creating RocksDB MultiGet workspace");
    }

    return nullptr;
}

extern "C" void geo_rocks_multi_get_destroy(GeoRocksMultiGet *multi_get)
{
    delete multi_get;
}

extern "C" void geo_rocks_multi_get_release(GeoRocksMultiGet *multi_get)
{
    if (multi_get == nullptr) {
        return;
    }

    multi_get_reset_values(multi_get, multi_get->result_count);
    multi_get->result_count = 0U;
}

extern "C" bool geo_rocks_multi_get_u64_be(GeoRocksDatabase *database,
                                             GeoRocksColumnFamily column_family,
                                             const GeoRocksSnapshot *snapshot,
                                             const uint64_t *keys,
                                             size_t key_count,
                                             GeoRocksMultiGet *multi_get,
                                             GeoRocksStatus *status)
{
    clear_status(status);

    if (database == nullptr || !valid_column_family(column_family) || keys == nullptr || key_count == 0U ||
        multi_get == nullptr || key_count > multi_get->keys.size() ||
        (snapshot != nullptr && snapshot->owner != database)) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB MultiGet argument");
        return false;
    }

    try {
        size_t reset_count = std::max(key_count, multi_get->result_count);
        bool sorted_input = true;

        multi_get_reset_values(multi_get, reset_count);

        for (size_t index = 0U; index < key_count; ++index) {
            uint64_t key = keys[index];

            for (size_t byte = 0U; byte < sizeof(key); ++byte) {
                size_t shift = (sizeof(key) - byte - 1U) * 8U;
                multi_get->encoded_keys[index][byte] = static_cast<char>((key >> shift) & UINT64_C(0xff));
            }

            multi_get->keys[index] = rocksdb::Slice(multi_get->encoded_keys[index].data(), sizeof(key));
            sorted_input = sorted_input && (index == 0U || keys[index - 1U] <= key);
        }

        rocksdb::ReadOptions options;
        options.snapshot = snapshot == nullptr ? nullptr : snapshot->snapshot;
        database->database->MultiGet(options,
                                     database->column_families[static_cast<size_t>(column_family)],
                                     key_count,
                                     multi_get->keys.data(),
                                     multi_get->values.data(),
                                     multi_get->statuses.data(),
                                     sorted_input);
        multi_get->result_count = key_count;
        return true;
    } catch (const std::bad_alloc &) {
        multi_get_reset_values(multi_get, key_count);
        multi_get->result_count = 0U;
        set_status(status, GEO_ROCKS_OUT_OF_MEMORY, "RocksDB MultiGet execution allocation failed");
    } catch (const std::exception &exception) {
        multi_get_reset_values(multi_get, key_count);
        multi_get->result_count = 0U;
        set_status(status, GEO_ROCKS_INTERNAL_ERROR, exception.what());
    } catch (...) {
        multi_get_reset_values(multi_get, key_count);
        multi_get->result_count = 0U;
        set_status(status, GEO_ROCKS_INTERNAL_ERROR, "Unknown exception during RocksDB MultiGet");
    }

    return false;
}

extern "C" GeoRocksStatusCode geo_rocks_multi_get_result(const GeoRocksMultiGet *multi_get,
                                                           size_t index,
                                                           const void **value,
                                                           size_t *value_size)
{
    if (multi_get == nullptr || index >= multi_get->result_count || value == nullptr || value_size == nullptr) {
        return GEO_ROCKS_INVALID_ARGUMENT;
    }

    const rocksdb::Status &status = multi_get->statuses[index];

    if (!status.ok()) {
        *value = nullptr;
        *value_size = 0U;
        return translate_status(status);
    }

    *value = multi_get->values[index].data();
    *value_size = multi_get->values[index].size();
    return GEO_ROCKS_OK;
}

extern "C" GeoRocksSnapshot *geo_rocks_snapshot_create(GeoRocksDatabase *database, GeoRocksStatus *status)
{
    clear_status(status);

    if (database == nullptr) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB snapshot database");
        return nullptr;
    }

    try {
        std::unique_ptr<GeoRocksSnapshot> snapshot(new GeoRocksSnapshot());
        snapshot->owner = database;
        snapshot->snapshot = database->database->GetSnapshot();

        if (snapshot->snapshot == nullptr) {
            set_status(status, GEO_ROCKS_INTERNAL_ERROR, "RocksDB did not create a snapshot");
            return nullptr;
        }

        return snapshot.release();
    } catch (const std::bad_alloc &) {
        set_status(status, GEO_ROCKS_OUT_OF_MEMORY, "RocksDB snapshot allocation failed");
    }

    return nullptr;
}

extern "C" void geo_rocks_snapshot_destroy(GeoRocksSnapshot *snapshot)
{
    if (snapshot == nullptr) {
        return;
    }

    snapshot->owner->database->ReleaseSnapshot(snapshot->snapshot);
    delete snapshot;
}

extern "C" uint64_t geo_rocks_snapshot_sequence(const GeoRocksSnapshot *snapshot)
{
    return snapshot == nullptr ? 0U : snapshot->snapshot->GetSequenceNumber();
}

extern "C" GeoRocksIterator *geo_rocks_iterator_create(GeoRocksDatabase *database,
                                                        GeoRocksColumnFamily column_family,
                                                        const GeoRocksSnapshot *snapshot,
                                                        GeoRocksStatus *status)
{
    clear_status(status);

    if (database == nullptr || !valid_column_family(column_family) || (snapshot != nullptr && snapshot->owner != database)) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB iterator argument");
        return nullptr;
    }

    try {
        rocksdb::ReadOptions options;
        options.snapshot = snapshot == nullptr ? nullptr : snapshot->snapshot;
        std::unique_ptr<GeoRocksIterator> iterator(new GeoRocksIterator());
        iterator->iterator.reset(database->database->NewIterator(options,
                                                                 database->column_families[static_cast<size_t>(column_family)]));

        if (!iterator->iterator) {
            set_status(status, GEO_ROCKS_INTERNAL_ERROR, "RocksDB did not create an iterator");
            return nullptr;
        }

        return iterator.release();
    } catch (const std::bad_alloc &) {
        set_status(status, GEO_ROCKS_OUT_OF_MEMORY, "RocksDB iterator allocation failed");
    }

    return nullptr;
}

extern "C" void geo_rocks_iterator_destroy(GeoRocksIterator *iterator)
{
    delete iterator;
}

extern "C" void geo_rocks_iterator_seek_first(GeoRocksIterator *iterator)
{
    if (iterator != nullptr) {
        iterator->iterator->SeekToFirst();
    }
}

extern "C" void geo_rocks_iterator_seek(GeoRocksIterator *iterator, const void *key, size_t key_size)
{
    if (iterator != nullptr && (key != nullptr || key_size == 0U)) {
        iterator->iterator->Seek(make_slice(key, key_size));
    }
}

extern "C" void geo_rocks_iterator_next(GeoRocksIterator *iterator)
{
    if (iterator != nullptr && iterator->iterator->Valid()) {
        iterator->iterator->Next();
    }
}

extern "C" bool geo_rocks_iterator_valid(const GeoRocksIterator *iterator)
{
    return iterator != nullptr && iterator->iterator->Valid();
}

extern "C" const void *geo_rocks_iterator_key(const GeoRocksIterator *iterator, size_t *size)
{
    if (!geo_rocks_iterator_valid(iterator) || size == nullptr) {
        return nullptr;
    }

    rocksdb::Slice key = iterator->iterator->key();
    *size = key.size();
    return key.data();
}

extern "C" const void *geo_rocks_iterator_value(const GeoRocksIterator *iterator, size_t *size)
{
    if (!geo_rocks_iterator_valid(iterator) || size == nullptr) {
        return nullptr;
    }

    rocksdb::Slice value = iterator->iterator->value();
    *size = value.size();
    return value.data();
}

extern "C" bool geo_rocks_iterator_status(const GeoRocksIterator *iterator, GeoRocksStatus *status)
{
    if (iterator == nullptr) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB iterator");
        return false;
    }

    return export_status(iterator->iterator->status(), status);
}

extern "C" bool geo_rocks_flush(GeoRocksDatabase *database, bool wait, GeoRocksStatus *status)
{
    clear_status(status);

    if (database == nullptr) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB flush database");
        return false;
    }

    rocksdb::FlushOptions options;
    options.wait = wait;
    std::vector<rocksdb::ColumnFamilyHandle *> handles(database->column_families.begin(), database->column_families.end());

    return export_status(database->database->Flush(options, handles), status);
}

extern "C" bool geo_rocks_get_property_u64(GeoRocksDatabase *database,
                                            GeoRocksColumnFamily column_family,
                                            const char *property,
                                            uint64_t *value,
                                            GeoRocksStatus *status)
{
    clear_status(status);

    if (database == nullptr || !valid_column_family(column_family) || property == nullptr || value == nullptr) {
        set_status(status, GEO_ROCKS_INVALID_ARGUMENT, "Invalid RocksDB property argument");
        return false;
    }

    bool found = database->database->GetIntProperty(database->column_families[static_cast<size_t>(column_family)], property, value);

    if (!found) {
        set_status(status, GEO_ROCKS_NOT_FOUND, "RocksDB integer property is unavailable");
    }

    return found;
}
