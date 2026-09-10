#include "geo_db_format.h"

#include "geobolt/geodoc.h"

#include <stdlib.h>
#include <string.h>

#define GEO_DB_OBJECT_MAGIC UINT32_C(0x314f4247)
#define GEO_DB_DELTA_MAGIC UINT32_C(0x31444247)
#define GEO_DB_IDEMPOTENCY_MAGIC UINT32_C(0x31494247)

static const char GEO_DB_CATALOG_FORMAT[] = "format";
static const char GEO_DB_CATALOG_NEXT_SEQUENCE[] = "next_sequence";
static const char GEO_DB_CATALOG_COMMITTED[] = "committed_operations";
static const char GEO_DB_CATALOG_APPLIED[] = "spatial_applied_sequence";
static const char GEO_DB_CATALOG_NEXT_SERVER_ID[] = "next_server_id";
static const char GEO_DB_CATALOG_NEXT_SECONDARY_INDEX_ID[] = "next_secondary_index_id";

static uint16_t db_load_u16_le(const unsigned char *input)
{
    return (uint16_t) input[0] | (uint16_t) ((uint16_t) input[1] << 8U);
}

static uint32_t db_load_u32_le(const unsigned char *input)
{
    return (uint32_t) input[0] |
           (uint32_t) input[1] << 8U |
           (uint32_t) input[2] << 16U |
           (uint32_t) input[3] << 24U;
}

static uint64_t db_load_u64_le(const unsigned char *input)
{
    return (uint64_t) db_load_u32_le(input) | (uint64_t) db_load_u32_le(input + 4U) << 32U;
}

static void db_store_u16_le(unsigned char *output, uint16_t value)
{
    output[0] = (unsigned char) value;
    output[1] = (unsigned char) (value >> 8U);
}

static void db_store_u32_le(unsigned char *output, uint32_t value)
{
    output[0] = (unsigned char) value;
    output[1] = (unsigned char) (value >> 8U);
    output[2] = (unsigned char) (value >> 16U);
    output[3] = (unsigned char) (value >> 24U);
}

static void db_store_u64_le(unsigned char *output, uint64_t value)
{
    db_store_u32_le(output, (uint32_t) value);
    db_store_u32_le(output + 4U, (uint32_t) (value >> 32U));
}

void geo_db_store_u64_be(unsigned char output[8], uint64_t value)
{
    for (size_t index = 0; index < 8U; ++index) {
        output[index] = (unsigned char) (value >> ((7U - index) * 8U));
    }
}

uint64_t geo_db_load_u64_be(const unsigned char input[8])
{
    uint64_t value = 0U;

    for (size_t index = 0; index < 8U; ++index) {
        value = value << 8U | input[index];
    }

    return value;
}

GeoDatabaseStatus geo_db_status_from_rocks(const GeoRocksStatus *status)
{
    if (!status || status->code == GEO_ROCKS_OK) {
        return GEO_DATABASE_OK;
    }

    switch (status->code) {
        case GEO_ROCKS_INVALID_ARGUMENT:
            return GEO_DATABASE_INVALID_ARGUMENT;
        case GEO_ROCKS_OUT_OF_MEMORY:
            return GEO_DATABASE_OUT_OF_MEMORY;
        case GEO_ROCKS_CORRUPTION:
            return GEO_DATABASE_CORRUPTION;
        case GEO_ROCKS_IO_ERROR:
            return GEO_DATABASE_IO_ERROR;
        case GEO_ROCKS_BUSY:
        case GEO_ROCKS_INTERNAL_ERROR:
        case GEO_ROCKS_NOT_FOUND:
        case GEO_ROCKS_OK:
            return GEO_DATABASE_FAILED_STATE;
    }

    return GEO_DATABASE_FAILED_STATE;
}

static bool catalog_get_u64(GeoRocksDatabase *rocks, const char *key, uint64_t *value, GeoRocksStatus *status)
{
    GeoRocksBuffer buffer = { 0 };
    bool found = geo_rocks_get(rocks, GEO_ROCKS_CF_CATALOG, NULL, key, strlen(key), &buffer, status);

    if (found && buffer.size == sizeof(uint64_t)) {
        *value = geo_db_load_u64_be(buffer.data);
    } else if (found) {
        status->code = GEO_ROCKS_CORRUPTION;
        found = false;
    }

    geo_rocks_buffer_release(&buffer);
    return found;
}

static bool catalog_put_u64(GeoRocksBatch *batch, const char *key, uint64_t value, GeoRocksStatus *status)
{
    unsigned char encoded[8];

    geo_db_store_u64_be(encoded, value);
    return geo_rocks_batch_put(batch, GEO_ROCKS_CF_CATALOG, key, strlen(key), encoded, sizeof(encoded), status);
}

bool geo_db_catalog_put(GeoRocksBatch *batch, const GeoDbCatalog *catalog, GeoDatabaseStatus *status)
{
    GeoRocksStatus rocks_status;
    bool succeeded = batch && catalog &&
                     catalog_put_u64(batch, GEO_DB_CATALOG_FORMAT, GEO_DB_FORMAT_VERSION, &rocks_status) &&
                     catalog_put_u64(batch, GEO_DB_CATALOG_NEXT_SEQUENCE, catalog->next_sequence, &rocks_status) &&
                     catalog_put_u64(batch, GEO_DB_CATALOG_COMMITTED, catalog->committed_operations, &rocks_status) &&
                     catalog_put_u64(batch, GEO_DB_CATALOG_APPLIED, catalog->spatial_applied_sequence, &rocks_status) &&
                     catalog_put_u64(batch, GEO_DB_CATALOG_NEXT_SERVER_ID, catalog->next_server_id, &rocks_status) &&
                     catalog_put_u64(batch,
                                     GEO_DB_CATALOG_NEXT_SECONDARY_INDEX_ID,
                                     catalog->next_secondary_index_id,
                                     &rocks_status);

    if (!succeeded && status) {
        *status = batch && catalog ? geo_db_status_from_rocks(&rocks_status) : GEO_DATABASE_INVALID_ARGUMENT;
    }

    return succeeded;
}

bool geo_db_catalog_load(GeoRocksDatabase *rocks,
                         bool create_if_missing,
                         GeoDbCatalog *catalog,
                         GeoDatabaseStatus *status)
{
    if (!rocks || !catalog) {
        if (status) {
            *status = GEO_DATABASE_INVALID_ARGUMENT;
        }
        return false;
    }

    GeoRocksStatus rocks_status;
    uint64_t format;

    if (catalog_get_u64(rocks, GEO_DB_CATALOG_FORMAT, &format, &rocks_status)) {
        bool succeeded = format == GEO_DB_FORMAT_VERSION &&
                         catalog_get_u64(rocks, GEO_DB_CATALOG_NEXT_SEQUENCE, &catalog->next_sequence, &rocks_status) &&
                         catalog_get_u64(rocks, GEO_DB_CATALOG_COMMITTED, &catalog->committed_operations, &rocks_status) &&
                         catalog_get_u64(rocks, GEO_DB_CATALOG_APPLIED, &catalog->spatial_applied_sequence, &rocks_status) &&
                         catalog_get_u64(rocks, GEO_DB_CATALOG_NEXT_SERVER_ID, &catalog->next_server_id, &rocks_status) &&
                         catalog_get_u64(rocks,
                                         GEO_DB_CATALOG_NEXT_SECONDARY_INDEX_ID,
                                         &catalog->next_secondary_index_id,
                                         &rocks_status) &&
                         catalog->next_sequence > 0U && catalog->spatial_applied_sequence > 0U &&
                         catalog->spatial_applied_sequence <= catalog->next_sequence && catalog->next_secondary_index_id > 0U;

        if (!succeeded && status) {
            *status = rocks_status.code == GEO_ROCKS_OK ? GEO_DATABASE_CORRUPTION : geo_db_status_from_rocks(&rocks_status);
        }
        return succeeded;
    }

    if (rocks_status.code != GEO_ROCKS_NOT_FOUND || !create_if_missing) {
        if (status) {
            *status = rocks_status.code == GEO_ROCKS_NOT_FOUND ? GEO_DATABASE_IO_ERROR : geo_db_status_from_rocks(&rocks_status);
        }
        return false;
    }

    *catalog = (GeoDbCatalog) {
        .next_sequence = 1U,
        .spatial_applied_sequence = 1U,
        .next_server_id = UINT64_C(1) << 63U,
        .next_secondary_index_id = 1U,
    };
    GeoRocksBatch *batch = geo_rocks_batch_create(rocks, 128U, &rocks_status);
    GeoDatabaseStatus database_status = GEO_DATABASE_OK;
    bool succeeded = batch && geo_db_catalog_put(batch, catalog, &database_status) &&
                     geo_rocks_write(rocks, batch, true, &rocks_status);

    geo_rocks_batch_destroy(batch);

    if (!succeeded && status) {
        *status = database_status != GEO_DATABASE_OK ? database_status : geo_db_status_from_rocks(&rocks_status);
    }

    return succeeded;
}

bool geo_db_object_put(GeoRocksBatch *batch,
                       uint64_t object_id,
                       uint64_t sequence,
                       uint64_t morton_code,
                       GeoDocView document,
                       GeoDatabaseStatus *status)
{
    if (!batch || !object_id || !sequence || (document.size && !document.data) || document.size > UINT32_MAX ||
        document.size > SIZE_MAX - GEO_DB_OBJECT_HEADER_SIZE) {
        if (status) {
            *status = GEO_DATABASE_INVALID_ARGUMENT;
        }
        return false;
    }

    size_t value_size = GEO_DB_OBJECT_HEADER_SIZE + document.size;
    unsigned char *value = calloc(1U, value_size);

    if (!value) {
        if (status) {
            *status = GEO_DATABASE_OUT_OF_MEMORY;
        }
        return false;
    }

    unsigned char key[8];
    db_store_u32_le(value, GEO_DB_OBJECT_MAGIC);
    db_store_u16_le(value + 4U, GEO_DB_FORMAT_VERSION);
    db_store_u16_le(value + 6U, GEO_DB_OBJECT_HEADER_SIZE);
    db_store_u32_le(value + 12U, (uint32_t) document.size);
    db_store_u64_le(value + 16U, sequence);
    db_store_u64_le(value + 24U, morton_code);
    if (document.size) {
        memcpy(value + GEO_DB_OBJECT_HEADER_SIZE, document.data, document.size);
    }
    geo_db_store_u64_be(key, object_id);

    GeoRocksStatus rocks_status;
    bool succeeded = geo_rocks_batch_put(batch, GEO_ROCKS_CF_OBJECTS, key, sizeof(key), value, value_size, &rocks_status);

    free(value);

    if (!succeeded && status) {
        *status = geo_db_status_from_rocks(&rocks_status);
    }

    return succeeded;
}

bool geo_db_object_delete(GeoRocksBatch *batch, uint64_t object_id, GeoDatabaseStatus *status)
{
    if (!batch || !object_id) {
        if (status) {
            *status = GEO_DATABASE_INVALID_ARGUMENT;
        }
        return false;
    }

    unsigned char key[8];
    GeoRocksStatus rocks_status;

    geo_db_store_u64_be(key, object_id);
    bool succeeded = geo_rocks_batch_delete(batch, GEO_ROCKS_CF_OBJECTS, key, sizeof(key), &rocks_status);

    if (!succeeded && status) {
        *status = geo_db_status_from_rocks(&rocks_status);
    }

    return succeeded;
}

bool geo_db_object_decode_view(const void *value,
                              size_t value_size,
                              uint64_t *sequence,
                              uint64_t *morton_code,
                              GeoDocView *document)
{
    if (!value || value_size < GEO_DB_OBJECT_HEADER_SIZE || !sequence || !morton_code || !document) {
        return false;
    }

    const unsigned char *bytes = value;
    uint32_t encoded_document_size = db_load_u32_le(bytes + 12U);

    if (db_load_u32_le(bytes) != GEO_DB_OBJECT_MAGIC || db_load_u16_le(bytes + 4U) != GEO_DB_FORMAT_VERSION ||
        db_load_u16_le(bytes + 6U) != GEO_DB_OBJECT_HEADER_SIZE || db_load_u32_le(bytes + 8U) != 0U ||
        encoded_document_size != value_size - GEO_DB_OBJECT_HEADER_SIZE || db_load_u64_le(bytes + 16U) == 0U) {
        return false;
    }

    GeoDocView decoded_document = { 0 };

    if (encoded_document_size &&
        geo_doc_open(bytes + GEO_DB_OBJECT_HEADER_SIZE, encoded_document_size, &decoded_document) != GEO_DOC_OK) {
            return false;
    }

    *sequence = db_load_u64_le(bytes + 16U);
    *morton_code = db_load_u64_le(bytes + 24U);
    *document = decoded_document;
    return true;
}

bool geo_db_object_decode(const void *value,
                          size_t value_size,
                          uint64_t *sequence,
                          uint64_t *morton_code,
                          const void **document,
                          size_t *document_size)
{
    if (!document || !document_size) {
        return false;
    }

    GeoDocView decoded_document;

    if (!geo_db_object_decode_view(value, value_size, sequence, morton_code, &decoded_document)) {
        return false;
    }

    *document = decoded_document.data;
    *document_size = decoded_document.size;
    return true;
}

bool geo_db_spatial_delta_put(GeoRocksBatch *batch,
                              const GeoDbSpatialDelta *delta,
                              GeoDatabaseStatus *status)
{
    if (!batch || !delta || !delta->object_id || !delta->sequence ||
        (delta->operation != GEO_DATABASE_UPSERT && delta->operation != GEO_DATABASE_DELETE)) {
        if (status) {
            *status = GEO_DATABASE_INVALID_ARGUMENT;
        }
        return false;
    }

    unsigned char key[8];
    unsigned char value[GEO_DB_SPATIAL_DELTA_SIZE] = { 0 };

    geo_db_store_u64_be(key, delta->sequence);
    db_store_u32_le(value, GEO_DB_DELTA_MAGIC);
    db_store_u16_le(value + 4U, GEO_DB_FORMAT_VERSION);
    db_store_u16_le(value + 6U, (uint16_t) delta->operation);
    db_store_u64_le(value + 8U, delta->sequence);
    db_store_u64_le(value + 16U, delta->object_id);
    db_store_u64_le(value + 24U, delta->morton_code);

    GeoRocksStatus rocks_status;
    bool succeeded = geo_rocks_batch_put(batch,
                                         GEO_ROCKS_CF_SPATIAL_DELTA,
                                         key,
                                         sizeof(key),
                                         value,
                                         sizeof(value),
                                         &rocks_status);

    if (!succeeded && status) {
        *status = geo_db_status_from_rocks(&rocks_status);
    }

    return succeeded;
}

bool geo_db_spatial_delta_decode(const void *key,
                                 size_t key_size,
                                 const void *value,
                                 size_t value_size,
                                 GeoDbSpatialDelta *delta)
{
    if (!key || key_size != sizeof(uint64_t) || !value || value_size != GEO_DB_SPATIAL_DELTA_SIZE || !delta) {
        return false;
    }

    const unsigned char *bytes = value;
    uint32_t operation = db_load_u16_le(bytes + 6U);
    uint64_t sequence = db_load_u64_le(bytes + 8U);

    if (db_load_u32_le(bytes) != GEO_DB_DELTA_MAGIC || db_load_u16_le(bytes + 4U) != GEO_DB_FORMAT_VERSION ||
        (operation != GEO_DATABASE_UPSERT && operation != GEO_DATABASE_DELETE) || sequence != geo_db_load_u64_be(key)) {
        return false;
    }

    *delta = (GeoDbSpatialDelta) {
        .object_id = db_load_u64_le(bytes + 16U),
        .sequence = sequence,
        .morton_code = db_load_u64_le(bytes + 24U),
        .operation = operation,
    };
    return delta->object_id != 0U && delta->sequence != 0U;
}

bool geo_db_idempotency_get(GeoRocksDatabase *rocks,
                            const void *key,
                            size_t key_size,
                            GeoDbIdempotencyRecord *record,
                            bool *found,
                            GeoDatabaseStatus *status)
{
    if (!rocks || !key || !key_size || !record || !found) {
        if (status) {
            *status = GEO_DATABASE_INVALID_ARGUMENT;
        }

        return false;
    }

    GeoRocksBuffer value = { 0 };
    GeoRocksStatus rocks_status;

    if (!geo_rocks_get(rocks, GEO_ROCKS_CF_IDEMPOTENCY, NULL, key, key_size, &value, &rocks_status)) {
        *found = false;

        if (rocks_status.code == GEO_ROCKS_NOT_FOUND) {
            return true;
        }

        if (status) {
            *status = geo_db_status_from_rocks(&rocks_status);
        }

        return false;
    }

    const unsigned char *bytes = value.data;
    bool valid = value.size == GEO_DB_IDEMPOTENCY_VALUE_SIZE &&
                 db_load_u32_le(bytes) == GEO_DB_IDEMPOTENCY_MAGIC &&
                 db_load_u16_le(bytes + 4U) == GEO_DB_FORMAT_VERSION &&
                 db_load_u16_le(bytes + 6U) == GEO_DB_IDEMPOTENCY_VALUE_SIZE &&
                 db_load_u32_le(bytes + 12U) == 0U;

    if (valid) {
        *record = (GeoDbIdempotencyRecord) {
            .operation_count = db_load_u32_le(bytes + 8U),
            .fingerprint_low = db_load_u64_le(bytes + 16U),
            .fingerprint_high = db_load_u64_le(bytes + 24U),
            .expires_at_seconds = db_load_u64_le(bytes + 32U),
            .first_sequence = db_load_u64_le(bytes + 40U),
        };
        valid = record->operation_count != 0U && record->expires_at_seconds != 0U && record->first_sequence != 0U;
    }

    geo_rocks_buffer_release(&value);
    *found = valid;

    if (!valid && status) {
        *status = GEO_DATABASE_CORRUPTION;
    }

    return valid;
}

bool geo_db_idempotency_put(GeoRocksBatch *batch,
                            const void *key,
                            size_t key_size,
                            const GeoDbIdempotencyRecord *record,
                            GeoDatabaseStatus *status)
{
    if (!batch || !key || !key_size || !record || !record->operation_count || !record->expires_at_seconds ||
        !record->first_sequence) {
        if (status) {
            *status = GEO_DATABASE_INVALID_ARGUMENT;
        }

        return false;
    }

    unsigned char value[GEO_DB_IDEMPOTENCY_VALUE_SIZE] = { 0 };

    db_store_u32_le(value, GEO_DB_IDEMPOTENCY_MAGIC);
    db_store_u16_le(value + 4U, GEO_DB_FORMAT_VERSION);
    db_store_u16_le(value + 6U, GEO_DB_IDEMPOTENCY_VALUE_SIZE);
    db_store_u32_le(value + 8U, record->operation_count);
    db_store_u64_le(value + 16U, record->fingerprint_low);
    db_store_u64_le(value + 24U, record->fingerprint_high);
    db_store_u64_le(value + 32U, record->expires_at_seconds);
    db_store_u64_le(value + 40U, record->first_sequence);

    GeoRocksStatus rocks_status;
    bool succeeded = geo_rocks_batch_put(batch,
                                         GEO_ROCKS_CF_IDEMPOTENCY,
                                         key,
                                         key_size,
                                         value,
                                         sizeof(value),
                                         &rocks_status);

    if (!succeeded && status) {
        *status = geo_db_status_from_rocks(&rocks_status);
    }

    return succeeded;
}
