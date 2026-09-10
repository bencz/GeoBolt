#ifndef GEOBOLT_GEODOC_H
#define GEOBOLT_GEODOC_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GeoDocBuilder GeoDocBuilder;

typedef enum {
    GEO_DOC_NULL = 0,
    GEO_DOC_BOOL,
    GEO_DOC_INT64,
    GEO_DOC_UINT64,
    GEO_DOC_DOUBLE,
    GEO_DOC_STRING,
    GEO_DOC_BYTES,
    GEO_DOC_ARRAY,
    GEO_DOC_OBJECT,
} GeoDocType;

typedef enum {
    GEO_DOC_OK = 0,
    GEO_DOC_INVALID_ARGUMENT,
    GEO_DOC_OUT_OF_MEMORY,
    GEO_DOC_INVALID_FORMAT,
    GEO_DOC_DUPLICATE_KEY,
    GEO_DOC_NOT_FOUND,
    GEO_DOC_TYPE_MISMATCH,
    GEO_DOC_LIMIT_EXCEEDED,
} GeoDocStatus;

typedef struct {
    void *data;
    size_t size;
} GeoDocBuffer;

typedef struct {
    const unsigned char *data;
    size_t size;
} GeoDocView;

typedef struct {
    const unsigned char *document;
    uint32_t node_index;
} GeoDocValue;

GeoDocBuilder *geo_doc_builder_create(void);
void geo_doc_builder_destroy(GeoDocBuilder *builder);
uint32_t geo_doc_builder_root(const GeoDocBuilder *builder);

GeoDocStatus geo_doc_builder_add_null(GeoDocBuilder *builder, uint32_t parent, const char *name, uint32_t *node);
GeoDocStatus geo_doc_builder_add_bool(GeoDocBuilder *builder, uint32_t parent, const char *name, bool value, uint32_t *node);
GeoDocStatus geo_doc_builder_add_int64(GeoDocBuilder *builder,
                                       uint32_t parent,
                                       const char *name,
                                       int64_t value,
                                       uint32_t *node);
GeoDocStatus geo_doc_builder_add_uint64(GeoDocBuilder *builder,
                                        uint32_t parent,
                                        const char *name,
                                        uint64_t value,
                                        uint32_t *node);
GeoDocStatus geo_doc_builder_add_double(GeoDocBuilder *builder,
                                        uint32_t parent,
                                        const char *name,
                                        double value,
                                        uint32_t *node);
GeoDocStatus geo_doc_builder_add_string(GeoDocBuilder *builder,
                                        uint32_t parent,
                                        const char *name,
                                        const char *value,
                                        size_t value_size,
                                        uint32_t *node);
GeoDocStatus geo_doc_builder_add_bytes(GeoDocBuilder *builder,
                                       uint32_t parent,
                                       const char *name,
                                       const void *value,
                                       size_t value_size,
                                       uint32_t *node);
GeoDocStatus geo_doc_builder_add_array(GeoDocBuilder *builder, uint32_t parent, const char *name, uint32_t *node);
GeoDocStatus geo_doc_builder_add_object(GeoDocBuilder *builder, uint32_t parent, const char *name, uint32_t *node);

GeoDocStatus geo_doc_builder_finish(const GeoDocBuilder *builder, GeoDocBuffer *document);
void geo_doc_buffer_release(GeoDocBuffer *buffer);

GeoDocStatus geo_doc_open(const void *data, size_t size, GeoDocView *view);
GeoDocValue geo_doc_root(GeoDocView view);
GeoDocStatus geo_doc_find_pointer(GeoDocView view, const char *json_pointer, GeoDocValue *value);
GeoDocType geo_doc_value_type(GeoDocValue value);
size_t geo_doc_value_count(GeoDocValue value);
GeoDocStatus geo_doc_value_at(GeoDocValue value, size_t index, GeoDocValue *child);
GeoDocStatus geo_doc_value_bool(GeoDocValue value, bool *result);
GeoDocStatus geo_doc_value_int64(GeoDocValue value, int64_t *result);
GeoDocStatus geo_doc_value_uint64(GeoDocValue value, uint64_t *result);
GeoDocStatus geo_doc_value_double(GeoDocValue value, double *result);
GeoDocStatus geo_doc_value_data(GeoDocValue value, const void **data, size_t *size);
GeoDocStatus geo_doc_value_name(GeoDocValue value, const char **name, size_t *size);

#ifdef __cplusplus
}
#endif

#endif
