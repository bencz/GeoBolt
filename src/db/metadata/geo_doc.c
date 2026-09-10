#include "geobolt/geodoc.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define GEO_DOC_MAGIC UINT32_C(0x314f4447)
#define GEO_DOC_VERSION 1U
#define GEO_DOC_HEADER_SIZE 32U
#define GEO_DOC_NODE_SIZE 32U
#define GEO_DOC_MAX_NODES UINT32_C(1048576)
#define GEO_DOC_MAX_DEPTH 256U

typedef struct {
    uint32_t magic;
    uint16_t version;
    uint16_t flags;
    uint32_t total_size;
    uint32_t node_count;
    uint32_t root_index;
    uint32_t payload_offset;
    uint64_t checksum;
} GeoDocHeader;

typedef struct {
    uint8_t type;
    uint8_t flags;
    uint16_t reserved;
    uint32_t name_offset;
    uint32_t name_size;
    uint32_t value_offset;
    uint32_t value_size;
    uint32_t child_count;
    uint32_t name_hash;
    uint32_t reserved32;
} GeoDocNode;

_Static_assert(sizeof(GeoDocHeader) == GEO_DOC_HEADER_SIZE, "GeoDoc header layout must remain explicit");
_Static_assert(sizeof(GeoDocNode) == GEO_DOC_NODE_SIZE, "GeoDoc node layout must remain explicit");

typedef struct {
    GeoDocType type;
    char *name;
    size_t name_size;
    unsigned char *value;
    size_t value_size;
    uint32_t *children;
    size_t child_count;
    size_t child_capacity;
} GeoDocBuilderNode;

struct GeoDocBuilder {
    GeoDocBuilderNode *nodes;
    size_t count;
    size_t capacity;
};

static uint16_t doc_load_u16(const unsigned char *data)
{
    return (uint16_t) data[0] | (uint16_t) ((uint16_t) data[1] << 8U);
}

static uint32_t doc_load_u32(const unsigned char *data)
{
    return (uint32_t) data[0] |
           (uint32_t) data[1] << 8U |
           (uint32_t) data[2] << 16U |
           (uint32_t) data[3] << 24U;
}

static uint64_t doc_load_u64(const unsigned char *data)
{
    return (uint64_t) doc_load_u32(data) | (uint64_t) doc_load_u32(data + 4U) << 32U;
}

static void doc_store_u16(unsigned char *data, uint16_t value)
{
    data[0] = (unsigned char) value;
    data[1] = (unsigned char) (value >> 8U);
}

static void doc_store_u32(unsigned char *data, uint32_t value)
{
    data[0] = (unsigned char) value;
    data[1] = (unsigned char) (value >> 8U);
    data[2] = (unsigned char) (value >> 16U);
    data[3] = (unsigned char) (value >> 24U);
}

static void doc_store_u64(unsigned char *data, uint64_t value)
{
    doc_store_u32(data, (uint32_t) value);
    doc_store_u32(data + 4U, (uint32_t) (value >> 32U));
}

static bool doc_size_add(size_t first, size_t second, size_t *result)
{
    if (first > SIZE_MAX - second) {
        return false;
    }

    *result = first + second;
    return true;
}

static bool doc_size_multiply(size_t first, size_t second, size_t *result)
{
    if (first != 0U && second > SIZE_MAX / first) {
        return false;
    }

    *result = first * second;
    return true;
}

static uint32_t doc_hash_name(const void *data, size_t size)
{
    const unsigned char *bytes = data;
    uint32_t hash = UINT32_C(2166136261);

    for (size_t index = 0; index < size; ++index) {
        hash ^= bytes[index];
        hash *= UINT32_C(16777619);
    }

    return hash;
}

static bool doc_utf8_valid(const void *data, size_t size)
{
    const unsigned char *bytes = data;
    size_t index = 0U;

    while (index < size) {
        unsigned char first = bytes[index++];

        if (first <= UINT8_C(0x7f)) {
            continue;
        }

        if (first >= UINT8_C(0xc2) && first <= UINT8_C(0xdf)) {
            if (index >= size || bytes[index] < UINT8_C(0x80) || bytes[index] > UINT8_C(0xbf)) {
                return false;
            }
            index++;
            continue;
        }

        if (first >= UINT8_C(0xe0) && first <= UINT8_C(0xef)) {
            if (index + 1U >= size) {
                return false;
            }

            unsigned char second = bytes[index];
            unsigned char third = bytes[index + 1U];
            bool second_valid = second >= UINT8_C(0x80) && second <= UINT8_C(0xbf);

            if (first == UINT8_C(0xe0)) {
                second_valid = second >= UINT8_C(0xa0) && second <= UINT8_C(0xbf);
            } else if (first == UINT8_C(0xed)) {
                second_valid = second >= UINT8_C(0x80) && second <= UINT8_C(0x9f);
            }

            if (!second_valid || third < UINT8_C(0x80) || third > UINT8_C(0xbf)) {
                return false;
            }
            index += 2U;
            continue;
        }

        if (first >= UINT8_C(0xf0) && first <= UINT8_C(0xf4)) {
            if (index + 2U >= size) {
                return false;
            }

            unsigned char second = bytes[index];
            unsigned char third = bytes[index + 1U];
            unsigned char fourth = bytes[index + 2U];
            bool second_valid = second >= UINT8_C(0x80) && second <= UINT8_C(0xbf);

            if (first == UINT8_C(0xf0)) {
                second_valid = second >= UINT8_C(0x90) && second <= UINT8_C(0xbf);
            } else if (first == UINT8_C(0xf4)) {
                second_valid = second >= UINT8_C(0x80) && second <= UINT8_C(0x8f);
            }

            if (!second_valid || third < UINT8_C(0x80) || third > UINT8_C(0xbf) ||
                fourth < UINT8_C(0x80) || fourth > UINT8_C(0xbf)) {
                return false;
            }
            index += 3U;
            continue;
        }

        return false;
    }

    return true;
}

static uint64_t doc_checksum(const unsigned char *data, size_t size)
{
    uint64_t checksum = UINT64_C(1469598103934665603);

    for (size_t index = 0; index < size; ++index) {
        unsigned char byte = index >= 24U && index < 32U ? 0U : data[index];

        checksum ^= byte;
        checksum *= UINT64_C(1099511628211);
    }

    return checksum;
}

static void builder_node_destroy(GeoDocBuilderNode *node)
{
    free(node->children);
    free(node->value);
    free(node->name);
    memset(node, 0, sizeof(*node));
}

GeoDocBuilder *geo_doc_builder_create(void)
{
    GeoDocBuilder *builder = calloc(1U, sizeof(*builder));

    if (!builder) {
        return NULL;
    }

    builder->nodes = calloc(1U, sizeof(*builder->nodes));

    if (!builder->nodes) {
        free(builder);
        return NULL;
    }

    builder->nodes[0].type = GEO_DOC_OBJECT;
    builder->count = 1U;
    builder->capacity = 1U;
    return builder;
}

void geo_doc_builder_destroy(GeoDocBuilder *builder)
{
    if (!builder) {
        return;
    }

    for (size_t index = 0; index < builder->count; ++index) {
        builder_node_destroy(builder->nodes + index);
    }

    free(builder->nodes);
    free(builder);
}

uint32_t geo_doc_builder_root(const GeoDocBuilder *builder)
{
    return builder ? 0U : UINT32_MAX;
}

static bool builder_reserve_nodes(GeoDocBuilder *builder, size_t required)
{
    if (required <= builder->capacity) {
        return true;
    }

    size_t capacity = builder->capacity ? builder->capacity : 8U;

    while (capacity < required) {
        if (capacity > SIZE_MAX / 2U) {
            return false;
        }

        capacity *= 2U;
    }

    if (capacity > GEO_DOC_MAX_NODES || capacity > SIZE_MAX / sizeof(*builder->nodes)) {
        return false;
    }

    GeoDocBuilderNode *nodes = realloc(builder->nodes, capacity * sizeof(*nodes));

    if (!nodes) {
        return false;
    }

    memset(nodes + builder->capacity, 0, (capacity - builder->capacity) * sizeof(*nodes));
    builder->nodes = nodes;
    builder->capacity = capacity;
    return true;
}

static bool builder_reserve_children(GeoDocBuilderNode *parent, size_t required)
{
    if (required <= parent->child_capacity) {
        return true;
    }

    size_t capacity = parent->child_capacity ? parent->child_capacity * 2U : 4U;

    while (capacity < required) {
        if (capacity > SIZE_MAX / 2U) {
            return false;
        }

        capacity *= 2U;
    }

    if (capacity > UINT32_MAX || capacity > SIZE_MAX / sizeof(*parent->children)) {
        return false;
    }

    uint32_t *children = realloc(parent->children, capacity * sizeof(*children));

    if (!children) {
        return false;
    }

    parent->children = children;
    parent->child_capacity = capacity;
    return true;
}

static bool builder_object_contains(const GeoDocBuilder *builder,
                                    const GeoDocBuilderNode *parent,
                                    const char *name,
                                    size_t name_size)
{
    for (size_t index = 0; index < parent->child_count; ++index) {
        const GeoDocBuilderNode *child = builder->nodes + parent->children[index];

        if (child->name_size == name_size && memcmp(child->name, name, name_size) == 0) {
            return true;
        }
    }

    return false;
}

static GeoDocStatus builder_add(GeoDocBuilder *builder,
                                uint32_t parent_index,
                                const char *name,
                                GeoDocType type,
                                const void *value,
                                size_t value_size,
                                uint32_t *node_index)
{
    if (!builder || parent_index >= builder->count || type > GEO_DOC_OBJECT || (value_size && !value)) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    GeoDocBuilderNode *parent = builder->nodes + parent_index;

    if (parent->type != GEO_DOC_OBJECT && parent->type != GEO_DOC_ARRAY) {
        return GEO_DOC_TYPE_MISMATCH;
    }

    size_t name_size = name ? strlen(name) : 0U;

    if ((parent->type == GEO_DOC_OBJECT && !name) || (parent->type == GEO_DOC_ARRAY && name)) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    if (name_size > UINT32_MAX || value_size > UINT32_MAX || builder->count >= GEO_DOC_MAX_NODES) {
        return GEO_DOC_LIMIT_EXCEEDED;
    }

    if ((name_size && !doc_utf8_valid(name, name_size)) ||
        (type == GEO_DOC_STRING && value_size && !doc_utf8_valid(value, value_size))) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    if (parent->type == GEO_DOC_OBJECT && builder_object_contains(builder, parent, name, name_size)) {
        return GEO_DOC_DUPLICATE_KEY;
    }

    char *name_copy = name_size ? malloc(name_size) : NULL;
    unsigned char *value_copy = value_size ? malloc(value_size) : NULL;

    if ((name_size && !name_copy) || (value_size && !value_copy)) {
        free(value_copy);
        free(name_copy);
        return GEO_DOC_OUT_OF_MEMORY;
    }

    if (name_size) {
        memcpy(name_copy, name, name_size);
    }

    if (value_size) {
        memcpy(value_copy, value, value_size);
    }

    if (!builder_reserve_nodes(builder, builder->count + 1U)) {
        free(value_copy);
        free(name_copy);
        return GEO_DOC_OUT_OF_MEMORY;
    }

    parent = builder->nodes + parent_index;

    if (!builder_reserve_children(parent, parent->child_count + 1U)) {
        free(value_copy);
        free(name_copy);
        return GEO_DOC_OUT_OF_MEMORY;
    }

    size_t new_index = builder->count++;
    builder->nodes[new_index] = (GeoDocBuilderNode) {
        .type = type,
        .name = name_copy,
        .name_size = name_size,
        .value = value_copy,
        .value_size = value_size,
    };
    parent = builder->nodes + parent_index;
    parent->children[parent->child_count++] = (uint32_t) new_index;

    if (node_index) {
        *node_index = (uint32_t) new_index;
    }

    return GEO_DOC_OK;
}

GeoDocStatus geo_doc_builder_add_null(GeoDocBuilder *builder, uint32_t parent, const char *name, uint32_t *node)
{
    return builder_add(builder, parent, name, GEO_DOC_NULL, NULL, 0U, node);
}

GeoDocStatus geo_doc_builder_add_bool(GeoDocBuilder *builder, uint32_t parent, const char *name, bool value, uint32_t *node)
{
    const unsigned char encoded = value ? 1U : 0U;

    return builder_add(builder, parent, name, GEO_DOC_BOOL, &encoded, sizeof(encoded), node);
}

GeoDocStatus geo_doc_builder_add_int64(GeoDocBuilder *builder,
                                       uint32_t parent,
                                       const char *name,
                                       int64_t value,
                                       uint32_t *node)
{
    uint64_t bits;
    unsigned char encoded[8];

    memcpy(&bits, &value, sizeof(bits));
    doc_store_u64(encoded, bits);
    return builder_add(builder, parent, name, GEO_DOC_INT64, encoded, sizeof(encoded), node);
}

GeoDocStatus geo_doc_builder_add_uint64(GeoDocBuilder *builder,
                                        uint32_t parent,
                                        const char *name,
                                        uint64_t value,
                                        uint32_t *node)
{
    unsigned char encoded[8];

    doc_store_u64(encoded, value);
    return builder_add(builder, parent, name, GEO_DOC_UINT64, encoded, sizeof(encoded), node);
}

GeoDocStatus geo_doc_builder_add_double(GeoDocBuilder *builder,
                                        uint32_t parent,
                                        const char *name,
                                        double value,
                                        uint32_t *node)
{
    if (!isfinite(value)) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    uint64_t bits;
    unsigned char encoded[8];

    if (value == 0.0) {
        value = 0.0;
    }

    memcpy(&bits, &value, sizeof(bits));
    doc_store_u64(encoded, bits);
    return builder_add(builder, parent, name, GEO_DOC_DOUBLE, encoded, sizeof(encoded), node);
}

GeoDocStatus geo_doc_builder_add_string(GeoDocBuilder *builder,
                                        uint32_t parent,
                                        const char *name,
                                        const char *value,
                                        size_t value_size,
                                        uint32_t *node)
{
    return builder_add(builder, parent, name, GEO_DOC_STRING, value, value_size, node);
}

GeoDocStatus geo_doc_builder_add_bytes(GeoDocBuilder *builder,
                                       uint32_t parent,
                                       const char *name,
                                       const void *value,
                                       size_t value_size,
                                       uint32_t *node)
{
    return builder_add(builder, parent, name, GEO_DOC_BYTES, value, value_size, node);
}

GeoDocStatus geo_doc_builder_add_array(GeoDocBuilder *builder, uint32_t parent, const char *name, uint32_t *node)
{
    return builder_add(builder, parent, name, GEO_DOC_ARRAY, NULL, 0U, node);
}

GeoDocStatus geo_doc_builder_add_object(GeoDocBuilder *builder, uint32_t parent, const char *name, uint32_t *node)
{
    return builder_add(builder, parent, name, GEO_DOC_OBJECT, NULL, 0U, node);
}

static int builder_compare_children(const GeoDocBuilder *builder, uint32_t first_index, uint32_t second_index)
{
    const GeoDocBuilderNode *first_node = builder->nodes + first_index;
    const GeoDocBuilderNode *second_node = builder->nodes + second_index;
    size_t common = first_node->name_size < second_node->name_size ? first_node->name_size : second_node->name_size;
    int comparison = memcmp(first_node->name, second_node->name, common);

    if (comparison != 0) {
        return comparison;
    }

    return (first_node->name_size > second_node->name_size) - (first_node->name_size < second_node->name_size);
}

static void builder_sift_children(const GeoDocBuilder *builder, uint32_t *children, size_t root, size_t count)
{
    for (;;) {
        size_t selected = root;
        size_t left = root * 2U + 1U;
        size_t right = left + 1U;

        if (left < count && builder_compare_children(builder, children[left], children[selected]) > 0) {
            selected = left;
        }
        if (right < count && builder_compare_children(builder, children[right], children[selected]) > 0) {
            selected = right;
        }
        if (selected == root) {
            return;
        }

        uint32_t temporary = children[root];

        children[root] = children[selected];
        children[selected] = temporary;
        root = selected;
    }
}

static void builder_sort_object_children(const GeoDocBuilder *builder, uint32_t *children, size_t count)
{
    for (size_t root = count / 2U; root > 0U; --root) {
        builder_sift_children(builder, children, root - 1U, count);
    }

    for (size_t remaining = count; remaining > 1U; --remaining) {
        uint32_t temporary = children[0];

        children[0] = children[remaining - 1U];
        children[remaining - 1U] = temporary;
        builder_sift_children(builder, children, 0U, remaining - 1U);
    }
}

static bool builder_measure(const GeoDocBuilder *builder, size_t *document_size)
{
    size_t node_bytes;
    size_t payload_size = 0U;

    if (!doc_size_multiply(builder->count, GEO_DOC_NODE_SIZE, &node_bytes) ||
        !doc_size_add(GEO_DOC_HEADER_SIZE, node_bytes, document_size)) {
        return false;
    }

    for (size_t index = 0; index < builder->count; ++index) {
        const GeoDocBuilderNode *node = builder->nodes + index;
        size_t child_bytes;

        if (!doc_size_multiply(node->child_count, sizeof(uint32_t), &child_bytes) ||
            !doc_size_add(payload_size, node->name_size, &payload_size) ||
            !doc_size_add(payload_size, node->value_size, &payload_size) ||
            !doc_size_add(payload_size, child_bytes, &payload_size)) {
            return false;
        }
    }

    return doc_size_add(*document_size, payload_size, document_size) && *document_size <= UINT32_MAX;
}

GeoDocStatus geo_doc_builder_finish(const GeoDocBuilder *builder, GeoDocBuffer *document)
{
    if (!builder || !document || builder->count == 0U || builder->count > GEO_DOC_MAX_NODES) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    document->data = NULL;
    document->size = 0U;
    size_t total_size;

    if (!builder_measure(builder, &total_size)) {
        return GEO_DOC_LIMIT_EXCEEDED;
    }

    unsigned char *output = calloc(1U, total_size);

    if (!output) {
        return GEO_DOC_OUT_OF_MEMORY;
    }

    size_t payload_offset = GEO_DOC_HEADER_SIZE + builder->count * GEO_DOC_NODE_SIZE;
    size_t cursor = payload_offset;
    doc_store_u32(output, GEO_DOC_MAGIC);
    doc_store_u16(output + 4U, GEO_DOC_VERSION);
    doc_store_u32(output + 8U, (uint32_t) total_size);
    doc_store_u32(output + 12U, (uint32_t) builder->count);
    doc_store_u32(output + 16U, 0U);
    doc_store_u32(output + 20U, (uint32_t) payload_offset);

    for (size_t index = 0; index < builder->count; ++index) {
        const GeoDocBuilderNode *node = builder->nodes + index;
        unsigned char *encoded = output + GEO_DOC_HEADER_SIZE + index * GEO_DOC_NODE_SIZE;
        encoded[0] = (unsigned char) node->type;

        if (node->name_size) {
            doc_store_u32(encoded + 4U, (uint32_t) cursor);
            doc_store_u32(encoded + 8U, (uint32_t) node->name_size);
            memcpy(output + cursor, node->name, node->name_size);
            cursor += node->name_size;
        }

        if (node->type == GEO_DOC_ARRAY || node->type == GEO_DOC_OBJECT) {
            size_t child_bytes = node->child_count * sizeof(uint32_t);
            uint32_t *children = node->child_count ? malloc(child_bytes) : NULL;

            if (node->child_count && !children) {
                free(output);
                return GEO_DOC_OUT_OF_MEMORY;
            }

            if (node->child_count) {
                memcpy(children, node->children, child_bytes);

                if (node->type == GEO_DOC_OBJECT) {
                    builder_sort_object_children(builder, children, node->child_count);
                }
            }

            doc_store_u32(encoded + 12U, (uint32_t) cursor);
            doc_store_u32(encoded + 16U, (uint32_t) child_bytes);
            doc_store_u32(encoded + 20U, (uint32_t) node->child_count);

            for (size_t child = 0; child < node->child_count; ++child) {
                doc_store_u32(output + cursor + child * sizeof(uint32_t), children[child]);
            }

            cursor += child_bytes;
            free(children);
        } else if (node->value_size) {
            doc_store_u32(encoded + 12U, (uint32_t) cursor);
            doc_store_u32(encoded + 16U, (uint32_t) node->value_size);
            memcpy(output + cursor, node->value, node->value_size);
            cursor += node->value_size;
        }

        doc_store_u32(encoded + 24U, doc_hash_name(node->name, node->name_size));
    }

    if (cursor != total_size) {
        free(output);
        return GEO_DOC_INVALID_FORMAT;
    }

    doc_store_u64(output + 24U, doc_checksum(output, total_size));
    document->data = output;
    document->size = total_size;
    return GEO_DOC_OK;
}

void geo_doc_buffer_release(GeoDocBuffer *buffer)
{
    if (!buffer) {
        return;
    }

    free(buffer->data);
    buffer->data = NULL;
    buffer->size = 0U;
}

static const unsigned char *view_node(const GeoDocView *view, uint32_t node_index)
{
    uint32_t count = doc_load_u32(view->data + 12U);

    if (node_index >= count) {
        return NULL;
    }

    return view->data + GEO_DOC_HEADER_SIZE + (size_t) node_index * GEO_DOC_NODE_SIZE;
}

static bool view_range_valid(const GeoDocView *view, uint32_t offset, uint32_t size)
{
    uint32_t payload_offset = doc_load_u32(view->data + 20U);

    return offset >= payload_offset && offset <= view->size && size <= view->size - offset;
}

static bool validate_node(const GeoDocView *view, uint32_t node_index)
{
    const unsigned char *node = view_node(view, node_index);
    GeoDocType type = node ? (GeoDocType) node[0] : GEO_DOC_OBJECT + 1;
    uint32_t name_offset = node ? doc_load_u32(node + 4U) : 0U;
    uint32_t name_size = node ? doc_load_u32(node + 8U) : 0U;
    uint32_t value_offset = node ? doc_load_u32(node + 12U) : 0U;
    uint32_t value_size = node ? doc_load_u32(node + 16U) : 0U;
    uint32_t child_count = node ? doc_load_u32(node + 20U) : 0U;

    if (!node || type > GEO_DOC_OBJECT || node[1] != 0U || doc_load_u16(node + 2U) != 0U || doc_load_u32(node + 28U) != 0U) {
        return false;
    }

    if ((name_size && !view_range_valid(view, name_offset, name_size)) ||
        (!name_size && name_offset != 0U) ||
        (name_size && !doc_utf8_valid(view->data + name_offset, name_size)) ||
        doc_load_u32(node + 24U) != doc_hash_name(name_size ? view->data + name_offset : NULL, name_size)) {
        return false;
    }

    if (type == GEO_DOC_ARRAY || type == GEO_DOC_OBJECT) {
        return value_size == child_count * sizeof(uint32_t) &&
               ((!value_size && value_offset == 0U) || view_range_valid(view, value_offset, value_size));
    }

    if (child_count != 0U || (value_size && !view_range_valid(view, value_offset, value_size)) ||
        (!value_size && value_offset != 0U)) {
        return false;
    }

    if (type == GEO_DOC_NULL) {
        return value_size == 0U;
    }
    if (type == GEO_DOC_BOOL) {
        return value_size == 1U && view->data[value_offset] <= 1U;
    }
    if (type == GEO_DOC_INT64 || type == GEO_DOC_UINT64 || type == GEO_DOC_DOUBLE) {
        if (value_size != sizeof(uint64_t)) {
            return false;
        }

        if (type == GEO_DOC_DOUBLE) {
            uint64_t bits = doc_load_u64(view->data + value_offset);
            double value;

            memcpy(&value, &bits, sizeof(value));
            return isfinite(value) && (value != 0.0 || bits == 0U);
        }
    }

    if (type == GEO_DOC_STRING) {
        return doc_utf8_valid(value_size ? view->data + value_offset : NULL, value_size);
    }

    return true;
}

static bool validate_tree(const GeoDocView *view, uint32_t node_index, uint8_t *states, unsigned depth)
{
    if (depth > GEO_DOC_MAX_DEPTH || states[node_index] != 0U) {
        return false;
    }

    states[node_index] = 1U;
    const unsigned char *node = view_node(view, node_index);
    GeoDocType type = (GeoDocType) node[0];

    if (type == GEO_DOC_ARRAY || type == GEO_DOC_OBJECT) {
        uint32_t value_offset = doc_load_u32(node + 12U);
        uint32_t child_count = doc_load_u32(node + 20U);
        const unsigned char *previous_name = NULL;
        uint32_t previous_name_size = 0U;

        for (uint32_t child_position = 0; child_position < child_count; ++child_position) {
            uint32_t child_index = doc_load_u32(view->data + value_offset + child_position * sizeof(uint32_t));
            const unsigned char *child = view_node(view, child_index);

            if (!child || (type == GEO_DOC_ARRAY && doc_load_u32(child + 8U) != 0U) ||
                !validate_tree(view, child_index, states, depth + 1U)) {
                return false;
            }

            if (type == GEO_DOC_OBJECT) {
                uint32_t name_offset = doc_load_u32(child + 4U);
                uint32_t name_size = doc_load_u32(child + 8U);
                const unsigned char *name = view->data + name_offset;

                if (previous_name) {
                    size_t common = previous_name_size < name_size ? previous_name_size : name_size;
                    int order = memcmp(previous_name, name, common);

                    if (order > 0 || (order == 0 && previous_name_size >= name_size)) {
                        return false;
                    }
                }

                previous_name = name;
                previous_name_size = name_size;
            }
        }
    }

    states[node_index] = 2U;
    return true;
}

GeoDocStatus geo_doc_open(const void *data, size_t size, GeoDocView *view)
{
    if (!data || !view || size < GEO_DOC_HEADER_SIZE || size > UINT32_MAX) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    const unsigned char *bytes = data;
    uint32_t node_count = doc_load_u32(bytes + 12U);
    size_t node_bytes;
    size_t expected_payload;

    if (doc_load_u32(bytes) != GEO_DOC_MAGIC || doc_load_u16(bytes + 4U) != GEO_DOC_VERSION ||
        doc_load_u16(bytes + 6U) != 0U || doc_load_u32(bytes + 8U) != size || node_count == 0U ||
        node_count > GEO_DOC_MAX_NODES || doc_load_u32(bytes + 16U) != 0U ||
        !doc_size_multiply(node_count, GEO_DOC_NODE_SIZE, &node_bytes) ||
        !doc_size_add(GEO_DOC_HEADER_SIZE, node_bytes, &expected_payload) ||
        doc_load_u32(bytes + 20U) != expected_payload || expected_payload > size ||
        doc_load_u64(bytes + 24U) != doc_checksum(bytes, size)) {
        return GEO_DOC_INVALID_FORMAT;
    }

    GeoDocView candidate = {
        .data = bytes,
        .size = size,
    };

    for (uint32_t index = 0; index < node_count; ++index) {
        if (!validate_node(&candidate, index)) {
            return GEO_DOC_INVALID_FORMAT;
        }
    }

    uint8_t *states = calloc(node_count, sizeof(*states));

    if (!states) {
        return GEO_DOC_OUT_OF_MEMORY;
    }

    bool valid = validate_tree(&candidate, 0U, states, 0U);

    for (uint32_t index = 0; valid && index < node_count; ++index) {
        valid = states[index] == 2U;
    }

    free(states);

    if (!valid || view_node(&candidate, 0U)[0] != GEO_DOC_OBJECT || doc_load_u32(view_node(&candidate, 0U) + 8U) != 0U) {
        return GEO_DOC_INVALID_FORMAT;
    }

    *view = candidate;
    return GEO_DOC_OK;
}

GeoDocValue geo_doc_root(GeoDocView view)
{
    return (GeoDocValue) {
        .document = view.data,
        .node_index = 0U,
    };
}

static GeoDocView value_view(GeoDocValue value)
{
    return (GeoDocView) {
        .data = value.document,
        .size = value.document ? doc_load_u32(value.document + 8U) : 0U,
    };
}

GeoDocType geo_doc_value_type(GeoDocValue value)
{
    GeoDocView view = value_view(value);
    const unsigned char *node = value.document ? view_node(&view, value.node_index) : NULL;

    return node ? (GeoDocType) node[0] : GEO_DOC_NULL;
}

size_t geo_doc_value_count(GeoDocValue value)
{
    GeoDocView view = value_view(value);
    const unsigned char *node = value.document ? view_node(&view, value.node_index) : NULL;

    return node ? doc_load_u32(node + 20U) : 0U;
}

GeoDocStatus geo_doc_value_at(GeoDocValue value, size_t index, GeoDocValue *child)
{
    if (!value.document || !child) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    GeoDocView view = value_view(value);
    const unsigned char *node = view_node(&view, value.node_index);
    GeoDocType type = node ? (GeoDocType) node[0] : GEO_DOC_NULL;
    uint32_t child_count = node ? doc_load_u32(node + 20U) : 0U;

    if (type != GEO_DOC_ARRAY && type != GEO_DOC_OBJECT) {
        return GEO_DOC_TYPE_MISMATCH;
    }
    if (index >= child_count) {
        return GEO_DOC_NOT_FOUND;
    }

    uint32_t child_index = doc_load_u32(view.data + doc_load_u32(node + 12U) + index * sizeof(uint32_t));
    *child = (GeoDocValue) {
        .document = value.document,
        .node_index = child_index,
    };
    return GEO_DOC_OK;
}

static GeoDocStatus value_scalar_data(GeoDocValue value,
                                      GeoDocType expected,
                                      const unsigned char **data,
                                      uint32_t *size)
{
    if (!value.document || !data || !size) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    GeoDocView view = value_view(value);
    const unsigned char *node = view_node(&view, value.node_index);

    if (!node || node[0] != expected) {
        return GEO_DOC_TYPE_MISMATCH;
    }

    *data = view.data + doc_load_u32(node + 12U);
    *size = doc_load_u32(node + 16U);
    return GEO_DOC_OK;
}

GeoDocStatus geo_doc_value_bool(GeoDocValue value, bool *result)
{
    const unsigned char *data;
    uint32_t size;
    GeoDocStatus status = value_scalar_data(value, GEO_DOC_BOOL, &data, &size);

    if (status == GEO_DOC_OK && result) {
        *result = data[0] != 0U;
    } else if (!result) {
        status = GEO_DOC_INVALID_ARGUMENT;
    }

    return status;
}

GeoDocStatus geo_doc_value_int64(GeoDocValue value, int64_t *result)
{
    const unsigned char *data;
    uint32_t size;
    GeoDocStatus status = value_scalar_data(value, GEO_DOC_INT64, &data, &size);

    if (status == GEO_DOC_OK && result) {
        uint64_t bits = doc_load_u64(data);

        memcpy(result, &bits, sizeof(*result));
    } else if (!result) {
        status = GEO_DOC_INVALID_ARGUMENT;
    }

    return status;
}

GeoDocStatus geo_doc_value_uint64(GeoDocValue value, uint64_t *result)
{
    const unsigned char *data;
    uint32_t size;
    GeoDocStatus status = value_scalar_data(value, GEO_DOC_UINT64, &data, &size);

    if (status == GEO_DOC_OK && result) {
        *result = doc_load_u64(data);
    } else if (!result) {
        status = GEO_DOC_INVALID_ARGUMENT;
    }

    return status;
}

GeoDocStatus geo_doc_value_double(GeoDocValue value, double *result)
{
    const unsigned char *data;
    uint32_t size;
    GeoDocStatus status = value_scalar_data(value, GEO_DOC_DOUBLE, &data, &size);

    if (status == GEO_DOC_OK && result) {
        uint64_t bits = doc_load_u64(data);

        memcpy(result, &bits, sizeof(*result));
    } else if (!result) {
        status = GEO_DOC_INVALID_ARGUMENT;
    }

    return status;
}

GeoDocStatus geo_doc_value_data(GeoDocValue value, const void **data, size_t *size)
{
    if (!value.document || !data || !size) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    GeoDocType type = geo_doc_value_type(value);

    if (type != GEO_DOC_STRING && type != GEO_DOC_BYTES) {
        return GEO_DOC_TYPE_MISMATCH;
    }

    GeoDocView view = value_view(value);
    const unsigned char *node = view_node(&view, value.node_index);
    *data = view.data + doc_load_u32(node + 12U);
    *size = doc_load_u32(node + 16U);
    return GEO_DOC_OK;
}

GeoDocStatus geo_doc_value_name(GeoDocValue value, const char **name, size_t *size)
{
    if (!value.document || !name || !size) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    GeoDocView view = value_view(value);
    const unsigned char *node = view_node(&view, value.node_index);

    if (!node) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    uint32_t name_size = doc_load_u32(node + 8U);
    *name = name_size ? (const char *) view.data + doc_load_u32(node + 4U) : NULL;
    *size = name_size;
    return GEO_DOC_OK;
}

#define GEO_DOC_POINTER_STACK_SEGMENT 256U

static bool pointer_decode_segment(const char **cursor,
                                   char stack_storage[GEO_DOC_POINTER_STACK_SEGMENT],
                                   char **decoded,
                                   size_t *decoded_size,
                                   char **heap_storage)
{
    const char *begin = *cursor;
    const char *end = begin;

    while (*end && *end != '/') {
        end++;
    }

    size_t maximum = (size_t) (end - begin);
    char *result = stack_storage;

    *heap_storage = NULL;

    if (maximum >= GEO_DOC_POINTER_STACK_SEGMENT) {
        *heap_storage = malloc(maximum + 1U);
        result = *heap_storage;
    }

    if (!result) {
        return false;
    }

    size_t written = 0U;

    for (const char *source = begin; source < end; ++source) {
        if (*source == '~') {
            source++;

            if (source >= end || (*source != '0' && *source != '1')) {
                free(*heap_storage);
                return false;
            }

            result[written++] = *source == '0' ? '~' : '/';
        } else {
            result[written++] = *source;
        }
    }

    result[written] = '\0';
    *cursor = end;
    *decoded = result;
    *decoded_size = written;
    return true;
}

static GeoDocStatus object_find_child(GeoDocValue object, const char *name, size_t name_size, GeoDocValue *child)
{
    uint32_t expected_hash = doc_hash_name(name, name_size);
    size_t count = geo_doc_value_count(object);

    for (size_t index = 0; index < count; ++index) {
        GeoDocValue candidate;
        const char *candidate_name;
        size_t candidate_size;

        if (geo_doc_value_at(object, index, &candidate) != GEO_DOC_OK ||
            geo_doc_value_name(candidate, &candidate_name, &candidate_size) != GEO_DOC_OK) {
            return GEO_DOC_INVALID_FORMAT;
        }

        GeoDocView view = value_view(candidate);
        const unsigned char *node = view_node(&view, candidate.node_index);

        bool metadata_matches = doc_load_u32(node + 24U) == expected_hash && candidate_size == name_size;

        if (metadata_matches && (name_size == 0U || memcmp(candidate_name, name, name_size) == 0)) {
            *child = candidate;
            return GEO_DOC_OK;
        }
    }

    return GEO_DOC_NOT_FOUND;
}

static bool parse_array_index(const char *text, size_t size, size_t *index)
{
    if (size == 0U || (size > 1U && text[0] == '0')) {
        return false;
    }

    size_t value = 0U;

    for (size_t position = 0; position < size; ++position) {
        unsigned digit = (unsigned) ((unsigned char) text[position] - (unsigned char) '0');

        if (digit > 9U || value > (SIZE_MAX - digit) / 10U) {
            return false;
        }

        value = value * 10U + digit;
    }

    *index = value;
    return true;
}

GeoDocStatus geo_doc_find_pointer(GeoDocView view, const char *json_pointer, GeoDocValue *value)
{
    if (!view.data || !json_pointer || !value) {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    GeoDocValue current = geo_doc_root(view);

    if (json_pointer[0] == '\0') {
        *value = current;
        return GEO_DOC_OK;
    }

    if (json_pointer[0] != '/') {
        return GEO_DOC_INVALID_ARGUMENT;
    }

    const char *cursor = json_pointer + 1U;

    for (;;) {
        char stack_segment[GEO_DOC_POINTER_STACK_SEGMENT];
        char *segment;
        char *heap_segment;
        size_t segment_size;

        if (!pointer_decode_segment(&cursor,
                                    stack_segment,
                                    &segment,
                                    &segment_size,
                                    &heap_segment)) {
            return GEO_DOC_INVALID_ARGUMENT;
        }

        GeoDocType type = geo_doc_value_type(current);
        GeoDocStatus status;

        if (type == GEO_DOC_OBJECT) {
            status = object_find_child(current, segment, segment_size, &current);
        } else if (type == GEO_DOC_ARRAY) {
            size_t index;

            status = parse_array_index(segment, segment_size, &index)
                         ? geo_doc_value_at(current, index, &current)
                         : GEO_DOC_INVALID_ARGUMENT;
        } else {
            status = GEO_DOC_TYPE_MISMATCH;
        }

        free(heap_segment);

        if (status != GEO_DOC_OK) {
            return status;
        }

        if (*cursor == '\0') {
            *value = current;
            return GEO_DOC_OK;
        }

        cursor++;
    }
}
