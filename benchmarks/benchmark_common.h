#ifndef GEOBOLT_BENCHMARK_COMMON_H
#define GEOBOLT_BENCHMARK_COMMON_H

#include "geo_index.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef GEO_BENCHMARK_ALLOCATOR
#define GEO_BENCHMARK_ALLOCATOR "unspecified"
#endif

#if defined(__clang__)
#define GEO_BENCHMARK_COMPILER "clang " __clang_version__
#elif defined(__GNUC__)
#define GEO_BENCHMARK_COMPILER "gcc " __VERSION__
#else
#define GEO_BENCHMARK_COMPILER "unknown"
#endif

static inline bool geo_benchmark_parse_size(const char *text, size_t *value)
{
    char *end = NULL;

    errno = 0;
    unsigned long long parsed = strtoull(text, &end, 10);

    if (errno || !text[0] || !end || *end || parsed > SIZE_MAX) {
        return false;
    }

    *value = (size_t) parsed;

    return true;
}

static inline void geo_benchmark_print_environment(void)
{
    printf("compiler=%s backend=%s allocator=%s\n",
           GEO_BENCHMARK_COMPILER,
           geo_simd_get_name(),
           GEO_BENCHMARK_ALLOCATOR);
}

#endif
