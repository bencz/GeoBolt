#include "../benchmarks/benchmark_common.h"

#include <dirent.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <sched.h>
#include <signal.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#define SOAK_LATENCY_BUCKETS 32U
#define SOAK_METRICS_ALIGNMENT 64U
#define SOAK_PATH_CAPACITY 512U
#define SOAK_DEFAULT_DURATION_SECONDS 300U
#define SOAK_DEFAULT_VEHICLES 100000U
#define SOAK_DEFAULT_PUBLISHERS 4U
#define SOAK_DEFAULT_QUERY_THREADS 8U
#define SOAK_DEFAULT_BATCH_SIZE 2048U
#define SOAK_DEFAULT_MIN_DELAY_MS 250U
#define SOAK_DEFAULT_MAX_DELAY_MS 3000U
#define SOAK_DEFAULT_MAX_SEGMENTS 32U
#define SOAK_DEFAULT_REPORT_SECONDS 10U
#define SOAK_DEFAULT_AUDIT_SECONDS 5U
#define SOAK_MICROBATCH_WINDOW_NS UINT64_C(50000000)

typedef struct {
    double latitude;
    double longitude;
} SoakHub;

static const SoakHub SOAK_HUBS[] = {
    { -23.5505, -46.6333 },
    { 40.7128, -74.0060 },
    { 51.5074, -0.1278 },
    { 48.8566, 2.3522 },
    { 19.4326, -99.1332 },
    { 35.6762, 139.6503 },
    { 1.3521, 103.8198 },
    { 28.6139, 77.2090 },
    { -1.2921, 36.8219 },
    { -33.8688, 151.2093 },
    { 37.7749, -122.4194 },
    { 25.2048, 55.2708 },
};

typedef struct {
    uint64_t id;
    uint64_t next_due_ns;
    double latitude;
    double longitude;
    uint8_t hub;
    bool active;
} SoakVehicle;

typedef struct {
    size_t vehicle_index;
    double latitude;
    double longitude;
    bool insertion;
} SoakLiveChange;

typedef struct {
    uint64_t operations;
    uint64_t updates;
    uint64_t inserts;
    uint64_t deletes;
    uint64_t queries;
    uint64_t audits;
    uint64_t failures;
    uint64_t results;
    uint64_t scanned;
    uint64_t segments;
    uint64_t latency_total_ns;
    uint64_t latency_max_ns;
    uint64_t latency_histogram[SOAK_LATENCY_BUCKETS];
    uint64_t thread_cpu_ns;
} SoakCounters;

typedef struct {
    _Atomic uint64_t operations;
    _Atomic uint64_t updates;
    _Atomic uint64_t inserts;
    _Atomic uint64_t deletes;
    _Atomic uint64_t queries;
    _Atomic uint64_t audits;
    _Atomic uint64_t failures;
    _Atomic uint64_t results;
    _Atomic uint64_t scanned;
    _Atomic uint64_t segments;
    _Atomic uint64_t latency_total_ns;
    _Atomic uint64_t latency_max_ns;
    _Atomic uint64_t latency_histogram[SOAK_LATENCY_BUCKETS];
    _Atomic uint64_t thread_cpu_ns;
    uint8_t padding[24];
} SoakMetrics;

typedef struct {
    size_t duration_seconds;
    size_t vehicle_count;
    size_t publisher_count;
    size_t query_thread_count;
    size_t batch_size;
    size_t minimum_delay_ms;
    size_t maximum_delay_ms;
    size_t maximum_segments;
    size_t report_seconds;
    size_t audit_seconds;
} SoakConfig;

typedef struct {
    GeoSegmentSet *set;
    SoakVehicle *vehicles;
    const SoakConfig *config;
    pthread_rwlock_t *audit_gate;
    _Atomic bool *start;
    _Atomic bool *stop;
    _Atomic bool *failed;
    _Atomic uint64_t *segment_sequence;
    SoakMetrics *metrics;
    const char *directory;
    size_t worker_index;
    size_t first_vehicle;
    size_t vehicle_count;
    uint64_t random_state;
} SoakPublisher;

typedef struct {
    GeoSegmentSet *set;
    const SoakConfig *config;
    _Atomic bool *start;
    _Atomic bool *stop;
    _Atomic bool *failed;
    SoakMetrics *metrics;
    size_t worker_index;
    uint64_t random_state;
} SoakQueryWorker;

typedef struct {
    GeoSegmentSet *set;
    SoakVehicle *vehicles;
    const SoakConfig *config;
    pthread_rwlock_t *audit_gate;
    _Atomic bool *start;
    _Atomic bool *stop;
    _Atomic bool *failed;
    _Atomic uint64_t *audited_active;
    SoakMetrics *metrics;
    uint64_t random_state;
} SoakAuditor;

typedef struct {
    uint64_t resident_bytes;
    uint64_t peak_resident_bytes;
    uint64_t bytes_read;
    uint64_t bytes_written;
    uint64_t minor_faults;
    uint64_t major_faults;
    uint64_t voluntary_switches;
    uint64_t involuntary_switches;
} SoakResources;

typedef struct {
    size_t count;
    double sum_x;
    double sum_y;
    double sum_xx;
    double sum_xy;
    uint64_t minimum_rss;
    uint64_t maximum_rss;
} SoakMemoryTrend;

_Static_assert(sizeof(SoakMetrics) % SOAK_METRICS_ALIGNMENT == 0,
               "worker metrics must occupy complete cache lines");

static volatile sig_atomic_t soak_interrupted = 0;

static uint64_t soak_monotonic_ns(void)
{
    struct timespec now;

    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0) {
        return 0;
    }

    return (uint64_t) now.tv_sec * UINT64_C(1000000000) + (uint64_t) now.tv_nsec;
}

static uint64_t soak_thread_cpu_ns(void)
{
    struct timespec now;

    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &now) != 0) {
        return 0;
    }

    return (uint64_t) now.tv_sec * UINT64_C(1000000000) + (uint64_t) now.tv_nsec;
}

static uint64_t soak_process_cpu_ns(void)
{
    struct timespec now;

    if (clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now) != 0) {
        return 0;
    }

    return (uint64_t) now.tv_sec * UINT64_C(1000000000) + (uint64_t) now.tv_nsec;
}

static void soak_sleep_ns(uint64_t duration_ns)
{
    struct timespec requested = {
        .tv_sec = (time_t) (duration_ns / UINT64_C(1000000000)),
        .tv_nsec = (long) (duration_ns % UINT64_C(1000000000)),
    };

    while (nanosleep(&requested, &requested) != 0 && errno == EINTR && !soak_interrupted) {
    }
}

static uint64_t soak_random_u64(uint64_t *state)
{
    uint64_t value = *state;

    value ^= value >> 12U;
    value ^= value << 25U;
    value ^= value >> 27U;
    *state = value;

    return value * UINT64_C(2685821657736338717);
}

static double soak_random_unit(uint64_t *state)
{
    return (double) (soak_random_u64(state) >> 11U) * (1.0 / 9007199254740992.0);
}

static uint64_t soak_random_range(uint64_t *state, uint64_t minimum, uint64_t maximum)
{
    uint64_t span = maximum - minimum + 1U;

    return minimum + soak_random_u64(state) % span;
}

static double soak_wrap_longitude(double longitude)
{
    longitude = fmod(longitude + 180.0, 360.0);

    if (longitude < 0.0) {
        longitude += 360.0;
    }

    return longitude - 180.0;
}

static void soak_quantize_location(double latitude, double longitude, double *quantized_latitude, double *quantized_longitude)
{
    GeoPoint decoded = geo_decode(geo_encode(latitude, longitude));

    *quantized_latitude = decoded.lat;
    *quantized_longitude = decoded.lng;
}

static void soak_initial_location(uint64_t *random_state, uint8_t *hub, double *latitude, double *longitude)
{
    size_t hub_count = sizeof(SOAK_HUBS) / sizeof(SOAK_HUBS[0]);
    bool global_vehicle = soak_random_u64(random_state) % 5U == 0;

    if (global_vehicle) {
        *hub = UINT8_MAX;
        *latitude = soak_random_unit(random_state) * 170.0 - 85.0;
        *longitude = soak_random_unit(random_state) * 360.0 - 180.0;
    } else {
        *hub = (uint8_t) (soak_random_u64(random_state) % hub_count);
        *latitude = SOAK_HUBS[*hub].latitude + (soak_random_unit(random_state) - 0.5) * 0.6;
        *longitude = SOAK_HUBS[*hub].longitude + (soak_random_unit(random_state) - 0.5) * 0.6;
    }

    soak_quantize_location(*latitude, *longitude, latitude, longitude);
}

static void soak_move_vehicle(const SoakVehicle *vehicle,
                              uint64_t *random_state,
                              double *latitude,
                              double *longitude)
{
    if (soak_random_u64(random_state) % 1000U == 0) {
        uint8_t new_hub;

        soak_initial_location(random_state, &new_hub, latitude, longitude);
        return;
    }

    double distance_meters = 2.0 + soak_random_unit(random_state) * 248.0;
    double angle = soak_random_unit(random_state) * 6.28318530717958647692;
    double latitude_delta = cos(angle) * distance_meters / 111320.0;
    double longitude_scale = cos(geo_to_radians(vehicle->latitude));

    if (fabs(longitude_scale) < 0.01) {
        longitude_scale = longitude_scale < 0.0 ? -0.01 : 0.01;
    }

    double longitude_delta = sin(angle) * distance_meters / (111320.0 * longitude_scale);

    *latitude = vehicle->latitude + latitude_delta;
    *longitude = soak_wrap_longitude(vehicle->longitude + longitude_delta);

    if (*latitude > 89.9 || *latitude < -89.9) {
        *latitude = *latitude > 0.0 ? 89.9 : -89.9;
    }

    soak_quantize_location(*latitude, *longitude, latitude, longitude);
}

static void soak_signal_handler(int signal_number)
{
    (void) signal_number;
    soak_interrupted = 1;
}

static unsigned soak_latency_bucket(uint64_t duration_ns)
{
    uint64_t microseconds = duration_ns / 1000U;

    if (!microseconds) {
        return 0;
    }

    unsigned bucket = 63U - (unsigned) __builtin_clzll(microseconds);

    return bucket < SOAK_LATENCY_BUCKETS ? bucket : SOAK_LATENCY_BUCKETS - 1U;
}

static void soak_counters_record_latency(SoakCounters *counters, uint64_t duration_ns)
{
    counters->operations++;
    counters->latency_total_ns += duration_ns;

    if (duration_ns > counters->latency_max_ns) {
        counters->latency_max_ns = duration_ns;
    }

    counters->latency_histogram[soak_latency_bucket(duration_ns)]++;
}

static void soak_metrics_publish(SoakMetrics *metrics, const SoakCounters *counters)
{
    atomic_store_explicit(&metrics->operations, counters->operations, memory_order_relaxed);
    atomic_store_explicit(&metrics->updates, counters->updates, memory_order_relaxed);
    atomic_store_explicit(&metrics->inserts, counters->inserts, memory_order_relaxed);
    atomic_store_explicit(&metrics->deletes, counters->deletes, memory_order_relaxed);
    atomic_store_explicit(&metrics->queries, counters->queries, memory_order_relaxed);
    atomic_store_explicit(&metrics->audits, counters->audits, memory_order_relaxed);
    atomic_store_explicit(&metrics->failures, counters->failures, memory_order_relaxed);
    atomic_store_explicit(&metrics->results, counters->results, memory_order_relaxed);
    atomic_store_explicit(&metrics->scanned, counters->scanned, memory_order_relaxed);
    atomic_store_explicit(&metrics->segments, counters->segments, memory_order_relaxed);
    atomic_store_explicit(&metrics->latency_total_ns, counters->latency_total_ns, memory_order_relaxed);
    atomic_store_explicit(&metrics->latency_max_ns, counters->latency_max_ns, memory_order_relaxed);
    atomic_store_explicit(&metrics->thread_cpu_ns, counters->thread_cpu_ns, memory_order_relaxed);

    for (unsigned i = 0; i < SOAK_LATENCY_BUCKETS; ++i) {
        atomic_store_explicit(&metrics->latency_histogram[i], counters->latency_histogram[i], memory_order_relaxed);
    }
}

static void soak_metrics_accumulate(SoakCounters *total, const SoakMetrics *metrics)
{
    total->operations += atomic_load_explicit(&metrics->operations, memory_order_relaxed);
    total->updates += atomic_load_explicit(&metrics->updates, memory_order_relaxed);
    total->inserts += atomic_load_explicit(&metrics->inserts, memory_order_relaxed);
    total->deletes += atomic_load_explicit(&metrics->deletes, memory_order_relaxed);
    total->queries += atomic_load_explicit(&metrics->queries, memory_order_relaxed);
    total->audits += atomic_load_explicit(&metrics->audits, memory_order_relaxed);
    total->failures += atomic_load_explicit(&metrics->failures, memory_order_relaxed);
    total->results += atomic_load_explicit(&metrics->results, memory_order_relaxed);
    total->scanned += atomic_load_explicit(&metrics->scanned, memory_order_relaxed);
    total->segments += atomic_load_explicit(&metrics->segments, memory_order_relaxed);
    total->latency_total_ns += atomic_load_explicit(&metrics->latency_total_ns, memory_order_relaxed);
    total->thread_cpu_ns += atomic_load_explicit(&metrics->thread_cpu_ns, memory_order_relaxed);

    uint64_t maximum = atomic_load_explicit(&metrics->latency_max_ns, memory_order_relaxed);

    if (maximum > total->latency_max_ns) {
        total->latency_max_ns = maximum;
    }

    for (unsigned i = 0; i < SOAK_LATENCY_BUCKETS; ++i) {
        total->latency_histogram[i] += atomic_load_explicit(&metrics->latency_histogram[i], memory_order_relaxed);
    }
}

static uint64_t soak_histogram_percentile(const SoakCounters *counters, uint64_t numerator, uint64_t denominator)
{
    if (!counters->operations) {
        return 0;
    }

    uint64_t target = (counters->operations * numerator + denominator - 1U) / denominator;
    uint64_t cumulative = 0;

    for (unsigned i = 0; i < SOAK_LATENCY_BUCKETS; ++i) {
        cumulative += counters->latency_histogram[i];

        if (cumulative >= target) {
            return UINT64_C(1) << i;
        }
    }

    return UINT64_C(1) << (SOAK_LATENCY_BUCKETS - 1U);
}

static bool soak_read_resources(SoakResources *resources)
{
    memset(resources, 0, sizeof(*resources));

    struct rusage usage;

    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return false;
    }

#if defined(__APPLE__)
    resources->peak_resident_bytes = (uint64_t) usage.ru_maxrss;
#else
    resources->peak_resident_bytes = (uint64_t) usage.ru_maxrss * 1024U;
#endif
    resources->minor_faults = (uint64_t) usage.ru_minflt;
    resources->major_faults = (uint64_t) usage.ru_majflt;
    resources->voluntary_switches = (uint64_t) usage.ru_nvcsw;
    resources->involuntary_switches = (uint64_t) usage.ru_nivcsw;

#if defined(__linux__)
    FILE *memory_file = fopen("/proc/self/statm", "r");
    unsigned long total_pages = 0;
    unsigned long resident_pages = 0;

    if (memory_file) {
        if (fscanf(memory_file, "%lu %lu", &total_pages, &resident_pages) == 2) {
            long page_size = sysconf(_SC_PAGESIZE);

            if (page_size > 0 && resident_pages <= UINT64_MAX / (uint64_t) page_size) {
                resources->resident_bytes = (uint64_t) resident_pages * (uint64_t) page_size;
            }
        }

        fclose(memory_file);
    }

    FILE *io_file = fopen("/proc/self/io", "r");

    if (io_file) {
        char label[64];
        uint64_t value;

        while (fscanf(io_file, "%63s %" SCNu64, label, &value) == 2) {
            if (strcmp(label, "read_bytes:") == 0) {
                resources->bytes_read = value;
            } else if (strcmp(label, "write_bytes:") == 0) {
                resources->bytes_written = value;
            }
        }

        fclose(io_file);
    }
#else
    resources->resident_bytes = resources->peak_resident_bytes;
#endif

    return true;
}

static void soak_memory_trend_add(SoakMemoryTrend *trend, double elapsed_minutes, uint64_t resident_bytes)
{
    double resident_mib = (double) resident_bytes / (1024.0 * 1024.0);

    trend->count++;
    trend->sum_x += elapsed_minutes;
    trend->sum_y += resident_mib;
    trend->sum_xx += elapsed_minutes * elapsed_minutes;
    trend->sum_xy += elapsed_minutes * resident_mib;

    if (!trend->minimum_rss || resident_bytes < trend->minimum_rss) {
        trend->minimum_rss = resident_bytes;
    }

    if (resident_bytes > trend->maximum_rss) {
        trend->maximum_rss = resident_bytes;
    }
}

static double soak_memory_trend_slope(const SoakMemoryTrend *trend)
{
    if (trend->count < 2) {
        return 0.0;
    }

    double count = (double) trend->count;
    double denominator = count * trend->sum_xx - trend->sum_x * trend->sum_x;

    if (fabs(denominator) < 1e-12) {
        return 0.0;
    }

    return (count * trend->sum_xy - trend->sum_x * trend->sum_y) / denominator;
}

static bool soak_parse_config(int argc, char **argv, SoakConfig *config)
{
    *config = (SoakConfig) {
        .duration_seconds = SOAK_DEFAULT_DURATION_SECONDS,
        .vehicle_count = SOAK_DEFAULT_VEHICLES,
        .publisher_count = SOAK_DEFAULT_PUBLISHERS,
        .query_thread_count = SOAK_DEFAULT_QUERY_THREADS,
        .batch_size = SOAK_DEFAULT_BATCH_SIZE,
        .minimum_delay_ms = SOAK_DEFAULT_MIN_DELAY_MS,
        .maximum_delay_ms = SOAK_DEFAULT_MAX_DELAY_MS,
        .maximum_segments = SOAK_DEFAULT_MAX_SEGMENTS,
        .report_seconds = SOAK_DEFAULT_REPORT_SECONDS,
        .audit_seconds = SOAK_DEFAULT_AUDIT_SECONDS,
    };

    size_t *values[] = {
        &config->duration_seconds,
        &config->vehicle_count,
        &config->publisher_count,
        &config->query_thread_count,
        &config->batch_size,
        &config->minimum_delay_ms,
        &config->maximum_delay_ms,
        &config->maximum_segments,
        &config->report_seconds,
        &config->audit_seconds,
    };

    if (argc > (int) (sizeof(values) / sizeof(values[0])) + 1) {
        return false;
    }

    for (int i = 1; i < argc; ++i) {
        if (!geo_benchmark_parse_size(argv[i], values[i - 1])) {
            return false;
        }
    }

    return config->duration_seconds &&
           config->vehicle_count &&
           config->publisher_count &&
           config->query_thread_count &&
           config->batch_size &&
           config->minimum_delay_ms &&
           config->maximum_delay_ms >= config->minimum_delay_ms &&
           config->maximum_segments &&
           config->report_seconds &&
           config->audit_seconds &&
           config->publisher_count <= config->vehicle_count &&
           config->batch_size <= config->vehicle_count &&
           config->duration_seconds <= UINT64_MAX / UINT64_C(1000000000) &&
           config->maximum_delay_ms <= UINT64_MAX / UINT64_C(1000000);
}

static bool soak_actor_due_less(const SoakVehicle *vehicles, size_t first, size_t second)
{
    if (vehicles[first].next_due_ns != vehicles[second].next_due_ns) {
        return vehicles[first].next_due_ns < vehicles[second].next_due_ns;
    }

    return vehicles[first].id < vehicles[second].id;
}

static void soak_actor_heap_sift_down(size_t *heap,
                                      size_t count,
                                      size_t position,
                                      const SoakVehicle *vehicles)
{
    size_t actor = heap[position];

    while (position * 2U + 1U < count) {
        size_t left = position * 2U + 1U;
        size_t right = left + 1U;
        size_t earliest = right < count && soak_actor_due_less(vehicles, heap[right], heap[left]) ? right : left;

        if (!soak_actor_due_less(vehicles, heap[earliest], actor)) {
            break;
        }

        heap[position] = heap[earliest];
        position = earliest;
    }

    heap[position] = actor;
}

static void soak_actor_heap_build(size_t *heap,
                                  size_t count,
                                  size_t first_vehicle,
                                  SoakVehicle *vehicles,
                                  uint64_t now_ns,
                                  uint64_t *random_state,
                                  const SoakConfig *config)
{
    uint64_t minimum_delay_ns = config->minimum_delay_ms * UINT64_C(1000000);
    uint64_t maximum_delay_ns = config->maximum_delay_ms * UINT64_C(1000000);

    for (size_t i = 0; i < count; ++i) {
        size_t vehicle_index = first_vehicle + i;

        vehicles[vehicle_index].next_due_ns = now_ns + soak_random_range(random_state, minimum_delay_ns, maximum_delay_ns);
        heap[i] = vehicle_index;
    }

    for (size_t position = count / 2U; position > 0; --position) {
        soak_actor_heap_sift_down(heap, count, position - 1U, vehicles);
    }
}

static void soak_actor_heap_reschedule_root(size_t *heap,
                                            size_t count,
                                            SoakVehicle *vehicles,
                                            uint64_t next_due_ns)
{
    vehicles[heap[0]].next_due_ns = next_due_ns;
    soak_actor_heap_sift_down(heap, count, 0, vehicles);
}

static bool soak_build_segment_file(const char *path, const GeoRecord *records, size_t count)
{
    GeoIndex *index = geo_index_create(count);
    bool succeeded = index &&
                     geo_index_add_records(index, records, count) &&
                     geo_index_build(index) &&
                     geo_index_save(index, path);

    geo_index_destroy(index);

    return succeeded;
}

static bool soak_publish_live_batch(SoakPublisher *publisher,
                                    const GeoRecord *records,
                                    const SoakLiveChange *changes,
                                    size_t count,
                                    SoakCounters *counters)
{
    uint64_t sequence = atomic_fetch_add_explicit(publisher->segment_sequence, 1, memory_order_relaxed);
    char path[SOAK_PATH_CAPACITY];
    int length = snprintf(path,
                          sizeof(path),
                          "%s/live-%zu-%020" PRIu64 ".geobolt",
                          publisher->directory,
                          publisher->worker_index,
                          sequence);

    if (length < 0 || (size_t) length >= sizeof(path)) {
        return false;
    }

    uint64_t start_ns = soak_monotonic_ns();
    bool built = soak_build_segment_file(path, records, count);
    bool succeeded = built && geo_segment_set_upsert_file(publisher->set, path);
    uint64_t duration_ns = soak_monotonic_ns() - start_ns;

    soak_counters_record_latency(counters, duration_ns);

    if (!succeeded) {
        fprintf(stderr,
                "publisher=%zu live_batch_failed stage=%s count=%zu path=%s errno=%d\n",
                publisher->worker_index,
                built ? "add_file" : "build_file",
                count,
                path,
                errno);
        unlink(path);
        return false;
    }

    for (size_t i = 0; i < count; ++i) {
        SoakVehicle *vehicle = publisher->vehicles + changes[i].vehicle_index;

        vehicle->latitude = changes[i].latitude;
        vehicle->longitude = changes[i].longitude;
        vehicle->active = true;

        if (changes[i].insertion) {
            counters->inserts++;
        } else {
            counters->updates++;
        }
    }

    counters->segments++;

    return true;
}

static bool soak_publish_delete_batch(SoakPublisher *publisher,
                                      const uint64_t *ids,
                                      const size_t *vehicle_indices,
                                      size_t count,
                                      SoakCounters *counters)
{
    uint64_t start_ns = soak_monotonic_ns();
    bool succeeded = geo_segment_set_remove_ids(publisher->set, ids, count);
    uint64_t duration_ns = soak_monotonic_ns() - start_ns;

    soak_counters_record_latency(counters, duration_ns);

    if (!succeeded) {
        fprintf(stderr,
                "publisher=%zu delete_batch_failed count=%zu errno=%d\n",
                publisher->worker_index,
                count,
                errno);
        return false;
    }

    for (size_t i = 0; i < count; ++i) {
        publisher->vehicles[vehicle_indices[i]].active = false;
    }

    counters->deletes += count;

    return true;
}

static void soak_publisher_fail(SoakPublisher *publisher, SoakCounters *counters)
{
    counters->failures++;
    soak_metrics_publish(publisher->metrics, counters);
    atomic_store_explicit(publisher->failed, true, memory_order_release);
    atomic_store_explicit(publisher->stop, true, memory_order_release);
}

static void *soak_publisher_main(void *argument)
{
    SoakPublisher *publisher = argument;
    size_t *actor_heap = malloc(publisher->vehicle_count * sizeof(*actor_heap));
    GeoRecord *live_records = malloc(publisher->config->batch_size * sizeof(*live_records));
    SoakLiveChange *live_changes = malloc(publisher->config->batch_size * sizeof(*live_changes));
    uint64_t *delete_ids = malloc(publisher->config->batch_size * sizeof(*delete_ids));
    size_t *delete_indices = malloc(publisher->config->batch_size * sizeof(*delete_indices));
    SoakCounters counters = { 0 };

    if (!actor_heap || !live_records || !live_changes || !delete_ids || !delete_indices) {
        soak_publisher_fail(publisher, &counters);
        free(delete_indices);
        free(delete_ids);
        free(live_changes);
        free(live_records);
        free(actor_heap);

        return NULL;
    }

    while (!atomic_load_explicit(publisher->start, memory_order_acquire) &&
           !atomic_load_explicit(publisher->stop, memory_order_acquire)) {
        sched_yield();
    }

    uint64_t thread_cpu_start = soak_thread_cpu_ns();
    uint64_t now_ns = soak_monotonic_ns();

    soak_actor_heap_build(actor_heap,
                          publisher->vehicle_count,
                          publisher->first_vehicle,
                          publisher->vehicles,
                          now_ns,
                          &publisher->random_state,
                          publisher->config);

    while (!atomic_load_explicit(publisher->stop, memory_order_acquire)) {
        now_ns = soak_monotonic_ns();

        if (publisher->vehicles[actor_heap[0]].next_due_ns > now_ns) {
            uint64_t wait_ns = publisher->vehicles[actor_heap[0]].next_due_ns - now_ns;

            soak_sleep_ns(wait_ns < UINT64_C(1000000) ? wait_ns : UINT64_C(1000000));
            continue;
        }

        size_t live_count = 0;
        size_t delete_count = 0;
        size_t event_count = 0;
        uint64_t minimum_delay_ns = publisher->config->minimum_delay_ms * UINT64_C(1000000);
        uint64_t maximum_delay_ns = publisher->config->maximum_delay_ms * UINT64_C(1000000);
        uint64_t collection_deadline_ns = now_ns + SOAK_MICROBATCH_WINDOW_NS;

        while (event_count < publisher->config->batch_size) {
            uint64_t next_due_ns = publisher->vehicles[actor_heap[0]].next_due_ns;

            if (next_due_ns > collection_deadline_ns) {
                break;
            }

            now_ns = soak_monotonic_ns();

            if (next_due_ns > now_ns) {
                soak_sleep_ns(next_due_ns - now_ns);
            }

            size_t vehicle_index = actor_heap[0];
            SoakVehicle *vehicle = publisher->vehicles + vehicle_index;
            uint64_t operation = soak_random_u64(&publisher->random_state) % 100U;
            uint64_t next_delay_ns = soak_random_range(&publisher->random_state, minimum_delay_ns, maximum_delay_ns);

            if (vehicle->active && operation < 3U) {
                delete_ids[delete_count] = vehicle->id;
                delete_indices[delete_count] = vehicle_index;
                delete_count++;
            } else {
                double latitude;
                double longitude;

                if (vehicle->active) {
                    soak_move_vehicle(vehicle, &publisher->random_state, &latitude, &longitude);
                } else {
                    soak_initial_location(&publisher->random_state, &vehicle->hub, &latitude, &longitude);
                }

                live_records[live_count] = (GeoRecord) {
                    .id = vehicle->id,
                    .z = geo_encode(latitude, longitude),
                };
                live_changes[live_count] = (SoakLiveChange) {
                    .vehicle_index = vehicle_index,
                    .latitude = latitude,
                    .longitude = longitude,
                    .insertion = !vehicle->active,
                };
                live_count++;
            }

            soak_actor_heap_reschedule_root(actor_heap,
                                            publisher->vehicle_count,
                                            publisher->vehicles,
                                            collection_deadline_ns + next_delay_ns);
            event_count++;
        }

        bool gate_locked = pthread_rwlock_rdlock(publisher->audit_gate) == 0;
        bool succeeded = gate_locked;

        if (succeeded && live_count) {
            succeeded = soak_publish_live_batch(publisher, live_records, live_changes, live_count, &counters);
        }

        if (succeeded && delete_count) {
            succeeded = soak_publish_delete_batch(publisher, delete_ids, delete_indices, delete_count, &counters);
        }

        if (gate_locked) {
            pthread_rwlock_unlock(publisher->audit_gate);
        }

        if (!succeeded) {
            soak_publisher_fail(publisher, &counters);
            break;
        }

        if ((counters.operations & 15U) == 0) {
            counters.thread_cpu_ns = soak_thread_cpu_ns() - thread_cpu_start;
            soak_metrics_publish(publisher->metrics, &counters);
        }
    }

    counters.thread_cpu_ns = soak_thread_cpu_ns() - thread_cpu_start;
    soak_metrics_publish(publisher->metrics, &counters);
    free(delete_indices);
    free(delete_ids);
    free(live_changes);
    free(live_records);
    free(actor_heap);

    return NULL;
}

static void soak_random_query_center(uint64_t *random_state, double *latitude, double *longitude)
{
    size_t hub_count = sizeof(SOAK_HUBS) / sizeof(SOAK_HUBS[0]);

    if (soak_random_u64(random_state) % 5U == 0) {
        *latitude = soak_random_unit(random_state) * 170.0 - 85.0;
        *longitude = soak_random_unit(random_state) * 360.0 - 180.0;
        return;
    }

    const SoakHub *hub = SOAK_HUBS + soak_random_u64(random_state) % hub_count;

    *latitude = hub->latitude + (soak_random_unit(random_state) - 0.5) * 0.4;
    *longitude = soak_wrap_longitude(hub->longitude + (soak_random_unit(random_state) - 0.5) * 0.4);
}

static double soak_random_query_radius(uint64_t *random_state)
{
    static const double radii[] = { 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0 };

    return radii[soak_random_u64(random_state) % (sizeof(radii) / sizeof(radii[0]))];
}

static bool soak_run_query(SoakQueryWorker *worker,
                           GeoIdResult *id_result,
                           uint64_t *result_count,
                           uint64_t *scanned_count)
{
    double latitude;
    double longitude;
    GeoSearchStats stats = { 0 };
    uint64_t operation = soak_random_u64(&worker->random_state) % 100U;
    const char *operation_name;
    bool succeeded;

    soak_random_query_center(&worker->random_state, &latitude, &longitude);

    if (operation < 70U) {
        operation_name = "radius_count";
        size_t count = 0;

        succeeded = geo_segment_set_search_radius_count(worker->set,
                                                        latitude,
                                                        longitude,
                                                        soak_random_query_radius(&worker->random_state),
                                                        &count,
                                                        &stats);
        *result_count = count;
    } else if (operation < 90U) {
        operation_name = "radius_ids";
        succeeded = geo_segment_set_search_radius_ids_reuse(worker->set,
                                                            latitude,
                                                            longitude,
                                                            soak_random_query_radius(&worker->random_state),
                                                            id_result,
                                                            &stats);
        *result_count = id_result->count;
    } else if (operation < 98U) {
        operation_name = "bbox_count";
        double latitude_span = 0.01 + soak_random_unit(&worker->random_state) * 0.5;
        double longitude_span = 0.01 + soak_random_unit(&worker->random_state) * 0.5;
        size_t count = 0;

        succeeded = geo_segment_set_search_bbox_count(worker->set,
                                                      fmax(-90.0, latitude - latitude_span),
                                                      fmin(90.0, latitude + latitude_span),
                                                      soak_wrap_longitude(longitude - longitude_span),
                                                      soak_wrap_longitude(longitude + longitude_span),
                                                      &count,
                                                      &stats);
        *result_count = count;
    } else {
        operation_name = "knn";
        size_t k = 1U + (size_t) (soak_random_u64(&worker->random_state) % 32U);
        GeoSearchResult *result = geo_segment_set_search_knn(worker->set, latitude, longitude, k, 100.0, &stats);

        succeeded = result != NULL;
        *result_count = result ? result->count : 0;
        geo_result_destroy(result);
    }

    *scanned_count = stats.records_scanned;

    if (!succeeded || *result_count > worker->config->vehicle_count) {
        fprintf(stderr,
                "query_worker=%zu operation=%s succeeded=%s results=%" PRIu64 " vehicle_limit=%zu\n",
                worker->worker_index,
                operation_name,
                succeeded ? "yes" : "no",
                *result_count,
                worker->config->vehicle_count);

        return false;
    }

    return true;
}

static void *soak_query_main(void *argument)
{
    SoakQueryWorker *worker = argument;
    GeoIdResult *id_result = geo_id_result_create(1024);
    SoakCounters counters = { 0 };

    if (!id_result) {
        counters.failures++;
        soak_metrics_publish(worker->metrics, &counters);
        atomic_store_explicit(worker->failed, true, memory_order_release);
        atomic_store_explicit(worker->stop, true, memory_order_release);

        return NULL;
    }

    while (!atomic_load_explicit(worker->start, memory_order_acquire) &&
           !atomic_load_explicit(worker->stop, memory_order_acquire)) {
        sched_yield();
    }

    uint64_t thread_cpu_start = soak_thread_cpu_ns();

    while (!atomic_load_explicit(worker->stop, memory_order_acquire)) {
        uint64_t result_count = 0;
        uint64_t scanned_count = 0;
        uint64_t start_ns = soak_monotonic_ns();
        bool succeeded = soak_run_query(worker, id_result, &result_count, &scanned_count);
        uint64_t duration_ns = soak_monotonic_ns() - start_ns;

        soak_counters_record_latency(&counters, duration_ns);
        counters.queries++;
        counters.results += result_count;
        counters.scanned += scanned_count;

        if (!succeeded) {
            counters.failures++;
            atomic_store_explicit(worker->failed, true, memory_order_release);
            atomic_store_explicit(worker->stop, true, memory_order_release);
            break;
        }

        if ((counters.operations & 4095U) == 0) {
            counters.thread_cpu_ns = soak_thread_cpu_ns() - thread_cpu_start;
            soak_metrics_publish(worker->metrics, &counters);
        }
    }

    counters.thread_cpu_ns = soak_thread_cpu_ns() - thread_cpu_start;
    soak_metrics_publish(worker->metrics, &counters);
    geo_id_result_destroy(id_result);

    return NULL;
}

static int soak_compare_u64(const void *first, const void *second)
{
    uint64_t first_value = *(const uint64_t *) first;
    uint64_t second_value = *(const uint64_t *) second;

    if (first_value < second_value) {
        return -1;
    }

    return first_value > second_value;
}

static bool soak_audit_once(SoakAuditor *auditor,
                            GeoIdResult *actual,
                            uint64_t *expected,
                            SoakCounters *counters)
{
    uint64_t start_ns = soak_monotonic_ns();

    if (pthread_rwlock_wrlock(auditor->audit_gate) != 0) {
        return false;
    }

    size_t active_count = 0;

    for (size_t i = 0; i < auditor->config->vehicle_count; ++i) {
        active_count += auditor->vehicles[i].active;
    }

    size_t global_count = 0;
    bool succeeded = geo_segment_set_search_radius_count(auditor->set,
                                                         0.0,
                                                         0.0,
                                                         25000.0,
                                                         &global_count,
                                                         NULL) &&
                     global_count == active_count;

    if (!succeeded) {
        fprintf(stderr, "audit_global_count_failed expected=%zu actual=%zu\n", active_count, global_count);
    }

    double latitude;
    double longitude;
    double radius_km = soak_random_query_radius(&auditor->random_state);

    soak_random_query_center(&auditor->random_state, &latitude, &longitude);

    if (succeeded) {
        succeeded = geo_segment_set_search_radius_ids_reuse(auditor->set,
                                                            latitude,
                                                            longitude,
                                                            radius_km,
                                                            actual,
                                                            NULL);
    }

    size_t expected_count = 0;

    for (size_t i = 0; succeeded && i < auditor->config->vehicle_count; ++i) {
        const SoakVehicle *vehicle = auditor->vehicles + i;

        if (vehicle->active &&
            geo_haversine_km(latitude, longitude, vehicle->latitude, vehicle->longitude) <= radius_km) {
            expected[expected_count++] = vehicle->id;
        }
    }

    if (succeeded) {
        qsort(actual->ids, actual->count, sizeof(*actual->ids), soak_compare_u64);
        qsort(expected, expected_count, sizeof(*expected), soak_compare_u64);
        succeeded = actual->count == expected_count &&
                    memcmp(actual->ids, expected, expected_count * sizeof(*expected)) == 0;

        if (!succeeded) {
            fprintf(stderr,
                    "audit_radius_ids_failed lat=%.8f lng=%.8f radius=%.3f expected=%zu actual=%zu\n",
                    latitude,
                    longitude,
                    radius_km,
                    expected_count,
                    actual->count);
        }
    }

    atomic_store_explicit(auditor->audited_active, active_count, memory_order_relaxed);
    pthread_rwlock_unlock(auditor->audit_gate);

    uint64_t duration_ns = soak_monotonic_ns() - start_ns;

    soak_counters_record_latency(counters, duration_ns);
    counters->audits++;
    counters->results += actual->count;

    return succeeded;
}

static void *soak_auditor_main(void *argument)
{
    SoakAuditor *auditor = argument;
    GeoIdResult *actual = geo_id_result_create(auditor->config->vehicle_count);
    uint64_t *expected = malloc(auditor->config->vehicle_count * sizeof(*expected));
    SoakCounters counters = { 0 };

    if (!actual || !expected) {
        counters.failures++;
        soak_metrics_publish(auditor->metrics, &counters);
        atomic_store_explicit(auditor->failed, true, memory_order_release);
        atomic_store_explicit(auditor->stop, true, memory_order_release);
        free(expected);
        geo_id_result_destroy(actual);

        return NULL;
    }

    while (!atomic_load_explicit(auditor->start, memory_order_acquire) &&
           !atomic_load_explicit(auditor->stop, memory_order_acquire)) {
        sched_yield();
    }

    uint64_t thread_cpu_start = soak_thread_cpu_ns();
    uint64_t interval_ns = auditor->config->audit_seconds * UINT64_C(1000000000);
    uint64_t next_audit_ns = soak_monotonic_ns();

    while (!atomic_load_explicit(auditor->stop, memory_order_acquire)) {
        uint64_t now_ns = soak_monotonic_ns();

        if (now_ns < next_audit_ns) {
            uint64_t wait_ns = next_audit_ns - now_ns;

            soak_sleep_ns(wait_ns < UINT64_C(100000000) ? wait_ns : UINT64_C(100000000));
            continue;
        }

        if (!soak_audit_once(auditor, actual, expected, &counters)) {
            counters.failures++;
            soak_metrics_publish(auditor->metrics, &counters);
            atomic_store_explicit(auditor->failed, true, memory_order_release);
            atomic_store_explicit(auditor->stop, true, memory_order_release);
            break;
        }

        counters.thread_cpu_ns = soak_thread_cpu_ns() - thread_cpu_start;
        soak_metrics_publish(auditor->metrics, &counters);
        next_audit_ns = soak_monotonic_ns() + interval_ns;
    }

    counters.thread_cpu_ns = soak_thread_cpu_ns() - thread_cpu_start;
    soak_metrics_publish(auditor->metrics, &counters);
    free(expected);
    geo_id_result_destroy(actual);

    return NULL;
}

static bool soak_initialize_vehicles(SoakVehicle *vehicles, size_t count, uint64_t *checksum)
{
    uint64_t random_state = UINT64_C(0x243f6a8885a308d3);

    *checksum = 0;

    for (size_t i = 0; i < count; ++i) {
        double latitude;
        double longitude;
        uint8_t hub;

        soak_initial_location(&random_state, &hub, &latitude, &longitude);

        vehicles[i] = (SoakVehicle) {
            .id = UINT64_C(1000000000) + i,
            .latitude = latitude,
            .longitude = longitude,
            .hub = hub,
            .active = true,
        };
        *checksum ^= vehicles[i].id + geo_encode(latitude, longitude) * UINT64_C(0x9e3779b97f4a7c15);
    }

    return true;
}

static bool soak_create_bootstrap(const char *path,
                                  const SoakVehicle *vehicles,
                                  size_t vehicle_count,
                                  size_t sort_threads)
{
    GeoRecord *records = malloc(vehicle_count * sizeof(*records));

    if (!records) {
        return false;
    }

    for (size_t i = 0; i < vehicle_count; ++i) {
        records[i] = (GeoRecord) {
            .id = vehicles[i].id,
            .z = geo_encode(vehicles[i].latitude, vehicles[i].longitude),
        };
    }

    GeoIndex *index = geo_index_create(vehicle_count);
    bool succeeded = index &&
                     geo_index_add_records(index, records, vehicle_count) &&
                     geo_index_build_parallel(index, sort_threads) &&
                     geo_index_save(index, path);

    geo_index_destroy(index);
    free(records);

    return succeeded;
}

static void soak_collect_metrics(const SoakMetrics *publisher_metrics,
                                 size_t publisher_count,
                                 const SoakMetrics *query_metrics,
                                 size_t query_count,
                                 const SoakMetrics *audit_metrics,
                                 SoakCounters *publishers,
                                 SoakCounters *queries,
                                 SoakCounters *audits)
{
    memset(publishers, 0, sizeof(*publishers));
    memset(queries, 0, sizeof(*queries));
    memset(audits, 0, sizeof(*audits));

    for (size_t i = 0; i < publisher_count; ++i) {
        soak_metrics_accumulate(publishers, publisher_metrics + i);
    }

    for (size_t i = 0; i < query_count; ++i) {
        soak_metrics_accumulate(queries, query_metrics + i);
    }

    soak_metrics_accumulate(audits, audit_metrics);
}

static void soak_print_progress(double elapsed_seconds,
                                double interval_seconds,
                                const SoakCounters *publishers,
                                const SoakCounters *queries,
                                uint64_t previous_events,
                                uint64_t previous_queries,
                                uint64_t active_vehicles,
                                const SoakResources *resources,
                                const GeoSegmentSet *set)
{
    uint64_t events = publishers->updates + publishers->inserts + publishers->deletes;
    double event_rate = interval_seconds > 0.0 ? (double) (events - previous_events) / interval_seconds : 0.0;
    double query_rate = interval_seconds > 0.0 ? (double) (queries->queries - previous_queries) / interval_seconds : 0.0;
    double query_p99_us = (double) soak_histogram_percentile(queries, 99, 100);
    double publish_p99_us = (double) soak_histogram_percentile(publishers, 99, 100);

    printf("progress elapsed=%.1fs events=%" PRIu64 " event_rate=%.0f/s queries=%" PRIu64
           " query_rate=%.0f/s active=%" PRIu64 " segments=%zu physical_records=%" PRIu64
           " compacting=%s rss=%.1fMiB peak_rss=%.1fMiB query_p99<=%.0fus publish_p99<=%.0fus failures=%" PRIu64 "\n",
           elapsed_seconds,
           events,
           event_rate,
           queries->queries,
           query_rate,
           active_vehicles,
           geo_segment_set_count(set),
           geo_segment_set_record_count(set),
           geo_segment_set_compaction_active(set) ? "yes" : "no",
           (double) resources->resident_bytes / (1024.0 * 1024.0),
           (double) resources->peak_resident_bytes / (1024.0 * 1024.0),
           query_p99_us,
           publish_p99_us,
           publishers->failures + queries->failures);
    fflush(stdout);
}

static bool soak_cleanup_directory(const char *directory)
{
    DIR *stream = opendir(directory);

    if (!stream) {
        return false;
    }

    bool succeeded = true;
    struct dirent *entry;

    while ((entry = readdir(stream)) != NULL) {
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }

        char path[SOAK_PATH_CAPACITY];
        int length = snprintf(path, sizeof(path), "%s/%s", directory, entry->d_name);

        if (length < 0 || (size_t) length >= sizeof(path) || unlink(path) != 0) {
            succeeded = false;
        }
    }

    if (closedir(stream) != 0) {
        succeeded = false;
    }

    if (succeeded && rmdir(directory) != 0) {
        succeeded = false;
    }

    return succeeded;
}

static bool soak_start_threads(const SoakConfig *config,
                               GeoSegmentSet *set,
                               SoakVehicle *vehicles,
                               const char *directory,
                               pthread_rwlock_t *audit_gate,
                               _Atomic bool *start,
                               _Atomic bool *stop,
                               _Atomic bool *failed,
                               _Atomic uint64_t *segment_sequence,
                               _Atomic uint64_t *audited_active,
                               pthread_t *publisher_threads,
                               SoakPublisher *publishers,
                               SoakMetrics *publisher_metrics,
                               pthread_t *query_threads,
                               SoakQueryWorker *query_workers,
                               SoakMetrics *query_metrics,
                               pthread_t *auditor_thread,
                               SoakAuditor *auditor,
                               SoakMetrics *audit_metrics,
                               size_t *publishers_started,
                               size_t *queries_started,
                               bool *auditor_started)
{
    size_t next_vehicle = 0;

    for (size_t i = 0; i < config->publisher_count; ++i) {
        size_t workers_left = config->publisher_count - i;
        size_t vehicles_left = config->vehicle_count - next_vehicle;
        size_t owned_vehicles = (vehicles_left + workers_left - 1U) / workers_left;

        publishers[i] = (SoakPublisher) {
            .set = set,
            .vehicles = vehicles,
            .config = config,
            .audit_gate = audit_gate,
            .start = start,
            .stop = stop,
            .failed = failed,
            .segment_sequence = segment_sequence,
            .metrics = publisher_metrics + i,
            .directory = directory,
            .worker_index = i,
            .first_vehicle = next_vehicle,
            .vehicle_count = owned_vehicles,
            .random_state = UINT64_C(0x9e3779b97f4a7c15) ^ (i + 1U) * UINT64_C(0xbf58476d1ce4e5b9),
        };
        next_vehicle += owned_vehicles;

        if (pthread_create(publisher_threads + i, NULL, soak_publisher_main, publishers + i) != 0) {
            return false;
        }

        (*publishers_started)++;
    }

    for (size_t i = 0; i < config->query_thread_count; ++i) {
        query_workers[i] = (SoakQueryWorker) {
            .set = set,
            .config = config,
            .start = start,
            .stop = stop,
            .failed = failed,
            .metrics = query_metrics + i,
            .worker_index = i,
            .random_state = UINT64_C(0x94d049bb133111eb) ^ (i + 1U) * UINT64_C(0x369dea0f31a53f85),
        };

        if (pthread_create(query_threads + i, NULL, soak_query_main, query_workers + i) != 0) {
            return false;
        }

        (*queries_started)++;
    }

    *auditor = (SoakAuditor) {
        .set = set,
        .vehicles = vehicles,
        .config = config,
        .audit_gate = audit_gate,
        .start = start,
        .stop = stop,
        .failed = failed,
        .audited_active = audited_active,
        .metrics = audit_metrics,
        .random_state = UINT64_C(0x3c6ef372fe94f82b),
    };

    if (pthread_create(auditor_thread, NULL, soak_auditor_main, auditor) != 0) {
        return false;
    }

    *auditor_started = true;

    return true;
}

static void soak_join_threads(pthread_t *publisher_threads,
                              size_t publisher_count,
                              pthread_t *query_threads,
                              size_t query_count,
                              pthread_t auditor_thread,
                              bool auditor_started)
{
    for (size_t i = 0; i < publisher_count; ++i) {
        pthread_join(publisher_threads[i], NULL);
    }

    for (size_t i = 0; i < query_count; ++i) {
        pthread_join(query_threads[i], NULL);
    }

    if (auditor_started) {
        pthread_join(auditor_thread, NULL);
    }
}

static void soak_print_final_report(const SoakConfig *config,
                                    double elapsed_seconds,
                                    uint64_t process_cpu_ns,
                                    const SoakCounters *publishers,
                                    const SoakCounters *queries,
                                    const SoakCounters *audits,
                                    const SoakResources *initial_resources,
                                    const SoakResources *final_resources,
                                    const SoakMemoryTrend *memory_trend,
                                    const GeoSegmentCompactionStats *compaction_stats,
                                    uint64_t compaction_runs,
                                    uint64_t compaction_failures)
{
    uint64_t events = publishers->updates + publishers->inserts + publishers->deletes;
    double query_average_us = queries->operations
                                  ? (double) queries->latency_total_ns / (double) queries->operations / 1000.0
                                  : 0.0;
    double publish_average_us = publishers->operations
                                    ? (double) publishers->latency_total_ns / (double) publishers->operations / 1000.0
                                    : 0.0;
    double cpu_percent = elapsed_seconds > 0.0
                             ? (double) process_cpu_ns / (elapsed_seconds * 1000000000.0) * 100.0
                             : 0.0;

    printf("\nSOAK SUMMARY\n");
    printf("duration=%.3fs vehicles=%zu publishers=%zu query_threads=%zu batch_size=%zu\n",
           elapsed_seconds,
           config->vehicle_count,
           config->publisher_count,
           config->query_thread_count,
           config->batch_size);
    printf("events=%" PRIu64 " updates=%" PRIu64 " inserts=%" PRIu64 " deletes=%" PRIu64
           " event_rate=%.0f/s segments_published=%" PRIu64 "\n",
           events,
           publishers->updates,
           publishers->inserts,
           publishers->deletes,
           elapsed_seconds > 0.0 ? (double) events / elapsed_seconds : 0.0,
           publishers->segments);
    printf("queries=%" PRIu64 " query_rate=%.0f/s results=%" PRIu64 " scanned=%" PRIu64
           " audits=%" PRIu64 " failures=%" PRIu64 "\n",
           queries->queries,
           elapsed_seconds > 0.0 ? (double) queries->queries / elapsed_seconds : 0.0,
           queries->results,
           queries->scanned,
           audits->audits,
           publishers->failures + queries->failures + audits->failures);
    printf("query_latency avg=%.2fus p50<=%" PRIu64 "us p95<=%" PRIu64 "us p99<=%" PRIu64
           "us max=%.2fms\n",
           query_average_us,
           soak_histogram_percentile(queries, 50, 100),
           soak_histogram_percentile(queries, 95, 100),
           soak_histogram_percentile(queries, 99, 100),
           (double) queries->latency_max_ns / 1000000.0);
    printf("publish_latency avg=%.2fus p50<=%" PRIu64 "us p95<=%" PRIu64 "us p99<=%" PRIu64
           "us max=%.2fms\n",
           publish_average_us,
           soak_histogram_percentile(publishers, 50, 100),
           soak_histogram_percentile(publishers, 95, 100),
           soak_histogram_percentile(publishers, 99, 100),
           (double) publishers->latency_max_ns / 1000000.0);
    printf("cpu_time=%.3fs aggregate_cpu=%.1f%% worker_cpu=%.3fs\n",
           (double) process_cpu_ns / 1000000000.0,
           cpu_percent,
           (double) (publishers->thread_cpu_ns + queries->thread_cpu_ns + audits->thread_cpu_ns) / 1000000000.0);
    printf("rss_initial=%.1fMiB rss_final=%.1fMiB rss_min=%.1fMiB rss_max=%.1fMiB peak_rss=%.1fMiB trend=%.3fMiB/min\n",
           (double) initial_resources->resident_bytes / (1024.0 * 1024.0),
           (double) final_resources->resident_bytes / (1024.0 * 1024.0),
           (double) memory_trend->minimum_rss / (1024.0 * 1024.0),
           (double) memory_trend->maximum_rss / (1024.0 * 1024.0),
           (double) final_resources->peak_resident_bytes / (1024.0 * 1024.0),
           soak_memory_trend_slope(memory_trend));
    printf("io_read=%.1fMiB io_written=%.1fMiB minor_faults=%" PRIu64 " major_faults=%" PRIu64
           " voluntary_cs=%" PRIu64 " involuntary_cs=%" PRIu64 "\n",
           (double) (final_resources->bytes_read - initial_resources->bytes_read) / (1024.0 * 1024.0),
           (double) (final_resources->bytes_written - initial_resources->bytes_written) / (1024.0 * 1024.0),
           final_resources->minor_faults - initial_resources->minor_faults,
           final_resources->major_faults - initial_resources->major_faults,
           final_resources->voluntary_switches - initial_resources->voluntary_switches,
           final_resources->involuntary_switches - initial_resources->involuntary_switches);
    printf("background_compaction runs=%" PRIu64 " failures=%" PRIu64 " input_segments=%zu records_written=%" PRIu64
           " max_workers=%zu partitions=%zu merge_ms=%.3f\n",
           compaction_runs,
           compaction_failures,
           compaction_stats->input_segments,
           compaction_stats->records_written,
           compaction_stats->worker_count,
           compaction_stats->partition_count,
           compaction_stats->merge_time_ms);
}

static SoakMetrics *soak_allocate_metrics(size_t count)
{
    if (!count || count > SIZE_MAX / sizeof(SoakMetrics)) {
        return NULL;
    }

    void *allocation = NULL;

    if (posix_memalign(&allocation, SOAK_METRICS_ALIGNMENT, count * sizeof(SoakMetrics)) != 0) {
        return NULL;
    }

    memset(allocation, 0, count * sizeof(SoakMetrics));

    return allocation;
}

int main(int argc, char **argv)
{
    SoakConfig config;

    if (!soak_parse_config(argc, argv, &config)) {
        fprintf(stderr,
                "usage: %s [seconds] [vehicles] [publishers] [query-threads] [batch] "
                "[min-delay-ms] [max-delay-ms] [max-segments] [report-seconds] [audit-seconds]\n",
                argv[0]);

        return 2;
    }

    if (config.vehicle_count > SIZE_MAX / sizeof(SoakVehicle) ||
        config.vehicle_count > UINT64_MAX - UINT64_C(1000000000) ||
        config.publisher_count > SIZE_MAX / sizeof(pthread_t) ||
        config.publisher_count > SIZE_MAX / sizeof(SoakPublisher) ||
        config.query_thread_count > SIZE_MAX / sizeof(pthread_t) ||
        config.query_thread_count > SIZE_MAX / sizeof(SoakQueryWorker) ||
        config.batch_size > SIZE_MAX / sizeof(GeoRecord) ||
        config.batch_size > SIZE_MAX / sizeof(SoakLiveChange) ||
        config.batch_size > SIZE_MAX / sizeof(uint64_t) ||
        config.batch_size > SIZE_MAX / sizeof(size_t)) {
        return 2;
    }

    char directory_template[] = "/tmp/geobolt-uber-soak-XXXXXX";
    char *directory = mkdtemp(directory_template);

    if (!directory) {
        perror("mkdtemp");
        return 1;
    }

    char manifest_path[SOAK_PATH_CAPACITY];
    char bootstrap_path[SOAK_PATH_CAPACITY];
    int manifest_length = snprintf(manifest_path, sizeof(manifest_path), "%s/active.manifest", directory);
    int bootstrap_length = snprintf(bootstrap_path, sizeof(bootstrap_path), "%s/bootstrap.geobolt", directory);

    if (manifest_length < 0 || (size_t) manifest_length >= sizeof(manifest_path) ||
        bootstrap_length < 0 || (size_t) bootstrap_length >= sizeof(bootstrap_path)) {
        soak_cleanup_directory(directory);
        return 1;
    }

    SoakVehicle *vehicles = calloc(config.vehicle_count, sizeof(*vehicles));
    pthread_t *publisher_threads = calloc(config.publisher_count, sizeof(*publisher_threads));
    SoakPublisher *publishers = calloc(config.publisher_count, sizeof(*publishers));
    SoakMetrics *publisher_metrics = soak_allocate_metrics(config.publisher_count);
    pthread_t *query_threads = calloc(config.query_thread_count, sizeof(*query_threads));
    SoakQueryWorker *query_workers = calloc(config.query_thread_count, sizeof(*query_workers));
    SoakMetrics *query_metrics = soak_allocate_metrics(config.query_thread_count);
    SoakMetrics *audit_metrics = soak_allocate_metrics(1);
    bool allocations_valid = vehicles &&
                             publisher_threads &&
                             publishers &&
                             publisher_metrics &&
                             query_threads &&
                             query_workers &&
                             query_metrics &&
                             audit_metrics;

    if (!allocations_valid) {
        free(audit_metrics);
        free(query_metrics);
        free(query_workers);
        free(query_threads);
        free(publisher_metrics);
        free(publishers);
        free(publisher_threads);
        free(vehicles);
        soak_cleanup_directory(directory);

        return 1;
    }

    uint64_t initial_checksum = 0;
    bool succeeded = soak_initialize_vehicles(vehicles, config.vehicle_count, &initial_checksum) &&
                     soak_create_bootstrap(bootstrap_path,
                                           vehicles,
                                           config.vehicle_count,
                                           config.publisher_count);
    GeoSegmentSet *set = succeeded ? geo_segment_set_create(manifest_path) : NULL;

    succeeded = succeeded &&
                set &&
                geo_segment_set_add_file(set, bootstrap_path) &&
                geo_segment_set_enable_background_compaction(set, directory, config.maximum_segments);

    if (!succeeded) {
        geo_segment_set_destroy(set);
        free(audit_metrics);
        free(query_metrics);
        free(query_workers);
        free(query_threads);
        free(publisher_metrics);
        free(publishers);
        free(publisher_threads);
        free(vehicles);
        fprintf(stderr, "failed to initialize soak data; artifacts retained at %s\n", directory);

        return 1;
    }

    geo_benchmark_print_environment();
    printf("workload=uber_actor_soak duration=%zus vehicles=%zu publishers=%zu query_threads=%zu batch=%zu "
           "delay_ms=%zu-%zu max_segments=%zu report=%zus audit=%zus initial_checksum=%" PRIu64 "\n",
           config.duration_seconds,
           config.vehicle_count,
           config.publisher_count,
           config.query_thread_count,
           config.batch_size,
           config.minimum_delay_ms,
           config.maximum_delay_ms,
           config.maximum_segments,
           config.report_seconds,
           config.audit_seconds,
           initial_checksum);
    printf("data_directory=%s\n", directory);
    fflush(stdout);

    signal(SIGINT, soak_signal_handler);
    signal(SIGTERM, soak_signal_handler);

    pthread_rwlock_t audit_gate;
    bool audit_gate_initialized = pthread_rwlock_init(&audit_gate, NULL) == 0;
    _Atomic bool start = false;
    _Atomic bool stop = false;
    _Atomic bool failed = false;
    _Atomic uint64_t segment_sequence = 1;
    _Atomic uint64_t audited_active = config.vehicle_count;
    pthread_t auditor_thread = { 0 };
    SoakAuditor auditor = { 0 };
    size_t publishers_started = 0;
    size_t queries_started = 0;
    bool auditor_started = false;

    succeeded = audit_gate_initialized &&
                soak_start_threads(&config,
                                   set,
                                   vehicles,
                                   directory,
                                   &audit_gate,
                                   &start,
                                   &stop,
                                   &failed,
                                   &segment_sequence,
                                   &audited_active,
                                   publisher_threads,
                                   publishers,
                                   publisher_metrics,
                                   query_threads,
                                   query_workers,
                                   query_metrics,
                                   &auditor_thread,
                                   &auditor,
                                   audit_metrics,
                                   &publishers_started,
                                   &queries_started,
                                   &auditor_started);

    SoakResources initial_resources = { 0 };
    SoakResources current_resources = { 0 };
    SoakMemoryTrend memory_trend = { 0 };
    uint64_t workload_start_ns = soak_monotonic_ns();
    uint64_t process_cpu_start_ns = soak_process_cpu_ns();

    soak_read_resources(&initial_resources);
    soak_memory_trend_add(&memory_trend, 0.0, initial_resources.resident_bytes);
    atomic_store_explicit(&start, true, memory_order_release);

    uint64_t duration_ns = config.duration_seconds * UINT64_C(1000000000);
    uint64_t deadline_ns = duration_ns <= UINT64_MAX - workload_start_ns ? workload_start_ns + duration_ns : UINT64_MAX;
    uint64_t next_sample_ns = workload_start_ns + UINT64_C(1000000000);
    uint64_t next_report_ns = workload_start_ns + config.report_seconds * UINT64_C(1000000000);
    uint64_t previous_events = 0;
    uint64_t previous_queries = 0;
    uint64_t previous_report_ns = workload_start_ns;

    while (succeeded &&
           !soak_interrupted &&
           !atomic_load_explicit(&stop, memory_order_acquire) &&
           soak_monotonic_ns() < deadline_ns) {
        soak_sleep_ns(UINT64_C(100000000));

        uint64_t now_ns = soak_monotonic_ns();

        if (now_ns >= next_sample_ns) {
            soak_read_resources(&current_resources);
            soak_memory_trend_add(&memory_trend,
                                  (double) (now_ns - workload_start_ns) / 60000000000.0,
                                  current_resources.resident_bytes);
            next_sample_ns = now_ns + UINT64_C(1000000000);
        }

        if (now_ns >= next_report_ns) {
            SoakCounters publisher_totals;
            SoakCounters query_totals;
            SoakCounters audit_totals;

            soak_collect_metrics(publisher_metrics,
                                 config.publisher_count,
                                 query_metrics,
                                 config.query_thread_count,
                                 audit_metrics,
                                 &publisher_totals,
                                 &query_totals,
                                 &audit_totals);
            soak_read_resources(&current_resources);
            soak_print_progress((double) (now_ns - workload_start_ns) / 1000000000.0,
                                (double) (now_ns - previous_report_ns) / 1000000000.0,
                                &publisher_totals,
                                &query_totals,
                                previous_events,
                                previous_queries,
                                atomic_load_explicit(&audited_active, memory_order_relaxed),
                                &current_resources,
                                set);
            previous_events = publisher_totals.updates + publisher_totals.inserts + publisher_totals.deletes;
            previous_queries = query_totals.queries;
            previous_report_ns = now_ns;
            next_report_ns = now_ns + config.report_seconds * UINT64_C(1000000000);
        }
    }

    uint64_t workload_stop_ns = soak_monotonic_ns();
    uint64_t process_cpu_stop_ns = soak_process_cpu_ns();

    atomic_store_explicit(&stop, true, memory_order_release);
    soak_join_threads(publisher_threads,
                      publishers_started,
                      query_threads,
                      queries_started,
                      auditor_thread,
                      auditor_started);

    GeoSegmentCompactionStats last_compaction_stats = { 0 };
    GeoSegmentCompactionStats compaction_stats = { 0 };
    uint64_t compaction_runs = 0;
    uint64_t compaction_failures = 0;
    bool compaction_succeeded = geo_segment_set_wait_for_background_compaction(set, &last_compaction_stats) &&
                                geo_segment_set_background_compaction_totals(set,
                                                                            &compaction_stats,
                                                                            &compaction_runs,
                                                                            &compaction_failures);
    SoakCounters publisher_totals;
    SoakCounters query_totals;
    SoakCounters audit_totals;

    soak_collect_metrics(publisher_metrics,
                         config.publisher_count,
                         query_metrics,
                         config.query_thread_count,
                         audit_metrics,
                         &publisher_totals,
                         &query_totals,
                         &audit_totals);
    soak_read_resources(&current_resources);
    soak_memory_trend_add(&memory_trend,
                          (double) (workload_stop_ns - workload_start_ns) / 60000000000.0,
                          current_resources.resident_bytes);

    double elapsed_seconds = (double) (workload_stop_ns - workload_start_ns) / 1000000000.0;
    uint64_t process_cpu_ns = process_cpu_stop_ns - process_cpu_start_ns;

    soak_print_final_report(&config,
                            elapsed_seconds,
                            process_cpu_ns,
                            &publisher_totals,
                            &query_totals,
                            &audit_totals,
                            &initial_resources,
                            &current_resources,
                            &memory_trend,
                            &compaction_stats,
                            compaction_runs,
                            compaction_failures);

    uint64_t completed_events = publisher_totals.updates + publisher_totals.inserts + publisher_totals.deletes;

    succeeded = succeeded &&
                compaction_succeeded &&
                compaction_failures == 0 &&
                !atomic_load_explicit(&failed, memory_order_acquire) &&
                !soak_interrupted &&
                completed_events > 0 &&
                query_totals.queries > 0 &&
                audit_totals.audits > 0;

    geo_segment_set_destroy(set);

    if (audit_gate_initialized) {
        pthread_rwlock_destroy(&audit_gate);
    }

    free(audit_metrics);
    free(query_metrics);
    free(query_workers);
    free(query_threads);
    free(publisher_metrics);
    free(publishers);
    free(publisher_threads);
    free(vehicles);

    if (succeeded) {
        succeeded = soak_cleanup_directory(directory);
    } else {
        fprintf(stderr, "soak failed; artifacts retained at %s\n", directory);
    }

    printf("status=%s\n", succeeded ? "PASS" : "FAIL");

    return succeeded ? 0 : 1;
}
