#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "geobolt/server.h"

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <pthread.h>
#include <signal.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#define GEOBOLTD_VERSION "0.1.0"
#define GEOBOLTD_MAX_TOKEN_SIZE 4096U

typedef struct {
    GeoServer *server;
    sigset_t signals;
    atomic_bool finished;
    atomic_bool failed;
} SignalThreadContext;

static void print_usage(FILE *stream, const char *program)
{
    fprintf(stream,
            "Usage: %s --data DIRECTORY [options]\n"
            "\n"
            "Options:\n"
            "  --bind ADDRESS          IPv4 address to bind (default: 127.0.0.1)\n"
            "  --port PORT             TCP port, or 0 for an ephemeral port (default: 7447)\n"
            "  --token-file PATH       Authentication token file (permissions must be 0600 or stricter)\n"
            "  --workers COUNT         Persistent database worker threads (default: 4)\n"
            "  --queue-capacity COUNT  Maximum queued commands (default: 1024)\n"
            "  --max-connections COUNT Maximum simultaneous clients (default: 4096)\n"
            "  --max-frame BYTES       Maximum request/response payload (default: 16777216)\n"
            "  --max-inflight BYTES    Global retained request-payload budget (default: 268435456)\n"
            "  --backlog COUNT         TCP listen backlog (default: 512)\n"
            "  --help                   Show this help\n"
            "  --version                Show the server version\n",
            program);
}

static const char *server_status_name(GeoServerStatus status)
{
    switch (status) {
        case GEO_SERVER_OK:
            return "ok";
        case GEO_SERVER_INVALID_ARGUMENT:
            return "invalid configuration";
        case GEO_SERVER_OUT_OF_MEMORY:
            return "out of memory";
        case GEO_SERVER_DATABASE_ERROR:
            return "database initialization failed";
        case GEO_SERVER_WORKER_ERROR:
            return "worker pool initialization failed";
        case GEO_SERVER_SOCKET_ERROR:
            return "socket bind/listen failed";
        case GEO_SERVER_SYNCHRONIZATION_ERROR:
            return "synchronization primitive initialization failed";
        case GEO_SERVER_DESCRIPTOR_ERROR:
            return "event or socket descriptor creation failed";
        case GEO_SERVER_EVENT_REGISTRATION_ERROR:
            return "event descriptor registration failed";
        case GEO_SERVER_EVENT_LOOP_ERROR:
            return "event loop failed";
        case GEO_SERVER_ALREADY_RUNNING:
            return "server is already running";
    }

    return "unknown error";
}

static bool parse_u64(const char *text, uint64_t maximum, uint64_t *value)
{
    if (!text || !text[0] || text[0] == '-') {
        return false;
    }

    errno = 0;
    char *end = NULL;
    unsigned long long parsed = strtoull(text, &end, 10);

    if (errno != 0 || !end || *end != '\0' || parsed > maximum) {
        return false;
    }

    *value = (uint64_t) parsed;

    return true;
}

static void erase_secret(char *secret, size_t size)
{
    volatile unsigned char *bytes = (volatile unsigned char *) secret;

    while (size) {
        *bytes++ = 0;
        size--;
    }
}

static char *read_token_file(const char *path, size_t *token_size)
{
    int descriptor = open(path, O_RDONLY | O_CLOEXEC);

    if (descriptor < 0) {
        fprintf(stderr, "geoboltd: cannot open token file '%s': %s\n", path, strerror(errno));
        return NULL;
    }

    struct stat metadata;

    if (fstat(descriptor, &metadata) != 0 || !S_ISREG(metadata.st_mode) || (metadata.st_mode & 0077) != 0) {
        fprintf(stderr, "geoboltd: token file must be a regular file with permissions 0600 or stricter\n");
        (void) close(descriptor);
        return NULL;
    }

    char *token = malloc(GEOBOLTD_MAX_TOKEN_SIZE + 3U);

    if (!token) {
        fprintf(stderr, "geoboltd: cannot allocate authentication token\n");
        (void) close(descriptor);
        return NULL;
    }

    size_t used = 0;

    while (used < GEOBOLTD_MAX_TOKEN_SIZE + 3U) {
        ssize_t received = read(descriptor, token + used, GEOBOLTD_MAX_TOKEN_SIZE + 3U - used);

        if (received > 0) {
            used += (size_t) received;
            continue;
        }

        if (received < 0 && errno == EINTR) {
            continue;
        }

        if (received < 0) {
            fprintf(stderr, "geoboltd: cannot read token file '%s': %s\n", path, strerror(errno));
            erase_secret(token, GEOBOLTD_MAX_TOKEN_SIZE + 3U);
            free(token);
            token = NULL;
        }

        break;
    }

    (void) close(descriptor);

    if (!token) {
        return NULL;
    }

    if (memchr(token, '\0', used)) {
        fprintf(stderr, "geoboltd: token must not contain NUL bytes\n");
        erase_secret(token, GEOBOLTD_MAX_TOKEN_SIZE + 3U);
        free(token);
        return NULL;
    }

    while (used && (token[used - 1U] == '\n' || token[used - 1U] == '\r')) {
        used--;
    }

    if (!used || used > GEOBOLTD_MAX_TOKEN_SIZE) {
        fprintf(stderr, "geoboltd: token must contain 1 to %u bytes\n", GEOBOLTD_MAX_TOKEN_SIZE);
        erase_secret(token, GEOBOLTD_MAX_TOKEN_SIZE + 3U);
        free(token);
        return NULL;
    }

    token[used] = '\0';
    *token_size = used;

    return token;
}

static void *signal_thread_main(void *argument)
{
    SignalThreadContext *context = argument;
    const struct timespec interval = {
        .tv_nsec = 100000000L,
    };

    while (!atomic_load_explicit(&context->finished, memory_order_acquire)) {
        int signal_number = sigtimedwait(&context->signals, NULL, &interval);

        if (signal_number == SIGINT || signal_number == SIGTERM) {
            geo_server_request_stop(context->server);
            break;
        }

        if (signal_number < 0 && errno != EAGAIN && errno != EINTR) {
            atomic_store_explicit(&context->failed, true, memory_order_release);
            geo_server_request_stop(context->server);
            break;
        }
    }

    return NULL;
}

int main(int argc, char **argv)
{
    const char *database_directory = NULL;
    const char *token_file = NULL;
    GeoServerConfig config = geo_server_default_config(NULL);

    for (int argument = 1; argument < argc; ++argument) {
        const char *option = argv[argument];

        if (strcmp(option, "--help") == 0) {
            print_usage(stdout, argv[0]);
            return EXIT_SUCCESS;
        }

        if (strcmp(option, "--version") == 0) {
            printf("geoboltd %s\n", GEOBOLTD_VERSION);
            return EXIT_SUCCESS;
        }

        if (argument + 1 >= argc) {
            fprintf(stderr, "geoboltd: option '%s' requires a value\n", option);
            print_usage(stderr, argv[0]);
            return EXIT_FAILURE;
        }

        const char *value = argv[++argument];
        uint64_t parsed;

        if (strcmp(option, "--data") == 0) {
            database_directory = value;
        } else if (strcmp(option, "--bind") == 0) {
            config.bind_address = value;
        } else if (strcmp(option, "--token-file") == 0) {
            token_file = value;
        } else if (strcmp(option, "--port") == 0 && parse_u64(value, UINT16_MAX, &parsed)) {
            config.port = (uint16_t) parsed;
        } else if (strcmp(option, "--workers") == 0 && parse_u64(value, SIZE_MAX, &parsed) && parsed) {
            config.worker_threads = (size_t) parsed;
        } else if (strcmp(option, "--queue-capacity") == 0 && parse_u64(value, SIZE_MAX, &parsed) && parsed) {
            config.work_queue_capacity = (size_t) parsed;
        } else if (strcmp(option, "--max-connections") == 0 && parse_u64(value, SIZE_MAX, &parsed) && parsed) {
            config.max_connections = (size_t) parsed;
        } else if (strcmp(option, "--max-frame") == 0 && parse_u64(value, SIZE_MAX, &parsed) && parsed) {
            config.max_frame_size = (size_t) parsed;
        } else if (strcmp(option, "--max-inflight") == 0 && parse_u64(value, SIZE_MAX, &parsed) && parsed) {
            config.max_inflight_payload_bytes = (size_t) parsed;
        } else if (strcmp(option, "--backlog") == 0 && parse_u64(value, INT_MAX, &parsed) && parsed) {
            config.listen_backlog = (int) parsed;
        } else {
            fprintf(stderr, "geoboltd: invalid option or value: %s %s\n", option, value);
            return EXIT_FAILURE;
        }
    }

    if (!database_directory) {
        fprintf(stderr, "geoboltd: --data is required\n");
        print_usage(stderr, argv[0]);
        return EXIT_FAILURE;
    }

    config.database_directory = database_directory;
    char *authentication_token = NULL;
    size_t authentication_token_size = 0;

    if (token_file) {
        authentication_token = read_token_file(token_file, &authentication_token_size);

        if (!authentication_token) {
            return EXIT_FAILURE;
        }

        config.authentication_token = authentication_token;
    }

    sigset_t signals;

    sigemptyset(&signals);
    sigaddset(&signals, SIGINT);
    sigaddset(&signals, SIGTERM);

    if (pthread_sigmask(SIG_BLOCK, &signals, NULL) != 0) {
        fprintf(stderr, "geoboltd: cannot configure signal handling\n");
        erase_secret(authentication_token, authentication_token_size);
        free(authentication_token);
        return EXIT_FAILURE;
    }

    GeoServerStatus server_status;
    GeoServer *server = geo_server_create(&config, &server_status);

    erase_secret(authentication_token, authentication_token_size);
    free(authentication_token);

    if (!server) {
        fprintf(stderr, "geoboltd: startup failed: %s\n", server_status_name(server_status));
        return EXIT_FAILURE;
    }

    SignalThreadContext signal_context = {
        .server = server,
        .signals = signals,
        .finished = ATOMIC_VAR_INIT(false),
        .failed = ATOMIC_VAR_INIT(false),
    };
    pthread_t signal_thread;

    if (pthread_create(&signal_thread, NULL, signal_thread_main, &signal_context) != 0) {
        fprintf(stderr, "geoboltd: cannot start signal monitor\n");
        geo_server_destroy(server);
        return EXIT_FAILURE;
    }

    printf("geoboltd: listening on %s:%u\n", config.bind_address, geo_server_port(server));
    fflush(stdout);
    server_status = geo_server_run(server);
    atomic_store_explicit(&signal_context.finished, true, memory_order_release);
    int signal_join_status = pthread_join(signal_thread, NULL);
    bool signal_monitor_failed = atomic_load_explicit(&signal_context.failed, memory_order_acquire);

    if (signal_join_status != 0) {
        signal_monitor_failed = true;
    }

    GeoServerStats stats;

    if (geo_server_get_stats(server, &stats)) {
        printf("geoboltd: stopped; connections=%llu requests=%llu rejected=%llu backpressure=%llu peak_payload_bytes=%llu\n",
               (unsigned long long) stats.accepted_connections,
               (unsigned long long) stats.completed_requests,
               (unsigned long long) stats.rejected_connections,
               (unsigned long long) stats.backpressure_rejections,
               (unsigned long long) stats.peak_inflight_payload_bytes);
    }

    geo_server_destroy(server);

    if (signal_monitor_failed) {
        fprintf(stderr, "geoboltd: signal monitor failed\n");
        return EXIT_FAILURE;
    }

    if (server_status != GEO_SERVER_OK) {
        fprintf(stderr, "geoboltd: event loop stopped with error: %s\n", server_status_name(server_status));
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
