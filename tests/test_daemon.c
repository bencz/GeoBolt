#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "geobolt/client.h"
#include "geobolt/geobolt.h"
#include "test_support.h"

#include <errno.h>
#include <fcntl.h>
#include <signal.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#define DAEMON_ASSERT(condition, message)                                                                                  \
    do {                                                                                                                   \
        if (!(condition)) {                                                                                                \
            fprintf(stderr, "FAIL: %s (%s:%d)\n", message, __FILE__, __LINE__);                                          \
            succeeded = false;                                                                                             \
            goto cleanup;                                                                                                  \
        }                                                                                                                  \
    } while (0)

static bool write_all(int descriptor, const void *data, size_t size)
{
    const unsigned char *bytes = data;

    while (size) {
        ssize_t written = write(descriptor, bytes, size);

        if (written < 0 && errno == EINTR) {
            continue;
        }

        if (written <= 0) {
            return false;
        }

        bytes += (size_t) written;
        size -= (size_t) written;
    }

    return true;
}

static bool read_daemon_port(int descriptor, uint16_t *port)
{
    char line[256];
    size_t used = 0;

    while (used + 1U < sizeof(line)) {
        ssize_t received = read(descriptor, line + used, 1U);

        if (received < 0 && errno == EINTR) {
            continue;
        }

        if (received != 1) {
            return false;
        }

        if (line[used++] == '\n') {
            break;
        }
    }

    line[used] = '\0';
    unsigned parsed_port;

    if (sscanf(line, "geoboltd: listening on 127.0.0.1:%u", &parsed_port) != 1 || parsed_port > UINT16_MAX) {
        return false;
    }

    *port = (uint16_t) parsed_port;

    return *port != 0;
}

static pid_t start_daemon(const char *database_directory, const char *token_file, int *output_descriptor)
{
    int output_pipe[2];

    if (pipe(output_pipe) != 0) {
        return -1;
    }

    pid_t child = fork();

    if (child == 0) {
        const char *daemon_binary = getenv("GEOBOLTD_TEST_BINARY");

        if (!daemon_binary || !daemon_binary[0]) {
            daemon_binary = "./bin/geoboltd";
        }

        (void) close(output_pipe[0]);

        if (dup2(output_pipe[1], STDOUT_FILENO) < 0) {
            _exit(126);
        }

        (void) close(output_pipe[1]);
        execl(daemon_binary,
              "geoboltd",
              "--data",
              database_directory,
              "--port",
              "0",
              "--token-file",
              token_file,
              "--workers",
              "2",
              "--queue-capacity",
              "16",
              (char *) NULL);
        _exit(127);
    }

    (void) close(output_pipe[1]);

    if (child < 0) {
        (void) close(output_pipe[0]);
        return -1;
    }

    *output_descriptor = output_pipe[0];

    return child;
}

int main(void)
{
    bool succeeded = true;
    char root_template[] = "/tmp/geobolt-daemon-test-XXXXXX";
    char *root = mkdtemp(root_template);
    char database_directory[512];
    char token_path[512];
    int token_descriptor = -1;
    int daemon_output = -1;
    pid_t daemon_process = -1;
    GeoClient *client = NULL;

    DAEMON_ASSERT(root != NULL, "daemon test root must be created");
    DAEMON_ASSERT(snprintf(database_directory, sizeof(database_directory), "%s/data", root) > 0,
                  "database path must fit");
    DAEMON_ASSERT(snprintf(token_path, sizeof(token_path), "%s/token", root) > 0, "token path must fit");

    token_descriptor = open(token_path, O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, S_IRUSR | S_IWUSR);
    DAEMON_ASSERT(token_descriptor >= 0, "secure token file must be created");

    static const char token[] = "daemon integration secret\n";

    DAEMON_ASSERT(write_all(token_descriptor, token, sizeof(token) - 1U), "token file must be written");
    DAEMON_ASSERT(close(token_descriptor) == 0, "token file must close");
    token_descriptor = -1;

    daemon_process = start_daemon(database_directory, token_path, &daemon_output);
    DAEMON_ASSERT(daemon_process > 0, "daemon process must start");

    uint16_t port;

    DAEMON_ASSERT(read_daemon_port(daemon_output, &port), "daemon must publish its bound port");

    GeoClientConfig client_config = geo_client_default_config("127.0.0.1", port);

    client_config.authentication_token = "daemon integration secret";
    GeoClientStatus client_status;

    client = geo_client_connect(&client_config, &client_status);
    DAEMON_ASSERT(client && client_status == GEO_CLIENT_OK, "driver must authenticate with daemon");

    GeoDatabaseMutation mutations[] = {
        { .object_id = 101U, .morton_code = geo_encode(-23.5505, -46.6333), .operation = GEO_DATABASE_UPSERT },
        { .object_id = 202U, .morton_code = geo_encode(40.7128, -74.0060), .operation = GEO_DATABASE_UPSERT },
    };

    DAEMON_ASSERT(geo_client_write(client, mutations, 2U) == GEO_CLIENT_OK, "daemon must durably commit a network batch");
    geo_client_close(client);
    client = NULL;

    DAEMON_ASSERT(kill(daemon_process, SIGTERM) == 0, "SIGTERM must request daemon shutdown");

    int child_status;

    DAEMON_ASSERT(waitpid(daemon_process, &child_status, 0) == daemon_process, "daemon process must be reaped");
    daemon_process = -1;
    DAEMON_ASSERT(WIFEXITED(child_status) && WEXITSTATUS(child_status) == EXIT_SUCCESS,
                  "daemon must exit successfully after SIGTERM");

    GeoDatabaseConfig database_config = geo_database_default_config();
    GeoDatabaseStatus database_status;
    GeoDatabase *database = geo_database_open(database_directory, &database_config, &database_status);

    DAEMON_ASSERT(database && database_status == GEO_DATABASE_OK, "database must reopen after daemon shutdown");

    size_t count = 0;

    DAEMON_ASSERT(geo_database_search_radius_count(database, 0.0, 0.0, 25000.0, &count, NULL) && count == 2U,
                  "daemon shutdown must preserve committed objects");
    geo_database_close(database);

cleanup:
    geo_client_close(client);

    if (daemon_process > 0) {
        (void) kill(daemon_process, SIGKILL);
        (void) waitpid(daemon_process, NULL, 0);
    }

    if (daemon_output >= 0) {
        (void) close(daemon_output);
    }

    if (token_descriptor >= 0) {
        (void) close(token_descriptor);
    }

    if (root && !geobolt_test_remove_tree(root)) {
        succeeded = false;
    }

    if (succeeded) {
        puts("[PASS] daemon configuration, authentication, durable write, and SIGTERM shutdown");
    }

    return succeeded ? EXIT_SUCCESS : EXIT_FAILURE;
}
