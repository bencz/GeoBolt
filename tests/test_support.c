#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#include "test_support.h"

#include <dirent.h>
#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

bool geobolt_test_remove_tree(const char *path)
{
    DIR *directory = opendir(path);

    if (!directory) {
        return errno == ENOENT;
    }

    bool succeeded = true;

    while (succeeded) {
        errno = 0;
        struct dirent *entry = readdir(directory);

        if (!entry) {
            succeeded = errno == 0;
            break;
        }

        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) {
            continue;
        }

        size_t path_length = strlen(path);
        size_t name_length = strlen(entry->d_name);

        if (path_length > SIZE_MAX - name_length - 2U) {
            succeeded = false;
            break;
        }

        char *child = malloc(path_length + name_length + 2U);

        if (!child) {
            succeeded = false;
            break;
        }

        snprintf(child, path_length + name_length + 2U, "%s/%s", path, entry->d_name);

        struct stat status;

        if (lstat(child, &status) != 0) {
            succeeded = false;
        } else if (S_ISDIR(status.st_mode)) {
            succeeded = geobolt_test_remove_tree(child);
        } else {
            succeeded = unlink(child) == 0;
        }

        free(child);
    }

    if (closedir(directory) != 0) {
        succeeded = false;
    }

    return succeeded && rmdir(path) == 0;
}
