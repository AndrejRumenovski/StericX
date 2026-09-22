#define _GNU_SOURCE
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
static long long ns(void) {
    struct timespec value;
    if (clock_gettime(CLOCK_MONOTONIC, &value)) { perror("clock_gettime"); exit(125); }
    return (long long)value.tv_sec * 1000000000LL + value.tv_nsec;
}
int main(int argc, char **argv) {
    if (argc < 3) return 125;
    long long started = ns();
    pid_t child = fork();
    if (child < 0) { perror("fork"); return 125; }
    if (child == 0) { execv(argv[2], &argv[2]); perror("execv"); _exit(127); }
    int status = 0;
    struct rusage usage;
    while (wait4(child, &status, 0, &usage) < 0) {
        if (errno != EINTR) { perror("wait4"); return 125; }
    }
    long long elapsed = ns() - started;
    int code = WIFEXITED(status) ? WEXITSTATUS(status) : -WTERMSIG(status);
    FILE *out = fopen(argv[1], "w");
    if (!out) { perror("fopen metrics"); return 125; }
    fprintf(out, "{\"pid\":%ld,\"wall_ns\":%lld,\"cpu_user_s\":%.6f,"
        "\"cpu_system_s\":%.6f,\"peak_rss_bytes\":%ld,\"minor_faults\":%ld,"
        "\"major_faults\":%ld,\"voluntary_context_switches\":%ld,"
        "\"involuntary_context_switches\":%ld,\"returncode\":%d}\n",
        (long)child, elapsed, usage.ru_utime.tv_sec + usage.ru_utime.tv_usec / 1e6,
        usage.ru_stime.tv_sec + usage.ru_stime.tv_usec / 1e6, usage.ru_maxrss * 1024,
        usage.ru_minflt, usage.ru_majflt, usage.ru_nvcsw, usage.ru_nivcsw, code);
    if (fclose(out)) { perror("fclose metrics"); return 125; }
    return code < 0 ? 128 - code : code;
}
