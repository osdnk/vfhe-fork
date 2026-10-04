// SPDX-FileCopyrightText: 2026 Antonio Guimarães <antonio.guimaraes@imdea.org>
// SPDX-License-Identifier: Apache-2.0
#include "util.h"

#include <pthread.h>
#include <stdatomic.h>
#include <stdlib.h>

// 0 until first read or after a reset: the default is decided lazily, so the
// environment is read when the library is first used rather than when it loads.
static atomic_uint_fast64_t thread_limit = 0;

// Set on every thread while it runs loop bodies, the caller's included: a loop
// started from inside a body runs on the thread it is called from instead of
// multiplying the threads already running.
static _Thread_local int inside_parallel_loop = 0;

static uint64_t default_thread_limit(void)
{
    const char *env = getenv("VFHE_NUM_THREADS");
    if (env != NULL && *env != '\0')
    {
        char *end;
        const unsigned long long value = strtoull(env, &end, 10);
        if (*end == '\0' && value > 0)
            return (uint64_t)value;
    }
    return 1;
}

uint64_t vfhe_num_threads(void)
{
    uint_fast64_t limit = atomic_load(&thread_limit);
    if (limit == 0)
    {
        uint_fast64_t unset = 0;
        atomic_compare_exchange_strong(&thread_limit, &unset, default_thread_limit());
        limit = atomic_load(&thread_limit);
    }
    return (uint64_t)limit;
}

void vfhe_set_num_threads(uint64_t n) { atomic_store(&thread_limit, n); }

uint64_t vfhe_threads_for(uint64_t requested, uint64_t n_items)
{
    if (inside_parallel_loop || n_items <= 1)
        return 1;
    const uint64_t limit = vfhe_num_threads();
    uint64_t n = (requested == 0 || requested > limit) ? limit : requested;
    return n < n_items ? n : n_items;
}

typedef struct
{
    void (*body)(void *ctx, uint64_t i);
    void *ctx;
    uint64_t n;
    atomic_uint_fast64_t next;
} ParallelLoop;

// Items are handed out one at a time, so threads that draw cheap items keep
// drawing while another finishes an expensive one.
static void run_loop(ParallelLoop *loop)
{
    inside_parallel_loop = 1;
    for (;;)
    {
        const uint64_t i = (uint64_t)atomic_fetch_add(&loop->next, 1);
        if (i >= loop->n)
            break;
        loop->body(loop->ctx, i);
    }
    inside_parallel_loop = 0;
}

static void *run_loop_thread(void *loop)
{
    run_loop((ParallelLoop *)loop);
    return NULL;
}

void vfhe_parallel_for(uint64_t n, uint64_t n_threads, void (*body)(void *ctx, uint64_t i),
                       void *ctx)
{
    const uint64_t threads = vfhe_threads_for(n_threads, n);
    if (threads <= 1)
    {
        for (uint64_t i = 0; i < n; i++)
            body(ctx, i);
        return;
    }
    ParallelLoop loop = {body, ctx, n, 0};
    pthread_t *helpers = (pthread_t *)safe_malloc((threads - 1) * sizeof(pthread_t));
    uint64_t started = 0;
    // A thread that cannot be created only means fewer of them: the caller
    // runs whatever the others do not take.
    while (started < threads - 1 &&
           pthread_create(&helpers[started], NULL, run_loop_thread, &loop) == 0)
        started++;
    run_loop(&loop);
    for (uint64_t t = 0; t < started; t++)
        pthread_join(helpers[t], NULL);
    free(helpers);
}
