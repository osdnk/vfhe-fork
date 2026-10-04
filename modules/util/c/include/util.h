// SPDX-FileCopyrightText: 2026 Antonio Guimarães <antonio.guimaraes@imdea.org>
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <stddef.h>
#include <stdint.h>

#include <engine.h>
#include <vfhe_cpu.h>

#ifdef __cplusplus
extern "C"
{
#endif

    // Index and modulus helpers
    void array_reduce_mod_N(uint64_t *out, uint64_t *in, uint64_t size, uint64_t p);
    void array_mod_switch_from_2k(uint64_t *out, uint64_t *in, uint64_t p, uint64_t q, uint64_t n);
    uint64_t double2int(double x);
    uint32_t int_rev(uint32_t b);
    void bit_rev(uint64_t *out, uint64_t *in, uint64_t n, uint64_t log_n);

    // Allocation that aborts rather than returning NULL
    void *safe_malloc(size_t size);
    void *safe_realloc(void *ptr, size_t size);
    // 64-byte aligned, as the SIMD kernels require; a large buffer also asks
    // for huge pages. Release with plain `free`.
    void *safe_aligned_malloc(size_t size);

    // Parallelism, under one library-wide limit on the threads any operation
    // uses: VFHE_NUM_THREADS if set, 1 otherwise, until vfhe_set_num_threads
    // changes it (0 restores that default).
    uint64_t vfhe_num_threads(void);
    void vfhe_set_num_threads(uint64_t n);
    // The threads an operation over `n_items` independent items should use
    // when its caller asked for `requested` (0: as many as the limit allows):
    // at most the limit and the items, and 1 inside a vfhe_parallel_for body,
    // so nested parallel operations do not multiply threads.
    uint64_t vfhe_threads_for(uint64_t requested, uint64_t n_items);
    // body(ctx, i) for every i < n, on vfhe_threads_for(n_threads, n) threads
    // including the caller's; returns when all are done. Bodies run
    // concurrently and in no particular order.
    void vfhe_parallel_for(uint64_t n, uint64_t n_threads, void (*body)(void *ctx, uint64_t i),
                           void *ctx);

    // Which engine this binary is (CPU capability lives in vfhe_cpu.h, which
    // this header includes).
    const char *vfhe_engine_active(void); // e.g. "portable", "avx512ifma"

#ifdef __cplusplus
}
#endif
