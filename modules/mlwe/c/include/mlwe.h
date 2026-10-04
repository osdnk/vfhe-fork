// SPDX-FileCopyrightText: 2026 Antonio Guimarães <antonio.guimaraes@imdea.org>
// SPDX-License-Identifier: Apache-2.0
#pragma once
#include <arith.h>
#include <arith_generic.h>

#ifdef __cplusplus
extern "C"
{
#endif
    void gen_sparse_ternary_array_modq(uint64_t *out, uint64_t size, uint64_t h, uint64_t q);
    /* LWE
     *
     * An LWE key or sample lives over the primes of `base` that `mask`
     * selects, one limb per prime: limb i is the residue modulo the prime at
     * the i-th set bit of `mask`, in ascending base index, and `l` is the
     * number of set bits. A key and the samples it touches share one mask.
     */

    typedef struct _LWE_Key
    {
        uint64_t **s;
        uint64_t n, l, mask;
        RNS_Base base;
        double sigma;
    } *LWE_Key;

    typedef struct _LWE
    {
        uint64_t **a;
        uint64_t *b;
        uint64_t n, l, mask;
        RNS_Base base;
    } *LWE;

    typedef struct _LWE_KS_Key
    {
        LWE ***s;
        uint64_t base_bit, t;
    } *LWE_KS_Key;

    // mlwe rns
    /* MLWE RNS */

    typedef struct _RNS_MLWE_Key
    {
        // The secret, one element per module component, held in the mul
        // domain where every use of it multiplies.
        ArithElement *s;
        ArithRing ring;
        uint64_t N, l, r;
        double sigma;
    } *RNS_MLWE_Key;

    // A module-LWE sample: `r` mask components and a body, over one ring.
    //
    // Every component of a sample is in the same domain, which the elements
    // themselves carry: `mlwe_domain` reads it, and the whole-sample
    // conversions move all of them together. A caller that puts components in
    // different domains breaks that invariant, and the arithmetic will refuse
    // them.
    typedef struct _MLWE
    {
        ArithElement *a, b;
        uint64_t r;
        ArithRing ring;
    } *MLWE;

    // Aliases that let a signature say which domain it expects its argument
    // in. They are the same type, so this is documentation, not enforcement:
    // check `mlwe_domain` where it matters.
    typedef MLWE RNS_MLWE;
    typedef MLWE RNSc_MLWE;

    ArithDomain mlwe_domain(MLWE c);
    MLWE mlwe_alloc_sample(ArithRing ring, uint64_t r);

    // A key-switch key: one gadget-decomposed key array per input component
    // (NULL marks a component that keeps the target key and passes through),
    // and the ring the key lives in. A key switch must accumulate in that
    // ring -- the caller's `out` is only guaranteed to be allocated for its
    // own, narrower one -- so it allocates its scratch there and copies the
    // finished result out. The key stays immutable and therefore shareable:
    // gp25 hands one key to every thread of a parallel bootstrap.
    //
    // `log_base` is the gadget the keys were generated for, and the key switch
    // decomposes against it: 0 for the RNS gadget (one key per prime), or the
    // radix base's log for the radix one (one key per prime and digit; see
    // `gadget_radix_digits`). A key generated for one and read as the other
    // decrypts to garbage, which is why it travels with the key rather than
    // with the call.
    typedef struct _RNS_MLWE_KS_Key
    {
        RNS_MLWE **s;
        uint64_t count;
        uint64_t mask;
        uint64_t log_base;
        ArithRing ring;
    } *RNS_MLWE_KS_Key;

    // mlwe rns
    RNS_MLWE_Key mlwe_alloc_RNS_key_special_primes(uint64_t N, uint64_t r, uint64_t l,
                                                   uint64_t special_primes, RNS_Base base,
                                                   double sigma);
    void free_polynomial_array(uint64_t size, IntPolynomial *p);
    RNS_MLWE_Key mlwe_get_RNS_key_from_array(uint64_t N, uint64_t r, uint64_t l, uint64_t *array,
                                             RNS_Base base, double sigma);
    RNS_MLWE_Key mlwe_alloc_key(ArithRing ring, uint64_t r, uint64_t l, double sigma);
    RNS_MLWE_Key mlwe_alloc_RNS_key(uint64_t N, uint64_t r, uint64_t l, RNS_Base base,
                                    double sigma);
    void free_RNS_mlwe_sample(RNS_MLWE c);
    void free_mlwe_RNS_key(RNS_MLWE_Key key);
    LWE mlwe_extract_LWE(RNSc_MLWE in, uint64_t idx);
    RNS_MLWE_Key mlwe_new_RNS_gaussian_key(uint64_t N, uint64_t r, uint64_t l, double key_sigma,
                                           RNS_Base base, double sigma);
    RNS_MLWE mlwe_alloc_RNS_sample(uint64_t N, uint64_t r, uint64_t mask, RNS_Base base);
    RNSc_MLWE mlwe_alloc_RNSc_sample(uint64_t N, uint64_t r, uint64_t mask, RNS_Base base);
    RNS_MLWE mlwe_new_RNS_sample(RNS_MLWE_Key key, uint64_t *m, uint64_t p);
    void mlwe_RNS_sample_of_zero(RNS_MLWE out, RNS_MLWE_Key key);
    void mlwe_RNSc_sample_of_zero(RNSc_MLWE out, RNS_MLWE_Key key);
    RNS_MLWE mlwe_new_RNS_sample_of_zero(RNS_MLWE_Key key);
    RNSc_MLWE mlwe_new_RNSc_sample_of_zero(RNS_MLWE_Key key);
    RNS_MLWE mlwe_new_RNS_trivial_sample_of_zero(uint64_t N, uint64_t r, uint64_t mask,
                                                 RNS_Base base);
    void mlwe_RNS_phase(ArithElement *out, RNS_MLWE in, RNS_MLWE_Key key);
    void mlwe_RNSc_to_RNS(RNS_MLWE out, RNSc_MLWE in);
    void mlwe_RNS_to_RNSc(RNSc_MLWE out, RNS_MLWE in);

    void mlwe_RNS_trivial_sample_of_zero(RNS_MLWE out);
    void mlwe_copy_RNS_sample(RNS_MLWE out, RNS_MLWE in);
    void mlwe_copy_RNSc_sample(RNSc_MLWE out, RNSc_MLWE in);
    void mlwe_RNSc_sample(RNSc_MLWE out, RNS_MLWE_Key key, const ArithElement *m);
    void mlwe_scale_RNS_mlwe_RNS(RNS_MLWE c, const uint64_t *per_component);
    void mlwe_add_RNSc_sample(RNSc_MLWE out, RNSc_MLWE in1, RNSc_MLWE in2);
    void mlwe_add_RNS_sample(RNS_MLWE out, RNS_MLWE in1, RNS_MLWE in2);
    void mlwe_sub_RNSc_sample(RNSc_MLWE out, RNSc_MLWE in1, RNSc_MLWE in2);
    void mlwe_RNSc_mul_by_xai(RNSc_MLWE out, RNSc_MLWE in, uint64_t a);
    void mlwe_RNSc_mul_by_xai_minus1(RNSc_MLWE out, RNSc_MLWE in, uint64_t a);
    void mlwe_RNS_mul_addto_by_poly(RNS_MLWE out, RNS_MLWE in, const ArithElement *poly);
    void mlwe_RNS_mul_subto_by_poly(RNS_MLWE out, RNS_MLWE in, const ArithElement *poly);
    // The linear combination out = sum_i coeff[i] * in[i] of samples with
    // plaintext coefficients, all in the mul domain over out's ring and rank,
    // `out` distinct from every in[i]. A NULL coeff[i].handle drops the term,
    // and no terms at all leave out = 0.
    void mlwe_RNS_linear_combination(RNS_MLWE out, RNS_MLWE *in, const ArithElement *coeff,
                                     uint64_t n);
    // out[j] = sum_i coeff[j * n_in + i] * in[i] for every j < n_out: several
    // linear combinations of the same samples (a plaintext matrix times a
    // vector of samples), on up to `n_threads` threads (0: the library's
    // limit; see vfhe_threads_for).
    void mlwe_RNS_linear_combinations(RNS_MLWE *out, RNS_MLWE *in, const ArithElement *coeff,
                                      uint64_t n_out, uint64_t n_in, uint64_t n_threads);
    void mlwe_automorphism_RNSc_GHS(RNSc_MLWE out, RNSc_MLWE in, uint64_t gen, RNS_MLWE_KS_Key ksk,
                                    uint64_t lvl);
    void mlwe_scale_RNSc_mlwe(RNSc_MLWE c, uint64_t scale);
    void mlwe_addto_RNSc_sample(RNSc_MLWE out, RNSc_MLWE in);
    RNS_MLWE *mlwe_alloc_RNS_sample_array(uint64_t size, uint64_t N, uint64_t r, uint64_t mask,
                                          RNS_Base base);
    RNS_MLWE *mlwe_alloc_RNS_sample_array2(uint64_t size, RNS_MLWE c);
    void free_RNS_mlwe_array(uint64_t size, RNS_MLWE *v);
    void free_mlwe_RNS_sample(void *p);
    void mlwe_scale_RNS_mlwe_addto(RNS_MLWE out, RNS_MLWE in, uint64_t scale);
    void mlwe_RNS_mul_by_poly(RNS_MLWE out, RNS_MLWE in, const ArithElement *poly);
    void mlwe_RNSc_extract_lwe(uint64_t *out, RNSc_MLWE in, uint64_t idx);
    void mlwe_add_RNSc_polynomial(RNSc_MLWE out, RNSc_MLWE in1, const ArithElement *in2);
    void mlwe_sub_RNSc_polynomial(RNSc_MLWE out, RNSc_MLWE in1, const ArithElement *in2);
    void mlwe_RNS_add_polynomial(RNS_MLWE out, RNS_MLWE in1, const ArithElement *in2);
    void mlwe_RNS_sub_polynomial(RNS_MLWE out, RNS_MLWE in1, const ArithElement *in2);

    // Rank of the extended (not-yet-relinearized) product of two rank-r ciphertexts:
    // r*(r+1)/2 quadratic components plus r linear ones.
    uint64_t mlwe_extended_rank(uint64_t r);
    void mlwe_tensor_product(ArithElement *out, RNS_MLWE in1, RNS_MLWE in2);
    void mlwe_multiply(RNS_MLWE out, RNS_MLWE in1, RNS_MLWE in2, RNS_MLWE_KS_Key ksk);
    // out[i] = in1[i] * in2[i] for every i < n, as mlwe_multiply with the one
    // `ksk` (NULL: extended products), on up to `n_threads` threads (0: the
    // library's limit). The inputs are in the mul domain; an input may appear
    // in several products.
    void mlwe_multiply_batch(RNS_MLWE *out, RNS_MLWE *in1, RNS_MLWE *in2, RNS_MLWE_KS_Key ksk,
                             uint64_t n, uint64_t n_threads);

    RNS_MLWE_Key mlwe_new_RNS_key_from_array(uint64_t *array, uint64_t N, uint64_t r, uint64_t l,
                                             RNS_Base base, double sigma);
    void mlwe_copy_array(RNS_MLWE *out, RNS_MLWE *in, uint64_t size);
    RNS_MLWE *mlwe_create_copy_array(RNS_MLWE *in, uint64_t size);

    // Wrap per-component gadget key arrays (borrowed, not deep-copied) into a
    // key-switch key, deriving the key's ring from its first real component.
    // `log_base` is the gadget the arrays were generated against (0 for the
    // RNS gadget), and each array must hold exactly as many keys as that
    // gadget decomposes an element of the key's ring into.
    RNS_MLWE_KS_Key mlwe_new_RNS_ks_key(RNS_MLWE **s, uint64_t count, uint64_t log_base);
    void free_mlwe_RNS_ks_key(RNS_MLWE_KS_Key key);
    void mlwe_RNSc_GHS_hybrid_keyswitch(RNSc_MLWE out, RNSc_MLWE in, RNS_MLWE_KS_Key ksk,
                                        uint64_t lvl);
    // Several automorphisms of one sample, key-switched with the input's
    // decomposition computed once (hoisting): `mlwe_hoist` copies `in` and
    // decomposes it against the gadget and ring of `ksk`; each
    // mlwe_automorphism_RNSc_GHS_hoisted then costs the key products alone,
    // and agrees with mlwe_automorphism_RNSc_GHS up to the noise. Any key with
    // that gadget, ring and pass-through pattern serves -- automorphism keys
    // for one level generated together all do -- and a key that does not is
    // refused with -1, leaving `out` untouched. The hoisted sample is read
    // only, so threads may share it.
    typedef struct _MLWE_Hoisted *MLWE_Hoisted;
    MLWE_Hoisted mlwe_hoist(RNSc_MLWE in, RNS_MLWE_KS_Key ksk);
    void free_mlwe_hoisted(MLWE_Hoisted h);
    int mlwe_automorphism_RNSc_GHS_hoisted(RNSc_MLWE out, MLWE_Hoisted h, uint64_t gen,
                                           RNS_MLWE_KS_Key ksk, uint64_t lvl);
    // out[i] = Aut_gens[i](h's sample) for every i < n, on up to `n_threads`
    // threads (0: the library's limit). Every key is checked first; -1 if one
    // does not fit, with no output written.
    int mlwe_automorphisms_RNSc_GHS_hoisted(RNSc_MLWE *out, MLWE_Hoisted h, const uint64_t *gens,
                                            RNS_MLWE_KS_Key *ksks, uint64_t n, uint64_t lvl,
                                            uint64_t n_threads);
    // out[i] = Aut_gens[i](in[i]) for every i < n, independently, on up to
    // `n_threads` threads (0: the library's limit). gens[i] == 1 copies, and
    // its key may be NULL.
    void mlwe_automorphism_RNSc_GHS_batch(RNSc_MLWE *out, RNSc_MLWE *in, const uint64_t *gens,
                                          RNS_MLWE_KS_Key *ksks, uint64_t n, uint64_t lvl,
                                          uint64_t n_threads);
    void mlwe_partial_trace(RNSc_MLWE out, RNSc_MLWE in, uint64_t *gens, RNS_MLWE_KS_Key *ksks,
                            uint64_t size, uint64_t lvl);
    void mlwe_trace(RNSc_MLWE out, RNSc_MLWE in, RNS_MLWE_KS_Key *ksks, uint64_t lvl);
    // The samples in `in` live over exactly the primes of `out`'s ring.
    void mlwe_full_packing_keyswitch(RNS_MLWE out, LWE *in, uint64_t size, RNS_MLWE_KS_Key ksk,
                                     uint64_t lvl);
    void mlwe_full_packing_keyswitch_scaled(RNSc_MLWE *vec, uint64_t ell, RNS_MLWE_KS_Key *ksks,
                                            uint64_t lvl);
    void mlwe_round_division(RNSc_MLWE out, ArithRing to);
    // Reduce every component into `to`, a quotient of the sample's ring, in
    // place and in either domain: the value is kept, not divided (compare
    // mlwe_round_division), so this is valid where the plaintext does not
    // depend on the modulus -- CKKS's level drop -- and not for BFV, whose
    // scaling is the modulus.
    void mlwe_mod_reduce(MLWE c, ArithRing to);
    // mlwe_round_division of each of the n distinct samples, in place, each
    // first moved to the canonical domain if it is not there, on up to
    // `n_threads` threads (0: the library's limit).
    void mlwe_round_division_batch(RNSc_MLWE *io, ArithRing to, uint64_t n, uint64_t n_threads);
    // Each of the n distinct samples moved to the mul domain in place (those
    // already there are left alone), on up to `n_threads` threads.
    void mlwe_RNSc_to_RNS_batch(RNS_MLWE *io, uint64_t n, uint64_t n_threads);

    // Gadget decomposition products: decompose `poly` against the gadget
    // `log_base` names -- the RNS one (one digit per prime) when it is 0, the
    // radix one (one digit per prime and power of 2^log_base) otherwise -- and
    // accumulate the products with the matching keys into `out`. `ksk` must
    // hold one key per digit, prime-major, generated against the same gadget.
    void gadget_mul_addto_polynomial(RNS_MLWE out, RNS_MLWE *ksk, const ArithElement *poly,
                                     uint64_t log_base);
    void gadget_mul_subto_polynomial(RNS_MLWE out, RNS_MLWE *ksk, const ArithElement *poly,
                                     uint64_t log_base);

    // How many base-2^log_base digits the radix gadget takes for one residue
    // modulo `prime`: enough to cover every value below it.
    uint64_t gadget_radix_digits(uint64_t prime, uint64_t log_base);

    // The gadget decomposition of one element, kept for several products: `n`
    // digits over the key's ring, in the order its keys are, each lifted and
    // ready to multiply (or, on a ring that is not fully split, canonical).
    typedef struct
    {
        ArithElement *digit;
        uint64_t n;
    } GadgetDigits;

    // Decompose `poly` as gadget_mul_*to_polynomial would against `ksk` and
    // `log_base`, keeping the digits; `ksk` fixes only the key ring and the
    // gadget, so any key array generated for both serves.
    void gadget_decompose(GadgetDigits *out, RNS_MLWE *ksk, const ArithElement *poly,
                          uint64_t log_base);
    void gadget_digits_free(GadgetDigits *digits);
    // out -= sum_i Aut_gen(digit_i) * ksk[i]. The image of a decomposition
    // under an automorphism is a decomposition of the image with the same
    // digit bounds, so this is the gadget product of Aut_gen(poly) -- up to
    // which representative each digit takes, not bit for bit. `gen` is odd and
    // below 2N; 1 is the identity.
    void gadget_mul_subto_automorphism(RNS_MLWE out, RNS_MLWE *ksk, const GadgetDigits *digits,
                                       uint64_t gen);

    // `ell` is the number of gadget keys per component of the MGSW key -- one
    // per prime for the RNS gadget (`log_base` 0), one per prime and digit for
    // the radix one -- and `log_base` is the gadget they were generated
    // against, as in `gadget_mul_addto_polynomial`.
    void mgsw_external_product(RNS_MLWE out, RNS_MLWE *mgsw, RNSc_MLWE in, uint64_t ell,
                               uint64_t special_primes, uint64_t log_base);
    void mgsw_CMUX(RNS_MLWE out, RNSc_MLWE in1, RNSc_MLWE in2, RNS_MLWE *mgsw, uint64_t ell,
                   uint64_t special_primes, uint64_t log_base);
    void mgsw_NCMUX(RNS_MLWE out, RNSc_MLWE in1, RNSc_MLWE in2, RNS_MLWE *mgsw, RNS_MLWE_KS_Key ksk,
                    uint64_t ell, uint64_t special_primes, uint64_t log_base);
    void mgsw_CMUX_to_coeff(RNS_MLWE out, RNSc_MLWE in1, RNSc_MLWE in2, RNS_MLWE *mgsw,
                            uint64_t ell, uint64_t special_primes, uint64_t log_base);
    void mgsw_NCMUX_to_coeff(RNS_MLWE out, RNSc_MLWE in1, RNSc_MLWE in2, RNS_MLWE *mgsw,
                             RNS_MLWE_KS_Key ksk, uint64_t ell, uint64_t special_primes,
                             uint64_t log_base);
    void gp25_RGSW_monomial_mul(RNS_MLWE *p0, uint64_t in_N, RNS_MLWE **e, uint64_t r_prec,
                                RNS_MLWE_KS_Key ksk, uint64_t ell, uint64_t special_primes);
    void gp25_RGSW_monomial_mul_mt(RNS_MLWE *p0, uint64_t in_N, RNS_MLWE **e, uint64_t r_prec,
                                   RNS_MLWE_KS_Key ksk, uint64_t ell, uint64_t special_primes,
                                   uint64_t num_threads);
    void gp25_sub_a_mt(RNS_MLWE *p0, uint64_t in_N, uint64_t *a, RNS_MLWE *s_sign, uint64_t ell,
                       uint64_t special_primes, uint64_t N, uint64_t num_threads);

    // lwe
    LWE_Key lwe_alloc_key(uint64_t n, uint64_t mask, RNS_Base base);
    LWE lwe_alloc_sample(uint64_t n, uint64_t mask, RNS_Base base);
    void free_lwe_sample(LWE c);
    LWE_Key lwe_new_key(uint64_t n, uint64_t mask, RNS_Base base, double sec_sigma,
                        double err_sigma);
    LWE_Key lwe_new_sparse_ternary_key(uint64_t n, uint64_t mask, RNS_Base base, uint64_t h,
                                       double err_sigma);
    void lwe_sample(LWE c, uint64_t *m, LWE_Key key);
    LWE lwe_new_sample(uint64_t *m, LWE_Key key);
    LWE lwe_new_trivial_sample(uint64_t *m, uint64_t n, uint64_t mask, RNS_Base base);
    void lwe_phase(uint64_t *out, LWE c, LWE_Key key);
    void lwe_subto(LWE out, LWE in);
    LWE_KS_Key lwe_new_KS_key(LWE_Key out_key, LWE_Key in_key, uint64_t t, uint64_t base_bit);
    void lwe_keyswitch(LWE out, LWE in, LWE_KS_Key ks_key);

#ifdef __cplusplus
}
#endif
