// Created by GPT-6 on 2026-09-11
// Modified by GPT-6 on 2026-09-12
// Local v197 implementation of the user's exact two-product square proposal.
// Project license: GNU GPL version 3; see ../../../../LICENSE.
// Source lineage and scope: ../../../../SOURCES_AND_REFERENCES.md.
#pragma once
#include <gmp.h>
#include <algorithm>
inline bool useTwoProductSquare(mpz_srcptr x, mpz_srcptr y) {
#ifdef STMS_LEGACY_INTEGER_SQUARE
    return false;
#else
    const auto a = mpz_size(x), b = mpz_size(y);
    return std::min(a, b) >= 8 && std::max(a, b) - std::min(a, b) <= std::min(a, b) / 4;
#endif
}
// Exact signed integers; callers retain their original truncation points.
// All outputs/scratch must be distinct from the inputs and each other.
// Returns re = x*x - y*y and im = x*y. The caller supplies the factor two
// by shifting im by b-1 instead of b; keep those toward-zero shifts unchanged.
inline void
integerSquare(mpz_ptr re, mpz_ptr im, mpz_srcptr x, mpz_srcptr y, mpz_ptr sum, mpz_ptr difference) {
    if (useTwoProductSquare(x, y)) {
        mpz_add(sum, x, y);
        mpz_sub(difference, x, y);
        mpz_mul(re, sum, difference);
    } else {
        mpz_mul(re, x, x);
        mpz_mul(sum, y, y);
        mpz_sub(re, re, sum);
    }
    mpz_mul(im, x, y);
}
