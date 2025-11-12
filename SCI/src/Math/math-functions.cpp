/*
Authors: Deevashwer Rathee
Copyright:
Copyright (c) 2021 Microsoft Research
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "Math/math-functions.h"

#include <complex>
#include <iomanip>

using namespace std;
using namespace sci;

#define KKOT_LIMIT 8
#define SQRT_LOOKUP_SCALE 2

MathFunctions::MathFunctions(int party, IOPack *iopack, OTPack *otpack) {
    this->party = party;
    this->iopack = iopack;
    this->otpack = otpack;
    this->aux = new AuxProtocols(party, iopack, otpack);
    this->equality = new Equality(party, iopack, otpack);
    this->xt = new XTProtocol(party, iopack, otpack);
    this->trunc = new Truncation(party, iopack, otpack);
    this->mult = new LinearOT(party, iopack, otpack);
}

MathFunctions::~MathFunctions() {
    delete aux;
    delete xt;
    delete trunc;
    delete equality;
    delete mult;
}

// A0 \in (1/4, 1)
uint64_t lookup_A0(uint64_t index, int m) {
    uint64_t k = 1ULL << m;
    double p = 1 + (double(index) / double(k));
    double A1 = 1.0 / (p * (p + 1.0 / double(k)));
    int32_t scale = m + 3;
    uint64_t mask = (1ULL << scale) - 1;
    uint64_t val = uint64_t(A1 * (1ULL << scale)) & mask;
    return val;
}

// A1 \in (1/2, 1)
uint64_t lookup_A1(uint64_t index, int m) {
    uint64_t k = 1ULL << m;
    double p = 1 + (double(index) / double(k));
    double z = (p * (p + (1.0 / double(k))));
    double A1 = ((1.0 / double(k * 2)) + sqrt(z)) / z;
    int32_t scale = 2 * m + 2;
    uint64_t mask = (1ULL << scale) - 1;
    uint64_t val = uint64_t(A1 * (1ULL << scale)) & mask;
    return val;
}

void MathFunctions::reciprocal_approximation(int32_t dim, int32_t m,
                                             uint64_t *dn, uint64_t *out,
                                             int32_t bw_dn, int32_t bw_out,
                                             int32_t s_dn, int32_t s_out) {
    assert(bw_out == m + s_dn + 4);
    assert(s_out == m + s_dn + 4);

    uint64_t s_dn_mask = (1ULL << s_dn) - 1;
    uint64_t m_mask = (1ULL << m) - 1;
    uint64_t s_min_m_mask = (1ULL << (s_dn - m)) - 1;

    uint64_t *tmp_1 = new uint64_t[dim];
    uint64_t *tmp_2 = new uint64_t[dim];

    for (int i = 0; i < dim; i++) {
        tmp_1[i] = dn[i] & s_dn_mask;
    }
    trunc->truncate_and_reduce(dim, tmp_1, tmp_2, s_dn - m, s_dn);

    int M = (1ULL << m);
    uint64_t c0_mask = (1ULL << (m + 4)) - 1;
    uint64_t c1_mask = (1ULL << (2 * m + 3)) - 1;
    uint64_t *c0 = new uint64_t[dim];
    uint64_t *c1 = new uint64_t[dim];
    if (party == ALICE) {
        uint64_t **spec;
        spec = new uint64_t *[dim];
        PRG128 prg;
        prg.random_data(c0, dim * sizeof(uint64_t));
        prg.random_data(c1, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            spec[i] = new uint64_t[M];
            c0[i] &= c0_mask;
            c1[i] &= c1_mask;
            for (int j = 0; j < M; j++) {
                int idx = (tmp_2[i] + j) & m_mask;
                spec[i][j] = (lookup_A0(idx, m) - c0[i]) & c0_mask;
                spec[i][j] <<= (2 * m + 3);
                spec[i][j] |= (lookup_A1(idx, m) - c1[i]) & c1_mask;
            }
        }
        aux->lookup_table<uint64_t>(spec, nullptr, nullptr, dim, m, 3 * m + 7);

        for (int i = 0; i < dim; i++)
            delete[] spec[i];
        delete[] spec;
    } else {
        aux->lookup_table<uint64_t>(nullptr, tmp_2, c1, dim, m, 3 * m + 7);

        for (int i = 0; i < dim; i++) {
            c0[i] = (c1[i] >> (2 * m + 3)) & c0_mask;
            c1[i] = c1[i] & c1_mask;
        }
    }
    for (int i = 0; i < dim; i++) {
        tmp_1[i] = dn[i] & s_min_m_mask;
    }
    uint8_t *zero_shares = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        zero_shares[i] = 0;
    }

    // Unsigned mult
    mult->hadamard_product(dim, c0, tmp_1, tmp_2, m + 4, s_dn - m, s_dn + 4,
                           false, false, MultMode::None, zero_shares, nullptr);

    xt->z_extend(dim, tmp_2, tmp_1, s_dn + 4, s_dn + m + 4, zero_shares);

    uint64_t out_mask = (1ULL << (s_dn + m + 4)) - 1;
    uint64_t scale_up = (1ULL << (s_dn - m + 1));
    for (int i = 0; i < dim; i++) {
        out[i] = ((c1[i] * scale_up) - tmp_1[i]) & out_mask;
    }

    delete[] tmp_1;
    delete[] tmp_2;
    delete[] c0;
    delete[] c1;
    delete[] zero_shares;
}

void MathFunctions::div(int32_t dim, uint64_t *nm, uint64_t *dn, uint64_t *out,
                        int32_t bw_nm, int32_t bw_dn, int32_t bw_out,
                        int32_t s_nm, int32_t s_dn, int32_t s_out,
                        bool signed_nm, bool compute_msnzb) {
    assert(s_out <= s_dn);

    // out_precision = iters * (2*m + 2)
    int32_t m, iters;
    m = (s_out <= 18
             ? ceil((s_out - 2) / 2.0)
             : ceil((ceil(s_out / 2.0) - 2) / 2.0));
    iters = (s_out <= 18 ? 0 : 1);

    int32_t s_tmp_dn;
    int32_t bw_adjust;
    int32_t s_adjust;
    uint64_t *tmp_dn;
    uint64_t *adjust;
    if (compute_msnzb) {
        s_tmp_dn = bw_dn - 1;
        bw_adjust = bw_dn + 1;
        s_adjust = bw_dn - 1 - s_dn;
        uint64_t mask_adjust = (bw_adjust == 64 ? -1 : ((1ULL << bw_adjust) - 1));
        // MSB is always 0, thus, not including it
        uint8_t *msnzb_vector_bool = new uint8_t[dim * bw_dn];
        uint64_t *msnzb_vector = new uint64_t[dim * bw_dn];
        aux->msnzb_one_hot(dn, msnzb_vector_bool, bw_dn, dim);
        aux->B2A(msnzb_vector_bool, msnzb_vector, dim * bw_dn, bw_adjust);
        // adjust: bw = bw_dn, scale = bw_dn - 1 - s_dn
        adjust = new uint64_t[dim];
        for (int i = 0; i < dim; i++) {
            adjust[i] = 0;
            for (int j = 0; j < bw_dn; j++) {
                adjust[i] += (1ULL << (bw_dn - 1 - j)) * msnzb_vector[i * bw_dn + j];
            }
            adjust[i] &= mask_adjust;
        }
        // tmp_dn: bw = bw_dn, scale = bw_dn - 1
        tmp_dn = new uint64_t[dim];
        mult->hadamard_product(dim, dn, adjust, tmp_dn, bw_dn, bw_dn + 1, bw_dn + 1,
                               false, false, MultMode::None);

        delete[] msnzb_vector_bool;
        delete[] msnzb_vector;
    } else {
        // tmp_dn: bw = s_dn + 1, scale = s_dn
        s_tmp_dn = s_dn;
        tmp_dn = dn;
    }

    uint64_t *tmp_1 = new uint64_t[dim];
    uint64_t *tmp_2 = new uint64_t[dim];
    // tmp_1: bw = s_tmp_dn + m + 4, scale = s_tmp_dn + m + 3
    reciprocal_approximation(dim, m, tmp_dn, tmp_1, bw_dn, s_tmp_dn + m + 4,
                             s_tmp_dn, s_tmp_dn + m + 4);

    uint64_t *w0 = new uint64_t[dim];
    // w0: bw = s_out + 1, scale = s_out
    trunc->truncate_and_reduce(dim, tmp_1, w0, s_tmp_dn + m + 3 - s_out,
                               s_tmp_dn + m + 4);

    uint8_t *msb_nm = new uint8_t[dim];
    aux->MSB(nm, msb_nm, dim, bw_nm);
    uint8_t *zero_shares = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        zero_shares[i] = 0;
    }

    // a0: bw = bw_out, scale = s_out
    uint64_t *a0 = new uint64_t[dim];
    // Mixed mult with w0 unsigned
    mult->hadamard_product(dim, nm, w0, tmp_1, bw_nm, s_out + 1, s_out + bw_nm,
                           signed_nm, false, MultMode::None, msb_nm, zero_shares);
    trunc->truncate_and_reduce(dim, tmp_1, tmp_2, s_nm, s_out + bw_nm);
    if ((bw_nm - s_nm) >= (bw_out - s_out)) {
        aux->reduce(dim, tmp_2, a0, bw_nm - s_nm + s_out, bw_out);
    } else {
        if (signed_nm) {
            xt->s_extend(dim, tmp_2, a0, s_out + bw_nm - s_nm, bw_out, msb_nm);
        } else {
            xt->z_extend(dim, tmp_2, a0, s_out + bw_nm - s_nm, bw_out, nullptr);
        }
    }

    if (compute_msnzb) {
        int32_t bw_tmp1 =
                (bw_out + s_adjust < bw_adjust ? bw_adjust : bw_out + s_adjust);
        // tmp_1: bw = bw_tmp1, scale = s_out + s_adjust
        mult->hadamard_product(dim, a0, adjust, tmp_1, bw_out, bw_adjust, bw_tmp1,
                               signed_nm, false, MultMode::None,
                               (signed_nm ? msb_nm : nullptr), zero_shares);
        // a0: bw = bw_out, scale = s_out
        trunc->truncate_and_reduce(dim, tmp_1, a0, s_adjust, bw_out + s_adjust);
    }

    // tmp_1: bw = s_tmp_dn + 2, scale = s_tmp_dn
    uint64_t s_plus_2_mask = (1ULL << (s_tmp_dn + 2)) - 1;
    for (int i = 0; i < dim; i++) {
        tmp_1[i] = tmp_dn[i] & s_plus_2_mask;
    }

    if (iters > 0) {
        // d0: bw = s_out + 2, scale = s_out
        uint64_t *d0 = new uint64_t[dim];
        mult->hadamard_product(dim, w0, tmp_1, tmp_2, s_out + 1, s_tmp_dn + 2,
                               s_out + s_tmp_dn + 2, false, false, MultMode::None,
                               zero_shares, zero_shares);
        trunc->truncate_and_reduce(dim, tmp_2, d0, s_tmp_dn, s_out + s_tmp_dn + 2);

        // e0: bw = s_out + 2, scale = s_out
        // e0 = 1 - d0
        uint64_t *e0 = new uint64_t[dim];
        for (int i = 0; i < dim; i++) {
            e0[i] = (party == ALICE ? (1ULL << (s_out)) : 0) - d0[i];
        }

        uint64_t e_mask = (1ULL << (s_out + 2)) - 1;
        uint64_t *a_curr = new uint64_t[dim];
        uint64_t *e_curr = new uint64_t[dim];
        uint64_t *a_prev = a0;
        uint64_t *e_prev = e0;
        for (int i = 0; i < iters - 1; i++) {
            // tmp_1: bw = 2*s_out+2, scale: 2*s_out
            mult->hadamard_product(dim, e_prev, e_prev, tmp_1, s_out + 2, s_out + 2,
                                   2 * s_out + 2, true, true, MultMode::None,
                                   zero_shares, zero_shares);
            // e_curr: bw = s_out + 2, scale: s_out
            trunc->truncate_and_reduce(dim, tmp_1, e_curr, s_out, 2 * s_out + 2);
            // e_prev = 1 + e_prev
            for (int j = 0; j < dim; j++) {
                e_prev[j] =
                        ((party == ALICE ? (1ULL << (s_out)) : 0) + e_prev[j]) & e_mask;
            }
            // tmp_1: bw = bw_out + s_out, scale: 2*s_out
            mult->hadamard_product(dim, a_prev, e_prev, tmp_1, bw_out, s_out + 2,
                                   bw_out + s_out, signed_nm, true, MultMode::None,
                                   (signed_nm ? msb_nm : nullptr), zero_shares);
            // a_curr: bw = bw_out, scale: s_out
            trunc->truncate_and_reduce(dim, tmp_1, a_curr, s_out, bw_out + s_out);
            memcpy(a_prev, a_curr, dim * sizeof(uint64_t));
            memcpy(e_prev, e_curr, dim * sizeof(uint64_t));
        }
        // e_prev = 1 + e_prev
        for (int j = 0; j < dim; j++) {
            e_prev[j] =
                    ((party == ALICE ? (1ULL << (s_out)) : 0) + e_prev[j]) & e_mask;
        }
        // out: bw = bw_out, scale: s_out
        // Mixed mult with e_prev unsigned
        mult->hadamard_product(dim, a_prev, e_prev, tmp_1, bw_out, s_out + 2,
                               bw_out + s_out, signed_nm, false, MultMode::None,
                               (signed_nm ? msb_nm : nullptr), zero_shares);
        trunc->truncate_and_reduce(dim, tmp_1, out, s_out, bw_out + s_out);

        delete[] d0;
        delete[] e0;
        delete[] a_curr;
        delete[] e_curr;
    } else {
        memcpy(out, a0, dim * sizeof(uint64_t));
    }

    delete[] tmp_1;
    delete[] tmp_2;
    delete[] w0;
    delete[] a0;
    delete[] msb_nm;
    delete[] zero_shares;
    if (compute_msnzb) {
        delete[] tmp_dn;
        delete[] adjust;
    }
}

uint64_t lookup_neg_exp(uint64_t val_in, int32_t s_in, int32_t s_out) {
    if (s_in < 0) {
        s_in *= -1;
        val_in *= (1ULL << (s_in));
        s_in = 0;
    }
    uint64_t res_val =
            exp(-1.0 * (val_in / double(1ULL << s_in))) * (1ULL << s_out);
    return res_val;
}

void MathFunctions::lookup_table_exp(int32_t dim, uint64_t *x, uint64_t *y,
                                     int32_t bw_x, int32_t bw_y, int32_t s_x,
                                     int32_t s_y) {
    assert(bw_y >= (s_y + 2));

    int LUT_size = KKOT_LIMIT;

    uint64_t bw_x_mask = (bw_x == 64 ? -1 : (1ULL << bw_x) - 1);
    uint64_t LUT_out_mask = ((s_y + 2) == 64 ? -1 : (1ULL << (s_y + 2)) - 1);

    uint64_t *tmp_1 = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        tmp_1[i] = (-1 * x[i]) & bw_x_mask;
    }
    int digit_size = LUT_size;
    int num_digits = ceil(double(bw_x) / digit_size);
    int last_digit_size = bw_x - (num_digits - 1) * digit_size;
    uint64_t *x_digits = new uint64_t[num_digits * dim];

    aux->digit_decomposition_sci(dim, tmp_1, x_digits, bw_x, digit_size);

    uint64_t digit_mask = (digit_size == 64 ? -1 : (1ULL << digit_size) - 1);
    uint64_t last_digit_mask =
            (last_digit_size == 64 ? -1 : (1ULL << last_digit_size) - 1);
    int N = (1ULL << digit_size);
    int last_N = (1ULL << last_digit_size);
    int N_digits = (digit_size == last_digit_size ? num_digits : num_digits - 1);
    uint64_t *digits_exp = new uint64_t[num_digits * dim];
    if (party == ALICE) {
        uint64_t **spec;
        spec = new uint64_t *[num_digits * dim];
        PRG128 prg;
        prg.random_data(digits_exp, num_digits * dim * sizeof(uint64_t));
        for (int i = 0; i < N_digits * dim; i++) {
            int digit_idx = i / dim;
            spec[i] = new uint64_t[N];
            digits_exp[i] &= LUT_out_mask;
            for (int j = 0; j < N; j++) {
                int idx = (x_digits[i] + j) & digit_mask;
                spec[i][j] = (lookup_neg_exp(idx, s_x - digit_size * digit_idx, s_y) -
                              digits_exp[i]) & LUT_out_mask;
            }
        }
        aux->lookup_table<uint64_t>(spec, nullptr, nullptr, N_digits * dim,
                                    digit_size, s_y + 2);

        if (digit_size != last_digit_size) {
            int offset = N_digits * dim;
            int digit_idx = N_digits;
            for (int i = offset; i < num_digits * dim; i++) {
                spec[i] = new uint64_t[last_N];
                digits_exp[i] &= LUT_out_mask;
                for (int j = 0; j < last_N; j++) {
                    int idx = (x_digits[i] + j) & last_digit_mask;
                    spec[i][j] = (lookup_neg_exp(idx, s_x - digit_size * digit_idx, s_y) -
                                  digits_exp[i]) &
                                 LUT_out_mask;
                }
            }
            aux->lookup_table<uint64_t>(spec + offset, nullptr, nullptr, dim,
                                        last_digit_size, s_y + 2);
        }

        for (int i = 0; i < num_digits * dim; i++)
            delete[] spec[i];
        delete[] spec;
    } else {
        aux->lookup_table<uint64_t>(nullptr, x_digits, digits_exp, N_digits * dim,
                                    digit_size, s_y + 2);
        if (digit_size != last_digit_size) {
            int offset = N_digits * dim;
            aux->lookup_table<uint64_t>(nullptr, x_digits + offset,
                                        digits_exp + offset, dim, last_digit_size,
                                        s_y + 2);
        }
        for (int i = 0; i < num_digits * dim; i++) {
            digits_exp[i] &= LUT_out_mask;
        }
    }

    uint8_t *zero_shares = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        zero_shares[i] = 0;
    }
    for (int i = 1; i < num_digits; i *= 2) {
        for (int j = 0; j < num_digits and j + i < num_digits; j += 2 * i) {
            mult->hadamard_product(dim, digits_exp + j * dim,
                                   digits_exp + (j + i) * dim, digits_exp + j * dim,
                                   s_y + 2, s_y + 2, 2 * s_y + 2, false, false,
                                   MultMode::None, zero_shares, zero_shares);
            trunc->truncate_and_reduce(dim, digits_exp + j * dim,
                                       digits_exp + j * dim, s_y, 2 * s_y + 2);
        }
    }
    xt->z_extend(dim, digits_exp, y, s_y + 2, bw_y, zero_shares);

    delete[] x_digits;
    delete[] tmp_1;
    delete[] digits_exp;
    delete[] zero_shares;
}

void MathFunctions::sigmoid(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                            int32_t bw_y, int32_t s_x, int32_t s_y) {
    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
    uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
    uint64_t mask_exp_out = ((s_y + 2) == 64 ? -1 : ((1ULL << (s_y + 2)) - 1));
    uint64_t mask_den = ((s_y + 2) == 64 ? -1 : ((1ULL << (s_y + 2)) - 1));
    uint8_t *zero_shares = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        zero_shares[i] = 0;
    }

    uint8_t *msb_x = new uint8_t[dim];
    aux->MSB(x, msb_x, dim, bw_x);

    // neg_x = -x + msb_x * (2x) (neg_x is always negative)
    uint64_t *tmp_1 = new uint64_t[dim];
    uint64_t *tmp_2 = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        tmp_1[i] = (-1 * x[i]) & mask_x;
        tmp_2[i] = (2 * x[i]) & mask_x;
    }
    uint64_t *neg_x = new uint64_t[dim];
    aux->multiplexer(msb_x, tmp_2, neg_x, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        neg_x[i] = (neg_x[i] + tmp_1[i]) & mask_x;
    }

    // den = tmp_1 = 1 + exp_neg_x
    uint64_t *exp_neg_x = new uint64_t[dim];
    lookup_table_exp(dim, neg_x, exp_neg_x, bw_x, s_y + 2, s_x, s_y);
    for (int i = 0; i < dim; i++) {
        tmp_1[i] =
                ((party == ALICE ? (1ULL << s_y) : 0) + exp_neg_x[i]) & mask_exp_out;
    }
    // den can't be 2^{s_y+1}, so 1 is subtracted if msb_den is 1
    uint8_t *msb_den = new uint8_t[dim];
    aux->MSB(tmp_1, msb_den, dim, s_y + 2);
    aux->B2A(msb_den, tmp_2, dim, s_y + 2);
    // den = tmp_1 = den - msb_den
    // tmp_2 = 1 (all 1 vector)
    for (int i = 0; i < dim; i++) {
        tmp_1[i] = (tmp_1[i] - tmp_2[i]) & mask_den;
        tmp_2[i] = (party == ALICE ? 1 : 0);
    }
    uint64_t *sig_neg_x = new uint64_t[dim];
    // sig_neg_x = 1/(1 + exp_neg_x)
    div(dim, tmp_2, tmp_1, sig_neg_x, 2, s_y + 2, s_y + 2, 0, s_y, s_y, true,
        false);

    // tmp_2 = num = 1 + msb_x * (exp_neg_x - 1)
    for (int i = 0; i < dim; i++) {
        tmp_1[i] = (exp_neg_x[i] - (party == ALICE ? 1ULL << s_y : 0)) & mask_den;
    }
    aux->multiplexer(msb_x, tmp_1, tmp_2, dim, s_y + 2, s_y + 2);
    for (int i = 0; i < dim; i++) {
        tmp_2[i] = (tmp_2[i] + (party == ALICE ? 1ULL << s_y : 0)) & mask_den;
    }

    mult->hadamard_product(dim, tmp_2, sig_neg_x, tmp_1, s_y + 2, s_y + 2,
                           2 * s_y + 2, false, false, MultMode::None, zero_shares,
                           zero_shares);
    trunc->truncate_and_reduce(dim, tmp_1, tmp_2, s_y, 2 * s_y + 2);

    if (bw_y <= (s_y + 2)) {
        for (int i = 0; i < dim; i++) {
            y[i] = tmp_2[i] & mask_y;
        }
    } else {
        xt->z_extend(dim, tmp_2, y, s_y + 2, bw_y, zero_shares);
    }

    delete[] msb_x;
    delete[] tmp_1;
    delete[] tmp_2;
    delete[] neg_x;
    delete[] exp_neg_x;
    delete[] msb_den;
    delete[] sig_neg_x;
}

void MathFunctions::tanh(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                         int32_t bw_y, int32_t s_x, int32_t s_y) {
    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
    uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));
    uint64_t mask_exp_out = ((s_y + 3) == 64 ? -1 : ((1ULL << (s_y + 3)) - 1));
    uint64_t mask_den = ((s_y + 2) == 64 ? -1 : ((1ULL << (s_y + 2)) - 1));

    uint8_t *msb_x = new uint8_t[dim];
    aux->MSB(x, msb_x, dim, bw_x);

    // neg_x = -x + msb_x * (2x) (neg_x is always negative)
    uint64_t *tmp_1 = new uint64_t[dim];
    uint64_t *tmp_2 = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        tmp_1[i] = (-1 * x[i]) & mask_x;
        tmp_2[i] = (2 * x[i]) & mask_x;
    }
    uint64_t *neg_x = new uint64_t[dim];
    aux->multiplexer(msb_x, tmp_2, neg_x, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        neg_x[i] = (neg_x[i] + tmp_1[i]) & mask_x;
    }

    uint64_t *exp_neg_2x = new uint64_t[dim];
    // scale of neg_x is reduced by 1 instead of mulitplying it with 2 to get
    // neg_2x
    lookup_table_exp(dim, neg_x, exp_neg_2x, bw_x, s_y + 2, s_x - 1, s_y);
    // den = tmp_1 = 1 + exp_neg_2x
    for (int i = 0; i < dim; i++) {
        tmp_1[i] =
                ((party == ALICE ? (1ULL << s_y) : 0) + exp_neg_2x[i]) & mask_exp_out;
    }
    // den can't be 2^{s_y+1}, so 1 is subtracted if msb_den is 1
    uint8_t *msb_den = new uint8_t[dim];
    aux->MSB(tmp_1, msb_den, dim, s_y + 2);
    aux->B2A(msb_den, tmp_2, dim, s_y + 2);
    // den = tmp_1 = den - msb_den
    // num = tmp_2 = 1 - exp_neg_2x
    for (int i = 0; i < dim; i++) {
        tmp_1[i] = (tmp_1[i] - tmp_2[i]) & mask_den;
        tmp_2[i] =
                ((party == ALICE ? (1ULL << s_y) : 0) - exp_neg_2x[i]) & mask_den;
    }
    uint64_t *tanh_neg_x = new uint64_t[dim];
    // tanh_neg_x = (1 - exp_neg_2x)/(1 + exp_neg_2x)
    div(dim, tmp_2, tmp_1, tanh_neg_x, s_y + 2, s_y + 2, s_y + 2, s_y, s_y, s_y,
        true, false);

    // tmp_2 = tanh_neg_x + msb_x * (-2 * tanh_neg_x)
    // tmp_1 = -2 * tanh_neg_x
    for (int i = 0; i < dim; i++) {
        tmp_1[i] = (-2 * tanh_neg_x[i]) & mask_exp_out;
    }
    aux->multiplexer(msb_x, tmp_1, tmp_2, dim, s_y + 2, s_y + 2);
    for (int i = 0; i < dim; i++) {
        tmp_2[i] = (tmp_2[i] + tanh_neg_x[i]) & mask_exp_out;
    }

    if (bw_y <= (s_y + 2)) {
        for (int i = 0; i < dim; i++) {
            y[i] = tmp_2[i] & mask_y;
        }
    } else {
        xt->s_extend(dim, tmp_2, y, s_y + 2, bw_y, msb_x);
    }

    delete[] msb_x;
    delete[] tmp_1;
    delete[] tmp_2;
    delete[] neg_x;
    delete[] exp_neg_2x;
    delete[] msb_den;
    delete[] tanh_neg_x;
}

uint64_t lookup_sqrt(int32_t index, int32_t m, int32_t exp_parity) {
    int32_t k = 1 << m;
    double u = (1.0 + (double(index) / double(k))) * (1 << exp_parity);
    double Y = 1.0 / sqrt(u);
    int32_t scale = m + SQRT_LOOKUP_SCALE;
    uint64_t val = (Y * (1ULL << scale));
    return val;
}

void MathFunctions::sqrt(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                         int32_t bw_y, int32_t s_x, int32_t s_y, bool inverse) {
    int32_t m, iters;
    if (s_y <= 14) {
        m = ceil(s_y / 2.0);
        iters = 1;
    } else {
        m = ceil((ceil(s_y / 2.0)) / 2.0);
        iters = 2;
    }
    assert(m <= KKOT_LIMIT);
    int32_t bw_adjust = bw_x - 1;
    uint64_t mask_adjust = (bw_adjust == 64 ? -1 : ((1ULL << bw_adjust) - 1));
    // MSB is always 0, thus, not including it
    uint8_t *msnzb_vector_bool = new uint8_t[dim * (bw_x - 1)];
    uint64_t *msnzb_vector = new uint64_t[dim * (bw_x - 1)];

    auto start = clock_start();
    uint64_t comm = iopack->get_comm();

    aux->msnzb_one_hot(x, msnzb_vector_bool, bw_x - 1, dim);

    comm = iopack->get_comm() - comm;
    long long t = time_from(start);

    cout << "MSNZB Time\t" << t / (1000.0) << " ms" << endl;
    cout << "MSNZB Bytes Sent\t" << comm << " bytes" << endl;


    // uint64_t *msnzb_index = new uint64_t[dim];
    // aux->msnzb_one_hot_index(x, msnzb_vector_bool, msnzb_index, bw_x - 1, dim);
    // cout << endl;
    // for (int i = 0; i < bw_x - 1; i++) {
    //     cout << (uint32_t)msnzb_vector_bool[i] << " ";
    // }
    // cout << endl;

    start = clock_start();
    comm = iopack->get_comm();

    aux->B2A(msnzb_vector_bool, msnzb_vector, dim * (bw_x - 1), bw_x - 1);
    uint64_t *adjust = new uint64_t[dim];
    uint8_t *exp_parity = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        adjust[i] = 0;
        exp_parity[i] = 0;
        for (int j = 0; j < (bw_x - 1); j++) {
            adjust[i] += (1ULL << (bw_x - 2 - j)) * msnzb_vector[i * (bw_x - 1) + j];
            if (((j - s_x) & 1)) {
                exp_parity[i] ^= msnzb_vector_bool[i * (bw_x - 1) + j];
            }
        }
        adjust[i] &= mask_adjust;
    }
    // adjusted_x: bw = bw_x - 1, scale = bw_x - 2, adjusted_x is the value in [1, 2).
    uint64_t *adjusted_x = new uint64_t[dim];
    mult->hadamard_product(dim, x, adjust, adjusted_x, bw_x - 1, bw_x - 1,
                           bw_x - 1, false, false, MultMode::None);
    // Top m bits of adjusted_x excluding MSB, which is always 1
    // adjusted_x_m: bw = m, scale = m
    uint64_t *adjusted_x_m = new uint64_t[dim];
    trunc->truncate_and_reduce(dim, adjusted_x, adjusted_x_m, bw_x - m - 2,
                               bw_x - 2);

    comm = iopack->get_comm() - comm;
    t = time_from(start);

    cout << "adjust x Time\t" << t / (1000.0) << " ms" << endl;
    cout << "adjust x Bytes Sent\t" << comm << " bytes" << endl;

    // m + 1 bits will be input to the lookup table
    int32_t M = (1LL << (m + 1));
    uint64_t Y_mask = (1ULL << (m + SQRT_LOOKUP_SCALE + 1)) - 1;
    uint64_t m_mask = (1ULL << m) - 1;
    // Y: bw = m + SQRT_LOOKUP_SCALE + 1, scale = m + SQRT_LOOKUP_SCALE
    uint64_t *Y = new uint64_t[dim];
    if (party == ALICE) {
        uint64_t **spec;
        spec = new uint64_t *[dim];
        PRG128 prg;
        prg.random_data(Y, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            spec[i] = new uint64_t[M];
            Y[i] &= Y_mask;
            for (int j = 0; j < M; j++) {
                // j = exp_parity || (adjusted_x_m) (LSB -> MSB)
                int32_t idx = (adjusted_x_m[i] + (j >> 1)) & m_mask;
                int32_t exp_parity_val = (exp_parity[i] ^ (j & 1));
                spec[i][j] = (lookup_sqrt(idx, m, exp_parity_val) - Y[i]) & Y_mask;
            }
        }
        aux->lookup_table<uint64_t>(spec, nullptr, nullptr, dim, m + 1,
                                    m + SQRT_LOOKUP_SCALE + 1);

        for (int i = 0; i < dim; i++)
            delete[] spec[i];
        delete[] spec;
    } else {
        // lut_in = exp_parity || adjusted_x_m
        uint64_t *lut_in = new uint64_t[dim];
        for (int i = 0; i < dim; i++) {
            lut_in[i] = ((adjusted_x_m[i] & m_mask) << 1) | (exp_parity[i] & 1);
        }
        aux->lookup_table<uint64_t>(nullptr, lut_in, Y, dim, m + 1,
                                    m + SQRT_LOOKUP_SCALE + 1);

        delete[] lut_in;
    }
    // X = (exp_parity ? 2 * adjusted_x : adjusted_x); scale = bw_x - 2
    // X: bw = bw_x
    uint64_t *X = new uint64_t[dim];
    uint64_t *tmp_1 = new uint64_t[dim];
    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
    xt->z_extend(dim, adjusted_x, X, bw_x - 1, bw_x);
    aux->multiplexer(exp_parity, X, tmp_1, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        X[i] = (X[i] + tmp_1[i]) & mask_x;
    }

    uint64_t *x_prev = new uint64_t[dim];
    if (inverse) {
        // x_prev : bw = m + SQRT_LOOKUP_SCALE + 1, scale = m + SQRT_LOOKUP_SCALE
        // x_prev \approx 0.5 < 1/sqrt(X) < 1
        memcpy(x_prev, Y, dim * sizeof(uint64_t));
    } else {
        // x_prev : bw = s_y + 1, scale = s_y
        // x_prev \approx 1 < sqrt(X) < 2
        mult->hadamard_product(dim, X, Y, tmp_1, bw_x, m + SQRT_LOOKUP_SCALE + 1,
                               bw_x - 1 + m + SQRT_LOOKUP_SCALE, false, false,
                               MultMode::None);
        assert((bw_x - 2 + m + SQRT_LOOKUP_SCALE) >= s_y);
        trunc->truncate_and_reduce(dim, tmp_1, x_prev,
                                   bw_x - 2 + m + SQRT_LOOKUP_SCALE - s_y,
                                   bw_x - 1 + m + SQRT_LOOKUP_SCALE);
    }
    // b_prev: bw = s_y + 2, scale = s_y
    uint64_t *b_prev = new uint64_t[dim];
    assert((bw_x - 2) >= s_y);
    trunc->truncate_and_reduce(dim, X, b_prev, bw_x - 2 - s_y, bw_x);
    // Y_prev: bw = m + SQRT_LOOKUP_SCALE + 1, scale = m + SQRT_LOOKUP_SCALE
    uint64_t *Y_prev = new uint64_t[dim];
    memcpy(Y_prev, Y, dim * sizeof(uint64_t));

    uint64_t b_mask = (1ULL << (s_y + 2)) - 1;
    uint64_t *x_curr = new uint64_t[dim];
    uint64_t *b_curr = new uint64_t[dim];
    uint64_t *Y_curr = new uint64_t[dim];
    uint64_t *Y_sq = new uint64_t[dim];

    start = clock_start();
    comm = iopack->get_comm();

    for (int i = 0; i < iters; i++) {
        if (i == 0) {
            // Y_sq: bw = 2m + 2SQRT_LOOKUP_SCALE + 1, scale = 2m + 2SQRT_LOOKUP_SCALE
            mult->hadamard_product(
                dim, Y_prev, Y_prev, Y_sq, m + SQRT_LOOKUP_SCALE + 1,
                m + SQRT_LOOKUP_SCALE + 1, 2 * m + 2 * SQRT_LOOKUP_SCALE + 1, false,
                false, MultMode::None);
            // tmp_1: bw = 2m+2SQRT_LOOKUP_SCALE+s_y+2, scale =
            // 2m+2SQRT_LOOKUP_SCALE+s_y
            mult->hadamard_product(dim, Y_sq, b_prev, tmp_1,
                                   2 * m + 2 * SQRT_LOOKUP_SCALE + 1, s_y + 2,
                                   2 * m + 2 * SQRT_LOOKUP_SCALE + s_y + 2, false,
                                   false, MultMode::None);
            // b_curr: bw = s_y + 2, scale = s_y
            trunc->truncate_and_reduce(dim, tmp_1, b_curr,
                                       2 * m + 2 * SQRT_LOOKUP_SCALE,
                                       2 * m + 2 * SQRT_LOOKUP_SCALE + s_y + 2);
        } else {
            // tmp_1: bw = 2*s_y + 3, scale = 2*s_y + 2
            mult->hadamard_product(dim, Y_prev, Y_prev, tmp_1, s_y + 2, s_y + 2,
                                   2 * s_y + 3, false, false, MultMode::None);
            // Y_sq: bw = s_y + 1, scale = s_y
            trunc->truncate_and_reduce(dim, tmp_1, Y_sq, s_y + 2, 2 * s_y + 3);
            // tmp_1: bw = 2s_y + 2, scale = 2s_y
            mult->hadamard_product(dim, Y_sq, b_prev, tmp_1, s_y + 1, s_y + 2,
                                   2 * s_y + 2, false, false, MultMode::None);
            // b_curr: bw = s_y + 2, scale = s_y
            trunc->truncate_and_reduce(dim, tmp_1, b_curr, s_y, 2 * s_y + 2);
        }
        for (int j = 0; j < dim; j++) {
            // Y_curr: bw = s_y + 2, scale = s_y + 1
            // Y_curr = (3 - b_curr)/2
            Y_curr[j] = ((party == ALICE ? (3ULL << s_y) : 0) - b_curr[j]) & b_mask;
        }
        if (inverse && (i == 0)) {
            // tmp_1: bw = s_y+m+SQRT_LOOKUP_SCALE+2, scale =
            // s_y+m+SQRT_LOOKUP_SCALE+1
            mult->hadamard_product(
                dim, x_prev, Y_curr, tmp_1, m + SQRT_LOOKUP_SCALE + 1, s_y + 2,
                s_y + m + SQRT_LOOKUP_SCALE + 2, false, false, MultMode::None);
            // x_curr: bw = s_y + 1, scale = s_y
            trunc->truncate_and_reduce(dim, tmp_1, x_curr, m + SQRT_LOOKUP_SCALE + 1,
                                       s_y + m + SQRT_LOOKUP_SCALE + 2);
        } else {
            // tmp_1: bw = 2*s_y + 2, scale = 2s_y + 1
            mult->hadamard_product(dim, x_prev, Y_curr, tmp_1, s_y + 1, s_y + 2,
                                   2 * s_y + 2, false, false, MultMode::None);
            // x_curr: bw = s_y + 1, scale = s_y
            trunc->truncate_and_reduce(dim, tmp_1, x_curr, s_y + 1, 2 * s_y + 2);
        }
        memcpy(x_prev, x_curr, dim * sizeof(uint64_t));
        memcpy(b_prev, b_curr, dim * sizeof(uint64_t));
        memcpy(Y_prev, Y_curr, dim * sizeof(uint64_t));
    }


    comm = iopack->get_comm() - comm;
    t = time_from(start);

    cout << "newton Time\t" << t / (1000.0) << " ms" << endl;
    cout << "newton Bytes Sent\t" << comm << " bytes" << endl;

    int32_t bw_sqrt_adjust = bw_x / 2;
    uint64_t mask_sqrt_adjust =
            (bw_sqrt_adjust == 64 ? -1 : ((1ULL << bw_sqrt_adjust) - 1));
    uint64_t *sqrt_adjust = new uint64_t[dim];
    int32_t sqrt_adjust_scale =
            (inverse ? floor((bw_x - 1 - s_x) / 2.0) : floor((s_x + 1) / 2.0));
    for (int i = 0; i < dim; i++) {
        sqrt_adjust[i] = 0;
        for (int j = 0; j < (bw_x - 1); j++) {
            if (inverse) {
                sqrt_adjust[i] +=
                        (1ULL << int(floor((s_x - j + 1) / 2.0) + sqrt_adjust_scale)) *
                        msnzb_vector[i * (bw_x - 1) + j];
            } else {
                sqrt_adjust[i] +=
                        (1ULL << int(floor((j - s_x) / 2.0) + sqrt_adjust_scale)) *
                        msnzb_vector[i * (bw_x - 1) + j];
            }
        }
        sqrt_adjust[i] &= mask_sqrt_adjust;
    }
    if (iters > 0 || (!inverse)) {
        // tmp_1: bw = s_y + 1 + bw_sqrt_adjust, scale = s_y + sqrt_adjust_scale
        mult->hadamard_product(dim, x_prev, sqrt_adjust, tmp_1, s_y + 1,
                               bw_sqrt_adjust, s_y + 1 + bw_sqrt_adjust, false,
                               false, MultMode::None);
        // x_curr: bw = s_y + floor(bw_x/2) + 1 - ceil(s_x/2), scale = s_y
        trunc->truncate_and_reduce(dim, tmp_1, x_prev, sqrt_adjust_scale,
                                   s_y + 1 + bw_sqrt_adjust);
        if (bw_y > (s_y + 1 + bw_sqrt_adjust - sqrt_adjust_scale)) {
            xt->z_extend(dim, x_prev, y, s_y + 1 + bw_sqrt_adjust - sqrt_adjust_scale,
                         bw_y);
        } else {
            aux->reduce(dim, x_prev, y, s_y + 1 + bw_sqrt_adjust - sqrt_adjust_scale,
                        bw_y);
        }
    } else {
        // tmp_1: bw = m + SQRT_LOOKUP_SCALE + 1 + bw_sqrt_adjust,
        //        scale = m + SQRT_LOOKUP_SCALE + sqrt_adjust_scale
        mult->hadamard_product(dim, x_prev, sqrt_adjust, tmp_1,
                               m + SQRT_LOOKUP_SCALE + 1, bw_sqrt_adjust,
                               m + SQRT_LOOKUP_SCALE + 1 + bw_sqrt_adjust, false,
                               false, MultMode::None);
        // x_curr: bw = m + floor(bw_x/2) + 1 - ceil(s_x/2), scale = m
        // If iters == 0, we know s_y = m
        trunc->truncate_and_reduce(dim, tmp_1, x_prev,
                                   sqrt_adjust_scale + SQRT_LOOKUP_SCALE,
                                   m + SQRT_LOOKUP_SCALE + 1 + bw_sqrt_adjust);
        if (bw_y > (m + 1 + bw_sqrt_adjust - sqrt_adjust_scale)) {
            xt->z_extend(dim, x_prev, y, m + 1 + bw_sqrt_adjust - sqrt_adjust_scale,
                         bw_y);
        } else {
            aux->reduce(dim, x_prev, y, m + 1 + bw_sqrt_adjust - sqrt_adjust_scale,
                        bw_y);
        }
    }

    delete[] msnzb_vector_bool;
    delete[] msnzb_vector;
    delete[] adjust;
    delete[] exp_parity;
    delete[] adjusted_x;
    delete[] X;
    delete[] tmp_1;
    delete[] x_prev;
    delete[] b_prev;
    delete[] Y_prev;
    delete[] x_curr;
    delete[] b_curr;
    delete[] Y_curr;
    delete[] Y_sq;
    delete[] sqrt_adjust;

    return;
}

void MathFunctions::ReLU(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                         uint64_t six) {
    bool six_comparison = false;
    if (six != 0)
        six_comparison = true;
    int32_t size = (six_comparison ? 2 * dim : dim);

    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

    uint64_t *tmp = new uint64_t[size];
    uint8_t *tmp_msb = new uint8_t[size];
    memcpy(tmp, x, dim * sizeof(uint64_t));
    if (six_comparison) {
        for (int i = 0; i < dim; i++) {
            tmp[dim + i] = (party == ALICE ? x[i] - six : x[i]) & mask_x;
        }
    }
    aux->MSB(tmp, tmp_msb, size, bw_x);
    for (int i = 0; i < size; i++) {
        if (party == ALICE) {
            tmp_msb[i] = tmp_msb[i] ^ 1;
        }
    }
    if (six_comparison)
        aux->AND(tmp_msb, tmp_msb + dim, tmp_msb + dim, dim);

    aux->multiplexer(tmp_msb, tmp, tmp, size, bw_x, bw_x);

    memcpy(y, tmp, dim * sizeof(uint64_t));
    if (six_comparison) {
        for (int i = 0; i < dim; i++) {
            y[i] = (y[i] - tmp[i + dim]) & mask_x;
        }
    }

    delete[] tmp;
    delete[] tmp_msb;
}

void MathFunctions::exp(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                        int32_t bw_y, int32_t s_x, int32_t s_y) {
    // use the limit formula to approximate exp
    // judge the x is positive or negative
    uint8_t *tmp = new uint8_t[dim];
    uint8_t *isPos = new uint8_t[dim];
    uint64_t *tmp1 = new uint64_t[dim];
    uint64_t *tmp2 = new uint64_t[dim];
    uint64_t *res = new uint64_t[dim];
    uint64_t *adjusted_x = new uint64_t[dim];
    aux->wrap_computation(x, tmp, dim, bw_x - 1);
    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            isPos[i] = 1 ^ tmp[i] ^ (uint8_t) (x[i] >> (bw_x - 1));
        }
    } else {
        for (int i = 0; i < dim; i++) {
            isPos[i] = tmp[i] ^ (uint8_t) (x[i] >> (bw_x - 1));
        }
    }

    // cout << "isPos: " << endl;
    // for (int i = 0; i < dim; i++) {
    //     cout << (uint32_t)isPos[i] << " ";
    // }
    // cout << endl;
    for (int i = 0; i < dim; i++) {
        tmp1[i] = x[i] * 2;
        tmp1[i] &= (1ULL << bw_x) - 1;
    }
    aux->multiplexer(isPos, tmp1, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        adjusted_x[i] = (1ULL << (bw_x)) + x[i] - tmp2[i];
        adjusted_x[i] &= (1ULL << bw_x) - 1;
    }
    // cout << "positive x: " << endl;
    // for (int i = 0; i < dim; i++) {
    //     cout << adjusted_x[i] << " ";
    // }
    // cout << endl;
    lookup_table_exp(dim, adjusted_x, res, bw_x, bw_y, s_x, s_y);
    if (party == ALICE) {
        std::fill(tmp1, tmp1 + dim, (1ULL << s_x));
        div(dim, tmp1, res, tmp2, bw_x, bw_x, bw_x, s_x, s_x, s_x, true, true);
    } else {
        std::fill(tmp1, tmp1 + dim, 0);
        div(dim, tmp1, res, tmp2, bw_x, bw_x, bw_x, s_x, s_x, s_x, true, true);
    }
    for (int i = 0; i < dim; i++) {
        tmp1[i] = tmp2[i] - res[i];
        tmp1[i] &= (1ULL << bw_x) - 1;
    }
    aux->multiplexer(isPos, tmp1, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        y[i] = res[i] + tmp2[i];
    }

    delete[] tmp;
    delete[] isPos;
    delete[] tmp1;
    delete[] tmp2;
    delete[] res;
    delete[] adjusted_x;
}


// void MathFunctions::exp(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
//                         int32_t bw_y, int32_t s_x, int32_t s_y) {
//     // use the limit formula to approximate exp
//     // judge the x is positive or negative
//     uint8_t *tmp = new uint8_t[dim];
//     uint8_t *isNeg = new uint8_t[dim];
//     uint64_t *tmp1 = new uint64_t[dim];
//     uint64_t *tmp2 = new uint64_t[dim];
//     uint64_t *adjusted_x = new uint64_t[dim];
//     aux->wrap_computation(x, tmp, dim, bw_x - 1);
//     for (int i = 0; i < dim; i++) {
//         isNeg[i] = tmp[i] ^ (uint8_t) (x[i] >> (bw_x - 1));
//     }
//     for (int i = 0; i < dim; i++) {
//         tmp1[i] = -x[i] * 2;
//         tmp1[i] &= (1ULL << bw_x) - 1;
//     }
//     aux->multiplexer(isNeg, tmp1, tmp2, dim, bw_x, bw_x);
//     for (int i = 0; i < dim; i++) {
//         adjusted_x[i] = x[i] + tmp1[i];
//         adjusted_x[i] &= (1ULL << bw_x) - 1;
//         // adjusted_x is x/(2**8) and change the scale to 30
//         adjusted_x[i] <<= (bw_x - s_x - 2 - 8);
//         // adjusted_x is (1 - x/(2**8)) with the scale to be 30
//         adjusted_x[i] = (1ULL << (bw_x - 3)) - adjusted_x[i];
//     }
//     // use the limit formula to approximate
//     for (int i = 0; i < 7; i++) {
//         mult->hadamard_product(dim, adjusted_x, adjusted_x, tmp1, bw_x - 2, bw_x - 2, 2 * (bw_x - 2));
//         trunc->truncate_and_reduce(dim, tmp1, adjusted_x, bw_x - 2, 2 * (bw_x - 2));
//     }
//     mult->hadamard_product(dim, adjusted_x, adjusted_x, tmp1, bw_x - 2, bw_x - 2, 2 * (bw_x - 2));
//     trunc->truncate_and_reduce(dim, tmp1, tmp2, bw_x - 2 + (bw_x - 2 - s_x), 2 * (bw_x - 2));
//     xt->z_extend(dim, tmp2, adjusted_x, s_x, bw_x);
//
//     delete[] tmp;
//     delete[] isNeg;
//     delete[] tmp1;
// }

// unused
// void MathFunctions::ln(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x, int32_t bw_y, int32_t s_x, int32_t s_y) {
//     int32_t bw_adjust = bw_x - 1;
//     uint64_t mask_adjust = (bw_adjust == 64 ? -1 : ((1ULL << bw_adjust) - 1));
//     // MSB is always 0, thus, not including it
//     uint8_t *msnzb_vector_bool = new uint8_t[dim * (bw_x - 1)];
//     uint64_t *msnzb_vector = new uint64_t[dim * (bw_x - 1)];
//     uint64_t *msnzb_index = new uint64_t[dim];
//     aux->msnzb_one_hot_index(x, msnzb_vector_bool, msnzb_index, bw_x - 1, dim);
//     aux->B2A(msnzb_vector_bool, msnzb_vector, dim * (bw_x - 1), bw_x - 1);
//     uint64_t *adjust = new uint64_t[dim];
//     uint8_t *exp_parity = new uint8_t[dim];
//     for (int i = 0; i < dim; i++) {
//         adjust[i] = 0;
//         exp_parity[i] = 0;
//         for (int j = 0; j < (bw_x - 1); j++) {
//             adjust[i] += (1ULL << (bw_x - 2 - j)) * msnzb_vector[i * (bw_x - 1) + j];
//             if (((j - s_x) & 1)) {
//                 exp_parity[i] ^= msnzb_vector_bool[i * (bw_x - 1) + j];
//             }
//         }
//         adjust[i] &= mask_adjust;
//     }
//     // adjusted_x: bw = bw_x - 1, scale = bw_x - 2, adjusted_x is the value in [1, 2).
//     uint64_t *adjusted_x = new uint64_t[dim];
//     mult->hadamard_product(dim, x, adjust, adjusted_x, bw_x - 1, bw_x - 1,
//                            bw_x - 1, false, false, MultMode::None);
//
//     // adjusted_x is now [1, 2) reduce into a more narrow range.
//     // compare with 1/2 get the comparison result which used for the following selection for ln x
//     uint64_t *tmp = new uint64_t[dim];
//     uint8_t *lethan2 = new uint8_t[dim];
//     uint8_t *lethan4 = new uint8_t[dim];
//     uint8_t *tmp_bool = new uint8_t[dim];
//     uint64_t *r_convert = new uint64_t[dim];
//     uint64_t *tmp2 = new uint64_t[dim]();
//     uint64_t *tmp4 = new uint64_t[dim]();
//     uint64_t *res = new uint64_t[dim]();
//
//
//     trunc->truncate_and_reduce(dim, adjusted_x, tmp, bw_x - 2 - 2, bw_x - 2);
//     // reduce the sign bit and the integer part 1
//     if (party == ALICE) {
//         for (int i = 0; i < dim; i++) {
//             tmp[i] = 4 - tmp[i];
//         }
//     }
//     equality->check_equality(lethan4, tmp, dim, 2);
//
//
//     // cout << "lethan4: " << endl;
//     // for (int i = 0; i < dim; i++) {
//     //     cout << (uint32_t)lethan4[i] << " ";
//     // }
//     // cout << endl;
//
//     cout << endl;
//     cout << "1-2: ";
//     for (int i = 0; i < dim; i++) {
//         cout << adjusted_x[i] << " ";
//     }
//     cout << endl;
//
//     if (party == ALICE) {
//         std::fill(tmp, tmp + dim, (0b11ULL << (bw_x - 3)));
//         div(dim, tmp, adjusted_x, tmp2, bw_x - 1, bw_x - 1, bw_x - 1, bw_x - 2, bw_x - 2, bw_x - 2, false, false);
//     } else {
//         std::fill(tmp, tmp + dim, 0);
//         div(dim, tmp, adjusted_x, tmp2, bw_x - 1, bw_x - 1, bw_x - 1, bw_x - 2, bw_x - 2, bw_x - 2, false, false);
//     }
//
//
//     cout << endl;
//     cout << "0.75-1.25: ";
//     for (int i = 0; i < dim; i++) {
//         cout << tmp2[i] << " ";
//     }
//     cout << endl;
//
//     for (int i = 0; i < dim; i++) {
//         r_convert[i] = adjusted_x[i];
//         r_convert[i] &= mask_adjust;
//     }
//
//
//     // Taylor series
//     uint64_t *resP2 = new uint64_t[dim]();
//     uint64_t *resP3 = new uint64_t[dim]();
//     uint64_t *resP4 = new uint64_t[dim]();
//     uint64_t *resP5 = new uint64_t[dim]();
//     uint64_t *approx_value = new uint64_t[dim]();
//     mult->hadamard_product(dim, res, res, tmp, bw_x - 2, bw_x - 2, 2 * (bw_x - 2));
//     trunc->truncate_and_reduce(dim, tmp, resP2, bw_x - 2, 2 * (bw_x - 2));
//
//     // cout << endl;
//     // for (int i = 0; i < dim; i++) {
//     //     cout << resP2[i] << " ";
//     // }
//     // cout << endl;
//
//     mult->hadamard_product(dim, res, resP2, resP3, bw_x - 2, bw_x - 2, 2 * (bw_x - 2));
//     trunc->truncate_and_reduce(dim, resP3, resP3, bw_x - 2, 2 * (bw_x - 2));
//
//     // cout << endl;
//     // for (int i = 0; i < dim; i++) {
//     //     cout << resP3[i] << " ";
//     // }
//     // cout << endl;
//
//     mult->hadamard_product(dim, resP2, resP2, resP4, bw_x - 2, bw_x - 2, 2 * (bw_x - 2));
//     trunc->truncate_and_reduce(dim, resP4, resP4, bw_x - 2, 2 * (bw_x - 2));
//
//     // cout << endl;
//     // for (int i = 0; i < dim; i++) {
//     //     cout << resP4[i] << " ";
//     // }
//     // cout << endl;
//
//     mult->hadamard_product(dim, resP2, resP3, resP5, bw_x - 2, bw_x - 2, 2 * (bw_x - 2));
//     trunc->truncate_and_reduce(dim, resP5, resP5, bw_x - 2, 2 * (bw_x - 2));
//
//     // cout << endl;
//     // for (int i = 0; i < dim; i++) {
//     //     cout << resP5[i] << " ";
//     // }
//     // cout << endl;
//
//
//     trunc->div_pow2(dim, resP2, tmp, 1, bw_x - 2);
//
//     aux->wrap_computation(resP3, tmp_bool, dim, bw_x - 2);
//     // here can replace by the multiplication of (1/3) to reduce the call of wrap and mux, but introduce one-bit error.
//     aux->multiplexer_two_plain(tmp_bool, ((1 << (bw_x - 2)) / 3), tmp4, dim, bw_x - 2, bw_x - 2);
//     for (int i = 0; i < dim; i++) {
//         resP3[i] /= 3;
//         resP3[i] -= tmp4[i];
//         resP3[i] &= (mask_adjust >> 1);
//     }
//
//     trunc->div_pow2(dim, resP4, tmp2, 2, bw_x - 2);
//
//     aux->wrap_computation(resP5, tmp_bool, dim, bw_x - 2);
//     aux->multiplexer_two_plain(tmp_bool, ((1 << (bw_x - 2)) / 5), tmp4, dim, bw_x - 2, bw_x - 2);
//     for (int i = 0; i < dim; i++) {
//         resP5[i] /= 5;
//         resP5[i] -= tmp4[i];
//         resP5[i] &= (mask_adjust >> 1);
//     }
//
//
//     // cout << endl;
//     // for (int i = 0; i < dim; i++) {
//     //     cout << res[i] << " ";
//     // }
//     // cout << endl;
//     // for (int i = 0; i < dim; i++) {
//     //     cout << tmp[i] << " ";
//     // }
//     // cout << endl;
//     // for (int i = 0; i < dim; i++) {
//     //     cout << resP3[i] << " ";
//     // }
//     // cout << endl;
//     // for (int i = 0; i < dim; i++) {
//     //     cout << tmp2[i] << " ";
//     // }
//     // cout << endl;
//     // for (int i = 0; i < dim; i++) {
//     //     cout << resP5[i] << " ";
//     // }
//     // cout << endl;
//
//
//     for (int i = 0; i < dim; i++) {
//         approx_value[i] = res[i] - tmp[i] + resP3[i] - tmp2[i] + resP5[i];
//         approx_value[i] &= (mask_adjust >> 1);
//     }
//
//
//     // combine the final result
//     // convert the scale back to s_x
//     trunc->truncate_and_reduce(dim, approx_value, approx_value, bw_x - 2 - 16, (bw_x - 2));
//     xt->z_extend(dim, approx_value, approx_value, 16, bw_x);
//
//
//     // lethan4
//     uint64_t *convert_res = new uint64_t[dim]();
//     if (party == ALICE) {
//         for (int i = 0; i < dim; i++) {
//             convert_res[i] = (uint64_t) (log(25.0 / 16) * (1 << (s_x))) - approx_value[i];
//             convert_res[i] &= ((1ULL << bw_x) - 1);
//         }
//     } else {
//         for (int i = 0; i < dim; i++) {
//             convert_res[i] = 0ULL - approx_value[i];
//             convert_res[i] &= ((1ULL << bw_x) - 1);
//         }
//     }
//
//     for (int i = 0; i < dim; i++) {
//         tmp[i] = approx_value[i] - convert_res[i];
//     }
//     aux->multiplexer(lethan4, tmp, tmp2, dim, bw_x, bw_x);
//     for (int i = 0; i < dim; i++) {
//         tmp2[i] += convert_res[i];
//     }
//
//     // tmp4 store the ln part
//
//
//     // the t
//     int msnzb_index_bits = ceil(log2(bw_x));
//     uint64_t *tmp5 = new uint64_t[dim]();
//     for (int i = 0; i < dim; i++) {
//         msnzb_index[i] *= (1 << s_x);
//     }
//     xt->z_extend(dim, msnzb_index, msnzb_index, s_x + msnzb_index_bits, bw_x);
//
//
//     for (int i = 0; i < dim; i++) {
//         tmp5[i] = msnzb_index[i] * (uint64_t) (log(2.0) * (1 << (s_x)));
//     }
//     trunc->truncate_and_reduce(dim, tmp5, tmp5, s_x, bw_x + s_x);
//     aux->wrap_computation(msnzb_index, tmp_bool, dim, bw_x);
//     aux->multiplexer_two_plain(
//         tmp_bool, (((uint64_t) (log(2.0) * (1ULL << (s_x)))) << (bw_x - s_x)) & ((1ULL << bw_x) - 1), tmp, dim, bw_x,
//         bw_x);
//     for (int i = 0; i < dim; i++) {
//         tmp5[i] -= tmp[i];
//     }
//
//     if (party == ALICE) {
//         for (int i = 0; i < dim; i++) {
//             tmp5[i] -= (uint64_t) (16 * log(2.0) * (1 << s_x));
//         }
//     }
//
//     for (int i = 0; i < dim; i++) {
//         tmp5[i] &= ((1ULL << bw_x) - 1);
//     }
//
//
//     // cout << endl;
//     // cout << "tmp4: ";
//     // for (int i = 0; i < dim; i++) {
//     //     cout << tmp4[i] << ", ";
//     // }
//     // cout << endl;
//     //
//     // cout << endl;
//     // cout << "tmp5: ";
//     // for (int i = 0; i < dim; i++) {
//     //     cout << tmp5[i] << ", ";
//     // }
//     // cout << endl;
//
//     for (int i = 0; i < dim; i++) {
//         y[i] = tmp5[i] + tmp4[i];
//         y[i] &= (1ULL << bw_x) - 1;
//     }
//
//     delete[] tmp;
//     delete[] tmp2;
//     delete[] tmp4;
//     delete[] tmp5;
//     delete[] lethan2;
//     delete[] lethan4;
//     delete[] res;
//     delete[] resP2;
//     delete[] resP3;
//     delete[] resP4;
//     delete[] resP5;
//     delete[] approx_value;
//     delete[] convert_res;
//     delete[] r_convert;
//     delete[] adjusted_x;
//     delete[] msnzb_vector_bool;
//     delete[] msnzb_vector;
//     delete[] msnzb_index;
//     delete[] exp_parity;
//     delete[] adjust;
//     delete[] tmp_bool;
// }

void MathFunctions::ln_v1(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x, int32_t bw_y, int32_t s_x, int32_t s_y) {
    int32_t bw_adjust = bw_x - 1;
    uint64_t mask_adjust = (bw_adjust == 64 ? -1 : ((1ULL << bw_adjust) - 1));
    // MSB is always 0, thus, not including it
    uint8_t *msnzb_vector_bool = new uint8_t[dim * (bw_x - 1)];
    uint64_t *msnzb_vector = new uint64_t[dim * (bw_x - 1)];
    uint64_t *msnzb_index = new uint64_t[dim];

    auto start = clock_start();
    uint64_t comm = iopack->get_comm();

    aux->msnzb_one_hot_index(x, msnzb_vector_bool, msnzb_index, bw_x - 1, dim);
    aux->B2A(msnzb_vector_bool, msnzb_vector, dim * (bw_x - 1), bw_x - 1);
    uint64_t *adjust = new uint64_t[dim];
    uint8_t *exp_parity = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        adjust[i] = 0;
        exp_parity[i] = 0;
        for (int j = 0; j < (bw_x - 1); j++) {
            adjust[i] += (1ULL << (bw_x - 2 - j)) * msnzb_vector[i * (bw_x - 1) + j];
            if (((j - s_x) & 1)) {
                exp_parity[i] ^= msnzb_vector_bool[i * (bw_x - 1) + j];
            }
        }
        adjust[i] &= mask_adjust;
    }
    // adjusted_x: bw = bw_x - 1, scale = bw_x - 2, adjusted_x is the value in [1, 2).
    uint64_t *adjusted_x = new uint64_t[dim];
    mult->hadamard_product(dim, x, adjust, adjusted_x, bw_x - 1, bw_x - 1,
                           bw_x - 1, false, false, MultMode::None);

    // adjusted_x is now [1, 2) reduce into a more narrow range.
    // compare with 1/2 get the comparison result which used for the following selection for ln x

    // cout << endl;
    // cout << "adjust_x [1-2): ";
    // for (int i = 0; i < dim; i++) {
    //     cout << adjusted_x[i] << ", ";
    // }
    // cout << endl;

    uint64_t *tmp = new uint64_t[dim];
    uint8_t *lethan2 = new uint8_t[dim];
    uint8_t *lethan4 = new uint8_t[dim];
    uint8_t *tmp_bool = new uint8_t[dim];
    uint64_t *r_convert = new uint64_t[dim];
    uint64_t *tmp2 = new uint64_t[dim]();
    uint64_t *tmp4 = new uint64_t[dim]();
    uint64_t *res = new uint64_t[dim]();


    trunc->truncate_and_reduce(dim, adjusted_x, tmp, bw_x - 2 - 1, bw_x - 2);
    // reduce the sign bit and the integer part 1
    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            tmp[i] = 2 - tmp[i];
        }
    }
    equality->check_equality(lethan2, tmp, dim, 1);

    // set the corresponding r, which is in [0, 1/2] that evaluated by the Taylor series
    for (int i = 0; i < dim; i++) {
        r_convert[i] = adjusted_x[i];
        r_convert[i] &= mask_adjust;
    }


    aux->wrap_computation(r_convert, tmp_bool, dim, bw_x - 1);
    aux->multiplexer_two_plain(tmp_bool, (1ULL << (bw_x - 1)) / 4, tmp4, dim, bw_x - 1, bw_x - 1);


    for (int i = 0; i < dim; i++) {
        r_convert[i] /= 4;
        r_convert[i] -= tmp4[i];
        r_convert[i] = adjusted_x[i] - r_convert[i];
        r_convert[i] &= mask_adjust;
    }

    for (int i = 0; i < dim; i++) {
        tmp4[i] = adjusted_x[i] - r_convert[i];
        tmp4[i] &= mask_adjust;
    }
    aux->multiplexer(lethan2, tmp4, res, dim, bw_x - 1, bw_x - 1);
    for (int i = 0; i < dim; i++) {
        res[i] += r_convert[i];
        res[i] &= mask_adjust;
    }

    // cout << endl;
    // cout << "res [1-1.5]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;

    // the res store x in [1-1.5], now reduce it into [1-1.25]
    trunc->truncate_and_reduce(dim, res, tmp, bw_x - 2 - 2, bw_x - 3); //the sign bit x - 1 < 0.5
    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            tmp[i] = 2 - tmp[i];
        }
    }
    equality->check_equality(lethan4, tmp, dim, 1);

    // set the corresponding r, which is in [0, 1/4] that evaluated by the Taylor series
    // sometimes the divide result is wrong

    for (int i = 0; i < dim; i++) {
        r_convert[i] = res[i];
        r_convert[i] &= mask_adjust;
    }


    aux->wrap_computation(r_convert, tmp_bool, dim, bw_x - 1);
    aux->multiplexer_two_plain(tmp_bool, 0x15555555ULL, tmp4, dim, bw_x - 1, bw_x - 1);


    for (int i = 0; i < dim; i++) {
        r_convert[i] /= 6;
        r_convert[i] -= tmp4[i];
        r_convert[i] = res[i] - r_convert[i];
        r_convert[i] &= mask_adjust;
    }


    for (int i = 0; i < dim; i++) {
        tmp4[i] = res[i] - r_convert[i];
        tmp4[i] &= mask_adjust;
    }
    aux->multiplexer(lethan4, tmp4, res, dim, bw_x - 1, bw_x - 1);
    for (int i = 0; i < dim; i++) {
        res[i] += r_convert[i];
        res[i] &= (mask_adjust >> 1); // mask_adjust is 31 bits, set the result be [0, 0.25] for Taylor series
    }


    comm = iopack->get_comm() - comm;
    long long t = time_from(start);

    cout << "MSNZB Time\t" << t / (1000.0) << " ms" << endl;
    cout << "MSNZB Bytes Sent\t" << comm << " bytes" << endl;


    // cout << endl;
    // cout << "res [0-0.25]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;

    // Taylor series
    uint64_t *resP2 = new uint64_t[dim]();
    uint64_t *resP3 = new uint64_t[dim]();
    uint64_t *resP4 = new uint64_t[dim]();
    uint64_t *resP5 = new uint64_t[dim]();
    // uint64_t *temp_constant = new uint64_t[dim]();
    uint64_t *approx_value = new uint64_t[dim]();
    mult->hadamard_product(dim, res, res, tmp, bw_x - 2, bw_x - 2, 2 * (bw_x - 2));
    trunc->truncate_and_reduce(dim, tmp, resP2, bw_x - 2, 2 * (bw_x - 2));
    // trunc->truncate_direct(dim, tmp, resP2, bw_x - 2, 2 * (bw_x - 2));

    mult->hadamard_product(dim, res, resP2, resP3, bw_x - 2, bw_x - 2, 2 * (bw_x - 2));
    trunc->truncate_and_reduce(dim, resP3, resP3, bw_x - 2, 2 * (bw_x - 2));
    // trunc->truncate_direct(dim, resP3, resP3, bw_x - 2, 2 * (bw_x - 2));

    mult->hadamard_product(dim, resP2, resP2, resP4, bw_x - 2, bw_x - 2, 2 * (bw_x - 2));
    trunc->truncate_and_reduce(dim, resP4, resP4, bw_x - 2, 2 * (bw_x - 2));
    // trunc->truncate_direct(dim, resP4, resP4, bw_x - 2, 2 * (bw_x - 2));


    mult->hadamard_product(dim, resP2, resP3, resP5, bw_x - 2, bw_x - 2, 2 * (bw_x - 2));
    trunc->truncate_and_reduce(dim, resP5, resP5, bw_x - 2, 2 * (bw_x - 2));
    // trunc->truncate_direct(dim, resP5, resP5, bw_x - 2, 2 * (bw_x - 2));


    trunc->div_pow2(dim, resP2, tmp, 1, bw_x - 2);

    // aux->wrap_computation(resP3, tmp_bool, dim, bw_x - 2);
    // aux->multiplexer_two_plain(tmp_bool, ((1 << (bw_x - 2)) / 3), tmp4, dim, bw_x - 2, bw_x - 2);
    // for (int i = 0; i < dim; i++) {
    //     resP3[i] /= 3;
    //     resP3[i] -= tmp4[i];
    //     resP3[i] &= (mask_adjust >> 1);
    // }

    aux->wrap_computation(resP3, tmp_bool, dim, bw_x - 2);
    aux->multiplexer_two_plain(tmp_bool, ((1ULL << (bw_x - 2)) / 3), tmp4, dim, bw_x - 2, bw_x - 2);
    for (int i = 0; i < dim; i++) {
        resP3[i] = resP3[i] * ((1ULL << (bw_x - 2)) / 3);
        resP3[i] -= tmp4[i] << (bw_x - 2);
    }
    trunc->truncate_and_reduce(dim, resP3, resP3, bw_x - 2, 2 * (bw_x - 2));
    // trunc->truncate_direct(dim, resP3, resP3, bw_x - 2, 2 * (bw_x - 2));


    trunc->div_pow2(dim, resP4, tmp2, 2, bw_x - 2);

    // aux->wrap_computation(resP5, tmp_bool, dim, bw_x - 2);
    // aux->multiplexer_two_plain(tmp_bool, ((1ULL << (bw_x - 2)) / 5), tmp4, dim, bw_x - 2, bw_x - 2);
    // for (int i = 0; i < dim; i++) {
    //     resP5[i] /= 5;  // divide will occur 0 to be 0.9999xx
    //     resP5[i] -= tmp4[i];
    //     resP5[i] &= (mask_adjust >> 1);
    // }

    aux->wrap_computation(resP5, tmp_bool, dim, bw_x - 2);
    aux->multiplexer_two_plain(tmp_bool, ((1ULL << (bw_x - 2)) / 5), tmp4, dim, bw_x - 2, bw_x - 2);
    for (int i = 0; i < dim; i++) {
        resP5[i] = resP5[i] * ((1ULL << (bw_x - 2)) / 5);
        resP5[i] -= tmp4[i] << (bw_x - 2);
    }
    trunc->truncate_and_reduce(dim, resP5, resP5, bw_x - 2, 2 * (bw_x - 2));
    // trunc->truncate_direct(dim, resP5, resP5, bw_x - 2, 2 * (bw_x - 2));

    // cout << endl;
    // cout << "res: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;
    //
    // cout << endl;
    // cout << "tmp: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << tmp[i] << ", ";
    // }
    // cout << endl;
    //
    // cout << endl;
    // cout << "resP3: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << resP3[i] << ", ";
    // }
    // cout << endl;
    //
    // cout << endl;
    // cout << "tmp2: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << tmp2[i] << ", ";
    // }
    // cout << endl;
    // cout << endl;
    // cout << "resP5: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << resP5[i] << ", ";
    // }
    // cout << endl;


    for (int i = 0; i < dim; i++) {
        approx_value[i] = res[i] - tmp[i] + resP3[i] - tmp2[i] + resP5[i];
        approx_value[i] &= (mask_adjust >> 1);
    }


    // combine the final result
    // convert the scale back to s_x
    trunc->truncate_and_reduce(dim, approx_value, approx_value, bw_x - 2 - 16, (bw_x - 2));
    // trunc->truncate_direct(dim, approx_value, approx_value, bw_x - 2 - 16, (bw_x - 2));
    xt->z_extend(dim, approx_value, approx_value, 16, bw_x);

    // cout << endl;
    // cout << "approx_value: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << approx_value[i] << ", ";
    // }
    // cout << endl;

    // lethan4
    uint64_t *convert_res = new uint64_t[dim]();
    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            convert_res[i] = approx_value[i] + (uint64_t) (log(6.0 / 5) * (1 << (s_x)));
            convert_res[i] &= ((1ULL << bw_x) - 1);
        }
    } else {
        for (int i = 0; i < dim; i++) {
            convert_res[i] = approx_value[i];
            convert_res[i] &= ((1ULL << bw_x) - 1);
        }
    }

    for (int i = 0; i < dim; i++) {
        tmp[i] = approx_value[i] - convert_res[i];
    }
    aux->multiplexer(lethan4, tmp, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        tmp2[i] += convert_res[i];
    }

    // lethan2
    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            convert_res[i] = tmp2[i] + (uint64_t) (log(4.0 / 3) * (1 << (s_x)));
            convert_res[i] &= ((1ULL << bw_x) - 1);
        }
    } else {
        for (int i = 0; i < dim; i++) {
            convert_res[i] = tmp2[i];
            convert_res[i] &= ((1ULL << bw_x) - 1);
        }
    }

    for (int i = 0; i < dim; i++) {
        tmp[i] = tmp2[i] - convert_res[i];
    }
    aux->multiplexer(lethan2, tmp, tmp4, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        tmp4[i] += convert_res[i];
    }
    // tmp4 store the ln part


    // the t
    int msnzb_index_bits = ceil(log2(bw_x));
    uint64_t *tmp5 = new uint64_t[dim]();
    for (int i = 0; i < dim; i++) {
        msnzb_index[i] *= (1 << s_x);
    }
    xt->z_extend(dim, msnzb_index, msnzb_index, s_x + msnzb_index_bits, bw_x);


    for (int i = 0; i < dim; i++) {
        tmp5[i] = msnzb_index[i] * (uint64_t) (log(2.0) * (1 << (s_x)));
    }
    trunc->truncate_and_reduce(dim, tmp5, tmp5, s_x, bw_x + s_x);
    aux->wrap_computation(msnzb_index, tmp_bool, dim, bw_x);
    aux->multiplexer_two_plain(
        tmp_bool, (((uint64_t) (log(2.0) * (1ULL << (s_x)))) << (bw_x - s_x)) & ((1ULL << bw_x) - 1), tmp, dim, bw_x,
        bw_x);
    for (int i = 0; i < dim; i++) {
        tmp5[i] -= tmp[i];
    }

    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            tmp5[i] -= (uint64_t) (16 * log(2.0) * (1 << s_x));
        }
    }

    for (int i = 0; i < dim; i++) {
        tmp5[i] &= ((1ULL << bw_x) - 1);
    }


    for (int i = 0; i < dim; i++) {
        y[i] = tmp5[i] + tmp4[i];
        y[i] &= (1ULL << bw_x) - 1;
    }

    delete[] tmp;
    delete[] tmp2;
    delete[] tmp4;
    delete[] tmp5;
    delete[] lethan2;
    delete[] lethan4;
    delete[] res;
    delete[] resP2;
    delete[] resP3;
    delete[] resP4;
    delete[] resP5;
    delete[] approx_value;
    delete[] convert_res;
    delete[] r_convert;
    delete[] adjusted_x;
    delete[] msnzb_vector_bool;
    delete[] msnzb_vector;
    delete[] msnzb_index;
    delete[] exp_parity;
    delete[] adjust;
    delete[] tmp_bool;
}

void MathFunctions::ln_v2(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x, int32_t bw_y, int32_t s_x, int32_t s_y) {
    int32_t bw_adjust = bw_x - 1;
    uint64_t mask_adjust = (bw_adjust == 64 ? -1 : ((1ULL << bw_adjust) - 1));
    // MSB is always 0, thus, not including it
    uint8_t *msnzb_vector_bool = new uint8_t[dim * (bw_x - 1)];
    uint64_t *msnzb_vector = new uint64_t[dim * (bw_x - 1)];
    uint64_t *msnzb_index = new uint64_t[dim];

    auto start = clock_start();
    uint64_t comm = iopack->get_comm();
    // aux->msnzb_one_hot(x, msnzb_vector_bool, bw_x - 1, dim);
    aux->msnzb_one_hot_index(x, msnzb_vector_bool, msnzb_index, bw_x - 1, dim);
    aux->B2A(msnzb_vector_bool, msnzb_vector, dim * (bw_x - 1), bw_x - 1);
    uint64_t *adjust = new uint64_t[dim];
    uint8_t *exp_parity = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        adjust[i] = 0;
        exp_parity[i] = 0;
        for (int j = 0; j < (bw_x - 1); j++) {
            adjust[i] += (1ULL << (bw_x - 2 - j)) * msnzb_vector[i * (bw_x - 1) + j];
            if (((j - s_x) & 1)) {
                exp_parity[i] ^= msnzb_vector_bool[i * (bw_x - 1) + j];
            }
        }
        adjust[i] &= mask_adjust;
    }
    // adjusted_x: bw = bw_x - 1, scale = bw_x - 2, adjusted_x is the value in [1, 2).
    uint64_t *adjusted_x = new uint64_t[dim];
    mult->hadamard_product(dim, x, adjust, adjusted_x, bw_x - 1, bw_x - 1,
                           bw_x - 1, false, false, MultMode::None);
    xt->z_extend(dim, adjusted_x, adjusted_x, bw_x - 1, bw_x);


    comm = iopack->get_comm() - comm;
    long long t = time_from(start);

    cout << "MSNZB Time\t" << t / (1000.0) << " ms" << endl;
    cout << "MSNZB Bytes Sent\t" << comm << " bytes" << endl;


    // adjusted_x is now [1, 2) reduce into a more narrow range.
    // compare with 1/2 get the comparison result which used for the following selection for ln x

    // cout << endl;
    // cout << "adjust_x [1-2): ";
    // for (int i = 0; i < dim; i++) {
    //     cout << adjusted_x[i] << ", ";
    // }
    // cout << endl;

    uint64_t *tmp = new uint64_t[dim];
    uint8_t *gethan1 = new uint8_t[dim];
    uint8_t *gethan2 = new uint8_t[dim];
    uint8_t *gethan3 = new uint8_t[dim];
    uint8_t *gethan4 = new uint8_t[dim];
    uint8_t *gethan5 = new uint8_t[dim];
    uint8_t *msbx = new uint8_t[dim]();
    uint8_t *tmp_bool = new uint8_t[dim];
    uint64_t *r_convert = new uint64_t[dim];
    uint64_t *tmp2 = new uint64_t[dim]();
    uint64_t *tmp4 = new uint64_t[dim]();
    uint64_t *res = new uint64_t[dim]();


    trunc->truncate_and_reduce(dim, adjusted_x, tmp, bw_x - 3, bw_x - 2);
    for (int i = 0; i < dim; i++) {
        gethan1[i] = tmp[i] & 1;
        r_convert[i] = adjusted_x[i];
    }

    // r_convert used to temp store the r*adjust_x r=3/4 for [1,2] to [1,1.5]
    // adjust_x: bw = bw_x, scale = bw_x - 2 bw = 32, scale = 30
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= (3ULL << (bw_x - 4));
    }
    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);
    // r_convert bw = bw_x, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - adjusted_x[i];
    }
    aux->multiplexer(gethan1, tmp4, res, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += adjusted_x[i];
    }

    // cout << endl;
    // cout << "res [1-1.5]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;


    // res: bw = bw_x - 1 scale = bw_x - 2, res, bw = 31, scale = 30, the 29-th (0.5, bw - 3) bit is zero
    trunc->truncate_and_reduce(dim, res, tmp, bw_x - 4, bw_x - 3);
    for (int i = 0; i < dim; i++) {
        gethan2[i] = tmp[i] & 1;
        r_convert[i] = res[i];
    }

    // r_convert used to temp store the r*res r=5/6 for [1,1.5] to [1,1.25]
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= ((5ULL << (bw_x - 2)) / 6);
    }

    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);

    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan2, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    // cout << endl;
    // cout << "res [1-1.25]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;


    // res: bw = bw_x - 1 scale = bw_x - 2, res, bw = 31, scale = 30, the 29-th (0.5, bw - 3) bit is zero
    trunc->truncate_and_reduce(dim, res, tmp, bw_x - 5, bw_x - 4);
    for (int i = 0; i < dim; i++) {
        gethan3[i] = tmp[i] & 1;
        r_convert[i] = res[i];
    }

    // r_convert used to temp store the r*res r=5/6 for [1,1.5] to [1,1.25]
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= ((9ULL << (bw_x - 2)) / 10);
    }

    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);

    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan3, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    // cout << endl;
    // cout << "res [1-1.125]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;


    // res: bw = bw_x - 1 scale = bw_x - 2, res, bw = 31, scale = 30, the 29-th (0.5, bw - 3) bit is zero
    trunc->truncate_and_reduce(dim, res, tmp, bw_x - 6, bw_x - 5);
    for (int i = 0; i < dim; i++) {
        gethan4[i] = tmp[i] & 1;
        r_convert[i] = res[i];
    }

    // r_convert used to temp store the r*res r=5/6 for [1,1.5] to [1,1.25]
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= ((17ULL << (bw_x - 2)) / 18);
    }

    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);

    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan4, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    // cout << endl;
    // cout << "res [1-1.0625]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;
    //
    // res: bw = bw_x - 1 scale = bw_x - 2, res, bw = 31, scale = 30, the 29-th (0.5, bw - 3) bit is zero
    trunc->truncate_and_reduce(dim, res, tmp, bw_x - 7, bw_x - 6);
    for (int i = 0; i < dim; i++) {
        gethan5[i] = tmp[i] & 1;
        r_convert[i] = res[i];
    }

    // r_convert used to temp store the r*res r=5/6 for [1,1.5] to [1,1.25]
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= ((33ULL << (bw_x - 2)) / 34);
    }

    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);

    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan5, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }
    // cout << endl;
    // cout << "res [1-1.03125]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;

    // - a + (b * (x-1))
    // a = 0.9698751292540353
    // b = 0.00030997309020363434
    uint64_t B = static_cast<uint64_t>(0.9845064 * (pow(2.0, bw_x - 2)));
    uint64_t A = static_cast<uint64_t>(0.000015139157 * (pow(2.0, bw_x - 2)));
    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            res[i] -= (1ULL << (bw_x - 2));
        }
    }

    xt->z_extend(dim, res, res, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        res[i] = B * res[i];
    }
    trunc->truncate_direct(dim, res, res, bw_x - 2, 2 * bw_x);

    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            res[i] = res[i] - A;
        }
    }

    // cout << endl;
    // cout << "ln(x): ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;


    // scale back
    uint64_t scale1 = static_cast<uint64_t>(0.28768207245178085 * (pow(2.0, bw_x - 2)));
    uint64_t scale2 = static_cast<uint64_t>(0.1823215567939546 * (pow(2.0, bw_x - 2)));
    uint64_t scale3 = static_cast<uint64_t>(0.10536051565782635 * (pow(2.0, bw_x - 2)));
    uint64_t scale4 = static_cast<uint64_t>(0.05715841383994862 * (pow(2.0, bw_x - 2)));
    uint64_t scale5 = static_cast<uint64_t>(0.02985296314968113 * (pow(2.0, bw_x - 2)));
    // uint64_t scale6 = static_cast<uint64_t>(0.01526747213078838 * (pow(2.0, bw_x - 2)));

    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            r_convert[i] = res[i] + scale5;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            r_convert[i] = res[i];
        }
    }
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan5, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            r_convert[i] = res[i] + scale4;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            r_convert[i] = res[i];
        }
    }
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan4, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            r_convert[i] = res[i] + scale3;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            r_convert[i] = res[i];
        }
    }
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan3, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            r_convert[i] = res[i] + scale2;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            r_convert[i] = res[i];
        }
    }
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan2, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            r_convert[i] = res[i] + scale1;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            r_convert[i] = res[i];
        }
    }
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan1, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    // cout << endl;
    // cout << "ln(x): ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;

    trunc->truncate_direct(dim, res, res, bw_x - 2 - s_x, bw_x);
    xt->z_extend(dim, res, res, s_x + 2, bw_x, msbx);

    // the t
    int msnzb_index_bits = ceil(log2(bw_x));
    uint64_t *tmp5 = new uint64_t[dim]();
    for (int i = 0; i < dim; i++) {
        msnzb_index[i] *= (1 << s_x);
    }
    xt->z_extend(dim, msnzb_index, msnzb_index, s_x + msnzb_index_bits, bw_x);

    for (int i = 0; i < dim; i++) {
        tmp5[i] = msnzb_index[i] * (uint64_t) (log(2.0) * (1 << (s_x)));
    }
    trunc->truncate_and_reduce(dim, tmp5, tmp5, s_x, bw_x + s_x);
    aux->wrap_computation(msnzb_index, tmp_bool, dim, bw_x);
    aux->multiplexer_two_plain(
        tmp_bool, (((uint64_t) (log(2.0) * (1ULL << (s_x)))) << (bw_x - s_x)) & ((1ULL << bw_x) - 1), tmp, dim, bw_x,
        bw_x);
    for (int i = 0; i < dim; i++) {
        tmp5[i] -= tmp[i];
    }

    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            tmp5[i] -= (uint64_t) (16 * log(2.0) * (1 << s_x));
        }
    }

    for (int i = 0; i < dim; i++) {
        tmp5[i] &= ((1ULL << bw_x) - 1);
    }


    for (int i = 0; i < dim; i++) {
        y[i] = tmp5[i] + res[i];
        y[i] &= (1ULL << bw_x) - 1;
    }

    delete[] tmp;
    delete[] tmp2;
    delete[] tmp4;
    delete[] gethan1;
    delete[] gethan2;
    delete[] gethan3;
    delete[] gethan4;
    delete[] gethan5;
    delete[] msbx;
    delete[] res;
    delete[] r_convert;
    delete[] adjusted_x;
    delete[] msnzb_vector_bool;
    delete[] msnzb_vector;
    delete[] msnzb_index;
    delete[] exp_parity;
    delete[] adjust;
    delete[] tmp_bool;
}


void MathFunctions::sqrt_v2(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                            int32_t bw_y, int32_t s_x, int32_t s_y, bool inverse) {
    int32_t bw_adjust = bw_x - 1;
    uint64_t mask_adjust = (bw_adjust == 64 ? -1 : ((1ULL << bw_adjust) - 1));
    // MSB is always 0, thus, not including it
    uint8_t *msnzb_vector_bool = new uint8_t[dim * (bw_x - 1)];
    uint64_t *msnzb_vector = new uint64_t[dim * (bw_x - 1)];
    uint64_t *msnzb_index = new uint64_t[dim];

    auto start = clock_start();
    uint64_t comm = iopack->get_comm();
    // aux->msnzb_one_hot(x, msnzb_vector_bool, bw_x - 1, dim);
    aux->msnzb_one_hot_index(x, msnzb_vector_bool, msnzb_index, bw_x - 1, dim);
    aux->B2A(msnzb_vector_bool, msnzb_vector, dim * (bw_x - 1), bw_x - 1);
    uint64_t *adjust = new uint64_t[dim];
    uint8_t *exp_parity = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        adjust[i] = 0;
        exp_parity[i] = 0;
        for (int j = 0; j < (bw_x - 1); j++) {
            adjust[i] += (1ULL << (bw_x - 2 - j)) * msnzb_vector[i * (bw_x - 1) + j];
            if (((j - s_x) & 1)) {
                exp_parity[i] ^= msnzb_vector_bool[i * (bw_x - 1) + j];
            }
        }
        adjust[i] &= mask_adjust;
    }
    // adjusted_x: bw = bw_x - 1, scale = bw_x - 2, adjusted_x is the value in [1, 2).
    uint64_t *adjusted_x = new uint64_t[dim];
    mult->hadamard_product(dim, x, adjust, adjusted_x, bw_x - 1, bw_x - 1,
                           bw_x - 1, false, false, MultMode::None);
    xt->z_extend(dim, adjusted_x, adjusted_x, bw_x - 1, bw_x);


    comm = iopack->get_comm() - comm;
    long long t = time_from(start);

    cout << "MSNZB Time\t" << t / (1000.0) << " ms" << endl;
    cout << "MSNZB Bytes Sent\t" << comm << " bytes" << endl;


    // adjusted_x is now [1, 2) reduce into a more narrow range.
    // compare with 1/2 get the comparison result which used for the following selection for ln x

    // cout << endl;
    // cout << "adjust_x [1-2): ";
    // for (int i = 0; i < dim; i++) {
    //     cout << adjusted_x[i] << ", ";
    // }
    // cout << endl;

    uint64_t *tmp = new uint64_t[dim];
    uint8_t *gethan1 = new uint8_t[dim];
    uint8_t *gethan2 = new uint8_t[dim];
    uint8_t *gethan3 = new uint8_t[dim];
    uint8_t *gethan4 = new uint8_t[dim];
    uint8_t *gethan5 = new uint8_t[dim];
    uint8_t *gethan6 = new uint8_t[dim];
    uint8_t *msbx = new uint8_t[dim]();
    uint8_t *tmp_bool = new uint8_t[dim];
    uint64_t *r_convert = new uint64_t[dim];
    uint64_t *tmp2 = new uint64_t[dim]();
    uint64_t *tmp4 = new uint64_t[dim]();
    uint64_t *res = new uint64_t[dim]();


    trunc->truncate_and_reduce(dim, adjusted_x, tmp, bw_x - 3, bw_x - 2);
    for (int i = 0; i < dim; i++) {
        gethan1[i] = tmp[i] & 1;
        r_convert[i] = adjusted_x[i];
    }

    // r_convert used to temp store the r*adjust_x r=3/4 for [1,2] to [1,1.5]
    // adjust_x: bw = bw_x, scale = bw_x - 2 bw = 32, scale = 30
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= (3ULL << (bw_x - 4));
    }
    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);
    // r_convert bw = bw_x, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - adjusted_x[i];
    }
    aux->multiplexer(gethan1, tmp4, res, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += adjusted_x[i];
    }

    // cout << endl;
    // cout << "res [1-1.5]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;


    // res: bw = bw_x - 1 scale = bw_x - 2, res, bw = 31, scale = 30, the 29-th (0.5, bw - 3) bit is zero
    trunc->truncate_and_reduce(dim, res, tmp, bw_x - 4, bw_x - 3);
    for (int i = 0; i < dim; i++) {
        gethan2[i] = tmp[i] & 1;
        r_convert[i] = res[i];
    }

    // r_convert used to temp store the r*res r=5/6 for [1,1.5] to [1,1.25]
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= ((5ULL << (bw_x - 2)) / 6);
    }

    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);

    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan2, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    // cout << endl;
    // cout << "res [1-1.25]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;


    // res: bw = bw_x - 1 scale = bw_x - 2, res, bw = 31, scale = 30, the 29-th (0.5, bw - 3) bit is zero
    trunc->truncate_and_reduce(dim, res, tmp, bw_x - 5, bw_x - 4);
    for (int i = 0; i < dim; i++) {
        gethan3[i] = tmp[i] & 1;
        r_convert[i] = res[i];
    }

    // r_convert used to temp store the r*res r=5/6 for [1,1.5] to [1,1.25]
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= ((9ULL << (bw_x - 2)) / 10);
    }

    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);

    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan3, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    // cout << endl;
    // cout << "res [1-1.125]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;


    // res: bw = bw_x - 1 scale = bw_x - 2, res, bw = 31, scale = 30, the 29-th (0.5, bw - 3) bit is zero
    trunc->truncate_and_reduce(dim, res, tmp, bw_x - 6, bw_x - 5);
    for (int i = 0; i < dim; i++) {
        gethan4[i] = tmp[i] & 1;
        r_convert[i] = res[i];
    }

    // r_convert used to temp store the r*res r=5/6 for [1,1.5] to [1,1.25]
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= ((17ULL << (bw_x - 2)) / 18);
    }

    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);

    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan4, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    // cout << endl;
    // cout << "res [1-1.0625]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;


    // res: bw = bw_x - 1 scale = bw_x - 2, res, bw = 31, scale = 30, the 29-th (0.5, bw - 3) bit is zero
    trunc->truncate_and_reduce(dim, res, tmp, bw_x - 7, bw_x - 6);
    for (int i = 0; i < dim; i++) {
        gethan5[i] = tmp[i] & 1;
        r_convert[i] = res[i];
    }

    // r_convert used to temp store the r*res r=5/6 for [1,1.5] to [1,1.25]
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= ((33ULL << (bw_x - 2)) / 34);
    }

    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);

    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan5, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    // cout << endl;
    // cout << "res [1-1.03125]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;

    trunc->truncate_and_reduce(dim, res, tmp, bw_x - 8, bw_x - 7);
    for (int i = 0; i < dim; i++) {
        gethan6[i] = tmp[i] & 1;
        r_convert[i] = res[i];
    }

    // r_convert used to temp store the r*res r=5/6 for [1,1.5] to [1,1.25]
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= ((65ULL << (bw_x - 2)) / 66);
    }

    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);

    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan6, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    // cout << endl;
    // cout << "res [1-1.015625]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;


    // a - (b * (x-1))
    // a = 0.9999849764044484
    // b = 0.4942084643407726
    uint64_t B = static_cast<uint64_t>(0.4942084643407726 * (pow(2.0, bw_x - 2)));
    uint64_t A = static_cast<uint64_t>(0.9999849764044484 * (pow(2.0, bw_x - 2)));
    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            res[i] -= (1ULL << (bw_x - 2));
        }
    }

    xt->z_extend(dim, res, res, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        res[i] = B * res[i];
    }
    trunc->truncate_direct(dim, res, res, bw_x - 2, 2 * bw_x);

    if (party == ALICE) {
        for (int i = 0; i < dim; i++) {
            res[i] = A - res[i];
        }
    } else {
        for (int i = 0; i < dim; i++) {
            res[i] = 0 - res[i];
        }
    }

    // cout << endl;
    // cout << "1/sqrt(x): ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;

    // scale back
    uint64_t scale1 = static_cast<uint64_t>(0.8660254037844386 * (pow(2.0, bw_x - 2)));
    uint64_t scale2 = static_cast<uint64_t>(0.9128709291752769 * (pow(2.0, bw_x - 2)));
    uint64_t scale3 = static_cast<uint64_t>(0.9486832980505138 * (pow(2.0, bw_x - 2)));
    uint64_t scale4 = static_cast<uint64_t>(0.97182531580755 * (pow(2.0, bw_x - 2)));
    uint64_t scale5 = static_cast<uint64_t>(0.985184366143778 * (pow(2.0, bw_x - 2)));
    uint64_t scale6 = static_cast<uint64_t>(0.9923953268977463 * (pow(2.0, bw_x - 2)));

    for (int i = 0; i < dim; i++) {
        r_convert[i] = res[i];
    }
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= scale6;
    }
    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);

    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan6, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }


    for (int i = 0; i < dim; i++) {
        r_convert[i] = res[i];
    }
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= scale5;
    }
    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);

    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan5, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }


    for (int i = 0; i < dim; i++) {
        r_convert[i] = res[i];
    }
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= scale4;
    }
    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);
    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan4, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }


    for (int i = 0; i < dim; i++) {
        r_convert[i] = res[i];
    }
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= scale3;
    }
    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);
    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan3, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    for (int i = 0; i < dim; i++) {
        r_convert[i] = res[i];
    }
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= scale2;
    }
    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);
    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan2, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }


    for (int i = 0; i < dim; i++) {
        r_convert[i] = res[i];
    }
    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= scale1;
    }
    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);
    // r_convert bw = bw_x - 1, scale = bw_x - 2
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan1, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    // cout << endl;
    // cout << "1/sqrt(x) [1, 2]: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << res[i] << ", ";
    // }
    // cout << endl;

    for (int i = 0; i < dim; i++) {
        r_convert[i] = res[i];
        gethan1[i] = msnzb_index[i] & 1;
    }

    xt->z_extend(dim, r_convert, r_convert, bw_x, 2 * bw_x, msbx);
    for (int i = 0; i < dim; i++) {
        r_convert[i] *= static_cast<uint64_t>(std::sqrt(2) * (pow(2.0, bw_x - 2)));;
    }
    trunc->truncate_direct(dim, r_convert, r_convert, bw_x - 2, 2 * bw_x);
    for (int i = 0; i < dim; i++) {
        tmp4[i] = r_convert[i] - res[i];
    }
    aux->multiplexer(gethan1, tmp4, tmp2, dim, bw_x, bw_x);
    for (int i = 0; i < dim; i++) {
        res[i] += tmp2[i];
    }

    int bw_sqrt_adjust = bw_x;
    int scale_sqrt_adjust = int(ceil(s_x / 2.0));
    uint64_t *sqrt_adjust = new uint64_t[dim]();
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < (bw_x - 1); j++) {
            sqrt_adjust[i] +=
                    (1ULL << (int(floor((s_x - j) / 2.0)) + scale_sqrt_adjust))
                    * msnzb_vector[i * (bw_x - 1) + j];
        }
    }

    // cout << endl;
    // cout << "sqrt_adjust: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << sqrt_adjust[i] << ", ";
    // }
    // cout << endl;

    mult->hadamard_product(dim, res, sqrt_adjust, tmp2, bw_x,
                           bw_sqrt_adjust, bw_x + bw_sqrt_adjust, false,
                           false, MultMode::None);
    // x_curr: bw = s_y + floor(bw_x/2) + 1 - ceil(s_x/2), scale = s_y
    trunc->truncate_and_reduce(dim, tmp2, y, scale_sqrt_adjust + (bw_x - 2 - s_x),
                               bw_x + bw_sqrt_adjust);


    // cout << endl;
    // cout << "final res: ";
    // for (int i = 0; i < dim; i++) {
    //     cout << y[i] << ", ";
    // }
    // cout << endl;

    delete[] tmp;
    delete[] tmp2;
    delete[] tmp4;
    delete[] gethan1;
    delete[] gethan2;
    delete[] gethan3;
    delete[] gethan4;
    delete[] gethan5;
    delete[] gethan6;
    delete[] msbx;
    delete[] res;
    delete[] r_convert;
    delete[] adjusted_x;
    delete[] sqrt_adjust;
    delete[] msnzb_vector_bool;
    delete[] msnzb_vector;
    delete[] msnzb_index;
    delete[] exp_parity;
    delete[] adjust;
    delete[] tmp_bool;
}


// uint64_t interval_endpoints[7] = {1 << 4, 1 << 8, 1 << 12, 1 << 16, 1 << 20, 1 << 24, 1 << 28};
// // the init values setting is not perfect if we only refer to the plaintext computation setting, in the MPC setting
// // the current setting gives no better approximate after iterations.
// uint64_t init_values[8] = {
//     (1 << 22) + (1 << 21), (1 << 20) + (1 << 19), (1 << 18) + (1 << 17), (1 << 16) + (1 << 15), (1 << 14) + (1 << 13),
//     (1 << 12) + (1 << 11), (1 << 10) + (1 << 9), (1 << 8) + (1 << 7)
// }; // the decimal under the scale = 16.
//
// void direct_trunc(int32_t dim, uint64_t *data, int32_t trunc_size) {
//     for (int i = 0; i < dim; i++) {
//         data[i] = data[i] >> trunc_size;
//     }
// }
//
//
// // the disadvantages are: imprecise initial value leads to more iterations to
// void MathFunctions::sqrtPSC(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x, int32_t bw_y, int32_t s_x, int32_t s_y,
//                             bool inverse) {
//     // first need the interval test
//     // confirm the interval we need test
//     // compare eq B2A compare eq B2A 32 / 4 = 8 intervals
//     auto *temp_data = new uint64_t[dim];
//     auto *temp_cmp = new uint8_t[dim];
//     auto *temp_B2A = new uint64_t[dim];
//     auto *t = new uint64_t[dim];
//     auto *init_guess = new uint64_t[dim]();
//     auto **interval_indict = new uint8_t *[8];
//     for (int i = 0; i < 8; i++) {
//         interval_indict[i] = new uint8_t[dim];
//     }
//     auto temp_init = new uint64_t *[8];
//     for (int i = 0; i < 8; i++) {
//         temp_init[i] = new uint64_t[dim];
//     }
//     memcpy(temp_data, x, dim * sizeof(uint64_t));
//
//     if (party == ALICE) {
//         // temp_cmp indicates whether Alice > Bob
//         for (int i = 0; i < 7; i++) {
//             for (int j = 0; j < dim; j++) {
//                 t[j] = temp_data[j] & ((1 << 4) - 1);
//             }
//             this->aux->mill->compare(temp_cmp, t, dim, 4, true);
//             this->aux->B2A(temp_cmp, temp_B2A, dim, bw_x - (i + 1) * 4);
//             direct_trunc(dim, temp_data, 4);
//             for (int j = 0; j < dim; j++) {
//                 temp_data[j] &= (1 << (bw_x - (i + 1) * 4)) - 1;
//                 temp_data[j] = temp_data[j] + temp_B2A[j];
//                 temp_data[j] &= (1 << (bw_x - (i + 1) * 4)) - 1;
//             }
//             for (int j = 0; j < dim; j++) {
//                 t[j] = (1 << (bw_x - (i + 1) * 4)) - temp_data[j];
//             }
//             this->equality->check_equality(temp_cmp, t, dim, bw_x - (i + 1) * 4);
//             memcpy(interval_indict[i], temp_cmp, dim * sizeof(uint8_t));
//         }
//         memcpy(interval_indict[7], temp_cmp, dim * sizeof(uint8_t));
//         // AND
//         for (int j = 6; j > 0; j--) {
//             this->aux->AND(interval_indict[j - 1], interval_indict[j], interval_indict[j], dim);
//         }
//
//         // for (int i = 0; i < 8; i++) {
//         //     for (int j = 0; j < dim; j++) {
//         //         cout << static_cast<int>(interval_indict[i][j]) << ", ";
//         //     }
//         //     cout << endl;
//         // }
//
//         for (int i = 0; i < 8; i++) {
//             this->aux->multiplexer_two_plain(interval_indict[i], init_values[i], temp_init[i], dim, bw_x, bw_y);
//         }
//
//         for (int i = 0; i < dim; i++) {
//             for (int j = 0; j < 8; j++) {
//                 init_guess[i] += temp_init[j][i];
//                 init_guess[i] &= ((1ULL << bw_y) - 1);
//             }
//         }
//
//         // for (int i = 0; i < 1; i++) {
//         //     std::cout << init_guess[i] << ", ";
//         // }
//         // std::cout << endl;
//
//         // newton iterative
//         newton_inverse_sqrt(dim, init_guess, x, bw_x, s_x);
//     } else {
//         for (int i = 0; i < 7; i++) {
//             for (int j = 0; j < dim; j++) {
//                 t[j] = (((1 << 4) - 1 - (temp_data[j] & ((1 << 4) - 1)))) & ((1 << 4) - 1);
//             }
//             this->aux->mill->compare(temp_cmp, t, dim, 4, true);
//             this->aux->B2A(temp_cmp, temp_B2A, dim, bw_x - (i + 1) * 4);
//             direct_trunc(dim, temp_data, 4);
//             for (int j = 0; j < dim; j++) {
//                 temp_data[j] &= (1 << (bw_x - (i + 1) * 4)) - 1;
//                 temp_data[j] = temp_data[j] + temp_B2A[j];
//                 temp_data[j] &= (1 << (bw_x - (i + 1) * 4)) - 1;
//             }
//             this->equality->check_equality(temp_cmp, temp_data, dim, bw_x - (i + 1) * 4);
//             memcpy(interval_indict[i], temp_cmp, dim * sizeof(uint8_t));
//         }
//         for (int i = 0; i < dim; i++) {
//             temp_cmp[i] = temp_cmp[i] ^ 1;
//         }
//         memcpy(interval_indict[7], temp_cmp, dim * sizeof(uint8_t));
//         // AND
//         for (int j = 6; j > 0; j--) {
//             for (int k = 0; k < dim; k++) {
//                 temp_cmp[k] = interval_indict[j - 1][k] ^ 1;
//             }
//             this->aux->AND(temp_cmp, interval_indict[j], interval_indict[j], dim);
//         }
//
//         // for (int i = 0; i < 8; i++) {
//         //     for (int j = 0; j < dim; j++) {
//         //         cout << static_cast<int>(interval_indict[i][j]) << ", ";
//         //     }
//         //     cout << endl;
//         // }
//
//         for (int i = 0; i < 8; i++) {
//             this->aux->multiplexer_two_plain(interval_indict[i], init_values[i], temp_init[i], dim, bw_x, bw_y);
//         }
//
//         for (int i = 0; i < dim; i++) {
//             for (int j = 0; j < 8; j++) {
//                 init_guess[i] += temp_init[j][i];
//                 init_guess[i] &= (1ULL << bw_y) - 1;
//             }
//         }
//
//         // for (int i = 0; i < 1; i++) {
//         //     std::cout << init_guess[i] << ", ";
//         // }
//         // std::cout << endl;
//
//         // newton iterative
//         newton_inverse_sqrt(dim, init_guess, x, bw_x, s_x);
//     }
//     memcpy(y, init_guess, dim * sizeof(uint64_t));
//     delete[] temp_data;
//     delete[] temp_cmp;
//     delete[] temp_B2A;
//     delete[] t;
//     delete[] init_guess;
//     for (int i = 0; i < 8; i++) {
//         delete[] interval_indict[i];
//     }
//     delete[] interval_indict;
//     for (int i = 0; i < 8; i++) {
//         delete[] temp_init[i];
//     }
//     delete[] temp_init;
// }
//
// // for the fixed key, the second number appears not to stable, for the reason of newton or the init values
// void MathFunctions::newton_inverse_sqrt(int32_t dim, uint64_t *init_guess, uint64_t *x,
//                                         int32_t bw, int32_t scale) {
//     int32_t iter_num = 5;
//     // auto *tmp = new uint64_t[dim];
//     auto *tmp1 = new uint64_t[dim];
//     auto *tmp11 = new uint64_t[dim];
//     auto *tmp2 = new uint64_t[dim];
//     auto *tmp22 = new uint64_t[dim];
//     for (int i = 0; i < iter_num; i++) {
//         this->mult->hadamard_product(dim, init_guess, x, tmp1, bw, bw, 2 * bw, true);
//         this->trunc->truncate_red_then_ext(dim, tmp1, tmp11, scale, 2 * bw);
//
//         // cout << "y: " << endl;
//         // cout << x[0] << " ";
//         // cout << endl;
//         //
//         // cout << "x*y: " << endl;
//         // cout << tmp11[0] << " ";
//         // cout << endl;
//
//         this->mult->hadamard_product(dim, init_guess, tmp11, tmp2, bw, bw, 2 * bw, true);
//         this->trunc->truncate_red_then_ext(dim, tmp2, tmp22, scale, 2 * bw);
//
//         // cout << "y * x^2: " << endl;
//         // cout << tmp22[0] << " ";
//         // cout << endl;
//
//         if (party == ALICE) {
//             for (int j = 0; j < dim; j++) {
//                 tmp22[j] = (3 << scale) - tmp22[j];
//                 tmp22[j] &= (1ULL << bw) - 1;
//             }
//         } else {
//             for (int j = 0; j < dim; j++) {
//                 tmp22[j] = 0 - tmp22[j];
//                 tmp22[j] &= (1ULL << bw) - 1;
//             }
//         }
//
//         // cout << "3 - y * x^2: " << endl;
//         // cout << tmp22[0] << " ";
//         // cout << endl;
//
//         this->mult->hadamard_product(dim, init_guess, tmp22, tmp1, bw, bw, 2 * bw, true);
//         this->trunc->truncate_red_then_ext(dim, tmp1, init_guess, scale, 2 * bw);
//         this->trunc->truncate_red_then_ext(dim, init_guess, init_guess, 1, bw);
//
//         // cout << "x_new: " << endl;
//         // cout << init_guess[0] << " ";
//         // cout << endl << endl;
//         // truncate
//     }
//     // for (int i = 0; i < 1; i++) {
//     //     cout << init_guess[i] << " ";
//     // }
//     // cout << endl;
//     delete[] tmp1;
//     delete[] tmp11;
//     delete[] tmp2;
//     delete[] tmp22;
//     // delete[] tmp;
// }
