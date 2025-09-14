/*
Authors: Deevashwer Rathee, Mayank Rathee
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

#include "BuildingBlocks/aux-protocols.h"

#include <csignal>
#include <iomanip>

#include "BuildingBlocks/truncation.h"
#include "BuildingBlocks/value-extension.h"
#include <omp.h>

using namespace std;
using namespace sci;

AuxProtocols::AuxProtocols(int party, IOPack *iopack, OTPack *otpack) {
    this->party = party;
    this->iopack = iopack;
    this->otpack = otpack;
    this->mill = new MillionaireProtocol(party, iopack, otpack);
    this->mill_and_eq = new MillionaireWithEquality(party, iopack, otpack);
}

AuxProtocols::~AuxProtocols() {
    delete mill;
    delete mill_and_eq;
}

void AuxProtocols::wrap_computation(uint64_t *x, uint8_t *y, int32_t size,
                                    int32_t bw_x) {
    assert(bw_x <= 64);
    uint64_t mask = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));

    uint64_t *tmp_x = new uint64_t[size];
    for (int i = 0; i < size; i++) {
        if (party == sci::ALICE)
            tmp_x[i] = x[i] & mask;
        else
            tmp_x[i] = (mask - x[i]) & mask; // 2^{bw_x} - 1 - x[i]
    }
    mill->compare(y, tmp_x, size, bw_x, true); // computing greater_than

    delete[] tmp_x;
}

void AuxProtocols::multiplexer(uint8_t *sel, uint64_t *x, uint64_t *y,
                               int32_t size, int32_t bw_x, int32_t bw_y) {
    assert(bw_x <= 64 && bw_y <= 64 && bw_y <= bw_x);
    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
    uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));

    uint64_t *corr_data = new uint64_t[size];
    uint64_t *data_S = new uint64_t[size];
    uint64_t *data_R = new uint64_t[size];

    // y = (sel_0 \xor sel_1) * (x_0 + x_1)
    // y = (sel_0 + sel_1 - 2*sel_0*sel_1)*x_0 + (sel_0 + sel_1 -
    // 2*sel_0*sel_1)*x_1 y = [sel_0*x_0 + sel_1*(x_0 - 2*sel_0*x_0)]
    //     + [sel_1*x_1 + sel_0*(x_1 - 2*sel_1*x_1)]
    for (int i = 0; i < size; i++) {
        corr_data[i] = (x[i] * (1 - 2 * uint64_t(sel[i]))) & mask_y;
    }
#pragma omp parallel num_threads(2)
    {
        if (omp_get_thread_num() == 1) {
            if (party == sci::ALICE) {
                otpack->iknp_reversed->recv_cot(data_R, (bool *) sel, size, bw_y);
            } else {
                // party == sci::BOB
                otpack->iknp_reversed->send_cot(data_S, corr_data, size, bw_y);
            }
        } else {
            if (party == sci::ALICE) {
                otpack->iknp_straight->send_cot(data_S, corr_data, size, bw_y);
            } else {
                // party == sci::BOB
                otpack->iknp_straight->recv_cot(data_R, (bool *) sel, size, bw_y);
            }
        }
    }
    for (int i = 0; i < size; i++) {
        y[i] = ((x[i] * uint64_t(sel[i]) + data_R[i] - data_S[i]) & mask_y);
    }

    delete[] corr_data;
    delete[] data_S;
    delete[] data_R;
}

void AuxProtocols::B2A(uint8_t *x, uint64_t *y, int32_t size, int32_t bw_y) {
    assert(bw_y <= 64 && bw_y >= 1);
    if (bw_y == 1) {
        for (int i = 0; i < size; i++) {
            y[i] = uint64_t(x[i]) & 1;
        }
        return;
    }
    uint64_t mask = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));

    if (party == sci::ALICE) {
        uint64_t *corr_data = new uint64_t[size];
        for (int i = 0; i < size; i++) {
            corr_data[i] = (-2 * uint64_t(x[i])) & mask;
        }
        otpack->iknp_straight->send_cot(y, corr_data, size, bw_y);
        for (int i = 0; i < size; i++) {
            y[i] = (uint64_t(x[i]) - y[i]) & mask;
        }
        delete[] corr_data;
    } else {
        // party == sci::BOB
        otpack->iknp_straight->recv_cot(y, (bool *) x, size, bw_y);
        for (int i = 0; i < size; i++) {
            y[i] = (uint64_t(x[i]) + y[i]) & mask;
        }
    }
}

void AuxProtocols::B2A_coprimes(uint8_t *x, uint64_t *y, int32_t size, vector<uint64_t> fragment_modulus) {
    uint64_t num_fragments = fragment_modulus.size();
    std::vector<int> msg_len(num_fragments); // the length I want is log fragment_modulus
    for (int i = 0; i < num_fragments; i++) {
        msg_len[i] = static_cast<int>(ceil(log2(fragment_modulus[i] + 1)));
    }
    if (party == sci::ALICE) {
        uint64_t *corr_data = new uint64_t[size * num_fragments];
        for (int i = 0; i < size * num_fragments; i += num_fragments) {
            for (int j = 0; j < num_fragments; j++) {
                corr_data[i + j] = (fragment_modulus[j] - 2 * uint64_t(x[i / num_fragments])) % fragment_modulus[j];
            }
        }
        otpack->iknp_straight->send_batched_cot_special(y, corr_data, msg_len, fragment_modulus, size, num_fragments);

        for (int i = 0; i < size * num_fragments; i += num_fragments) {
            for (int j = 0; j < num_fragments; j++) {
                // cout<< "A set corr " << corr_data[i + j] << " ";
                // cout<< "A rece H0 " << fragment_modulus[j] - (y[i + j] % fragment_modulus[j]) << " ";
                y[i + j] = (uint64_t(x[i / num_fragments]) + fragment_modulus[j] - (y[i + j] % fragment_modulus[j])) %
                           fragment_modulus[j];
            }
            // cout << endl;
        }
        delete[] corr_data;
    } else {
        // party == sci::BOB
        otpack->iknp_straight->recv_batched_cot_special(y, (bool *) x, msg_len, size, num_fragments);
        for (int i = 0; i < size * num_fragments; i += num_fragments) {
            for (int j = 0; j < num_fragments; j++) {
                // cout << "choice: " << uint64_t(x[i / num_fragments]) << " ";
                // cout << "B rece y : " << (y[i + j] % fragment_modulus[j]) << " ";
                y[i + j] = (uint64_t(x[i / num_fragments]) + (y[i + j] % fragment_modulus[j]) + fragment_modulus[j]) %
                           fragment_modulus[j];
            }
            // cout << endl;
        }
    }
}

template<typename T>
void AuxProtocols::lookup_table(T **spec, T *x, T *y, int32_t size,
                                int32_t bw_x, int32_t bw_y) {
    if (party == sci::ALICE) {
        assert(x == nullptr);
        assert(y == nullptr);
    } else {
        // party == sci::BOB
        assert(spec == nullptr);
    }
    assert(bw_x <= 8 && bw_x >= 1);
    int32_t T_size = sizeof(T) * 8;
    assert(bw_y <= T_size);

    T mask_x = (bw_x == T_size ? -1 : ((1ULL << bw_x) - 1));
    T mask_y = (bw_y == T_size ? -1 : ((1ULL << bw_y) - 1));
    uint64_t N = 1 << bw_x;

    if (party == sci::ALICE) {
        PRG128 prg;
        T **data = new T *[size];
        for (int i = 0; i < size; i++) {
            data[i] = new T[N];
            for (uint64_t j = 0; j < N; j++) {
                data[i][j] = spec[i][j];
            }
        }

        otpack->kkot[bw_x - 1]->send(data, size, bw_y);

        for (int i = 0; i < size; i++)
            delete[] data[i];
        delete[] data;
    } else {
        // party == sci::BOB
        uint8_t *choice = new uint8_t[size];
        for (int i = 0; i < size; i++) {
            choice[i] = x[i] & mask_x;
        }
        otpack->kkot[bw_x - 1]->recv(y, choice, size, bw_y);

        delete[] choice;
    }
}

void AuxProtocols::MSB(uint64_t *x, uint8_t *msb_x, int32_t size,
                       int32_t bw_x) {
    assert(bw_x <= 64);
    int32_t shift = bw_x - 1;
    uint64_t shift_mask = (shift == 64 ? -1 : ((1ULL << shift) - 1));

    uint64_t *tmp_x = new uint64_t[size];
    uint8_t *msb_xb = new uint8_t[size];
    for (int i = 0; i < size; i++) {
        tmp_x[i] = x[i] & shift_mask;
        msb_xb[i] = (x[i] >> shift) & 1;
        if (party == sci::BOB)
            tmp_x[i] = (shift_mask - tmp_x[i]) & shift_mask;
    }

    mill->compare(msb_x, tmp_x, size, bw_x - 1, true); // computing greater_than

    for (int i = 0; i < size; i++) {
        msb_x[i] = msb_x[i] ^ msb_xb[i];
    }

    delete[] tmp_x;
    delete[] msb_xb;
}

void AuxProtocols::MSB_to_Wrap(uint64_t *x, uint8_t *msb_x, uint8_t *wrap_x,
                               int32_t size, int32_t bw_x) {
    assert(bw_x <= 64);
    if (party == sci::ALICE) {
        PRG128 prg;
        prg.random_bool((bool *) wrap_x, size);
        uint8_t **spec = new uint8_t *[size];
        for (int i = 0; i < size; i++) {
            spec[i] = new uint8_t[4];
            uint8_t msb_xb = (x[i] >> (bw_x - 1)) & 1;
            for (int j = 0; j < 4; j++) {
                uint8_t bits_j[2]; // j0 || j1 (LSB to MSB)
                uint8_to_bool(bits_j, j, 2);
                spec[i][j] = (((1 ^ msb_x[i] ^ bits_j[0]) * (msb_xb ^ bits_j[1])) ^
                              (msb_xb * bits_j[1]) ^ wrap_x[i]) &
                             1;
            }
        }
        lookup_table<uint8_t>(spec, nullptr, nullptr, size, 2, 1);

        for (int i = 0; i < size; i++)
            delete[] spec[i];
        delete[] spec;
    } else {
        // party == sci::BOB
        uint8_t *lut_in = new uint8_t[size];
        for (int i = 0; i < size; i++) {
            lut_in[i] = (((x[i] >> (bw_x - 1)) & 1) << 1) | msb_x[i];
        }
        lookup_table<uint8_t>(nullptr, lut_in, wrap_x, size, 2, 1);

        delete[] lut_in;
    }
}

void AuxProtocols::AND(uint8_t *x, uint8_t *y, uint8_t *z, int32_t size) {
    int32_t old_size = size;
    size = ((old_size + 1) / 2) * 2;
    uint8_t *a = new uint8_t[size];
    uint8_t *b = new uint8_t[size];
    uint8_t *c = new uint8_t[size];
    a[size - 1] = 0;
    b[size - 1] = 0; // if size is odd, last element should be set
    memcpy(a, x, sizeof(uint8_t) * old_size);
    memcpy(b, y, sizeof(uint8_t) * old_size);
    switch (party) {
        case sci::ALICE: {
            PRG128 prg;
            prg.random_bool((bool *) c, size);
            uint8_t **ot_messages; // (size/2) X 16
            ot_messages = new uint8_t *[size / 2];
            for (int i = 0; i < size; i += 2)
                ot_messages[i / 2] = new uint8_t[16];
            for (int j = 0; j < 16; j++) {
                uint8_t bits_j[4]; // a01 || b01 || a11 || b11 (LSB->MSB)
                sci::uint8_to_bool(bits_j, j, 4);
                for (int i = 0; i < size; i += 2) {
                    ot_messages[i / 2][j] =
                            ((((a[i + 1] ^ bits_j[2]) & (b[i + 1] ^ bits_j[3])) ^ c[i + 1])
                             << 1) |
                            (((a[i] ^ bits_j[0]) & (b[i] ^ bits_j[1])) ^ c[i]);
                }
            }
            // otpack->kkot_16->send(ot_messages, size/2, 2);
            otpack->kkot[3]->send(ot_messages, size / 2, 2);
            for (int i = 0; i < size; i += 2)
                delete[] ot_messages[i / 2];
            delete[] ot_messages;
            break;
        }
        case sci::BOB: {
            uint8_t *ot_selection = new uint8_t[(size_t) size / 2];
            uint8_t *ot_result = new uint8_t[(size_t) size / 2];
            for (int i = 0; i < size; i += 2) {
                ot_selection[i / 2] =
                        (b[i + 1] << 3) | (a[i + 1] << 2) | (b[i] << 1) | a[i];
            }
            // otpack->kkot_16->recv(ot_result, ot_selection, size/2, 2);
            otpack->kkot[3]->recv(ot_result, ot_selection, size / 2, 2);
            for (int i = 0; i < size; i += 2) {
                c[i] = ot_result[i / 2] & 1;
                c[i + 1] = ot_result[i / 2] >> 1;
            }
            delete[] ot_selection;
            delete[] ot_result;
            break;
        }
    }
    memcpy(z, c, sizeof(uint8_t) * old_size);
    delete[] a;
    delete[] b;
    delete[] c;
    return;
}

void AuxProtocols::reduce(int32_t dim, uint64_t *x, uint64_t *y, int32_t bw_x,
                          int32_t bw_y) {
    assert(bw_y <= bw_x);
    uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));

    for (int i = 0; i < dim; i++) {
        y[i] = x[i] & mask_y;
    }
}

void AuxProtocols::digit_decomposition(int32_t dim, uint64_t *x,
                                       uint64_t *x_digits, int32_t bw_x,
                                       int32_t digit_size) {
    assert(false && "Inefficient version of digit decomposition called");
    int num_digits = ceil(double(bw_x) / digit_size);
    int last_digit_size = bw_x - (num_digits - 1) * digit_size;
    uint64_t digit_mask = (digit_size == 64 ? -1 : (1ULL << digit_size) - 1);
    uint64_t last_digit_mask =
            (last_digit_size == 64 ? -1 : (1ULL << last_digit_size) - 1);

    Truncation trunc(this->party, this->iopack, this->otpack);
    for (int i = 0; i < num_digits; i++) {
        trunc.truncate_and_reduce(dim, x, x_digits + i * dim, i * digit_size, bw_x);
        uint64_t mask = (i == (num_digits - 1) ? last_digit_mask : digit_mask);
        for (int j = 0; j < dim; j++) {
            x_digits[i * dim + j] &= mask;
        }
    }
}

void AuxProtocols::digit_decomposition_sci(int32_t dim, uint64_t *x,
                                           uint64_t *x_digits, int32_t bw_x,
                                           int32_t digit_size,
                                           bool all_digit_size) {
    int num_digits = ceil(double(bw_x) / digit_size);
    int last_digit_size = bw_x - (num_digits - 1) * digit_size;
    uint64_t digit_mask = (digit_size == 64 ? -1 : (1ULL << digit_size) - 1);
    uint64_t last_digit_mask =
            (last_digit_size == 64 ? -1 : (1ULL << last_digit_size) - 1);
    for (int i = 0; i < num_digits; i++) {
        for (int j = 0; j < dim; j++) {
            x_digits[i * dim + j] = (x[j] >> (i * digit_size));
            x_digits[i * dim + j] &=
                    (i == (num_digits - 1)) ? last_digit_mask : digit_mask;
        }
    }
    uint8_t *wrap_ = new uint8_t[dim * (num_digits - 1)];
    uint8_t *ones_ = new uint8_t[dim * (num_digits - 1)];
    uint8_t *dp_wrap_entering = new uint8_t[dim * num_digits];
    uint8_t *dp_temp = new uint8_t[dim * num_digits];
    uint64_t *dp_wrap_arith = new uint64_t[dim * num_digits];
    // Fill wrap_ and ones_
    uint64_t *temp_x_digits = new uint64_t[dim * (num_digits - 1)];

    for (int i = 0; i < (num_digits - 1); i++) {
        for (int j = 0; j < dim; j++) {
            if (party == sci::ALICE)
                temp_x_digits[i * dim + j] = x_digits[i * dim + j] & digit_mask;
            else
                temp_x_digits[i * dim + j] =
                        (digit_mask - x_digits[i * dim + j]) & digit_mask;
        }
    }
    this->mill_and_eq->compare_with_eq(wrap_, ones_, temp_x_digits,
                                       (dim * (num_digits - 1)), digit_size);

    // DP steps proceed
    for (int i = 0; i < num_digits; i++) {
        if (i > 0) {
            this->AND(ones_ + (i - 1) * dim, dp_wrap_entering + (i - 1) * dim,
                      dp_temp + (i - 1) * dim, dim);
        }
        for (int j = 0; j < dim; j++) {
            if (i == 0) {
                dp_wrap_entering[i * dim + j] = 0;
            } else {
                dp_wrap_entering[i * dim + j] =
                        wrap_[(i - 1) * dim + j] ^ dp_temp[(i - 1) * dim + j];
            }
        }
    }
    this->B2A(dp_wrap_entering, dp_wrap_arith, num_digits * dim, digit_size);
    for (int i = 0; i < num_digits; i++) {
        for (int j = 0; j < dim; j++) {
            x_digits[i * dim + j] += dp_wrap_arith[i * dim + j];
            uint64_t temp_mask =
                    (i == (num_digits - 1)) ? last_digit_mask : digit_mask;
            x_digits[i * dim + j] &= temp_mask;
        }
        if (all_digit_size) {
            if (i == (num_digits - 1)) {
                XTProtocol *xt =
                        new XTProtocol(this->party, this->iopack, this->otpack);
                uint64_t *temp_last_digs = new uint64_t[dim];
                xt->z_extend(dim, x_digits + (num_digits - 1) * dim, temp_last_digs,
                             last_digit_size, digit_size);
                for (int j = 0; j < dim; j++) {
                    x_digits[i * dim + j] = temp_last_digs[j];
                    x_digits[i * dim + j] &= digit_mask;
                }
                delete xt;
                delete[] temp_last_digs;
            }
        }
    }

    delete[] wrap_;
    delete[] ones_;
    delete[] dp_wrap_entering;
    delete[] dp_temp;
    delete[] dp_wrap_arith;
    delete[] temp_x_digits;
}

uint64_t lookup_msnzb(uint64_t index) {
    uint64_t ret = 0ULL;
    ret = floor(log2(index));
    if (index == 0) {
        ret = 0ULL;
    }
    // In the above step only at max log(64) = 6 bits are filled.
    ret <<= 1;
    // Last bit stores 1 if index is 0, else 0.
    if (index == 0) {
        ret ^= 1ULL;
    }
    return ret;
}

void AuxProtocols::msnzb_sci(uint64_t *x, uint64_t *msnzb_index, int32_t bw_x,
                             int32_t size, int32_t digit_size) {
    // The protocol only works when digit_size divides bw_x.
    int32_t last_digit_size = bw_x % digit_size;
    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
    uint64_t digit_mask = (digit_size == 64 ? -1 : ((1ULL << digit_size) - 1));
    uint64_t last_digit_mask =
            (last_digit_size == 64 ? -1 : ((1ULL << last_digit_size) - 1));
    if (last_digit_size == 0) {
        last_digit_mask = digit_mask;
        last_digit_size = digit_size;
    }
    int32_t num_digits = ceil((bw_x * 1.0) / digit_size);
    uint64_t *x_digits = new uint64_t[num_digits * size];

    XTProtocol *xt = new XTProtocol(this->party, this->iopack, this->otpack);

    // Extract digits
    this->digit_decomposition_sci(size, x, x_digits, bw_x, digit_size);

    // Use LUTs for MSNZB on digits
    int D = (1 << digit_size);
    int DLast = (1 << last_digit_size);
    uint8_t *z_ = new uint8_t[num_digits * size];
    uint64_t *msnzb_ = new uint64_t[num_digits * size];
    uint64_t *msnzb_extended = new uint64_t[num_digits * size];
    int lookup_output_bits = (ceil(log2(digit_size))) + 1;
    int mux_bits = ceil(log2(bw_x));
    uint64_t msnzb_mask = (1ULL << (lookup_output_bits - 1)) - 1;
    uint64_t mux_mask = (1ULL << mux_bits) - 1;
    if (party == ALICE) {
        uint64_t **spec;
        spec = new uint64_t *[num_digits * size];
        PRG128 prg;
        prg.random_data(z_, size * sizeof(uint8_t));
        prg.random_data(msnzb_, size * sizeof(uint64_t));
        for (int i = 0; i < (num_digits - 1) * size; i++) {
            spec[i] = new uint64_t[D];
            z_[i] &= 1;
            msnzb_[i] &= msnzb_mask;
            for (int j = 0; j < D; j++) {
                int idx = (x_digits[i] + j) & digit_mask;
                uint64_t lookup_val = lookup_msnzb(idx);
                spec[i][j] = ((lookup_val >> 1) - msnzb_[i]) & msnzb_mask;
                spec[i][j] <<= 1;
                spec[i][j] |=
                        ((uint64_t) (((uint8_t) (lookup_val & 1ULL)) ^ z_[i]) & 1ULL);
            }
        }
        for (int i = (num_digits - 1) * size; i < num_digits * size; i++) {
            spec[i] = new uint64_t[DLast];
            z_[i] &= 1;
            msnzb_[i] &= msnzb_mask;
            for (int j = 0; j < DLast; j++) {
                int idx = (x_digits[i] + j) & last_digit_mask;
                uint64_t lookup_val = lookup_msnzb(idx);
                spec[i][j] = ((lookup_val >> 1) - msnzb_[i]) & msnzb_mask;
                spec[i][j] <<= 1;
                spec[i][j] |=
                        ((uint64_t) (((uint8_t) (lookup_val & 1ULL)) ^ z_[i]) & 1ULL);
            }
        }
        if (last_digit_size == digit_size) {
            this->lookup_table<uint64_t>(spec, nullptr, nullptr, num_digits * size,
                                         digit_size, lookup_output_bits);
        } else {
            this->lookup_table<uint64_t>(spec, nullptr, nullptr,
                                         (num_digits - 1) * size, digit_size,
                                         lookup_output_bits);
            this->lookup_table<uint64_t>(spec + (num_digits - 1) * size, nullptr,
                                         nullptr, size, last_digit_size,
                                         lookup_output_bits);
        }

        // Zero extend to mux_bits
        xt->z_extend(num_digits * size, msnzb_, msnzb_extended,
                     lookup_output_bits - 1, mux_bits);

        for (int i = 0; i < num_digits * size; i++) {
            delete[] spec[i];
        }
        delete[] spec;
    } else {
        // BOB
        if (last_digit_size == digit_size) {
            this->lookup_table<uint64_t>(nullptr, x_digits, msnzb_, num_digits * size,
                                         digit_size, lookup_output_bits);
        } else {
            this->lookup_table<uint64_t>(nullptr, x_digits, msnzb_,
                                         (num_digits - 1) * size, digit_size,
                                         lookup_output_bits);
            this->lookup_table<uint64_t>(nullptr, x_digits + (num_digits - 1) * size,
                                         msnzb_ + (num_digits - 1) * size, size,
                                         last_digit_size, lookup_output_bits);
        }

        for (int i = 0; i < (num_digits * size); i++) {
            z_[i] = (uint8_t) (msnzb_[i] & 1ULL);
            msnzb_[i] >>= 1;
        }

        // Zero extend to mux_bits
        xt->z_extend(num_digits * size, msnzb_, msnzb_extended,
                     lookup_output_bits - 1, mux_bits);

        for (int i = 0; i < num_digits; i++) {
            for (int j = 0; j < size; j++) {
                msnzb_extended[i * size + j] += (i * digit_size);
                msnzb_extended[i * size + j] &= mux_mask;
            }
        }
    }

    // Combine MSNZB of digits
    uint8_t *dp_zeros_ = new uint8_t[(num_digits - 1) * size];
    uint8_t *one_xor_zeros_ = new uint8_t[(num_digits - 1) * size];
    uint8_t *dp_zeros_final = new uint8_t[num_digits * size];

    if (party == ALICE) {
        for (int i = 0; i < size; i++) {
            dp_zeros_final[(num_digits - 1) * size + i] =
                    z_[(num_digits - 1) * size + i];
        }
        for (int i = 0; i < (num_digits - 1); i++) {
            for (int j = 0; j < size; j++) {
                one_xor_zeros_[i * size + j] = z_[i * size + j];
            }
        }
    } else {
        for (int i = 0; i < size; i++) {
            dp_zeros_final[(num_digits - 1) * size + i] =
                    (1 ^ z_[(num_digits - 1) * size + i]);
        }
        for (int i = 0; i < (num_digits - 1); i++) {
            for (int j = 0; j < size; j++) {
                one_xor_zeros_[i * size + j] = (1 ^ z_[i * size + j]);
            }
        }
    }
    for (int i = (num_digits - 2); i >= 0; i--) {
        if (i == (num_digits - 2)) {
            for (int j = 0; j < size; j++) {
                dp_zeros_[i * size + j] = z_[(i + 1) * size + j];
            }
        } else {
            this->AND(dp_zeros_ + (i + 1) * size, z_ + (i + 1) * size,
                      dp_zeros_ + i * size, size);
        }
    }
    this->AND(dp_zeros_, one_xor_zeros_, dp_zeros_final, (num_digits - 1) * size);

    uint64_t *msnzb_muxed = new uint64_t[num_digits * size];
    this->multiplexer(dp_zeros_final, msnzb_extended, msnzb_muxed,
                      num_digits * size, mux_bits, mux_bits);

    for (int i = 0; i < size; i++) {
        msnzb_index[i] = 0ULL;
        for (int j = 0; j < num_digits; j++) {
            msnzb_index[i] += msnzb_muxed[j * size + i];
            msnzb_index[i] &= mux_mask;
        }
    }

    delete xt;
    delete[] x_digits;
    delete[] z_;
    delete[] msnzb_;
    delete[] msnzb_extended;
    delete[] dp_zeros_;
    delete[] one_xor_zeros_;
    delete[] dp_zeros_final;
    delete[] msnzb_muxed;
    return;
}

void AuxProtocols::msnzb_one_hot(uint64_t *x, uint8_t *one_hot_vector,
                                 int32_t bw_x, int32_t size,
                                 int32_t digit_size) {
    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
    int msnzb_index_bits = ceil(log2(bw_x));
    uint64_t msnzb_index_mask = (1ULL << msnzb_index_bits) - 1;

    uint64_t *msnzb_index = new uint64_t[size];

    this->msnzb_sci(x, msnzb_index, bw_x, size, digit_size);

    // cout << endl;
    // cout << "msnzb_index " << endl;
    // for (int i = 0; i < size; i++) {
    //     cout << msnzb_index[i] << " ";
    // }
    // cout << endl;

    // use LUT to get the one-hot representation
    int D = 1 << msnzb_index_bits;
    uint64_t *xor_mask = new uint64_t[size];
    if (party == ALICE) {
        uint64_t **spec;
        spec = new uint64_t *[size];
        PRG128 prg;
        prg.random_data(one_hot_vector, size * bw_x * sizeof(uint8_t));
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < bw_x; j++) {
                one_hot_vector[i * bw_x + j] &= 1;
            }
            xor_mask[i] = 0ULL;
            for (int j = (bw_x - 1); j >= 0; j--) {
                xor_mask[i] <<= 1;
                xor_mask[i] ^= (uint64_t) one_hot_vector[i * bw_x + j];
            }
        }
        for (int i = 0; i < size; i++) {
            spec[i] = new uint64_t[D];
            for (int j = 0; j < D; j++) {
                int idx = (msnzb_index[i] + j) & msnzb_index_mask;
                uint64_t lookup_val = (1ULL << idx);
                lookup_val ^= xor_mask[i];
                spec[i][j] = lookup_val;
            }
        }
        this->lookup_table<uint64_t>(spec, nullptr, nullptr, size, msnzb_index_bits,
                                     bw_x);

        for (int i = 0; i < size; i++) {
            delete[] spec[i];
        }
        delete[] spec;
    } else {
        // BOB
        uint64_t *temp = new uint64_t[size];
        this->lookup_table<uint64_t>(nullptr, msnzb_index, temp, size,
                                     msnzb_index_bits, bw_x);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < bw_x; j++) {
                one_hot_vector[i * bw_x + j] = (uint8_t) (temp[i] & 1ULL);
                temp[i] >>= 1;
            }
        }
        delete[] temp;
    }
    delete[] xor_mask;
    delete[] msnzb_index;
}

// add the new function that output both the index and the one_hot vector(use the uniShare for the one-hot)
// single instance is better.
// void AuxProtocols::msnzb_one_hot_index(uint64_t *x, uint8_t *one_hot_vector, uint64_t *msnzb_index,
//                                        int32_t bw_x, int32_t size, int32_t digit_size) {
//     this->msnzb_sci(x, msnzb_index, bw_x, size, digit_size);
//     int t = 1ULL << static_cast<int>(ceil(log2(bw_x)));
//     uint8_t *tmp = new uint8_t[t];
//     for (int i = 0; i < size; i++) {
//         uniShare_naive_bool(tmp, t, msnzb_index[i]);
//         memcpy(one_hot_vector + i * (bw_x - 1), tmp, (bw_x - 1) * sizeof(uint8_t));
//     }
// }

void AuxProtocols::msnzb_one_hot_index(uint64_t *x, uint8_t *one_hot_vector, uint64_t *msnzb_index,
                                       int32_t bw_x, int32_t size, int32_t digit_size) {
    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
    int msnzb_index_bits = ceil(log2(bw_x));
    uint64_t msnzb_index_mask = (1ULL << msnzb_index_bits) - 1;

    this->msnzb_sci(x, msnzb_index, bw_x, size, digit_size);

    // use LUT to get the one-hot representation
    int D = 1 << msnzb_index_bits;
    uint64_t *xor_mask = new uint64_t[size];
    if (party == ALICE) {
        uint64_t **spec;
        spec = new uint64_t *[size];
        PRG128 prg;
        prg.random_data(one_hot_vector, size * bw_x * sizeof(uint8_t));
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < bw_x; j++) {
                one_hot_vector[i * bw_x + j] &= 1;
            }
            xor_mask[i] = 0ULL;
            for (int j = (bw_x - 1); j >= 0; j--) {
                xor_mask[i] <<= 1;
                xor_mask[i] ^= (uint64_t) one_hot_vector[i * bw_x + j];
            }
        }
        for (int i = 0; i < size; i++) {
            spec[i] = new uint64_t[D];
            for (int j = 0; j < D; j++) {
                int idx = (msnzb_index[i] + j) & msnzb_index_mask;
                uint64_t lookup_val = (1ULL << idx);
                lookup_val ^= xor_mask[i];
                spec[i][j] = lookup_val;
            }
        }
        this->lookup_table<uint64_t>(spec, nullptr, nullptr, size, msnzb_index_bits,
                                     bw_x);

        for (int i = 0; i < size; i++) {
            delete[] spec[i];
        }
        delete[] spec;
    } else {
        // BOB
        uint64_t *temp = new uint64_t[size];
        this->lookup_table<uint64_t>(nullptr, msnzb_index, temp, size,
                                     msnzb_index_bits, bw_x);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < bw_x; j++) {
                one_hot_vector[i * bw_x + j] = (uint8_t) (temp[i] & 1ULL);
                temp[i] >>= 1;
            }
        }
        delete[] temp;
    }
    delete[] xor_mask;
}

void AuxProtocols::msnzb_GC(uint64_t *x, uint8_t *one_hot_vector, int32_t bw_x,
                            int32_t size) {
    int batch_size = (bw_x * size < (1 << 20)) ? bw_x * size : (1 << 20);
    SemiHonestParty<NetIO> *prot_exec =
            setup_semi_honest<NetIO>(iopack->io_GC, party, batch_size);
    SplitIKNP<NetIO> *iknp_s_base = otpack->iknp_straight;
    if (party == ALICE) {
        static_cast<SemiHonestGen<NetIO> *>(prot_exec)->setup_keys(iknp_s_base->k0,
                                                                   iknp_s_base->s);
    } else {
        static_cast<SemiHonestEva<NetIO> *>(prot_exec)->setup_keys(iknp_s_base->k0,
                                                                   iknp_s_base->k1);
    }

    Integer *x0_bits = new Integer[size];
    Integer *x1_bits = new Integer[size];
    Bit *msnzb_vector = new Bit[size * bw_x];
    Bit *out_mask0 = new Bit[size * bw_x];
    Bit *out = new Bit[size * bw_x];
    if (party == ALICE) {
        PRG128 prg;
        prg.random_bool((bool *) one_hot_vector, size * bw_x);
    }

    for (int i = 0; i < size; ++i) {
        if (party == ALICE) {
            x0_bits[i] = Integer(bw_x, x[i], ALICE);
            x1_bits[i] = Integer(bw_x, 0, BOB);
        } else {
            // party == BOB
            x0_bits[i] = Integer(bw_x, 0, ALICE);
            x1_bits[i] = Integer(bw_x, x[i], BOB);
        }
    }
    for (int i = 0; i < size * bw_x; ++i) {
        if (party == ALICE) {
            out_mask0[i] = Bit(one_hot_vector[i], ALICE);
        } else {
            // party == BOB
            out_mask0[i] = Bit(0, ALICE);
        }
    }

    for (int i = 0; i < size; i++) {
        Integer x0 = x0_bits[i];
        Integer x1 = x1_bits[i];
        Integer X = x0 + x1;

        msnzb_vector[i * bw_x] = X[bw_x - 1];
        for (int j = 1; j < bw_x; j++) {
            msnzb_vector[i * bw_x + j] =
                    (msnzb_vector[i * bw_x + (j - 1)] | X[bw_x - j - 1]);
        }
        for (int j = 1; j < bw_x; j++) {
            msnzb_vector[i * bw_x + bw_x - j] =
            (msnzb_vector[i * bw_x + bw_x - j] ^
             msnzb_vector[i * bw_x + bw_x - j - 1]);
        }
        for (int j = 0; j < bw_x; j++) {
            out[i * bw_x + j] =
                    (msnzb_vector[i * bw_x + bw_x - j - 1] ^ out_mask0[i * bw_x + j]);
        }
    }

    for (int i = 0; i < bw_x * size; i++) {
        if (party == ALICE) {
            out[i].reveal<bool>(BOB);
        } else {
            // party == BOB
            one_hot_vector[i] = out[i].reveal<bool>(BOB);
        }
    }
    iopack->io_GC->flush();

    delete[] x0_bits;
    delete[] x1_bits;
    delete[] msnzb_vector;
    delete[] out_mask0;
    delete[] out;

    delete circ_exec;
    delete prot_exec;

    return;
}

template void AuxProtocols::lookup_table(uint64_t **spec, uint64_t *x,
                                         uint64_t *y, int32_t size,
                                         int32_t bw_x, int32_t bw_y);

template void AuxProtocols::lookup_table(uint8_t **spec, uint8_t *x, uint8_t *y,
                                         int32_t size, int32_t bw_x,
                                         int32_t bw_y);


// added functions

void print_block(const block128 &b) {
    uint8_t bytes[16];
    _mm_storeu_si128((__m128i *) bytes, b); // copy to byte array
    for (int i = 0; i < 16; ++i)
        std::cout << std::hex << std::setw(2) << std::setfill('0') << (int) bytes[i];
    std::cout << std::endl;
}

// GGM_leaves: the result of (n-1, n)-OT
// length: the length of GGM_leaves
// offset: the punctured index of GGM_leaves, which is set by one of the parties.
void AuxProtocols::nMinus1OUTNOT(block128 *seed, uint64_t length, uint64_t offset, uint8_t thread_used = 0) {
    const int tree_depth = static_cast<int>(std::ceil(std::log2(length)));
    const int leaves_size = static_cast<int>(std::pow(2, tree_depth));
    auto *GGM_leaves = new block128[leaves_size]();
    PRG128 prg;
    // std::cout << "leaves_size: " << leaves_size << std::endl;
    if (party == sci::ALICE) {
        block128 my_seed;
        prg.random_block(&my_seed, 1);
        auto **layer_OT_sum_seed = new block128 *[2];
        layer_OT_sum_seed[0] = new block128[tree_depth];
        layer_OT_sum_seed[1] = new block128[tree_depth];
        // build tree
        GGM_leaves[0] = my_seed;
        auto *temp_data = new block128[2];
        // parameters for layer0
        uint32_t points_this_layer = 1;
        uint32_t step_size = leaves_size / points_this_layer;
        uint32_t step_size_last_layer = 2 * step_size;
        uint32_t points_last_layer = points_this_layer / 2;
        for (int i = 1; i <= tree_depth; i++) {
            points_this_layer *= 2;
            step_size = leaves_size / points_this_layer;
            step_size_last_layer = 2 * step_size;
            points_last_layer = points_this_layer / 2;
            for (int j = 0; j < points_last_layer; j++) {
                // see through the nodes in each layer
                prg.reseed(&GGM_leaves[j * step_size_last_layer]);
                prg.random_block(temp_data, 2);
                GGM_leaves[j * step_size_last_layer] = temp_data[0];
                GGM_leaves[j * step_size_last_layer + step_size] = temp_data[1];
            }
            layer_OT_sum_seed[0][i - 1] = GGM_leaves[0];
            layer_OT_sum_seed[1][i - 1] = GGM_leaves[step_size];
            for (int k = 2; k < points_this_layer; k += 2) {
                layer_OT_sum_seed[0][i - 1] = _mm_xor_si128(layer_OT_sum_seed[0][i - 1],
                                                            GGM_leaves[k * step_size]);
                layer_OT_sum_seed[1][i - 1] = _mm_xor_si128(layer_OT_sum_seed[1][i - 1],
                                                            GGM_leaves[(k + 1) * step_size]);
            }
            // #ifndef NDEBUG
            //             std::cout << std::endl;
            //             for (int ii = 0; ii < leaves_size; ii++) {
            //                 print_block(GGM_leaves[ii]);
            //             }
            //             std::cout << std::endl;
            // #endif
        }

        delete[] temp_data;
        if (thread_used == 0) {
            otpack->iknp_straight->send(layer_OT_sum_seed[0], layer_OT_sum_seed[1], tree_depth);
        } else {
            otpack->iknp_reversed->send(layer_OT_sum_seed[0], layer_OT_sum_seed[1], tree_depth);
        }
        // #ifndef NDEBUG
        //         std::cout << std::endl;
        //         for (int i = 0; i < leaves_size; i++) {
        //             print_block(GGM_leaves[i]);
        //         }
        //         std::cout << std::endl;
        // #endif
        for (int i = 0; i < 2; i++) {
            delete[] layer_OT_sum_seed[i];
        }
        delete[] layer_OT_sum_seed;
    } else {
        auto r = new bool[tree_depth];
        auto *layer_OT_sum_seed = new block128[tree_depth];
        for (int i = 1; i <= tree_depth; i++) {
            r[i - 1] = (offset >> (tree_depth - i) & 1) ^ 1; // extract i-th bit
        }
        if (thread_used == 0) {
            otpack->iknp_straight->recv(layer_OT_sum_seed, r, tree_depth);
        } else {
            otpack->iknp_reversed->recv(layer_OT_sum_seed, r, tree_depth);
        }

        // build the punctured tree
        // the parameters of layer1
        uint32_t skip_index = r[0] ^ 1;
        uint32_t points_this_layer = 2;
        uint32_t step_size = leaves_size / points_this_layer;
        uint32_t step_size_last_layer = 2 * step_size;
        uint32_t points_last_layer = points_this_layer / 2;
        GGM_leaves[r[0] * step_size] = layer_OT_sum_seed[0];
        auto *temp_data = new block128[2];
        for (int i = 2; i <= tree_depth; i++) {
            points_this_layer *= 2;
            step_size = leaves_size / points_this_layer;
            step_size_last_layer = 2 * step_size;
            points_last_layer = points_this_layer / 2;
            // first build the tree
            for (int j = 0; j < points_last_layer; j++) {
                if (j == skip_index) {
                    continue;
                }
                prg.reseed(&GGM_leaves[j * step_size_last_layer]);
                prg.random_block(temp_data, 2);
                GGM_leaves[j * step_size_last_layer] = temp_data[0];
                GGM_leaves[j * step_size_last_layer + step_size] = temp_data[1];
            }
            // #ifndef NDEBUG
            //             std::cout << std::endl;
            //             for (int ii = 0; ii < leaves_size; ii++) {
            //                 print_block(GGM_leaves[ii]);
            //             }
            //             std::cout << std::endl;
            // #endif
            // then resume the value
            block128 temp = layer_OT_sum_seed[i - 1];
            for (int j = 0; j < points_this_layer; j += 2) {
                temp = _mm_xor_si128(temp, GGM_leaves[(j + r[i - 1]) * step_size]);
            }
            GGM_leaves[skip_index * step_size_last_layer + step_size * r[i - 1]] = temp;
            skip_index = (skip_index << 1) + (r[i - 1] ^ 1);
            // #ifndef NDEBUG
            //             std::cout << std::endl;
            //             for (int ii = 0; ii < leaves_size; ii++) {
            //                 print_block(GGM_leaves[ii]);
            //             }
            //             std::cout << std::endl;
            // #endif
        }
        // #ifndef NDEBUG
        //         std::cout << std::endl;
        //         for (int i = 0; i < leaves_size; i++) {
        //             print_block(GGM_leaves[i]);
        //         }
        //         std::cout << std::endl;
        // #endif
        delete[] temp_data;
        delete[] r;
        delete[] layer_OT_sum_seed;
    }
    memcpy(seed, GGM_leaves, length * 16 * sizeof(uint8_t));
    delete[] GGM_leaves;
}


// batch aims to make n times m length to nm length
void AuxProtocols::nMinus1OUTNOT_batch(block128 *seed, uint64_t batch_size, uint64_t length, uint64_t *offset) {
    const int tree_depth = static_cast<int>(std::ceil(std::log2(length)));
    const int leaves_size = static_cast<int>(std::pow(2, tree_depth));
    auto *batch_GGM_leaves = new block128[batch_size * leaves_size]();
    PRG128 prg;
    if (party == sci::ALICE) {
        block128 *my_seed = new block128[batch_size];
        prg.random_block(my_seed, batch_size);
        auto **layer_OT_sum_seed = new block128 *[2];
        layer_OT_sum_seed[0] = new block128[tree_depth * batch_size];
        layer_OT_sum_seed[1] = new block128[tree_depth * batch_size];
        // build tree
        for (int i = 0; i < batch_size; i++) {
            batch_GGM_leaves[i * leaves_size] = my_seed[i];
        }
        auto *temp_data = new block128[2];

        // parameters for layer0
        uint32_t points_this_layer = 1;
        uint32_t step_size = leaves_size / points_this_layer;
        uint32_t step_size_last_layer = 2 * step_size;
        uint32_t points_last_layer = points_this_layer / 2;
        for (int i = 1; i <= tree_depth; i++) {
            points_this_layer *= 2;
            step_size = leaves_size / points_this_layer;
            step_size_last_layer = 2 * step_size;
            points_last_layer = points_this_layer / 2;
            for (int j = 0; j < points_last_layer; j++) {
                // see through the nodes in each layer
                for (int k = 0; k < batch_size; k++) {
                    prg.reseed(&batch_GGM_leaves[k * leaves_size + j * step_size_last_layer]);
                    prg.random_block(temp_data, 2);
                    batch_GGM_leaves[k * leaves_size + j * step_size_last_layer] = temp_data[0];
                    batch_GGM_leaves[k * leaves_size + j * step_size_last_layer + step_size] = temp_data[1];
                }
            }
            for (int k = 0; k < batch_size; k++) {
                layer_OT_sum_seed[0][k * tree_depth + i - 1] = batch_GGM_leaves[k * leaves_size + 0];
                layer_OT_sum_seed[1][k * tree_depth + i - 1] = batch_GGM_leaves[k * leaves_size + step_size];
                for (int j = 2; j < points_this_layer; j += 2) {
                    layer_OT_sum_seed[0][k * tree_depth + i - 1] = _mm_xor_si128(
                        layer_OT_sum_seed[0][k * tree_depth + i - 1],
                        batch_GGM_leaves[k * leaves_size + j * step_size]);
                    layer_OT_sum_seed[1][k * tree_depth + i - 1] = _mm_xor_si128(
                        layer_OT_sum_seed[1][k * tree_depth + i - 1],
                        batch_GGM_leaves[k * leaves_size + (j + 1) * step_size]);
                }
            }
            // std::cout << std::endl;
            // for (int ii = 0; ii < leaves_size; ii++) {
            //     print_block(batch_GGM_leaves[ii]);
            // }
            // std::cout << std::endl;
        }
        delete[] temp_data;
        otpack->iknp_straight->send(layer_OT_sum_seed[0], layer_OT_sum_seed[1], (tree_depth * batch_size));
        // #ifndef NDEBUG
        // std::cout << std::endl;
        // for (int i = 0; i < leaves_size; i++) {
        //     print_block(batch_GGM_leaves[i]);
        // }
        // std::cout << std::endl;
        // #endif
        for (int i = 0; i < 2; i++) {
            delete[] layer_OT_sum_seed[i];
        }
        delete[] layer_OT_sum_seed;
    } else {
        auto *r = new bool[tree_depth * batch_size];
        auto *layer_OT_sum_seed = new block128[tree_depth * batch_size];
        for (int i = 1; i <= tree_depth; i++) {
            for (int k = 0; k < batch_size; k++) {
                r[k * tree_depth + i - 1] = (offset[k] >> (tree_depth - i) & 1) ^ 1; // extract i-th bit
            }
        }
        otpack->iknp_straight->recv(layer_OT_sum_seed, r, tree_depth * batch_size);
        // build the punctured tree
        // the parameters of layer1
        uint32_t points_this_layer = 2;
        uint32_t step_size = leaves_size / points_this_layer;
        uint32_t step_size_last_layer = 2 * step_size;
        uint32_t points_last_layer = points_this_layer / 2;
        uint32_t skip_index = 0;
        auto *temp_data = new block128[2];
        for (int k = 0; k < batch_size; k++) {
            points_this_layer = 2;
            step_size = leaves_size / points_this_layer;
            step_size_last_layer = 2 * step_size;
            points_last_layer = points_this_layer / 2;
            skip_index = r[0 + k * tree_depth] ^ 1;
            batch_GGM_leaves[k * leaves_size + r[0 + k * tree_depth] * step_size] = layer_OT_sum_seed[0 + k * tree_depth];
            for (int i = 2; i <= tree_depth; i++) {
                points_this_layer *= 2;
                step_size = leaves_size / points_this_layer;
                step_size_last_layer = 2 * step_size;
                points_last_layer = points_this_layer / 2;
                // first build the tree
                for (int j = 0; j < points_last_layer; j++) {
                    if (j == skip_index) {
                        continue;
                    }
                    prg.reseed(&batch_GGM_leaves[k * leaves_size + j * step_size_last_layer]);
                    prg.random_block(temp_data, 2);
                    batch_GGM_leaves[k * leaves_size + j * step_size_last_layer] = temp_data[0];
                    batch_GGM_leaves[k * leaves_size + j * step_size_last_layer + step_size] = temp_data[1];
                }

                // std::cout << std::endl;
                // for (int ii = 0; ii < leaves_size; ii++) {
                //     print_block(batch_GGM_leaves[ii]);
                // }
                // std::cout << std::endl;

                block128 temp = layer_OT_sum_seed[k * tree_depth + i - 1];
                for (int j = 0; j < points_this_layer; j += 2) {
                    temp = _mm_xor_si128(
                        temp, batch_GGM_leaves[k * leaves_size + (j + r[k * tree_depth + i - 1]) * step_size]);
                }
                batch_GGM_leaves[k * leaves_size + skip_index * step_size_last_layer + step_size * r[k * tree_depth + i - 1]]
                        = temp;
                skip_index = (skip_index << 1) + (r[k * tree_depth + i - 1] ^ 1);

                // std::cout << std::endl;
                // for (int ii = 0; ii < leaves_size; ii++) {
                //     print_block(batch_GGM_leaves[ii]);
                // }
                // std::cout << std::endl;
            }
        }
        // #ifndef NDEBUG
        // std::cout << std::endl;
        // for (int i = 0; i < leaves_size; i++) {
        //     print_block(batch_GGM_leaves[i]);
        // }
        // std::cout << std::endl;
        // #endif
        delete[] temp_data;
        delete[] r;
        delete[] layer_OT_sum_seed;
    }
    for (int i = 0; i < batch_size; i++) {
        memcpy(seed + i * length, batch_GGM_leaves + i * leaves_size, length * 16 * sizeof(uint8_t));
    }
    // memcpy(seed, batch_GGM_leaves, batch_size * length * 16 * sizeof(uint8_t));
    delete[] batch_GGM_leaves;
}


void AuxProtocols::unpack_bits(const uint8_t *x_packed, uint8_t *x, int length) {
    for (int i = 0; i < length; ++i) {
        x[i] = (x_packed[i / 8] >> (i % 8)) & 1;
    }
}

void AuxProtocols::pack_bits(const uint8_t *x, uint8_t *x_packed, int length) {
    for (int i = 0; i < length; ++i) {
        x_packed[i / 8] |= (x[i] & 1) << (i % 8);
    }
}

void AuxProtocols::pack_bits(const uint8_t *x, uint64_t *x_packed, int length) {
    for (int i = 0; i < length; ++i) {
        x_packed[i / 64] |= (x[i] & 1) << (i % 64);
    }
}

void AuxProtocols::unpack_bits(const uint64_t *x_packed, uint8_t *x, int length) {
    for (int i = 0; i < length; ++i) {
        x[i] = (x_packed[i / 64] >> (i % 64)) & 1;
    }
}


// offset: the arithmetic share of the unit index
// length: size of the array
// uniShr: initialize the uniShr outside (for memory management) and all zeros with new T[length]()
void AuxProtocols::uniShare_naive_bool(uint8_t *uniShr, int length, const uint64_t offset, uint8_t thread_used = 0) {
    auto *seeds = new block128[length]();
    PRG128 prg;
    auto *x = new uint8_t[length]();
    if (party == sci::ALICE) {
        // set the offset directly in the index
        x[offset] = 1;
        nMinus1OUTNOT(seeds, length, 0, thread_used);
        // #ifndef NDEBUG
        //         std::cout << std::endl;
        //         for (int ii = 0; ii < length; ii++) {
        //             print_block(seeds[ii]);
        //         }
        //         std::cout << std::endl;
        // #endif
        // generate the shift translation shares
        auto **shift_translate = new uint8_t *[length];
        for (int i = 0; i < length; i++) {
            shift_translate[i] = new uint8_t[length]();
            prg.reseed(&seeds[i]);
            prg.random_bool((bool *) shift_translate[i], length);
        }
        auto *a = new uint8_t[length]();
        auto *b = new uint8_t[length]();
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                a[j] ^= shift_translate[i][j];
                b[i] ^= shift_translate[(i - j + length) % length][j];
            }
        }
        memcpy(uniShr, b, length * sizeof(uint8_t));
        for (int i = 0; i < length; i++) {
            x[i] = x[i] ^ a[i];
        }
        uint8_t *x_packed = new uint8_t[(length + 7) / 8]();
        pack_bits(x, x_packed, length);
        iopack->io->send_data(x_packed, ((length + 7) / 8) * sizeof(uint8_t));

        // #ifndef NDEBUG
        //         std::cout << std::endl;
        //         for (int i = 0; i < length; i++) {
        //             std::cout << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(a[i]) << " ";
        //         }
        //         std::cout << std::endl;
        //         for (int i = 0; i < length; i++) {
        //             std::cout << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b[i]) << " ";
        //         }
        //         std::cout << std::endl;
        // #endif
        for (int i = 0; i < length; i++) {
            delete[] shift_translate[i];
        }
        delete[] shift_translate;
        delete[] x_packed;
        delete[] a;
        delete[] b;
    } else {
        nMinus1OUTNOT(seeds, length, offset, thread_used);
        // #ifndef NDEBUG
        //         std::cout << std::endl;
        //         for (int ii = 0; ii < length; ii++) {
        //             print_block(seeds[ii]);
        //         }
        //         std::cout << std::endl;
        // #endif
        // generate the shift translation shares
        auto **shift_translate = new uint8_t *[length]();
        for (int i = 0; i < length; i++) {
            shift_translate[i] = new uint8_t[length]();
            if (i == offset) {
                continue;
            }
            prg.reseed(&seeds[i]);
            prg.random_bool(reinterpret_cast<bool *>(shift_translate[i]), length);
        }

        auto *c = new uint8_t[length]();
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                c[i] ^= shift_translate[(i - j + length) % length][j];
                c[(j + offset + length) % length] ^= shift_translate[i][j];
            }
        }

        // #ifndef NDEBUG
        //         std::cout << std::endl;
        //         for (int i = 0; i < length; i++) {
        //             std::cout << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(c[i]) << " ";
        //         }
        //         std::cout << std::endl;
        // #endif

        auto *xMa = new uint8_t[length]();
        auto *shift_xMa = new uint8_t[length]();
        auto *xMa_packed = new uint8_t[(length + 7) / 8];
        iopack->io->recv_data(xMa_packed, ((length + 7) / 8) * sizeof(uint8_t));
        unpack_bits(xMa_packed, xMa, length);
        delete[] xMa_packed;

        std::memcpy(shift_xMa, xMa + length - offset, offset);
        std::memcpy(shift_xMa + offset, xMa, length - offset);

        // #ifndef NDEBUG
        //         std::cout << std::endl;
        //         for (int i = 0; i < length; i++) {
        //             std::cout << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(xMa[i]) << " ";
        //         }
        //         std::cout << std::endl;
        // #endif
        // #ifndef NDEBUG
        //         std::cout << std::endl;
        //         std::cout << "the shifted vector" << std::endl;
        //         for (int i = 0; i < length; i++) {
        //             std::cout << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(shift_xMa[i]) << " ";
        //         }
        //         std::cout << std::endl;
        // #endif

        for (int i = 0; i < length; i++) {
            uniShr[i] = c[i] ^ shift_xMa[i];
        }

        for (int i = 0; i < length; i++) {
            delete[] shift_translate[i];
        }
        delete[] shift_translate;
        delete[] c;
        delete[] xMa;
        delete[] shift_xMa;
    }
    delete[] seeds;
    delete[] x;
}


void AuxProtocols::uniShare_naive_bool_batch(uint8_t *uniShr, int batch_size, int length, uint64_t *offset) {
    auto *seeds = new block128[length * batch_size]();
    PRG128 prg;
    auto *x = new uint8_t[length * batch_size]();
    if (party == sci::ALICE) {
        // set the offset directly in the index
        for (int i = 0; i < batch_size; i++) {
            x[offset[i] + i * length] = 1;
        }
        nMinus1OUTNOT_batch(seeds, batch_size, length, nullptr);
        // std::cout << std::endl;
        // for (int ii = 0; ii < length; ii++) {
        //     print_block(seeds[ii]);
        // }
        // std::cout << std::endl;

        // generate the shift translation shares
        auto **shift_translate = new uint8_t *[length * batch_size];
        auto *a = new uint8_t[length * batch_size]();
        auto *b = new uint8_t[length * batch_size]();
        uint8_t *x_packed = new uint8_t[(length * batch_size + 7) / 8]();
        for (int k = 0; k < batch_size; k++) {
            for (int i = 0; i < length; i++) {
                shift_translate[i + k * length] = new uint8_t[length]();
                prg.reseed(&seeds[i + k * length]);
                prg.random_bool((bool *) shift_translate[i + k * length], length);
            }
            for (int i = 0; i < length; i++) {
                for (int j = 0; j < length; j++) {
                    a[j + k * length] ^= shift_translate[i + k * length][j];
                    b[i + k * length] ^= shift_translate[(i - j + length) % length + k * length][j];
                }
            }
        }
        memcpy(uniShr, b, batch_size * length * sizeof(uint8_t));
        for (int i = 0; i < length * batch_size; i++) {
            x[i] = x[i] ^ a[i];
        }

        pack_bits(x, x_packed, length * batch_size);
        iopack->io->send_data(x_packed, ((length * batch_size + 7) / 8) * sizeof(uint8_t));

        // std::cout << std::endl;
        // for (int i = 0; i < length; i++) {
        //     std::cout << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b[length + i]) << " ";
        // }
        // std::cout << std::endl;

        for (int i = 0; i < length * batch_size; i++) {
            delete[] shift_translate[i];
        }
        delete[] shift_translate;
        delete[] x_packed;
        delete[] a;
        delete[] b;
    }
    else {
        nMinus1OUTNOT_batch(seeds, batch_size, length, offset);
        // #ifndef NDEBUG
        // std::cout << std::endl;
        // for (int ii = 0; ii < length; ii++) {
        //     print_block(seeds[ii]);
        // }
        // std::cout << std::endl;
        // #endif
        // generate the shift translation shares
        auto **shift_translate = new uint8_t *[length * batch_size]();
        auto *c = new uint8_t[length * batch_size]();
        for (int k = 0; k < batch_size; k++) {
            for (int i = 0; i < length; i++) {
                shift_translate[i + k * length] = new uint8_t[length]();
                if (i == offset[k]) {
                    continue;
                }
                prg.reseed(&seeds[i + k * length]);
                prg.random_bool(reinterpret_cast<bool *>(shift_translate[i + k * length]), length);
            }
            for (int i = 0; i < length; i++) {
                for (int j = 0; j < length; j++) {
                    c[i + k * length] ^= shift_translate[(i - j + length) % length + k * length][j];
                    c[(j + offset[k] + length) % length + k * length] ^= shift_translate[i + k * length][j];
                }
            }
        }

        auto *xMa = new uint8_t[length * batch_size]();
        auto *shift_xMa = new uint8_t[length * batch_size]();
        auto *xMa_packed = new uint8_t[(length * batch_size + 7) / 8];
        iopack->io->recv_data(xMa_packed, ((length * batch_size + 7) / 8) * sizeof(uint8_t));
        unpack_bits(xMa_packed, xMa, length * batch_size);
        delete[] xMa_packed;

        for (int k = 0; k < batch_size; k++) {
            std::memcpy(k * length + shift_xMa, k * length + xMa + length - offset[k], offset[k]);
            std::memcpy(k * length + shift_xMa + offset[k], k * length + xMa, length - offset[k]);
        }
        for (int k = 0; k < batch_size; k++) {
            for (int i = 0; i < length; i++) {
                uniShr[k * length + i] = c[k * length + i] ^ shift_xMa[k * length + i];
            }
        }


        // std::cout << std::endl;
        // std::cout << "the shifted vector" << std::endl;
        // for (int i = 0; i < length; i++) {
        //     std::cout << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(uniShr[length + i]) << " ";
        // }
        // std::cout << std::endl;

        for (int i = 0; i < length * batch_size; i++) {
            delete[] shift_translate[i];
        }
        delete[] shift_translate;
        delete[] c;
        delete[] xMa;
        delete[] shift_xMa;
    }
    delete[] seeds;
    delete[] x;
}


void AuxProtocols::multiplexer_two_plain(uint8_t *sel, uint64_t *x, uint64_t *y,
                                         int32_t size, int32_t bw_x, int32_t bw_y) {
    assert(bw_x <= 64 && bw_y <= 64 && bw_y <= bw_x);
    uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));

    uint64_t *corr_data = new uint64_t[size];
    uint64_t *data = new uint64_t[size];

    // y = (sel_0 \xor sel_1) * (x_0)
    // y = (sel_0 + sel_1 - 2*sel_0*sel_1)*x_0
    // y = [sel_0*x_0 + sel_1*(x_0 - 2*sel_0*x_0)]
    if (party == sci::ALICE) {
        // sel_0
        for (int i = 0; i < size; i++) {
            corr_data[i] = (x[i] * (1 - 2 * uint64_t(sel[i]))) & mask_y;
        }
        otpack->iknp_straight->send_cot(data, corr_data, size, bw_y);
        for (int i = 0; i < size; i++) {
            y[i] = (x[i] * uint64_t(sel[i]) + data[i]) & mask_y;
        }
    } else {
        // sel_1
        // party == sci::BOB
        otpack->iknp_straight->recv_cot(data, (bool *) sel, size, bw_y);
        for (int i = 0; i < size; i++) {
            y[i] = data[i] & mask_y;
        }
    }
    delete[] corr_data;
    delete[] data;
}

void AuxProtocols::multiplexer_two_plain(uint8_t *sel, uint64_t x, uint64_t *y,
                                         int32_t size, int32_t bw_x, int32_t bw_y) {
    assert(bw_x <= 64 && bw_y <= 64 && bw_y <= bw_x);
    uint64_t mask_y = (bw_y == 64 ? -1 : ((1ULL << bw_y) - 1));

    uint64_t *corr_data = new uint64_t[size];
    uint64_t *data = new uint64_t[size];

    // y = (sel_0 \xor sel_1) * (x_0)
    // y = (sel_0 + sel_1 - 2*sel_0*sel_1)*x_0
    // y = [sel_0*x_0 + sel_1*(x_0 - 2*sel_0*x_0)]

    if (party == sci::ALICE) {
        // sel_0
        for (int i = 0; i < size; i++) {
            corr_data[i] = (x * (1 - 2 * uint64_t(sel[i]))) & mask_y;
        }
        otpack->iknp_straight->send_cot(data, corr_data, size, bw_y);
        for (int i = 0; i < size; i++) {
            y[i] = (x * uint64_t(sel[i]) - data[i]) & mask_y;
        }
    } else {
        // sel_1
        // party == sci::BOB
        otpack->iknp_straight->recv_cot(data, (bool *) sel, size, bw_y);
        for (int i = 0; i < size; i++) {
            y[i] = data[i] & mask_y;
        }
    }

    delete[] corr_data;
    delete[] data;
}

void AuxProtocols::multiplexer_bShr(uint8_t *sel, uint8_t **x, uint8_t **y, int32_t sel_size, int32_t length) {
    int pack_len = (length + 63) / 64;
    uint64_t *corr_data = new uint64_t[sel_size * pack_len];
    uint64_t *x_packed = new uint64_t[sel_size * pack_len]();
    uint64_t *data_S = new uint64_t[sel_size * pack_len];
    uint64_t *data_R = new uint64_t[sel_size * pack_len];
    uint8_t *sel_ext = new uint8_t[sel_size * pack_len];
    for (int i = 0; i < sel_size; i++) {
        for (int j = 0; j < pack_len; j++) {
            sel_ext[i * pack_len + j] = sel[i];
        }
    }


    // y = (sel_0 \xor sel_1) * (x_0 \xor x_1)
    // y = (sel_0 \xor sel_1)*x_0 + (sel_0 \xor sel_1*x_1
    // y = [sel_0*x_0 + sel_1*(x_0 - 2*sel_0*x_0)]
    //     + [sel_1*x_1 + sel_0*(x_1 - 2*sel_1*x_1)]
    for (int i = 0; i < sel_size; i++) {
        pack_bits(x[i], &x_packed[i * pack_len], length);
    }
    for (int i = 0; i < sel_size * pack_len; i++) {
        corr_data[i] = x_packed[i];
    }

#pragma omp parallel num_threads(2)
    {
        if (omp_get_thread_num() == 1) {
            if (party == sci::ALICE) {
                otpack->iknp_reversed->recv_cot_BShr(data_R, (bool *) sel_ext, sel_size * pack_len, 64);
            } else {
                // party == sci::BOB
                otpack->iknp_reversed->send_cot_BShr(data_S, corr_data, sel_size * pack_len, 64);
            }
        } else {
            if (party == sci::ALICE) {
                otpack->iknp_straight->send_cot_BShr(data_S, corr_data, sel_size * pack_len, 64);
            } else {
                // party == sci::BOB
                otpack->iknp_straight->recv_cot_BShr(data_R, (bool *) sel_ext, sel_size * pack_len, 64);
            }
        }
    }

    uint64_t *sel_mask = new uint64_t[sel_size * pack_len];
    for (int i = 0; i < sel_size * pack_len; i++) {
        sel_mask[i] = sel_ext[i] ? ~0ULL : 0ULL;
    }

    for (int i = 0; i < sel_size * pack_len; i++) {
        corr_data[i] = (x_packed[i] & sel_mask[i]) ^ data_R[i] ^ data_S[i];
    }

    for (int i = 0; i < sel_size; i++) {
        unpack_bits(&corr_data[i * pack_len], y[i], length);
        // corr_data[i] = x[i];
    }
    delete[] corr_data;
    delete[] data_S;
    delete[] data_R;
    delete[] sel_ext;
    delete[] x_packed;
    delete[] sel_mask;
}

void AuxProtocols::msnzb_sci_tree(uint64_t *x, uint64_t *msnzb_index, int32_t bw_x,
                                  int32_t size, int32_t digit_size) {
    // The protocol only works when num_digits = ceil((bw_x * 1.0) / digit_size) is a power of 2.
    int32_t last_digit_size = bw_x % digit_size;
    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
    uint64_t digit_mask = (digit_size == 64 ? -1 : ((1ULL << digit_size) - 1));
    uint64_t last_digit_mask =
            (last_digit_size == 64 ? -1 : ((1ULL << last_digit_size) - 1));
    if (last_digit_size == 0) {
        last_digit_mask = digit_mask;
        last_digit_size = digit_size;
    }
    int32_t num_digits = ceil((bw_x * 1.0) / digit_size);
    int32_t depth = (int32_t) log2(num_digits);
    uint64_t *x_digits = new uint64_t[size];
    uint64_t *temp = new uint64_t[size];
    uint64_t *temp1 = new uint64_t[size];
    uint64_t *temp2 = new uint64_t[size];
    uint8_t **branch_choices = new uint8_t *[size];
    for (int i = 0; i < size; i++) {
        branch_choices[i] = new uint8_t[size];
    }

    XTProtocol *xt = new XTProtocol(this->party, this->iopack, this->otpack);

    // use the tree to search out the digits, and record the branches
    Truncation trunc(this->party, this->iopack, this->otpack);
    Equality eq(this->party, this->iopack, this->otpack);

    trunc.truncate_and_reduce(size, x, temp, bw_x / 2, bw_x);
    if (this->party == sci::ALICE) {
        for (int i = 0; i < size; i++) {
            temp2[i] = (1ULL << (bw_x / 2)) - temp[i];
            temp2[i] &= (1ULL << (bw_x / 2)) - 1;
        }
        eq.check_equality(branch_choices[0], temp2, size, bw_x / 2);
        for (int i = 0; i < size; i++) {
            temp2[i] = x[i] - temp[i];
            temp2[i] &= (1ULL << (bw_x / 2)) - 1;
        }
        multiplexer(branch_choices[0], temp2, temp1, size, bw_x / 2, bw_x / 2);
        for (int i = 0; i < size; i++) {
            temp[i] = temp1[i] + temp[i];
            temp[i] &= (1ULL << (bw_x / 2)) - 1;
        }
    } else {
        eq.check_equality(branch_choices[0], temp, size, bw_x / 2);
        for (int i = 0; i < size; i++) {
            temp2[i] = x[i] - temp[i];
            temp2[i] &= (1ULL << (bw_x / 2)) - 1;
        }
        multiplexer(branch_choices[0], temp2, temp1, size, bw_x / 2, bw_x / 2);
        for (int i = 0; i < size; i++) {
            temp[i] = temp1[i] + temp[i];
            temp[i] &= (1ULL << (bw_x / 2)) - 1;
        }
    }

    int bw_temp = bw_x / 2;
    for (int i = 1; i < depth; i++) {
        trunc.truncate_and_reduce(size, temp, temp1, bw_temp / 2, bw_temp);
        if (this->party == sci::ALICE) {
            for (int j = 0; j < size; j++) {
                temp2[j] = (1ULL << (bw_temp / 2)) - temp1[j];
                temp2[j] &= (1ULL << (bw_temp / 2)) - 1;
            }
            eq.check_equality(branch_choices[i], temp2, size, bw_temp / 2);
            for (int j = 0; j < size; j++) {
                temp2[j] = temp[j] - temp1[j];
                temp2[j] &= (1ULL << (bw_temp / 2)) - 1;
            }
            multiplexer(branch_choices[i], temp2, temp, size, bw_temp / 2, bw_temp / 2);
            for (int j = 0; j < size; j++) {
                temp[j] = temp[j] + temp1[j];
                temp[j] &= (1ULL << (bw_temp / 2)) - 1;
            }
        } else {
            eq.check_equality(branch_choices[i], temp1, size, bw_temp / 2);
            for (int j = 0; j < size; j++) {
                temp2[j] = temp[j] - temp1[j];
                temp2[j] &= (1ULL << (bw_temp / 2)) - 1;
            }
            multiplexer(branch_choices[i], temp2, temp, size, bw_temp / 2, bw_temp / 2);
            for (int j = 0; j < size; j++) {
                temp[j] = temp[j] + temp1[j];
                temp[j] &= (1ULL << (bw_temp / 2)) - 1;
            }
        }
        bw_temp /= 2;
    }

    // in here provide the partial digits.
    // Use LUTs for MSNZB on digits

    int D = (1 << digit_size);
    int DLast = (1 << last_digit_size);
    uint8_t *z_ = new uint8_t[size];
    uint64_t *msnzb_ = new uint64_t[size];
    uint64_t *msnzb_extended = new uint64_t[size];
    int lookup_output_bits = (ceil(log2(digit_size))) + 1;
    int mux_bits = ceil(log2(bw_x));
    uint64_t msnzb_mask = (1ULL << (lookup_output_bits - 1)) - 1;
    uint64_t mux_mask = (1ULL << mux_bits) - 1;
    if (party == ALICE) {
        uint64_t **spec;
        spec = new uint64_t *[size];
        PRG128 prg;
        prg.random_data(z_, size * sizeof(uint8_t));
        prg.random_data(msnzb_, size * sizeof(uint64_t));
        for (int i = 0; i < size; i++) {
            spec[i] = new uint64_t[D];
            z_[i] &= 1;
            msnzb_[i] &= msnzb_mask;
            for (int j = 0; j < D; j++) {
                // int idx = (x_digits[i] + j) & digit_mask;
                int idx = (temp[i] + j) & digit_mask;
                uint64_t lookup_val = lookup_msnzb(idx);
                spec[i][j] = ((lookup_val >> 1) - msnzb_[i]) & msnzb_mask;
                spec[i][j] <<= 1;
                spec[i][j] |=
                        ((uint64_t) (((uint8_t) (lookup_val & 1ULL)) ^ z_[i]) & 1ULL);
            }
        }
        this->lookup_table<uint64_t>(spec, nullptr, nullptr, size,
                                     digit_size, lookup_output_bits);

        // Zero extend to mux_bits
        xt->z_extend(size, msnzb_, msnzb_index,
                     lookup_output_bits - 1, mux_bits);

        for (int i = 0; i < size; i++) {
            delete[] spec[i];
        }
        delete[] spec;
    } else {
        // BOB
        this->lookup_table<uint64_t>(nullptr, temp, msnzb_, size,
                                     digit_size, lookup_output_bits);
        // this->lookup_table<uint64_t>(nullptr, x_digits, msnzb_, size,
        // digit_size, lookup_output_bits);
        for (int i = 0; i < size; i++) {
            z_[i] = (uint8_t) (msnzb_[i] & 1ULL);
            msnzb_[i] >>= 1;
        }
        // Zero extend to mux_bits
        xt->z_extend(size, msnzb_, msnzb_index,
                     lookup_output_bits - 1, mux_bits);

        for (int j = 0; j < size; j++) {
            msnzb_index[j] &= mux_mask;
        }
    }

    // Combine MSNZB of digits
    int t = bw_x;
    for (int i = 0; i < depth; i++) {
        t /= 2;
        if (party == ALICE) {
            for (int j = 0; j < size; j++) {
                branch_choices[i][j] ^= 1;
            }
        }
        multiplexer_two_plain(branch_choices[i], t, temp, size, mux_bits, mux_bits);
        for (int j = 0; j < size; j++) {
            msnzb_index[j] += temp[j];
            msnzb_index[j] &= mux_mask;
        }
    }

    delete xt;
    delete[] x_digits;
    delete[] z_;
    delete[] temp;
    delete[] temp1;
    delete[] temp2;
    delete[] msnzb_;
    delete[] msnzb_extended;
    for (int i = 0; i < size; i++) {
        delete[] branch_choices[i];
    }
    delete[] branch_choices;
    return;
}


void AuxProtocols::msnzb_one_hot_tree(uint64_t *x, uint8_t *one_hot_vector,
                                      int32_t bw_x, int32_t size,
                                      int32_t digit_size) {
    uint64_t mask_x = (bw_x == 64 ? -1 : ((1ULL << bw_x) - 1));
    int msnzb_index_bits = ceil(log2(bw_x));
    uint64_t msnzb_index_mask = (1ULL << msnzb_index_bits) - 1;

    uint64_t *msnzb_index = new uint64_t[size];

    // return the msnzb_index of x
    this->msnzb_sci_tree(x, msnzb_index, bw_x, size, digit_size);

    // cout << endl;
    // cout <<  "msnzb_index " << endl;
    // for (int i = 0; i < size; i++) {
    //     cout << msnzb_index[i] << " ";
    // }
    // cout << endl;

    // use LUT to get the one-hot representation
    int D = 1 << msnzb_index_bits;
    uint64_t *xor_mask = new uint64_t[size];
    if (party == ALICE) {
        uint64_t **spec;
        spec = new uint64_t *[size];
        PRG128 prg;
        prg.random_data(one_hot_vector, size * bw_x * sizeof(uint8_t));
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < bw_x; j++) {
                one_hot_vector[i * bw_x + j] &= 1;
            }
            xor_mask[i] = 0ULL;
            for (int j = (bw_x - 1); j >= 0; j--) {
                xor_mask[i] <<= 1;
                xor_mask[i] ^= (uint64_t) one_hot_vector[i * bw_x + j];
            }
        }
        for (int i = 0; i < size; i++) {
            spec[i] = new uint64_t[D];
            for (int j = 0; j < D; j++) {
                int idx = (msnzb_index[i] + j) & msnzb_index_mask;
                uint64_t lookup_val = (1ULL << idx);
                lookup_val ^= xor_mask[i];
                spec[i][j] = lookup_val;
            }
        }
        this->lookup_table<uint64_t>(spec, nullptr, nullptr, size, msnzb_index_bits,
                                     bw_x);

        for (int i = 0; i < size; i++) {
            delete[] spec[i];
        }
        delete[] spec;
    } else {
        // BOB
        uint64_t *temp = new uint64_t[size];
        this->lookup_table<uint64_t>(nullptr, msnzb_index, temp, size,
                                     msnzb_index_bits, bw_x);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < bw_x; j++) {
                one_hot_vector[i * bw_x + j] = (uint8_t) (temp[i] & 1ULL);
                temp[i] >>= 1;
            }
        }
        delete[] temp;
    }
    delete[] xor_mask;
    delete[] msnzb_index;
}
