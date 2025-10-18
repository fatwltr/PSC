//
// Created by a1141 on 25-7-15.
//
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


#include <iomanip>

#include "frequency.h"

#include "BuildingBlocks/truncation.h"
#include "BuildingBlocks/value-extension.h"

using namespace std;
using namespace sci;

Frequency::Frequency(int party, IOPack *iopack, OTPack *otpack) {
    this->party = party;
    this->iopack = iopack;
    this->otpack = otpack;
    this->aux = new AuxProtocols(party, iopack, otpack);
    this->xt = new XTProtocol(party, iopack, otpack);
    this->trunc = new Truncation(party, iopack, otpack);
    this->mult = new LinearOT(party, iopack, otpack);
    this->equality = new Equality(party, iopack, otpack);
}

Frequency::~Frequency() {
    delete this->aux;
    delete this->trunc;
    delete this->mult;
    delete this->equality;
    delete this->xt;
}

// to count the frequency of specific number.
void Frequency::count_eq(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand, int32_t bw_data,
                         int32_t bw_res) {
    uint8_t *temp = new uint8_t[num_data];
    uint64_t mask_data = (bw_data == 64 ? -1 : ((1ULL << bw_data) - 1));
    uint64_t mask_res = (bw_res == 64 ? -1 : ((1ULL << bw_res) - 1));
    uint64_t *t = new uint64_t[num_data];
    for (int i = 0; i < num_stand; i++) {
        if (party == ALICE) {
            for (int j = 0; j < num_data; j++) {
                t[j] = (stand[i] - data[j] + num_stand) % num_stand;
            }
            equality->check_equality(temp, t, num_data, bw_data);
            aux->B2A(temp, t, num_data, bw_res);
            for (int j = 0; j < num_data; j++) {
                res[i] += t[j];
                res[i] &= mask_res;
            }
            res[i] &= mask_res;
        } else {
            equality->check_equality(temp, data, num_data, bw_data);
            aux->B2A(temp, t, num_data, bw_res);
            for (int j = 0; j < num_data; j++) {
                res[i] += t[j];
                res[i] &= mask_res;
            }
            res[i] &= mask_res;
        }
    }
    delete[] temp;
    delete[] t;
}


void Frequency::count_eq_batch(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                               int32_t bw_data,
                               int32_t bw_res) {
    uint8_t *temp = new uint8_t[num_data * num_stand]();
    uint64_t mask_data = (bw_data == 64 ? -1 : ((1ULL << bw_data) - 1));
    uint64_t mask_res = (bw_res == 64 ? -1 : ((1ULL << bw_res) - 1));
    uint64_t *t = new uint64_t[num_data * num_stand]();
    if (party == ALICE) {
        for (int i = 0; i < num_stand; i++) {
            for (int j = 0; j < num_data; j++) {
                t[i * num_data + j] = (stand[i] - data[j] + num_stand) % num_stand;
            }
        }
        equality->check_equality(temp, t, num_data * num_stand, bw_data);
        aux->B2A(temp, t, num_data * num_stand, bw_res);
        for (int i = 0; i < num_stand; i++) {
            for (int j = 0; j < num_data; j++) {
                res[i] += t[i * num_data + j];
                res[i] &= mask_res;
            }
            res[i] &= mask_res;
        }
    } else {
        for (int i = 0; i < num_stand; i++) {
            for (int j = 0; j < num_data; j++) {
                t[i * num_data + j] = (data[j] + num_stand) % num_stand;
            }
        }
        equality->check_equality(temp, t, num_data * num_stand, bw_data);
        aux->B2A(temp, t, num_data * num_stand, bw_res);
        for (int i = 0; i < num_stand; i++) {
            for (int j = 0; j < num_data; j++) {
                res[i] += t[i * num_data + j];
                res[i] &= mask_res;
            }
            res[i] &= mask_res;
        }
    }

    delete[] temp;
    delete[] t;
}


// to count the frequency of specific number.
// the input is in Z/{num_stand}Z and the result is in Z/{2^{bw_res}}Z.
void Frequency::count_shift(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                            int32_t bw_data, int32_t bw_res) {
    uint8_t *t_shifted = new uint8_t[num_stand];
    uint64_t mask_data = (bw_data == 64 ? -1 : ((1ULL << bw_data) - 1));
    uint64_t mask_res = (bw_res == 64 ? -1 : ((1ULL << bw_res) - 1));
    uint64_t *t = new uint64_t[num_stand];
    for (int i = 0; i < num_data; i++) {
        // too many communication rounds, batch and packet some num_data
        aux->uniShare_naive_bool(t_shifted, num_stand, data[i] % num_stand, 0);
        aux->B2A(t_shifted, t, num_stand, bw_res); // B2A need batch more message
        for (int j = 0; j < num_stand; j++) {
            res[j] += t[j];
            res[j] &= mask_res;
        }
    }
    delete[] t_shifted;
    delete[] t;
}


void Frequency::count_shift_batch(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                                  int32_t bw_data, int32_t bw_res) {
    uint8_t *t_shifted = new uint8_t[num_stand * num_data];
    uint64_t mask_res = (bw_res == 64 ? -1 : ((1ULL << bw_res) - 1));
    uint64_t *t = new uint64_t[num_stand * num_data];
    for (int i = 0; i < num_data; i++) {
        data[i] = data[i] % num_stand;
    }
    aux->uniShare_naive_bool_batch(t_shifted, num_data, num_stand, data);
    aux->B2A(t_shifted, t, num_data * num_stand, bw_res);
    for (int i = 0; i < num_data; i++) {
        for (int j = 0; j < num_stand; j++) {
            res[j] += t[i * num_stand + j];
            res[j] &= mask_res;
        }
    }
    delete[] t;
    delete[] t_shifted;
}


// used for the mode function
void Frequency::count_eq_inner_self(uint64_t *res, uint64_t *data, int num_data, int num_stand, int32_t bw_data,
                                    int32_t bw_res) {
    uint8_t *temp = new uint8_t[num_data];
    uint64_t mask_data = (bw_data == 64 ? -1 : ((1ULL << bw_data) - 1));
    uint64_t mask_res = (bw_res == 64 ? -1 : ((1ULL << bw_res) - 1));
    uint64_t *t = new uint64_t[num_data];
    for (int i = 0; i < num_data; i++) {
        if (party == ALICE) {
            for (int j = 0; j < num_data; j++) {
                t[j] = (data[i] - data[j] + num_stand) % num_stand;
            }
            equality->check_equality(temp, t, num_data, bw_data);
            aux->B2A(temp, t, num_data, bw_res);
            cout << "t: ";
            for (int j = 0; j < num_data; j++) {
                cout << t[j] << ", ";
                res[i] += t[j];
                res[i] &= mask_res;
            }
            res[i] &= mask_res;
            cout << res[i] << endl;
        } else {
            for (int j = 0; j < num_data; j++) {
                t[j] = (data[j] - data[i] + num_stand) % num_stand;
            }
            equality->check_equality(temp, t, num_data, bw_data);
            aux->B2A(temp, t, num_data, bw_res);
            cout << "t: ";
            for (int j = 0; j < num_data; j++) {
                cout << t[j] << ", ";
                res[i] += t[j];
                res[i] &= mask_res;
            }
            res[i] &= mask_res;
            cout << res[i] << endl;
        }
    }
    delete[] temp;
    delete[] t;
}


// without use of CRT, compare with the origin numbers
// not complete
void Frequency::mode_naive(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                           int32_t bw_data, int32_t bw_res, uint8_t eq = 0) {
    uint64_t *count;
    int count_len;
    if (eq == 1) {
        // if (num_data > num_stand) {
        count_len = num_stand;
        count = new uint64_t[count_len]();
        count_eq_batch(count, data, num_data, stand, num_stand, bw_data, bw_res + 1);
        // } else {
        //     count_len = num_data;
        //     count = new uint64_t[count_len]();
        //     count_eq_inner_self(count, data, num_data, num_stand, bw_data, bw_res + 1);
        //
        //     // the count function is ok, but the mode paradigm is different from the two others.
        //     delete[] count;
        //     return;
        // }
    } else {
        count_len = num_stand;
        count = new uint64_t[count_len]();
        count_shift_batch(count, data, num_data, stand, num_stand, bw_data, bw_res + 1);
    }
    // cout << endl << bw_res + 1 << " " << "count: ";
    // for (int i = 0; i < count_len; i++) {
    //     cout << count[i] << ", ";
    // }
    // cout << endl;
    auto *maxpool_oracle = new MaxPoolProtocol<uint64_t>(
        party, RING, iopack, bw_res + 1, 4, 0, otpack);
    uint64_t tt;
    maxpool_oracle->funcMaxMPC_spcial(count_len, count, &tt, res, bw_res);
    delete[] count;
    delete maxpool_oracle;
}


void Frequency::mode_Mor(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand, int32_t bw_data,
                         int32_t bw_res) {
    res[0] = data[0];
    uint8_t *temp = new uint8_t[2]();
    uint64_t count_shr = 0;
    uint64_t *temp64_2 = new uint64_t[2]();
    if (party == ALICE) {
        count_shr = 1;
    }
    uint64_t temp64 = 0;
    for (int i = 1; i < num_data; i++) {
        if (party == ALICE) {
            temp64 = (1ULL << bw_data) - (res[0] - data[i]);
        } else {
            temp64 = (res[0] - data[i]);
        }
        temp64 &= (1ULL << bw_data) - 1;
        // cout << temp64 << " ";
        equality->check_equality(temp, &temp64, 1, bw_data);
        aux->B2A(temp, &temp64, 1, bw_res);
        temp64 *= 2;
        if (party == ALICE) {
            temp64 -= 1;
        }
        count_shr += temp64;
        count_shr &= (1ULL << bw_res) - 1;
        // cout << endl << "count share: " << count_shr << endl;
        if (party == ALICE) {
            temp64 = count_shr;
        } else {
            temp64 = (1ULL << bw_res) - count_shr;
        }
        // check whether the count equals to zero
        equality->check_equality(temp, &temp64, 1, bw_res);

        temp64_2[0] = data[i] - res[0];
        // if count == 0, in the next iteration need to be set to 1;
        if (party == ALICE) {
            temp64_2[1] = 1 - count_shr;
        } else {
            temp64_2[1] = 0 - count_shr;
        }
        temp[1] = temp[0];
        aux->multiplexer(temp, temp64_2, temp64_2, 2, bw_res, bw_res);
        // aux->multiplexer(temp, res, res, 2, bw_res, bw_res);
        res[0] = temp64_2[0] + res[0];
        count_shr += temp64_2[1];
    }
    res[0] &= (1ULL << bw_res) - 1;
    delete[] temp;
    delete[] temp64_2;
}

void Frequency::mode_Mor_batch(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                               int32_t bw_data,
                               int32_t bw_res) {
    int batches = 1;
    if (num_data == 1024) {
        batches = 32;
    } else if (num_data == 2048) {
        batches = 32;
    } else if (num_data == 4096) {
        batches = 64;
    } else if (num_data == 8192) {
        batches = 64;
    } else if (num_data == 16384) {
        batches = 128;
    }
    int batch_size = num_data / batches;
    uint64_t *internel_mode = new uint64_t[batches]();
    for (int i = 0; i < batches; i++) {
        internel_mode[i] = data[i * batch_size];
    }
    uint64_t *count_shr = new uint64_t[batches]();
    uint64_t *temp64 = new uint64_t[batches]();
    uint64_t *temp64_1 = new uint64_t[batches]();
    uint64_t *temp64_2 = new uint64_t[batches]();
    uint8_t *temp = new uint8_t[batches]();
    if (party == ALICE) {
        std::fill(count_shr, count_shr + batches, 1);
    }
    for (int i = 1; i < batch_size; i++) {
        for (int j = 0; j < batches; j++) {
            if (party == ALICE) {
                temp64[j] = (1ULL << bw_data) - (internel_mode[j] - data[i + j * batch_size]);
            } else {
                temp64[j] = (internel_mode[j] - data[i + j * batch_size]);
            }
        }
        // cout << temp64 << " ";
        equality->check_equality(temp, temp64, batches, bw_data);
        aux->B2A(temp, temp64, batches, bw_res);
        for (int j = 0; j < batches; j++) {
            temp64[j] *= 2;
            if (party == ALICE) {
                temp64[j] -= 1;
            }
            count_shr[j] += temp64[j];
            if (party == ALICE) {
                temp64[j] = count_shr[j];
            } else {
                temp64[j] = (1ULL << bw_res) - count_shr[j];
            }
        }

        // check whether the count equals to zero
        equality->check_equality(temp, temp64, batches, bw_res);

        for (int j = 0; j < batches; j++) {
            temp64_1[j] = data[i + j * batch_size] - internel_mode[j];
            // if count == 0, in the next iteration need to be set to 1;
            if (party == ALICE) {
                temp64_2[j] = 1 - count_shr[j];
            } else {
                temp64_2[j] = 0 - count_shr[j];
            }
        }
        aux->multiplexer(temp, temp64_1, temp64_1, batches, bw_res, bw_res);
        aux->multiplexer(temp, temp64_2, temp64_2, batches, bw_res, bw_res);
        for (int j = 0; j < batches; j++) {
            internel_mode[j] = temp64_2[j] + internel_mode[j];
            count_shr[j] += temp64_1[j];
        }
    }
    if (batches == 1) {
        res[0] = internel_mode[0];
    } else {
        xt->z_extend(batches, count_shr, count_shr, bw_res, bw_res + 1);
        int l = bw_res + 1;
        uint64_t max_count = count_shr[0];
        uint64_t max_temp = count_shr[0];
        uint64_t max_temp2 = count_shr[0];
        uint64_t t_trans;
        uint64_t compare_with;
        uint8_t equal = 0;
        uint8_t count_bigger = 0;
        res[0] = internel_mode[0];
        for (int i = 1; i < batches; i++) {
            if (party == ALICE) {
                t_trans = (1ULL << bw_data) - (internel_mode[i] - res[0]);
            } else {
                t_trans = (internel_mode[i] - res[0]);
            }
            equality->check_equality(&equal, &t_trans, 1, bw_data);
            max_temp2 = max_count + count_shr[i];

            if (party == sci::ALICE) {
                uint64_t Za = (max_count - count_shr[i]) & (1ULL << l) - 1;
                compare_with = Za & ((1 << (l - 1)) - 1);
                aux->mill->compare(&count_bigger, &compare_with, 1, l);
                uint8_t ba = Za < (1 << (l - 1));
                count_bigger ^= ba;
                max_temp = (max_count - count_shr[i]) & (1ULL << l) - 1;
                aux->multiplexer(&count_bigger, &max_temp, &max_temp, 1, l, l);
                max_temp += count_shr[i];
                max_temp &= (1ULL << l) - 1;
            } else {
                uint64_t Zb = (max_count - count_shr[i]) & (1ULL << l) - 1;
                compare_with = (1 << (l - 1)) - (Zb & ((1 << (l - 1)) - 1));
                aux->mill->compare(&count_bigger, &compare_with, 1, l);
                max_temp = (max_count - count_shr[i]) & (1ULL << l) - 1;
                uint8_t bb = Zb < (1 << (l - 1));
                count_bigger ^= bb;
                count_bigger ^= 1;
                aux->multiplexer(&count_bigger, &max_temp, &max_temp, 1, l, l);
                max_temp += count_shr[i];
                max_temp &= (1ULL << l) - 1;
            }
            max_temp = 2 * max_temp - max_count - count_shr[i]; // a - b
            max_count = max_temp - max_temp2; // max_temp = a + b
            aux->multiplexer(&equal, &max_count, &max_count, 1, l, l);
            max_count += max_temp2;
            max_count = res[0] - internel_mode[i];
            aux->multiplexer(&count_bigger, &max_count, &max_count, 1, l, l);
            res[0] = max_count + internel_mode[i];
        }
    }
    res[0] &= (1ULL << bw_res) - 1;
    delete[] temp;
    delete[] internel_mode;
    delete[] count_shr;
    delete[] temp64;
    delete[] temp64_1;
    delete[] temp64_2;
}


// the 20th primes start from 2
const vector<uint64_t> primes = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71
};

// 保存最优结果
struct Result {
    uint64_t N = 0;
    uint64_t minSum = numeric_limits<uint64_t>::max();
    vector<int> bestExponents;
};

// DFS 搜索
void dfs(int i, const vector<uint64_t> &primes, int maxExp, uint64_t product,
         uint64_t sum, vector<int> &exponents, uint64_t m,
         Result &result) {
    if (product > m) {
        if (sum < result.minSum) {
            result.minSum = sum;
            result.N = product;
            result.bestExponents = exponents;
        }
        return;
    }
    if (i >= primes.size())
        return;

    for (int a = 0; a <= maxExp; ++a) {
        uint64_t powVal = 1;
        for (int j = 0; j < a; ++j) {
            if (powVal > (numeric_limits<uint64_t>::max() / primes[i])) return; // 防止溢出
            powVal *= primes[i];
        }
        uint64_t newProduct = product * powVal;
        if (newProduct < 0) break;
        uint64_t newSum = sum + (a > 0 ? powVal : 0);
        if (newSum >= result.minSum) break;

        exponents.push_back(a);
        dfs(i + 1, primes, maxExp, newProduct, newSum, exponents, m, result);
        exponents.pop_back();
    }
}

vector<uint64_t> findMinSumOverM(uint64_t m, int maxExp = 10) {
    Result result;
    vector<int> exponents;

    dfs(0, primes, maxExp, 1, 0, exponents, m, result);

    // cout << "m = " << m << endl;
    // cout << "N = " << result.N << endl;
    // cout << "最小的 ∑ p^a = " << result.minSum << endl;

    map<int, int> factorMap;
    vector<uint64_t> powers;
    for (size_t i = 0; i < result.bestExponents.size(); ++i) {
        int exp = result.bestExponents[i];
        if (exp > 0) {
            uint64_t value = 1;
            for (int j = 0; j < exp; ++j) {
                value *= primes[i];
            }
            powers.push_back(value);
        }
    }
    return powers;
}

int64_t mod_inverse(int64_t a, int64_t m) {
    int64_t m0 = m, t, q;
    int64_t x0 = 0, x1 = 1;

    if (m == 1) return 0;

    while (a > 1) {
        q = a / m;
        t = m;

        m = a % m;
        a = t;
        t = x0;

        x0 = x1 - q * x0;
        x1 = t;
    }

    if (x1 < 0)
        x1 += m0;

    return x1;
}


void Frequency::mode_CRT_eq(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                            int32_t bw_data, int32_t bw_res) {
    vector<uint64_t> fragment_modulus = findMinSumOverM(num_stand);
    uint64_t num_fragments = fragment_modulus.size();
    uint64_t prod = 1;
    for (int i = 0; i < num_fragments; i++) {
        prod = prod * fragment_modulus[i];
    }
    // for (int i = 0; i < num_fragments; i++) {
    //     cout << fragment_modulus[i] << " ";
    // }
    // cout << endl;
    auto *t_ge_stand = new uint8_t[num_data];
    auto *conversion_fragments = new uint64_t[num_data * num_fragments];
    auto **data_fragments = new uint64_t *[num_fragments];
    auto res_t = new uint64_t *[num_fragments];
    for (int i = 0; i < num_fragments; ++i) {
        data_fragments[i] = new uint64_t[num_data];
        res_t[i] = new uint64_t[fragment_modulus[i]]();
    }
    auto *maxpool_oracle = new MaxPoolProtocol<uint64_t>(
        party, RING, iopack, bw_res + 1, 4, 0, otpack);
    auto *max_fragments = new uint64_t[num_fragments]();
    if (party == sci::ALICE) {
        uint64_t *temp = new uint64_t[num_data];
        for (int i = 0; i < num_data; ++i) {
            temp[i] = num_stand - data[i];
        }
        aux->mill->compare(t_ge_stand, temp, num_data, bw_data);
        aux->B2A_coprimes(t_ge_stand, conversion_fragments, num_data, fragment_modulus);
        for (int i = 0; i < num_fragments; i++) {
            for (int j = 0; j < num_data; j++) {
                data_fragments[i][j] = (data[j] % fragment_modulus[i] +
                                        (fragment_modulus[i] - num_stand % fragment_modulus[i]) * conversion_fragments
                                        [j * num_fragments + i]) % fragment_modulus[i];
            }
        }
        // for (int i = 0; i < num_fragments; i++) {
        //     for (int j = 0; j < num_data; j++) {
        //         cout << data_fragments[i][j] << ", ";
        //     }
        //     cout << endl;
        // }
        uint64_t tt;
        uint8_t t_ge;
        uint64_t N_trans;
        for (int i = 0; i < num_fragments; i++) {
            int bw_data_t = ceil(log2(fragment_modulus[i] + 1));
            // the bw_res + 1 aims to locally transfer the share compare to two numbers for mill->compare()
            // count_shift_batch(res_t[i], data_fragments[i], num_data, stand, fragment_modulus[i], bw_data_t, bw_res + 1);
            count_eq_batch(res_t[i], data_fragments[i], num_data, stand, fragment_modulus[i], bw_data_t, bw_res + 1);
            // cout << "bw_res + 1: " << bw_res + 1 << endl;
            // for (int j = 0; j < fragment_modulus[i]; j++) {
            //     cout << res_t[i][j] << ", ";
            // }
            // cout << endl;
            maxpool_oracle->funcMaxMPC_spcial(fragment_modulus[i], res_t[i], &tt, &max_fragments[i], bw_data_t);
            // cout << max_fragments[i] << ", ";
            // convert the share to mod $prod$
            tt = (1 << bw_data_t) - max_fragments[i];
            aux->mill->compare(&t_ge, &tt, 1, bw_data_t);
            t_ge ^= 1;
            aux->B2A_coprimes(&t_ge, &N_trans, 1, vector<uint64_t>{prod});
            // cout << "N_trans: " << N_trans << ", ";
            max_fragments[i] = ((max_fragments[i] + prod - (1 << bw_data_t) * N_trans % prod)) % prod;
            // cout << max_fragments[i] << endl;
        }
        // cout << endl;
        // for (int i = 0; i < num_fragments * num_data; i += num_fragments) {
        //     for (int j = 0; j < num_fragments; j++) {
        //         cout << conversion_fragments[i + j] << ", ";
        //     }
        //     cout << endl;
        // }
        //
        delete[] temp;
    } else {
        aux->mill->compare(t_ge_stand, data, num_data, bw_data);
        for (int i = 0; i < num_data; i++) {
            t_ge_stand[i] = t_ge_stand[i] ^ 1;
        }
        aux->B2A_coprimes(t_ge_stand, conversion_fragments, num_data, fragment_modulus);
        for (int i = 0; i < num_fragments; i++) {
            for (int j = 0; j < num_data; j++) {
                data_fragments[i][j] = (data[j] % fragment_modulus[i] +
                                        (fragment_modulus[i] - num_stand % fragment_modulus[i]) * conversion_fragments
                                        [j * num_fragments + i]) % fragment_modulus[i];
            }
        }
        // for (int i = 0; i < num_fragments; i++) {
        //     for (int j = 0; j < num_data; j++) {
        //         cout << data_fragments[i][j] << ", ";
        //     }
        //     cout << endl;
        // }
        uint64_t tt;
        uint8_t t_ge;
        uint64_t N_trans;
        for (int i = 0; i < num_fragments; i++) {
            int bw_data_t = ceil(log2(fragment_modulus[i] + 1));
            // count_shift_batch(res_t[i], data_fragments[i], num_data, stand, fragment_modulus[i], bw_data_t, bw_res + 1);
            count_eq_batch(res_t[i], data_fragments[i], num_data, stand, fragment_modulus[i], bw_data_t, bw_res + 1);
            // cout << "bw_res + 1: " << bw_res + 1 << endl;
            // for (int j = 0; j < fragment_modulus[i]; j++) {
            //     cout << res_t[i][j] << ", ";
            // }
            // cout << endl;
            maxpool_oracle->funcMaxMPC_spcial(fragment_modulus[i], res_t[i], &tt, &max_fragments[i], bw_data_t);
            // cout << max_fragments[i] << ", ";
            // convert to mod $prod$
            aux->mill->compare(&t_ge, &max_fragments[i], 1, bw_data_t);
            aux->B2A_coprimes(&t_ge, &N_trans, 1, vector<uint64_t>{prod});
            // cout << "N_trans: " << N_trans << ", ";
            max_fragments[i] = ((max_fragments[i] + prod - (1 << bw_data_t) * N_trans % prod)) % prod;
            // cout << max_fragments[i] << endl;
        }
        // cout << endl;

        // for (int i = 0; i < num_fragments * num_data; i += num_fragments) {
        //     for (int j = 0; j < num_fragments; j++) {
        //         cout << conversion_fragments[i + j] << ", ";
        //     }
        //     cout << endl;
        // }
    }

    uint64_t Mi = 1;
    uint64_t yi = 1;
    cout << endl << "prod: " << prod << endl;
    res[0] = 0;
    for (int i = 0; i < num_fragments; i++) {
        Mi = prod / fragment_modulus[i];
        yi = mod_inverse(Mi, fragment_modulus[i]);
        // cout << endl << "yi: " << yi << " Mi: " << Mi << endl;
        res[0] += (max_fragments[i] * Mi * yi) % prod;
        res[0] %= prod;
    }

    delete[] t_ge_stand;
    delete[] conversion_fragments;
    delete maxpool_oracle;
    for (int i = 0; i < num_fragments; i++) {
        delete[] data_fragments[i];
        delete[] res_t[i];
    }
    delete[] max_fragments;
    delete[] data_fragments;
    delete[] res_t;
}

void Frequency::mode_CRT_shift(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                               int32_t bw_data, int32_t bw_res) {
    vector<uint64_t> fragment_modulus = findMinSumOverM(num_stand);
    uint64_t num_fragments = fragment_modulus.size();
    uint64_t prod = 1;
    for (int i = 0; i < num_fragments; i++) {
        prod = prod * fragment_modulus[i];
    }
    // for (int i = 0; i < num_fragments; i++) {
    //     cout << fragment_modulus[i] << " ";
    // }
    // cout << endl;
    auto *t_ge_stand = new uint8_t[num_data];
    auto *conversion_fragments = new uint64_t[num_data * num_fragments];
    auto **data_fragments = new uint64_t *[num_fragments];
    auto res_t = new uint64_t *[num_fragments];
    for (int i = 0; i < num_fragments; ++i) {
        data_fragments[i] = new uint64_t[num_data];
        res_t[i] = new uint64_t[fragment_modulus[i]]();
    }
    auto *maxpool_oracle = new MaxPoolProtocol<uint64_t>(
        party, RING, iopack, bw_res + 1, 4, 0, otpack);
    auto *max_fragments = new uint64_t[num_fragments]();
    if (party == sci::ALICE) {
        uint64_t *temp = new uint64_t[num_data];
        for (int i = 0; i < num_data; ++i) {
            temp[i] = num_stand - data[i];
        }
        aux->mill->compare(t_ge_stand, temp, num_data, bw_data);
        aux->B2A_coprimes(t_ge_stand, conversion_fragments, num_data, fragment_modulus);
        for (int i = 0; i < num_fragments; i++) {
            for (int j = 0; j < num_data; j++) {
                data_fragments[i][j] = (data[j] % fragment_modulus[i] +
                                        (fragment_modulus[i] - num_stand % fragment_modulus[i]) * conversion_fragments
                                        [j * num_fragments + i]) % fragment_modulus[i];
            }
        }
        // for (int i = 0; i < num_fragments; i++) {
        //     for (int j = 0; j < num_data; j++) {
        //         cout << data_fragments[i][j] << ", ";
        //     }
        //     cout << endl;
        // }
        uint64_t tt;
        uint8_t t_ge;
        uint64_t N_trans;
        for (int i = 0; i < num_fragments; i++) {
            int bw_data_t = ceil(log2(fragment_modulus[i] + 1));
            // the bw_res + 1 aims to locally transfer the share compare to two numbers for mill->compare()
            count_shift_batch(res_t[i], data_fragments[i], num_data, stand, fragment_modulus[i], bw_data_t, bw_res + 1);
            // cout << "bw_res + 1: " << bw_res + 1 << endl;
            // for (int j = 0; j < fragment_modulus[i]; j++) {
            //     cout << res_t[i][j] << ", ";
            // }
            // cout << endl;
            maxpool_oracle->funcMaxMPC_spcial(fragment_modulus[i], res_t[i], &tt, &max_fragments[i], bw_data_t);
            // cout << max_fragments[i] << ", ";
            // convert the share to mod $prod$
            tt = (1 << bw_data_t) - max_fragments[i];
            aux->mill->compare(&t_ge, &tt, 1, bw_data_t);
            t_ge ^= 1;
            aux->B2A_coprimes(&t_ge, &N_trans, 1, vector<uint64_t>{prod});
            // cout << "N_trans: " << N_trans << ", ";
            max_fragments[i] = ((max_fragments[i] + prod - (1 << bw_data_t) * N_trans % prod)) % prod;
            // cout << max_fragments[i] << endl;
        }
        // cout << endl;
        // for (int i = 0; i < num_fragments * num_data; i += num_fragments) {
        //     for (int j = 0; j < num_fragments; j++) {
        //         cout << conversion_fragments[i + j] << ", ";
        //     }
        //     cout << endl;
        // }
        //
        delete[] temp;
    } else {
        aux->mill->compare(t_ge_stand, data, num_data, bw_data);
        for (int i = 0; i < num_data; i++) {
            t_ge_stand[i] = t_ge_stand[i] ^ 1;
        }
        aux->B2A_coprimes(t_ge_stand, conversion_fragments, num_data, fragment_modulus);
        for (int i = 0; i < num_fragments; i++) {
            for (int j = 0; j < num_data; j++) {
                data_fragments[i][j] = (data[j] % fragment_modulus[i] +
                                        (fragment_modulus[i] - num_stand % fragment_modulus[i]) * conversion_fragments
                                        [j * num_fragments + i]) % fragment_modulus[i];
            }
        }
        // for (int i = 0; i < num_fragments; i++) {
        //     for (int j = 0; j < num_data; j++) {
        //         cout << data_fragments[i][j] << ", ";
        //     }
        //     cout << endl;
        // }
        uint64_t tt;
        uint8_t t_ge;
        uint64_t N_trans;
        for (int i = 0; i < num_fragments; i++) {
            int bw_data_t = ceil(log2(fragment_modulus[i] + 1));
            count_shift_batch(res_t[i], data_fragments[i], num_data, stand, fragment_modulus[i], bw_data_t, bw_res + 1);
            // cout << "bw_res + 1: " << bw_res + 1 << endl;
            // for (int j = 0; j < fragment_modulus[i]; j++) {
            //     cout << res_t[i][j] << ", ";
            // }
            // cout << endl;
            maxpool_oracle->funcMaxMPC_spcial(fragment_modulus[i], res_t[i], &tt, &max_fragments[i], bw_data_t);
            // cout << max_fragments[i] << ", ";
            // convert to mod $prod$
            aux->mill->compare(&t_ge, &max_fragments[i], 1, bw_data_t);
            aux->B2A_coprimes(&t_ge, &N_trans, 1, vector<uint64_t>{prod});
            // cout << "N_trans: " << N_trans << ", ";
            max_fragments[i] = ((max_fragments[i] + prod - (1 << bw_data_t) * N_trans % prod)) % prod;
            // cout << max_fragments[i] << endl;
        }
        // cout << endl;

        // for (int i = 0; i < num_fragments * num_data; i += num_fragments) {
        //     for (int j = 0; j < num_fragments; j++) {
        //         cout << conversion_fragments[i + j] << ", ";
        //     }
        //     cout << endl;
        // }
    }

    uint64_t Mi = 1;
    uint64_t yi = 1;
    cout << endl << "prod: " << prod << endl;
    res[0] = 0;
    for (int i = 0; i < num_fragments; i++) {
        Mi = prod / fragment_modulus[i];
        yi = mod_inverse(Mi, fragment_modulus[i]);
        // cout << endl << "yi: " << yi << " Mi: " << Mi << endl;
        res[0] += (max_fragments[i] * Mi * yi) % prod;
        res[0] %= prod;
    }

    delete[] t_ge_stand;
    delete[] conversion_fragments;
    delete maxpool_oracle;
    for (int i = 0; i < num_fragments; i++) {
        delete[] data_fragments[i];
        delete[] res_t[i];
    }
    delete[] max_fragments;
    delete[] data_fragments;
    delete[] res_t;
}


// void Frequency::mode_CRT_shift(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
//                                int32_t bw_data, int32_t bw_res) {
//     vector<uint64_t> fragment_modulus = findMinSumOverM(num_stand);
//     uint64_t num_fragments = fragment_modulus.size();
//     uint64_t prod = 1;
//     for (int i = 0; i < num_fragments; i++) {
//         prod = prod * fragment_modulus[i];
//     }
//     uint64_t *res_count = new uint64_t[num_stand]();
//     auto **data_fragments = new uint64_t *[num_fragments];
//     auto res_t = new uint64_t *[num_fragments];
//     for (int i = 0; i < num_fragments; ++i) {
//         data_fragments[i] = new uint64_t[num_data];
//         res_t[i] = new uint64_t[fragment_modulus[i]]();
//     }
//     auto *maxpool_oracle = new MaxPoolProtocol<uint64_t>(
//         party, RING, iopack, bw_res + 1, 4, 0, otpack);
//     auto *max_fragments = new uint64_t[num_fragments]();
//     count_shift_batch(res_count, data, num_data, stand, num_stand, bw_data, bw_res + 1);
//     if (party == sci::ALICE) {
//         uint64_t tt;
//         uint8_t t_ge;
//         uint64_t N_trans;
//         for (int i = 0; i < num_fragments; i++) {
//             int bw_data_t = ceil(log2(fragment_modulus[i] + 1));
//             for (int j = 0; j < num_stand; j++) {
//                 res_t[i][j % fragment_modulus[i]] += res_count[j];
//                 res_t[i][j % fragment_modulus[i]] &= (1ULL << (bw_res + 1)) - 1;
//             }
//             maxpool_oracle->funcMaxMPC_spcial(fragment_modulus[i], res_t[i], &tt, &max_fragments[i], bw_data_t);
//             tt = (1 << bw_data_t) - max_fragments[i];
//             aux->mill->compare(&t_ge, &tt, 1, bw_data_t);
//             t_ge ^= 1;
//             aux->B2A_coprimes(&t_ge, &N_trans, 1, vector<uint64_t>{prod});
//             max_fragments[i] = ((max_fragments[i] + prod - (1 << bw_data_t) * N_trans % prod)) % prod;
//         }
//     } else {
//         uint64_t tt;
//         uint8_t t_ge;
//         uint64_t N_trans;
//         for (int i = 0; i < num_fragments; i++) {
//             int bw_data_t = ceil(log2(fragment_modulus[i] + 1));
//             for (int j = 0; j < num_stand; j++) {
//                 res_t[i][j % fragment_modulus[i]] += res_count[j];
//                 res_t[i][j % fragment_modulus[i]] &= (1ULL << (bw_res + 1)) - 1;
//             }
//             maxpool_oracle->funcMaxMPC_spcial(fragment_modulus[i], res_t[i], &tt, &max_fragments[i], bw_data_t);
//             aux->mill->compare(&t_ge, &max_fragments[i], 1, bw_data_t);
//             aux->B2A_coprimes(&t_ge, &N_trans, 1, vector<uint64_t>{prod});
//             max_fragments[i] = ((max_fragments[i] + prod - (1 << bw_data_t) * N_trans % prod)) % prod;
//         }
//     }
//
//     uint64_t Mi = 1;
//     uint64_t yi = 1;
//     cout << endl << "prod: " << prod << endl;
//     res[0] = 0;
//     for (int i = 0; i < num_fragments; i++) {
//         Mi = prod / fragment_modulus[i];
//         yi = mod_inverse(Mi, fragment_modulus[i]);
//         res[0] += (max_fragments[i] * Mi * yi) % prod;
//         res[0] %= prod;
//     }
//
//     delete[] res_count;
//     delete maxpool_oracle;
//     for (int i = 0; i < num_fragments; i++) {
//         delete[] data_fragments[i];
//         delete[] res_t[i];
//     }
//     delete[] max_fragments;
//     delete[] data_fragments;
//     delete[] res_t;
// }

// num_data: the number of value range. num_samples: the number of the collected set.

void Frequency::pack_bits(uint8_t **x, uint8_t *x_packed, int rows, int cols) {
    int total_bits = rows * cols;
    int packed_len = (total_bits + 7) / 8; // ceil division
    std::fill(x_packed, x_packed + packed_len, 0); // clear output

    int bit_index = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++, bit_index++) {
            if (x[i][j]) {
                int byte_index = bit_index / 8;
                int bit_offset = bit_index % 8;
                x_packed[byte_index] |= (1 << bit_offset);
            }
        }
    }
}

void Frequency::unpack_bits(uint8_t *x_packed, uint8_t **x, int rows, int cols) {
    int total_bits = rows * cols;
    int bit_index = 0;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++, bit_index++) {
            int byte_index = bit_index / 8;
            int bit_offset = bit_index % 8;
            x[i][j] = (x_packed[byte_index] >> bit_offset) & 1;
        }
    }
}

void Frequency::unpack_bits(const uint8_t *x_packed, uint8_t *x, int length) {
    for (int i = 0; i < length; ++i) {
        x[i] = (x_packed[i / 8] >> (i % 8)) & 1;
    }
}

void Frequency::pack_bits(const uint8_t *x, uint8_t *x_packed, int length) {
    for (int i = 0; i < length; ++i) {
        x_packed[i / 8] |= (x[i] & 1) << (i % 8);
    }
}


void Frequency::pack_data(uint64_t *x, uint64_t *x_packed, int length, int bw) {
    int total_bits = length * bw;
    int num_words = (total_bits + 63) / 64;

    uint64_t bit_pos = 0;

    for (int i = 0; i < length; i++) {
        uint64_t val = x[i] & ((bw == 64) ? ~0ULL : ((1ULL << bw) - 1));
        int word_idx = bit_pos / 64;
        int bit_off = bit_pos % 64;

        x_packed[word_idx] |= (val << bit_off);

        if (bit_off + bw > 64 && word_idx + 1 < num_words) {
            x_packed[word_idx + 1] |= (val >> (64 - bit_off));
        }

        bit_pos += bw;
    }
}

void Frequency::unpack_data(uint64_t *x_packed, uint64_t *x, int length, int bw) {
    int total_bits = length * bw;
    int num_words = (total_bits + 63) / 64;

    uint64_t bit_pos = 0;
    uint64_t mask = (bw == 64) ? ~0ULL : ((1ULL << bw) - 1);

    for (int i = 0; i < length; i++) {
        int word_idx = bit_pos / 64;
        int bit_off = bit_pos % 64;

        uint64_t val = (x_packed[word_idx] >> bit_off);
        if (bit_off + bw > 64 && word_idx + 1 < num_words) {
            int shift = 64 - bit_off;
            if (shift < 64) {
                val |= x_packed[word_idx + 1] << shift;
            }
        }

        x[i] = val & mask;
        bit_pos += bw;
    }
}

void print_block_(const block128 &b) {
    uint8_t bytes[16];
    _mm_storeu_si128((__m128i *) bytes, b); // copy to byte array
    for (int i = 0; i < 16; ++i)
        std::cout << std::hex << std::setw(2) << std::setfill('0') << (int) bytes[i];
    std::cout << std::endl;
}

void Frequency::count_sort(uint64_t *res, uint64_t *frequency, int num_stand, int num_data,
                           int32_t bw_data, int32_t bw_res) {
    int n_th = omp_get_num_procs() / 2;
    cout << "n_th = " << n_th << endl;
    omp_set_num_threads(n_th);
    uint64_t *endpoints = new uint64_t[num_stand]();
    endpoints[0] = frequency[0];
    for (int i = 1; i < num_stand; i++) {
        endpoints[i] = endpoints[i - 1] + frequency[i];
        endpoints[i] &= (1ULL << bw_res) - 1;
    }

    // cout << endl;
    // for (int i = 0; i < num_stand; i++) {
    //     cout << endpoints[i] << " ";
    // }
    // cout << endl;

    num_data += 1;

    // share conversion 2^bw_data -> num_data
    uint8_t *t = new uint8_t[num_stand]();
    uint64_t *appendix = new uint64_t[num_stand]();
    this->aux->wrap_computation(endpoints, t, num_stand, bw_res);
    std::vector<uint64_t> v(1, num_data);
    this->aux->B2A_coprimes(t, appendix, num_stand, v);
    // cout << endl;
    // for (int i = 0; i < num_stand; i++) {
    //     cout << (int) appendix[i] << " ";
    // }
    // cout << endl;
#pragma omp parallel for
    for (int i = 0; i < num_stand; i++) {
        endpoints[i] = (endpoints[i] % num_data) + (num_data - (1ULL << bw_res) % num_data) * (appendix[i]);
        endpoints[i] %= num_data;
    }

    // cout << num_data << " ";
    // cout << endl;
    // for (int i = 0; i < num_stand; i++) {
    //     cout << endpoints[i] << " ";
    // }
    // cout << endl;

    uint64_t mask_data = (bw_data == 64 ? -1 : ((1ULL << bw_data) - 1));
    uint64_t mask_res = (bw_res == 64 ? -1 : ((1ULL << bw_res) - 1));
    auto *seeds = new block128[num_stand * num_data * 2]();
    // auto *prg = new PRG128[num_stand * num_data * 2]();
    PRG128 prg;

    auto **uniShr = new uint8_t *[num_stand];
    for (int i = 0; i < num_stand; i++) {
        uniShr[i] = new uint8_t[num_data]();
    }
    auto **interval_indicate = new uint8_t *[num_stand];
    for (int i = 0; i < num_stand; i++) {
        interval_indicate[i] = new uint8_t[num_data]();
    }
    // int save_memory = num_data * 2 > 2048 ? 2048 : num_data * 2;

    // int save_memory = num_data * 2 * num_data * 2 > 11ULL * 1024 * 1024 * 1024 ? 11ULL * 1024 * 1024 * 1024 / (num_data * 2) : num_data * 2;
    // save_memory = save_memory > num_data * 2 ? num_data * 2 : save_memory;

    uint64_t mem_limit = 11ULL * 1024 * 1024 * 1024; // 11 GiB
    uint64_t total = static_cast<uint64_t>(num_data) * 2 * num_data * 2;

    uint64_t save_memory64 = total > mem_limit
                                 ? mem_limit / (static_cast<uint64_t>(num_data) * 2)
                                 : static_cast<uint64_t>(num_data) * 2;

    uint64_t save_memory = save_memory64;

    cout << "save_memory = " << save_memory << endl;

    if (party == sci::ALICE) {
        // set the offset directly in the index


        // do not need to generate a 2*num_data seed, we can set the second half seed directly as offset < num_data
        this->aux->nMinus1OUTNOT_batch(seeds, num_stand, 2 * num_data, nullptr);
        // generate the shift translation shares

        // cout << "------" << endl;
        // cout << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         print_block_(seeds[i * num_data * 2 + j]);
        //     }
        // }
        // cout << endl;

        auto **shift_translate = new uint8_t *[num_data * 2];
        for (int i = 0; i < num_data * 2; i++) {
            shift_translate[i] = new uint8_t[save_memory]();
        }
        auto *a = new uint8_t [num_data * 2 * num_stand]();
        auto *b = new uint8_t [num_data * 2 * num_stand]();

        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         prg[i * num_data * 2 + j].reseed(&seeds[i * num_data * 2 + j]);
        //     }
        // }

        for (int i = 0; i < num_stand; i++) {
            for (int k = 0; k < (num_data * 2 / save_memory); k++) {
                for (int j = 0; j < num_data * 2; j++) {
                    prg.reseed(&seeds[i * num_data * 2 + j]);
                    prg.random_bool((bool *) shift_translate[j], save_memory);
                }

                // cout << "shift_translate" << endl;
                // for (int j = 0; j < num_data * 2; j++) {
                //     for (int l = 0; l < save_memory; l++) {
                //         cout << (int) shift_translate[j][l] << " ";
                //     }
                //     cout << endl;
                // }

                for (int j = 0; j < num_data * 2; j++) {
                    for (int l = 0; l < save_memory; l++) {
                        a[i * num_data * 2 + l + k * save_memory] ^= shift_translate[j][l];
                        b[i * num_data * 2 + j] ^= shift_translate[
                            (j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
                    }
                }
            }
            int last = num_data * 2 % save_memory;
            if (last == 0) {
                continue;
            }
            int k = (num_data * 2 / save_memory);
            for (int j = 0; j < num_data * 2; j++) {
                prg.reseed(&seeds[i * num_data * 2 + j]);
                prg.random_bool((bool *) shift_translate[j], last);
            }
            // cout << "shift_translate" << endl;
            // for (int j = 0; j < num_data * 2; j++) {
            //     for (int l = 0; l < last; l++) {
            //         cout << (int) shift_translate[j][l] << " ";
            //     }
            //     cout << endl;
            // }
            for (int j = 0; j < num_data * 2; j++) {
                for (int l = 0; l < last; l++) {
                    a[i * num_data * 2 + l + k * save_memory] ^= shift_translate[j][l];
                    b[i * num_data * 2 + j] ^= shift_translate[
                        (j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
                }
            }
        }

        // cout << endl;
        // cout << "----------x---------------" << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         cout << (int) x[i][j] << ", ";
        //     }
        //     cout << endl;
        // }

        // cout << endl;
        // cout << "----------a---------------" << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         cout << (int) a[i][j] << ", ";
        //     }
        //     cout << endl;
        // }

        auto *x = new uint8_t [num_data * 2 * num_stand]();
        for (int i = 0; i < num_stand; i++) {
            std::fill(x + i * num_data * 2, x + i * num_data * 2 + endpoints[i], 1);
            std::fill(x + i * num_data * 2 + num_data + endpoints[i], x + i * num_data * 2 + 2 * num_data, 1);
        }

#pragma omp parallel for
        for (int j = 0; j < num_stand * num_data * 2; j++) {
            x[j] = x[j] ^ a[j];
        }
        // cout << endl;
        // cout << "----------xMa---------------" << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         cout << (int) x[i][j] << ", ";
        //     }
        //     cout << endl;
        // }


        // cout << endl;
        // cout << "------------b-------------" << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         cout << (int) b[i][j] << ", ";
        //     }
        //     cout << endl;
        // }


        uint8_t *x_packed = new uint8_t[(2 * num_data * num_stand + 7) / 8]();
        pack_bits(x, x_packed, num_stand * 2 * num_data);
        iopack->io->send_data(x_packed, ((2 * num_data * num_stand + 7) / 8) * sizeof(uint8_t));

        // allocate the uniShr, and the
        auto **mux_choice = new uint8_t *[num_stand];
        for (int i = 0; i < num_stand; i++) {
            mux_choice[i] = new uint8_t[num_data]();
#pragma omp parallel for
            for (int j = 0; j < num_data; j++) {
                mux_choice[i][j] = b[i * num_data * 2 + j + num_data] ^ b[i * num_data * 2 + j];
            }
        }
        uint64_t *temp = new uint64_t[num_stand];
#pragma omp parallel for
        for (int i = 0; i < num_stand; i++) {
            temp[i] = num_data - endpoints[i];
        }

        this->aux->mill->compare(t, temp, num_stand, ceil(log2(num_data + 1)));

        // the result is >=
        // cout << "t: " << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     cout << (int)t[i] << " ";
        // }
        // cout << endl;


        // mux_choice padding a 0-string
        this->aux->multiplexer_bShr(t, mux_choice, uniShr, num_stand, num_data);
#pragma omp parallel for
        for (int i = 0; i < num_stand; i++) {
            for (int j = 0; j < num_data; j++) {
                uniShr[i][j] ^= b[i * num_data * 2 + j];
            }
        }

        for (int i = 0; i < num_data * 2; i++) {
            delete[] shift_translate[i];
        }
        delete[] shift_translate;
        for (int i = 0; i < num_stand; i++) {
            delete[] mux_choice[i];
        }
        delete[] x_packed;
        delete[] x;
        delete[] temp;
        delete[] mux_choice;
        delete[] a;
        delete[] b;
    } else {
        this->aux->nMinus1OUTNOT_batch(seeds, num_stand, 2 * num_data, endpoints);

        // cout << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     cout << endpoints[i] << endl;
        //     for (int j = 0; j < num_data * 2; j++) {
        //         print_block_(seeds[i * num_data * 2 + j]);
        //     }
        // }
        // cout << endl;

        // generate the shift translation shares
        auto **shift_translate = new uint8_t *[num_data * 2];
        for (int i = 0; i < num_data * 2; i++) {
            shift_translate[i] = new uint8_t[save_memory]();
        }
        auto *c = new uint8_t [num_stand * num_data * 2]();

        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         prg.reseed(&seeds[i * num_data * 2 + j]);
        //     }
        // }

        for (int i = 0; i < num_stand; i++) {
            for (int k = 0; k < (num_data * 2 / save_memory); k++) {
                for (int j = 0; j < num_data * 2; j++) {
                    prg.reseed(&seeds[i * num_data * 2 + j]);
                    prg.random_bool((bool *) shift_translate[j], save_memory);
                }
                // cout << "shift_translate" << endl;
                // for (int j = 0; j < num_data * 2; j++) {
                //     for (int l = 0; l < save_memory; l++) {
                //         cout << (int) shift_translate[j][l] << " ";
                //     }
                //     cout << endl;
                // }
                for (int j = 0; j < num_data * 2; j++) {
                    for (int l = 0; l < save_memory; l++) {
                        c[i * num_data * 2 + j] ^= shift_translate[
                            (j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
                        c[i * num_data * 2 + ((l + k * save_memory + endpoints[i] + 2 * num_data) % (2 * num_data))] ^=
                                shift_translate[j][l];
                    }
                }
            }
            int last = num_data * 2 % save_memory;
            if (last == 0) {
                continue;
            }
            int k = (num_data * 2 / save_memory);
            for (int j = 0; j < num_data * 2; j++) {
                if (j == endpoints[i]) {
                    std::fill(shift_translate[j], shift_translate[j] + save_memory, 0);
                    continue;
                }
                prg.reseed(&seeds[i * num_data * 2 + j]);
                prg.random_bool((bool *) shift_translate[j], last);
            }
            // cout << "shift_translate" << endl;
            // for (int j = 0; j < num_data * 2; j++) {
            //     for (int l = 0; l < last; l++) {
            //         cout << (int) shift_translate[j][l] << " ";
            //     }
            //     cout << endl;
            // }
            for (int j = 0; j < num_data * 2; j++) {
                for (int l = 0; l < last; l++) {
                    c[i * num_data * 2 + j] ^= shift_translate[
                        (j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
                    c[i * num_data * 2 + (l + k * save_memory + endpoints[i] + 2 * num_data) % (2 * num_data)] ^=
                            shift_translate[j][l];
                }
            }
        }

        // cout << endl;
        // cout << "------------c-------------" << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         cout << (int) c[i][j] << ", ";
        //     }
        //     cout << endl;
        // }

        auto *xMa = new uint8_t [num_stand * num_data * 2]();

        auto *shift_xMa = new uint8_t [num_stand * num_data * 2];
        auto *xMa_packed = new uint8_t[(2 * num_data * num_stand + 7) / 8];
        iopack->io->recv_data(xMa_packed, ((2 * num_data * num_stand + 7) / 8) * sizeof(uint8_t));
        unpack_bits(xMa_packed, xMa, num_stand * 2 * num_data);

        // cout << endl;
        // cout << "------------xMa-------------" << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         cout << (int) xMa[i][j] << ", ";
        //     }
        //     cout << endl;
        // }

        for (int i = 0; i < num_stand; i++) {
            std::memcpy(shift_xMa + i * num_data * 2, xMa + i * num_data * 2 + (2 * num_data - endpoints[i]),
                        endpoints[i]);
            std::memcpy(shift_xMa + i * num_data * 2 + endpoints[i], xMa + i * num_data * 2,
                        2 * num_data - endpoints[i]);
        }

        // cout << endl;
        // cout << "------------shift_xMa-------------" << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         cout << (int) shift_xMa[i][j] << ", ";
        //     }
        //     cout << endl;
        // }

#pragma omp parallel for
        for (int j = 0; j < num_stand * num_data * 2; j++) {
            xMa[j] = c[j] ^ shift_xMa[j];
        }


        // cout << endl;
        // cout << "---------c xor shift_xMa-----------" << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         cout << (int) xMa[i][j] << ", ";
        //     }
        //     cout << endl;
        // }


        auto **mux_choice = new uint8_t *[num_stand];
        for (int i = 0; i < num_stand; i++) {
            mux_choice[i] = new uint8_t[num_data]();
#pragma omp parallel for
            for (int j = 0; j < num_data; j++) {
                mux_choice[i][j] = xMa[i * num_data * 2 + j] ^ xMa[i * num_data * 2 + j + num_data];
            }
        }

        this->aux->mill->compare(t, endpoints, num_stand, ceil(log2(num_data + 1)));

        // cout << "t: " << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     cout << (int)t[i] << " ";
        // }
        // cout << endl;


        this->aux->multiplexer_bShr(t, mux_choice, uniShr, num_stand, num_data);
        for (int i = 0; i < num_stand; i++) {
#pragma omp parallel for
            for (int j = 0; j < num_data; j++) {
                uniShr[i][j] ^= xMa[i * num_data * 2 + j];
            }
        }

        for (int i = 0; i < 2 * num_data; i++) {
            delete[] shift_translate[i];
        }
        delete[] shift_translate;
        for (int i = 0; i < num_stand; i++) {
            delete[] mux_choice[i];
        }
        delete[] c;
        delete[] xMa;
        delete[] xMa_packed;
        delete[] shift_xMa;
        // delete[] temp;
        delete[] mux_choice;
    }


    // cout << "temp checkpoint. " << endl;
    //
    // cout << "uniShr: " << endl;
    // for (int i = 0; i < num_stand; i++) {
    //     cout << endl;
    //     for (int j = 0; j < num_data; j++) {
    //         cout << (int) uniShr[i][j] << ", ";
    //     }
    // }
    // cout << endl;

    // xor unishr
    if (party == ALICE) {
#pragma omp parallel for
        for (int j = 0; j < num_data; j++) {
            interval_indicate[0][j] = 1 ^ uniShr[0][j];
        }
    } else {
#pragma omp parallel for
        for (int j = 0; j < num_data; j++) {
            interval_indicate[0][j] = uniShr[0][j];
        }
    }
    for (int i = 1; i < num_stand; i++) {
        for (int j = 0; j < num_data; j++) {
            interval_indicate[i][j] = uniShr[i][j] ^ uniShr[i - 1][j];
        }
    }

    // cout << "interval_indicate" << endl;
    // for (int i = 0; i < num_stand; i++) {
    //     cout << endl;
    //     for (int j = 0; j < num_data; j++) {
    //         cout << (int) interval_indicate[i][j] << ", ";
    //     }
    // }
    // cout << endl;


    uint64_t *temp_res = new uint64_t[num_data - 1]();
    // MUX, can batch   little gain as we set num_stand not too large
    for (int i = 0; i < num_stand; i++) {
        aux->multiplexer_two_plain(interval_indicate[i], (i + 1), temp_res, num_data - 1, bw_data, bw_data);
#pragma omp parallel for
        for (int j = 0; j < num_data - 1; j++) {
            res[j] += temp_res[j];
            // res[j] &= (1ULL << bw_data) - 1;
        }
    }


    // cout << "sorted_res" << endl;
    // for (int j = 0; j < num_data - 1; j++) {
    //     cout << (int) res[j] << ", ";
    // }
    //
    // cout << endl;

    delete[] endpoints;
    delete[] t;
    delete[] appendix;
    delete[] seeds;
    // delete[] prg;
    for (int i = 0; i < num_stand; i++) {
        // delete[] x[i];
        delete[] uniShr[i];
        delete[] interval_indicate[i];
    }
    // delete[] x;
    delete[] uniShr;
    delete[] interval_indicate;
}



// slower version
// void Frequency::count_sort(uint64_t *res, uint64_t *frequency, int num_stand, int num_data,
//                            int32_t bw_data, int32_t bw_res) {
//     int n_th = omp_get_num_procs() / 4;
//     cout << "n_th = " << n_th << endl;
//     omp_set_num_threads(n_th);
//     uint64_t *endpoints = new uint64_t[num_stand]();
//     endpoints[0] = frequency[0];
//     for (int i = 1; i < num_stand; i++) {
//         endpoints[i] = endpoints[i - 1] + frequency[i];
//         endpoints[i] &= (1ULL << bw_res) - 1;
//     }
//
//     // cout << endl;
//     // for (int i = 0; i < num_stand; i++) {
//     //     cout << endpoints[i] << " ";
//     // }
//     // cout << endl;
//
//     num_data += 1;
//
//     // share conversion 2^bw_data -> num_data
//     uint8_t *t = new uint8_t[num_stand]();
//     uint64_t *appendix = new uint64_t[num_stand]();
//     this->aux->wrap_computation(endpoints, t, num_stand, bw_res);
//     std::vector<uint64_t> v(1, num_data);
//     this->aux->B2A_coprimes(t, appendix, num_stand, v);
//     // cout << endl;
//     // for (int i = 0; i < num_stand; i++) {
//     //     cout << (int) appendix[i] << " ";
//     // }
//     // cout << endl;
// #pragma omp parallel for
//     for (int i = 0; i < num_stand; i++) {
//         endpoints[i] = (endpoints[i] % num_data) + (num_data - (1ULL << bw_res) % num_data) * (appendix[i]);
//         endpoints[i] %= num_data;
//     }
//
//     // cout << num_data << " ";
//     // cout << endl;
//     // for (int i = 0; i < num_stand; i++) {
//     //     cout << endpoints[i] << " ";
//     // }
//     // cout << endl;
//
//     uint64_t mask_data = (bw_data == 64 ? -1 : ((1ULL << bw_data) - 1));
//     uint64_t mask_res = (bw_res == 64 ? -1 : ((1ULL << bw_res) - 1));
//     auto *seeds = new block128[num_stand * num_data * 2]();
//     // auto *prg = new PRG128[num_stand * num_data * 2]();
//     PRG128 prg;
//
//     auto **uniShr = new uint8_t *[num_stand];
//     for (int i = 0; i < num_stand; i++) {
//         uniShr[i] = new uint8_t[num_data]();
//     }
//     auto **interval_indicate = new uint8_t *[num_stand];
//     for (int i = 0; i < num_stand; i++) {
//         interval_indicate[i] = new uint8_t[num_data]();
//     }
//     // int save_memory = num_data * 2 > 2048 ? 2048 : num_data * 2;
//
//     // int save_memory = num_data * 2 * num_data * 2 > 11ULL * 1024 * 1024 * 1024 ? 11ULL * 1024 * 1024 * 1024 / (num_data * 2) : num_data * 2;
//     // save_memory = save_memory > num_data * 2 ? num_data * 2 : save_memory;
//
//     uint64_t mem_limit = 11ULL * 1024 * 1024 * 1024;  // 11 GiB
//     uint64_t total = static_cast<uint64_t>(num_data) * 2 * num_data * 2;
//
//     uint64_t save_memory64 = total > mem_limit
//         ? mem_limit / (static_cast<uint64_t>(num_data) * 2)
//         : static_cast<uint64_t>(num_data) * 2;
//
//     uint64_t save_memory = save_memory64;
//
//     cout << "save_memory = " << save_memory << endl;
//
//     if (party == sci::ALICE) {
//         // set the offset directly in the index
//
//
//         // do not need to generate a 2*num_data seed, we can set the second half seed directly as offset < num_data
//         this->aux->nMinus1OUTNOT_batch(seeds, num_stand, 2 * num_data, nullptr);
//         // generate the shift translation shares
//
//         // cout << "------" << endl;
//         // cout << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         print_block_(seeds[i * num_data * 2 + j]);
//         //     }
//         // }
//         // cout << endl;
//
//         auto **shift_translate = new uint8_t *[num_data * 2];
//         for (int i = 0; i < num_data * 2; i++) {
//             shift_translate[i] = new uint8_t[save_memory]();
//         }
//         auto **a = new uint8_t *[num_stand];
//         for (int i = 0; i < num_stand; i++) {
//             a[i] = new uint8_t[num_data * 2]();
//         }
//         auto **b = new uint8_t *[num_stand];
//         for (int i = 0; i < num_stand; i++) {
//             b[i] = new uint8_t[num_data * 2]();
//         }
//
//         // for (int i = 0; i < num_stand; i++) {
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         prg[i * num_data * 2 + j].reseed(&seeds[i * num_data * 2 + j]);
//         //     }
//         // }
//
//         for (int i = 0; i < num_stand; i++) {
//             for (int k = 0; k < (num_data * 2 / save_memory); k++) {
//                 for (int j = 0; j < num_data * 2; j++) {
//                     prg.reseed(&seeds[i * num_data * 2 + j]);
//                     prg.random_bool((bool *) shift_translate[j], save_memory);
//                 }
//
//                 // cout << "shift_translate" << endl;
//                 // for (int j = 0; j < num_data * 2; j++) {
//                 //     for (int l = 0; l < save_memory; l++) {
//                 //         cout << (int) shift_translate[j][l] << " ";
//                 //     }
//                 //     cout << endl;
//                 // }
//
//                 for (int j = 0; j < num_data * 2; j++) {
//                     for (int l = 0; l < save_memory; l++) {
//                         a[i][l + k * save_memory] ^= shift_translate[j][l];
//                         b[i][j] ^= shift_translate[(j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
//                     }
//                 }
//             }
//             int last = num_data * 2 % save_memory;
//             if (last == 0) {
//                 continue;
//             }
//             int k = (num_data * 2 / save_memory);
//             for (int j = 0; j < num_data * 2; j++) {
//                 prg.reseed(&seeds[i * num_data * 2 + j]);
//                 prg.random_bool((bool *) shift_translate[j], last);
//             }
//             // cout << "shift_translate" << endl;
//             // for (int j = 0; j < num_data * 2; j++) {
//             //     for (int l = 0; l < last; l++) {
//             //         cout << (int) shift_translate[j][l] << " ";
//             //     }
//             //     cout << endl;
//             // }
//             for (int j = 0; j < num_data * 2; j++) {
//                 for (int l = 0; l < last; l++) {
//                     a[i][l + k * save_memory] ^= shift_translate[j][l];
//                     b[i][j] ^= shift_translate[(j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
//                 }
//             }
//         }
//
//         // cout << endl;
//         // cout << "----------x---------------" << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         cout << (int) x[i][j] << ", ";
//         //     }
//         //     cout << endl;
//         // }
//
//         // cout << endl;
//         // cout << "----------a---------------" << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         cout << (int) a[i][j] << ", ";
//         //     }
//         //     cout << endl;
//         // }
//
//         auto **x = new uint8_t *[num_stand];
//         for (int i = 0; i < num_stand; i++) {
//             x[i] = new uint8_t[num_data * 2]();
//         }
//         for (int i = 0; i < num_stand; i++) {
//             std::fill(x[i], x[i] + endpoints[i], 1);
//             std::fill(x[i] + num_data + endpoints[i], x[i] + 2 * num_data, 1);
//         }
//
//         for (int i = 0; i < num_stand; i++) {
// #pragma omp parallel for
//             for (int j = 0; j < num_data * 2; j++) {
//                 x[i][j] = x[i][j] ^ a[i][j];
//             }
//         }
//
//         // cout << endl;
//         // cout << "----------xMa---------------" << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         cout << (int) x[i][j] << ", ";
//         //     }
//         //     cout << endl;
//         // }
//
//
//         // cout << endl;
//         // cout << "------------b-------------" << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         cout << (int) b[i][j] << ", ";
//         //     }
//         //     cout << endl;
//         // }
//
//
//         uint8_t *x_packed = new uint8_t[(2 * num_data * num_stand + 7) / 8]();
//         pack_bits(x, x_packed, num_stand, 2 * num_data);
//         iopack->io->send_data(x_packed, ((2 * num_data * num_stand + 7) / 8) * sizeof(uint8_t));
//
//         // allocate the uniShr, and the
//         auto **mux_choice = new uint8_t *[num_stand];
//         for (int i = 0; i < num_stand; i++) {
//             mux_choice[i] = new uint8_t[num_data]();
// #pragma omp parallel for
//             for (int j = 0; j < num_data; j++) {
//                 mux_choice[i][j] = b[i][j + num_data] ^ b[i][j];
//             }
//         }
//         uint64_t *temp = new uint64_t[num_stand];
// #pragma omp parallel for
//         for (int i = 0; i < num_stand; i++) {
//             temp[i] = num_data - endpoints[i];
//         }
//
//         this->aux->mill->compare(t, temp, num_stand, ceil(log2(num_data + 1)));
//
//         // the result is >=
//         // cout << "t: " << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     cout << (int)t[i] << " ";
//         // }
//         // cout << endl;
//
//
//         // mux_choice padding a 0-string
//         this->aux->multiplexer_bShr(t, mux_choice, uniShr, num_stand, num_data);
// #pragma omp parallel for
//         for (int i = 0; i < num_stand; i++) {
//             for (int j = 0; j < num_data; j++) {
//                 uniShr[i][j] ^= b[i][j];
//             }
//         }
//
//         for (int i = 0; i < num_data * 2; i++) {
//             delete[] shift_translate[i];
//         }
//         delete[] shift_translate;
//         for (int i = 0; i < num_stand; i++) {
//             delete[] a[i];
//             delete[] b[i];
//             delete[] x[i];
//             delete[] mux_choice[i];
//         }
//         delete[] x_packed;
//         delete[] x;
//         delete[] temp;
//         delete[] mux_choice;
//         delete[] a;
//         delete[] b;
//     }
//     else {
//         this->aux->nMinus1OUTNOT_batch(seeds, num_stand, 2 * num_data, endpoints);
//
//         // cout << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     cout << endpoints[i] << endl;
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         print_block_(seeds[i * num_data * 2 + j]);
//         //     }
//         // }
//         // cout << endl;
//
//         // generate the shift translation shares
//         auto **shift_translate = new uint8_t *[num_data * 2];
//         for (int i = 0; i < num_data * 2; i++) {
//             shift_translate[i] = new uint8_t[save_memory]();
//         }
//         auto **c = new uint8_t *[num_stand];
//         for (int i = 0; i < num_stand; i++) {
//             c[i] = new uint8_t[num_data * 2]();
//         }
//
//         // for (int i = 0; i < num_stand; i++) {
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         prg.reseed(&seeds[i * num_data * 2 + j]);
//         //     }
//         // }
//
//         for (int i = 0; i < num_stand; i++) {
//             for (int k = 0; k < (num_data * 2 / save_memory); k++) {
//                 for (int j = 0; j < num_data * 2; j++) {
//                     prg.reseed(&seeds[i * num_data * 2 + j]);
//                     prg.random_bool((bool *) shift_translate[j], save_memory);
//                 }
//                 // cout << "shift_translate" << endl;
//                 // for (int j = 0; j < num_data * 2; j++) {
//                 //     for (int l = 0; l < save_memory; l++) {
//                 //         cout << (int) shift_translate[j][l] << " ";
//                 //     }
//                 //     cout << endl;
//                 // }
//                 for (int j = 0; j < num_data * 2; j++) {
//                     for (int l = 0; l < save_memory; l++) {
//                         c[i][j] ^= shift_translate[(j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
//                         c[i][(l + k * save_memory + endpoints[i] + 2 * num_data) % (2 * num_data)] ^= shift_translate[
//                             j][l];
//                     }
//                 }
//             }
//             int last = num_data * 2 % save_memory;
//             if (last == 0) {
//                 continue;
//             }
//             int k = (num_data * 2 / save_memory);
//             for (int j = 0; j < num_data * 2; j++) {
//                 if (j == endpoints[i]) {
//                     std::fill(shift_translate[j], shift_translate[j] + save_memory, 0);
//                     continue;
//                 }
//                 prg.reseed(&seeds[i * num_data * 2 + j]);
//                 prg.random_bool((bool *) shift_translate[j], last);
//             }
//             // cout << "shift_translate" << endl;
//             // for (int j = 0; j < num_data * 2; j++) {
//             //     for (int l = 0; l < last; l++) {
//             //         cout << (int) shift_translate[j][l] << " ";
//             //     }
//             //     cout << endl;
//             // }
//             for (int j = 0; j < num_data * 2; j++) {
//                 for (int l = 0; l < last; l++) {
//                     c[i][j] ^= shift_translate[(j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
//                     c[i][(l + k * save_memory + endpoints[i] + 2 * num_data) % (2 * num_data)] ^= shift_translate[j][l];
//                 }
//             }
//         }
//
//         // cout << endl;
//         // cout << "------------c-------------" << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         cout << (int) c[i][j] << ", ";
//         //     }
//         //     cout << endl;
//         // }
//
//         auto **xMa = new uint8_t *[num_stand];
//         for (int i = 0; i < num_stand; i++) {
//             xMa[i] = new uint8_t[num_data * 2]();
//         }
//
//         auto **shift_xMa = new uint8_t *[num_stand];
//         for (int i = 0; i < num_stand; i++) {
//             shift_xMa[i] = new uint8_t[num_data * 2]();
//         }
//         auto *xMa_packed = new uint8_t[(2 * num_data * num_stand + 7) / 8];
//         iopack->io->recv_data(xMa_packed, ((2 * num_data * num_stand + 7) / 8) * sizeof(uint8_t));
//         unpack_bits(xMa_packed, xMa, num_stand, 2 * num_data);
//
//         // cout << endl;
//         // cout << "------------xMa-------------" << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         cout << (int) xMa[i][j] << ", ";
//         //     }
//         //     cout << endl;
//         // }
//
//         for (int i = 0; i < num_stand; i++) {
//             std::memcpy(shift_xMa[i], xMa[i] + (2 * num_data - endpoints[i]), endpoints[i]);
//             std::memcpy(shift_xMa[i] + endpoints[i], xMa[i], 2 * num_data - endpoints[i]);
//         }
//
//         // cout << endl;
//         // cout << "------------shift_xMa-------------" << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         cout << (int) shift_xMa[i][j] << ", ";
//         //     }
//         //     cout << endl;
//         // }
//
//         for (int i = 0; i < num_stand; i++) {
// #pragma omp parallel for
//             for (int j = 0; j < num_data * 2; j++) {
//                 xMa[i][j] = c[i][j] ^ shift_xMa[i][j];
//             }
//         }
//
//
//         // cout << endl;
//         // cout << "---------c xor shift_xMa-----------" << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     for (int j = 0; j < num_data * 2; j++) {
//         //         cout << (int) xMa[i][j] << ", ";
//         //     }
//         //     cout << endl;
//         // }
//
//
//         auto **mux_choice = new uint8_t *[num_stand];
//         for (int i = 0; i < num_stand; i++) {
//             mux_choice[i] = new uint8_t[num_data]();
// #pragma omp parallel for
//             for (int j = 0; j < num_data; j++) {
//                 mux_choice[i][j] = xMa[i][j] ^ xMa[i][j + num_data];
//             }
//         }
//
//         this->aux->mill->compare(t, endpoints, num_stand, ceil(log2(num_data + 1)));
//
//         // cout << "t: " << endl;
//         // for (int i = 0; i < num_stand; i++) {
//         //     cout << (int)t[i] << " ";
//         // }
//         // cout << endl;
//
//
//         this->aux->multiplexer_bShr(t, mux_choice, uniShr, num_stand, num_data);
//         for (int i = 0; i < num_stand; i++) {
// #pragma omp parallel for
//             for (int j = 0; j < num_data; j++) {
//                 uniShr[i][j] ^= xMa[i][j];
//             }
//         }
//
//         for (int i = 0; i < 2 * num_data; i++) {
//             delete[] shift_translate[i];
//         }
//         delete[] shift_translate;
//         for (int i = 0; i < num_stand; i++) {
//             delete[] c[i];
//             delete[] mux_choice[i];
//             delete[] xMa[i];
//             delete[] shift_xMa[i];
//         }
//         delete[] c;
//         delete[] xMa;
//         delete[] xMa_packed;
//         delete[] shift_xMa;
//         // delete[] temp;
//         delete[] mux_choice;
//     }
//
//
//     // cout << "temp checkpoint. " << endl;
//     //
//     // cout << "uniShr: " << endl;
//     // for (int i = 0; i < num_stand; i++) {
//     //     cout << endl;
//     //     for (int j = 0; j < num_data; j++) {
//     //         cout << (int) uniShr[i][j] << ", ";
//     //     }
//     // }
//     // cout << endl;
//
//     // xor unishr
//     if (party == ALICE) {
// #pragma omp parallel for
//         for (int j = 0; j < num_data; j++) {
//             interval_indicate[0][j] = 1 ^ uniShr[0][j];
//         }
//     } else {
// #pragma omp parallel for
//         for (int j = 0; j < num_data; j++) {
//             interval_indicate[0][j] = uniShr[0][j];
//         }
//     }
//     for (int i = 1; i < num_stand; i++) {
//         for (int j = 0; j < num_data; j++) {
//             interval_indicate[i][j] = uniShr[i][j] ^ uniShr[i - 1][j];
//         }
//     }
//
//     // cout << "interval_indicate" << endl;
//     // for (int i = 0; i < num_stand; i++) {
//     //     cout << endl;
//     //     for (int j = 0; j < num_data; j++) {
//     //         cout << (int) interval_indicate[i][j] << ", ";
//     //     }
//     // }
//     // cout << endl;
//
//
//     uint64_t *temp_res = new uint64_t[num_data - 1]();
//     // MUX, can batch   little gain as we set num_stand not too large
//     for (int i = 0; i < num_stand; i++) {
//         aux->multiplexer_two_plain(interval_indicate[i], (i + 1), temp_res, num_data - 1, bw_data, bw_data);
// #pragma omp parallel for
//         for (int j = 0; j < num_data - 1; j++) {
//             res[j] += temp_res[j];
//             // res[j] &= (1ULL << bw_data) - 1;
//         }
//     }
//
//
//     // cout << "sorted_res" << endl;
//     // for (int j = 0; j < num_data - 1; j++) {
//     //     cout << (int) res[j] << ", ";
//     // }
//     //
//     // cout << endl;
//
//     delete[] endpoints;
//     delete[] t;
//     delete[] appendix;
//     delete[] seeds;
//     // delete[] prg;
//     for (int i = 0; i < num_stand; i++) {
//         // delete[] x[i];
//         delete[] uniShr[i];
//         delete[] interval_indicate[i];
//     }
//     // delete[] x;
//     delete[] uniShr;
//     delete[] interval_indicate;
// }



void random_permutation(uint64_t *perm, int n) {
    // initialize with 0, 1, 2, ..., n-1
    for (int i = 0; i < n; ++i) {
        perm[i] = i;
    }

    // random device + generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Fisher–Yates shuffle
    for (int i = n - 1; i > 0; --i) {
        std::uniform_int_distribution<> dis(0, i);
        int j = dis(gen);
        std::swap(perm[i], perm[j]);
    }
}

void plain_shuffle(uint64_t *perm, uint64_t *data, int length) {
    std::vector<uint64_t> tmp(length);

    for (int i = 0; i < length; i++) {
        // cout << "perm[i]= " << perm[i] << endl;
        tmp[perm[i]] = data[i];
    }

    for (int i = 0; i < length; i++) {
        data[i] = tmp[i];
    }
}

void Frequency::oblivious_shuffle(uint64_t *res, uint64_t *perm, uint64_t *data, int num_data, int32_t bw_data) {
    // opv
    if (party == ALICE) {
        block128 *seeds = new block128[num_data * num_data]();
        uint64_t *shuffle_matrix = new uint64_t[num_data * num_data]();
        PRG128 prg;
        this->aux->nMinus1OUTNOT_batch(seeds, num_data, num_data, nullptr);

        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     print_block_(seeds[i]);
        // }
        // cout << endl;
        for (int i = 0; i < num_data; i++) {
            for (int j = 0; j < num_data; j++) {
                prg.reseed(&seeds[i * num_data + j]);
                prg.random_data(&shuffle_matrix[i * num_data + j], sizeof(uint64_t));
            }
        }
        delete[] seeds;

        // cout << endl;
        // cout << "shuffle_matrix: " << endl;
        // for (int i = 0; i < num_data; i++) {
        //     for (int j = 0; j < num_data; j++) {
        //         cout << (shuffle_matrix[i * num_data + j] & ((1ULL << 3) - 1)) << " ";
        //     }
        //     cout << endl;
        // }
        // cout << endl;

        uint64_t *a = new uint64_t[num_data]();
        uint64_t *b = new uint64_t[num_data](); // related to the -b in the protocol description
        for (int i = 0; i < num_data; i++) {
            for (int j = 0; j < num_data; j++) {
                a[i] += shuffle_matrix[i * num_data + j];
                b[j] -= shuffle_matrix[i * num_data + j];
            }
        }

        // cout << "a: " << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (a[i] & ((1ULL << 3) - 1)) << " ";
        // }
        // cout << endl;

        // cout << "b: " << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (b[i] & ((1ULL << 3) - 1)) << " ";
        // }
        // cout << endl;

        uint64_t *xMa = new uint64_t[num_data]();
        int packed_length = (num_data * bw_data + 63) / 64;
        uint64_t *xMa_packed = new uint64_t[packed_length]();
        for (int i = 0; i < num_data; i++) {
            xMa[i] = data[i] + a[i];
        }

        pack_data(xMa, xMa_packed, num_data, bw_data);
        iopack->io->send_data(xMa_packed, packed_length * sizeof(uint64_t));

        for (int i = 0; i < num_data; i++) {
            res[i] = b[i];
        }
        delete[] xMa;
        delete[] xMa_packed;
        delete[] shuffle_matrix;
        delete[] a;
        delete[] b;
    } else {
        block128 *seeds = new block128[num_data * num_data]();
        uint64_t *shuffle_matrix = new uint64_t[num_data * num_data]();
        this->aux->nMinus1OUTNOT_batch(seeds, num_data, num_data, perm);

        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     print_block_(seeds[i]);
        // }
        // cout << endl;
        PRG128 prg;
        for (int i = 0; i < num_data; i++) {
            for (int j = 0; j < num_data; j++) {
                prg.reseed(&seeds[i * num_data + j]);
                prg.random_data(&shuffle_matrix[i * num_data + j], sizeof(uint64_t));
            }
        }
        delete[] seeds;

        // cout << endl;
        // cout << "shuffle_matrix: " << endl;
        // for (int i = 0; i < num_data; i++) {
        //     for (int j = 0; j < num_data; j++) {
        //         cout << (shuffle_matrix[i * num_data + j] & ((1ULL << 3) - 1)) << " ";
        //     }
        //     cout << endl;
        // }
        // cout << endl;

        uint64_t *c = new uint64_t[num_data](); // c is b-\pi(a)
        for (int i = 0; i < num_data; i++) {
            for (int j = 0; j < num_data; j++) {
                c[perm[i]] -= shuffle_matrix[i * num_data + j];
                c[j] += shuffle_matrix[i * num_data + j];
            }
        }


        // cout << "c: " << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (c[i] & ((1ULL << 3) - 1)) << " ";
        // }
        // cout << endl;

        uint64_t *xMa = new uint64_t[num_data]();
        int packed_length = (num_data * bw_data + 63) / 64;
        uint64_t *xMa_packed = new uint64_t[packed_length]();
        iopack->io->recv_data(xMa_packed, packed_length * sizeof(uint64_t));
        unpack_data(xMa_packed, xMa, num_data, bw_data);
        plain_shuffle(perm, xMa, num_data);
        for (int i = 0; i < num_data; i++) {
            res[i] = xMa[i] + c[i];
        }

        delete[] xMa;
        delete[] xMa_packed;
        delete[] shuffle_matrix;
        delete[] c;
    }
}

void Frequency::oblivious_shuffle_reverse(uint64_t *res, uint64_t *perm, uint64_t *data, int num_data,
                                          int32_t bw_data) {
    // opv
    if (party == BOB) {
        block128 *seeds = new block128[num_data * num_data]();
        uint64_t *shuffle_matrix = new uint64_t[num_data * num_data]();
        PRG128 prg;
        this->aux->nMinus1OUTNOT_batch_reverse(seeds, num_data, num_data, nullptr);
        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     print_block_(seeds[i]);
        // }
        // cout << endl;
        for (int i = 0; i < num_data; i++) {
            for (int j = 0; j < num_data; j++) {
                prg.reseed(&seeds[i * num_data + j]);
                prg.random_data(&shuffle_matrix[i * num_data + j], sizeof(uint64_t));
            }
        }
        delete[] seeds;

        // cout << endl;
        // cout << "shuffle_matrix: " << endl;
        // for (int i = 0; i < num_data; i++) {
        //     for (int j = 0; j < num_data; j++) {
        //         cout << (shuffle_matrix[i * num_data + j] & ((1ULL << 3) - 1)) << " ";
        //     }
        //     cout << endl;
        // }
        // cout << endl;

        uint64_t *a = new uint64_t[num_data]();
        uint64_t *b = new uint64_t[num_data](); // related to the -b in the protocol description
        for (int i = 0; i < num_data; i++) {
            for (int j = 0; j < num_data; j++) {
                a[i] += shuffle_matrix[i * num_data + j];
                b[j] -= shuffle_matrix[i * num_data + j];
            }
        }

        // cout << "a: " << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (a[i] & ((1ULL << 3) - 1)) << " ";
        // }
        // cout << endl;

        // cout << "b: " << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (b[i] & ((1ULL << 3) - 1)) << " ";
        // }
        // cout << endl;

        uint64_t *xMa = new uint64_t[num_data]();
        int packed_length = (num_data * bw_data + 63) / 64;
        uint64_t *xMa_packed = new uint64_t[packed_length]();
        for (int i = 0; i < num_data; i++) {
            xMa[i] = data[i] + a[i];
        }

        pack_data(xMa, xMa_packed, num_data, bw_data);
        iopack->io->send_data(xMa_packed, packed_length * sizeof(uint64_t));


        for (int i = 0; i < num_data; i++) {
            res[i] = b[i];
        }
        delete[] xMa;
        delete[] xMa_packed;
        delete[] shuffle_matrix;
        delete[] a;
        delete[] b;
    } else {
        block128 *seeds = new block128[num_data * num_data]();
        uint64_t *shuffle_matrix = new uint64_t[num_data * num_data]();
        this->aux->nMinus1OUTNOT_batch_reverse(seeds, num_data, num_data, perm);
        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     print_block_(seeds[i]);
        // }
        // cout << endl;
        PRG128 prg;
        for (int i = 0; i < num_data; i++) {
            for (int j = 0; j < num_data; j++) {
                prg.reseed(&seeds[i * num_data + j]);
                prg.random_data(&shuffle_matrix[i * num_data + j], sizeof(uint64_t));
            }
        }
        delete[] seeds;

        // cout << endl;
        // cout << "shuffle_matrix: " << endl;
        // for (int i = 0; i < num_data; i++) {
        //     for (int j = 0; j < num_data; j++) {
        //         cout << (shuffle_matrix[i * num_data + j] & ((1ULL << 3) - 1)) << " ";
        //     }
        //     cout << endl;
        // }
        // cout << endl;

        uint64_t *c = new uint64_t[num_data](); // c is b-\pi(a)
        for (int i = 0; i < num_data; i++) {
            for (int j = 0; j < num_data; j++) {
                c[perm[i]] -= shuffle_matrix[i * num_data + j];
                c[j] += shuffle_matrix[i * num_data + j];
            }
        }

        // cout << "c: " << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (c[i] & ((1ULL << 3) - 1)) << " ";
        // }
        // cout << endl;

        uint64_t *xMa = new uint64_t[num_data]();
        int packed_length = (num_data * bw_data + 63) / 64;
        uint64_t *xMa_packed = new uint64_t[packed_length]();
        iopack->io->recv_data(xMa_packed, packed_length * sizeof(uint64_t));
        unpack_data(xMa_packed, xMa, num_data, bw_data);
        plain_shuffle(perm, xMa, num_data);
        for (int i = 0; i < num_data; i++) {
            res[i] = xMa[i] + c[i];
        }
        delete[] xMa;
        delete[] xMa_packed;
        delete[] shuffle_matrix;
        delete[] c;
    }
}


// Function to partition the array
int Frequency::partition(uint64_t *arr, int low, int high, int bw) {
    int pivot = arr[high]; // choose the last element as pivot
    int i = low - 1; // index of smaller element
    uint64_t mask_l = (1ULL << bw) - 1;
    uint8_t res;
    uint8_t res_t;
    uint64_t Za;
    uint64_t compare_with;
    for (int j = low; j < high; j++) {
        // cout << pivot << " " << arr[j] << endl;
        // conduct a private cmp
        if (party == sci::ALICE) {
            Za = (pivot - arr[j]) & mask_l;
            compare_with = Za & ((1 << (bw - 1)) - 1);
            // std::cout << "wa: " << compare_with << std::endl;
            aux->mill->compare(&res, &compare_with, 1, bw);
            // std::cout << "res: " << static_cast<uint32_t>(res) << std::endl;
            uint8_t ba = Za < (1 << (bw - 1));
            // std::cout << "ba: " << static_cast<uint32_t>(ba) << std::endl;
            res ^= ba;
            // std::cout << static_cast<uint32_t>(res) << " ";
            iopack->io->send_data(&res, sizeof(uint8_t));
            iopack->io->recv_data(&res_t, sizeof(uint8_t));
        } else {
            uint64_t Zb = (pivot - arr[j]) & mask_l;
            compare_with = (1 << (bw - 1)) - (Zb & ((1 << (bw - 1)) - 1));
            // std::cout << "t/2 - wb: " << compare_with << std::endl;
            aux->mill->compare(&res, &compare_with, 1, bw);
            // std::cout << "res: " << static_cast<uint32_t>(res) << std::endl;
            uint8_t bb = Zb < (1 << (bw - 1));
            // std::cout << "bb: " << static_cast<uint32_t>(bb) << std::endl;
            res ^= bb;
            res ^= 1;
            // std::cout << static_cast<uint32_t>(res) << " ";
            iopack->io->recv_data(&res_t, sizeof(uint8_t));
            iopack->io->send_data(&res, sizeof(uint8_t));
        }
        // cout << (int)(res_t ^ res) << endl;
        if (res_t ^ res) {
            i++;
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return i + 1;
}

// Quicksort function
void Frequency::quickSort(uint64_t *arr, int low, int high, int bw) {
    if (low < high) {
        int pi = partition(arr, low, high, bw); // partition index

        // Recursively sort elements before and after partition
        quickSort(arr, low, pi - 1, bw);
        quickSort(arr, pi + 1, high, bw);
    }
}


void Frequency::shuffle_sort(uint64_t *res, uint64_t *data, int num_stand, int num_data,
                             int32_t bw_data, int32_t bw_res) {
    // for (int i = 0; i < num_data; i++) {
    //     cout << data[i] << " ";
    // }
    // cout << endl;
    uint64_t *perm = new uint64_t[num_data]();
    int bw_append = ceil(log2(num_data + 1));
    int bw_tot = bw_data + bw_append;
    random_permutation(perm, num_data);
    // cout << "bw_tot = " << bw_tot << endl;
    // cout << "bw_data = " << bw_data << endl;
    // for (int i = 0; i < num_data; i++) {
    //     cout << perm[i] << " ";
    // }
    // cout << endl;
    if (party == ALICE) {
        // element extend to make each element distinct.
        for (int i = 0; i < num_data; i++) {
            data[i] = (data[i] << bw_append);
            data[i] += i;
        }
        // for (int i = 0; i < num_data; i++) {
        //     cout << perm[i] << " ";
        // }
        // cout << endl;
        oblivious_shuffle(res, nullptr, data, num_data, bw_tot);
        for (int i = 0; i < num_data; i++) {
            data[i] = res[i];
            // data[i] &= (1ULL << bw_tot) - 1;
        }

        // cout << "----------------" << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << perm[i] << " ";
        // }
        // cout << endl;

        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (data[i] & ((1ULL << bw_append) - 1)) << " ";
        // }
        // cout << endl;

        oblivious_shuffle_reverse(res, perm, nullptr, num_data, bw_tot);
        plain_shuffle(perm, data, num_data);
        for (int i = 0; i < num_data; i++) {
            data[i] = data[i] + res[i];
            // data[i] &= (1ULL << bw_tot) - 1;
        }

        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (data[i] & ((1ULL << bw_append) - 1)) << " ";
        // }
        // cout << endl;
    } else {
        for (int i = 0; i < num_data; i++) {
            data[i] = data[i] << bw_append;
        }
        // for (int i = 0; i < num_data; i++) {
        //     cout << perm[i] << " ";
        // }
        // cout << endl;
        oblivious_shuffle(res, perm, nullptr, num_data, bw_tot);
        // for (int i = 0; i < num_data; i++) {
        //     cout << perm[i] << " ";
        // }
        // cout << endl;
        // cout << "----------------" << endl;
        plain_shuffle(perm, data, num_data);
        for (int i = 0; i < num_data; i++) {
            data[i] = data[i] + res[i];
            // data[i] &= (1ULL << bw_tot) - 1;
        }


        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (res[i] & ((1ULL << bw_append) - 1)) << " ";
        // }
        // cout << endl;

        oblivious_shuffle_reverse(res, nullptr, data, num_data, bw_tot);
        for (int i = 0; i < num_data; i++) {
            data[i] = res[i];
            // data[i] &= (1ULL << bw_tot) - 1;
        }

        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (data[i] & ((1ULL << bw_append) - 1)) << " ";
        // }
        // cout << endl;
    }

    xt->z_extend(num_data, data, data, bw_tot, bw_tot + 1); // for compare, we need one more bit
    // cout << endl;
    // for (int i = 0; i < num_data; i++) {
    //     cout << (data[i] & ((1ULL << (bw_tot + 1)) - 1)) << " ";
    // }
    // cout << endl;
    // cout << endl;
    // for (int i = 0; i < num_data; i++) {
    //     cout << (data[i] & ((1ULL << bw_tot) - 1)) << " ";
    // }
    // cout << endl;
    // quick sort
    quickSort(data, 0, num_data - 1, bw_tot + 1);
    // cout << endl;
    // for (int i = 0; i < num_data; i++) {
    //     cout << (data[i] & ((1ULL << bw_tot) - 1)) << " ";
    // }
    // cout << endl;
    trunc->truncate_and_reduce(num_data, data, res, bw_append, bw_tot);
    // cout << endl;
    // for (int i = 0; i < num_data; i++) {
    //     cout << (res[i] & ((1ULL << bw_data) - 1)) << " ";
    // }
    // cout << endl;
    delete[] perm;
}


void Frequency::compare_and_swap(uint64_t *arr1, uint64_t *arr2, int length, int32_t bw_data, int32_t bw_res) {
    uint64_t *Za = new uint64_t[length]();
    uint64_t *Zb = new uint64_t[length]();

    uint64_t *compare_with = new uint64_t[length]();
    uint64_t *temp_mux = new uint64_t[length]();
    uint8_t *res = new uint8_t[length]();
    uint8_t *ba = new uint8_t[length]();
    uint8_t *bb = new uint8_t[length]();

    if (party == sci::ALICE) {
        for (int i = 0; i < length; i++) {
            Za[i] = (arr1[i] - arr2[i]) & (1 << (bw_data)) - 1;
            compare_with[i] = Za[i] & ((1 << (bw_data - 1)) - 1);
        }
        aux->mill->compare(res, compare_with, length, bw_data);
        for (int i = 0; i < length; i++) {
            ba[i] = Za[i] < (1ULL << (bw_data - 1));
            res[i] ^= ba[i];
            temp_mux[i] = (arr1[i] - arr2[i]) & (1ULL << (bw_data)) - 1;
        }
        aux->multiplexer(res, temp_mux, temp_mux, length, bw_data, bw_data);
        for (int i = 0; i < length; i++) {
            temp_mux[i] = temp_mux[i] + arr2[i];
            arr2[i] = arr1[i] + arr2[i] - temp_mux[i];
            arr1[i] = temp_mux[i];
        }
    } else {
        for (int i = 0; i < length; i++) {
            Zb[i] = (arr1[i] - arr2[i]) & (1ULL << (bw_data)) - 1;
            compare_with[i] = (1 << (bw_data - 1)) - (Zb[i] & ((1 << (bw_data - 1)) - 1));
        }
        aux->mill->compare(res, compare_with, length, bw_data);
        for (int i = 0; i < length; i++) {
            temp_mux[i] = (arr1[i] - arr2[i]) & (1ULL << (bw_data)) - 1;
            bb[i] = Zb[i] < (1ULL << (bw_data - 1));
            res[i] ^= bb[i];
            res[i] ^= 1;
        }
        aux->multiplexer(res, temp_mux, temp_mux, length, bw_data, bw_data);
        for (int i = 0; i < length; i++) {
            temp_mux[i] = temp_mux[i] + arr2[i];
            arr2[i] = arr1[i] + arr2[i] - temp_mux[i];
            arr1[i] = temp_mux[i];
        }
    }
}


// Merge two sorted halves using your compare_and_swap(arr1, arr2, len)
void Frequency::batcher_merge(uint64_t *data, int num_stand, int num_data,
                              int32_t bw_data, int32_t bw_res) {
    if (num_data <= 1) return;

    int mid = num_data / 2;

    // Pairwise compare and swap between the two halves
    compare_and_swap(data, data + mid, mid, bw_data, bw_res);

    // Recursively merge inside each half
    batcher_merge(data, num_stand, mid, bw_data, bw_res);
    batcher_merge(data + mid, num_stand, mid, bw_data, bw_res);
}

// before call this function, we should extend the bit length of the data.
void Frequency::batcher_sort(uint64_t *data, int num_stand, int num_data,
                             int32_t bw_data, int32_t bw_res) {
    if (num_data <= 1) return;

    int mid = num_data / 2;

    // Recursively sort both halves
    batcher_sort(data, num_stand, mid, bw_data, bw_res);
    batcher_sort(data + mid, num_stand, mid, bw_data, bw_res);

    // Merge the two sorted halves
    batcher_merge(data, num_stand, num_data, bw_data, bw_res);
}

void Frequency::batcher_network_sort(uint64_t *data, int num_stand, int num_data,
                                     int32_t bw_data, int32_t bw_res) {
    if (num_data <= 1) return;
    xt->z_extend(num_data, data, data, bw_data, bw_data);
    batcher_sort(data, num_stand, num_data, bw_data + 1, bw_res);
}


// Partition function for top-k selection
int Frequency::partitionTopK(uint64_t *arr, int low, int high, int bw) {
    uint64_t pivot = arr[high]; // choose last element as pivot
    int i = low - 1;
    uint64_t mask_l = (1ULL << bw) - 1;
    uint8_t res, res_t;
    uint64_t compare_with;
    uint64_t Za;

    for (int j = low; j < high; j++) {
        if (party == sci::ALICE) {
            Za = (arr[j] - pivot) & mask_l; // flip comparison to get top-k
            compare_with = Za & ((1ULL << (bw - 1)) - 1);
            aux->mill->compare(&res, &compare_with, 1, bw);
            uint8_t ba = Za < (1ULL << (bw - 1));
            res ^= ba;
            iopack->io->send_data(&res, sizeof(uint8_t));
            iopack->io->recv_data(&res_t, sizeof(uint8_t));
        } else {
            uint64_t Zb = (arr[j] - pivot) & mask_l;
            compare_with = (1ULL << (bw - 1)) - (Zb & ((1ULL << (bw - 1)) - 1));
            aux->mill->compare(&res, &compare_with, 1, bw);
            uint8_t bb = Zb < (1ULL << (bw - 1));
            res ^= bb;
            res ^= 1;
            iopack->io->recv_data(&res_t, sizeof(uint8_t));
            iopack->io->send_data(&res, sizeof(uint8_t));
        }

        if (res_t ^ res) {
            i++;
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return i + 1;
}

// Top-k select function: modifies arr in-place so top k elements are at arr[0..k-1]
void Frequency::topKSelect(uint64_t *arr, int low, int high, int k, int bw) {
    if (low <= high) {
        int pi = partitionTopK(arr, low, high, bw);
        int left_count = pi - low + 1; // number of elements >= pivot

        if (left_count == k) {
            return; // top-k elements are arr[low..pi]
        } else if (left_count > k) {
            topKSelect(arr, low, pi - 1, k, bw); // top k are in left partition
        } else {
            topKSelect(arr, pi + 1, high, k - left_count, bw); // need more from right partition
        }
    }
}


void Frequency::shuffle_topk(uint64_t *res, uint64_t *data, int k, int num_data,
                             int32_t bw_data, int32_t bw_res) {
    // for (int i = 0; i < num_data; i++) {
    //     cout << data[i] << " ";
    // }
    // cout << endl;
    uint64_t *perm = new uint64_t[num_data]();
    int bw_append = ceil(log2(num_data + 1));
    int bw_tot = bw_data + bw_append;
    random_permutation(perm, num_data);
    // cout << "bw_tot = " << bw_tot << endl;
    // cout << "bw_data = " << bw_data << endl;
    // for (int i = 0; i < num_data; i++) {
    //     cout << perm[i] << " ";
    // }
    // cout << endl;
    if (party == ALICE) {
        // element extend to make each element distinct.
        for (int i = 0; i < num_data; i++) {
            data[i] = (data[i] << bw_append);
            data[i] += i;
        }
        // for (int i = 0; i < num_data; i++) {
        //     cout << perm[i] << " ";
        // }
        // cout << endl;
        oblivious_shuffle(res, nullptr, data, num_data, bw_tot);
        for (int i = 0; i < num_data; i++) {
            data[i] = res[i];
            // data[i] &= (1ULL << bw_tot) - 1;
        }

        // cout << "----------------" << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << perm[i] << " ";
        // }
        // cout << endl;

        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (data[i] & ((1ULL << bw_append) - 1)) << " ";
        // }
        // cout << endl;

        oblivious_shuffle_reverse(res, perm, nullptr, num_data, bw_tot);
        plain_shuffle(perm, data, num_data);
        for (int i = 0; i < num_data; i++) {
            data[i] = data[i] + res[i];
            // data[i] &= (1ULL << bw_tot) - 1;
        }

        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (data[i] & ((1ULL << bw_append) - 1)) << " ";
        // }
        // cout << endl;
    } else {
        for (int i = 0; i < num_data; i++) {
            data[i] = data[i] << bw_append;
        }
        // for (int i = 0; i < num_data; i++) {
        //     cout << perm[i] << " ";
        // }
        // cout << endl;
        oblivious_shuffle(res, perm, nullptr, num_data, bw_tot);
        // for (int i = 0; i < num_data; i++) {
        //     cout << perm[i] << " ";
        // }
        // cout << endl;
        // cout << "----------------" << endl;
        plain_shuffle(perm, data, num_data);
        for (int i = 0; i < num_data; i++) {
            data[i] = data[i] + res[i];
            // data[i] &= (1ULL << bw_tot) - 1;
        }


        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (res[i] & ((1ULL << bw_append) - 1)) << " ";
        // }
        // cout << endl;

        oblivious_shuffle_reverse(res, nullptr, data, num_data, bw_tot);
        for (int i = 0; i < num_data; i++) {
            data[i] = res[i];
            // data[i] &= (1ULL << bw_tot) - 1;
        }

        // cout << endl;
        // for (int i = 0; i < num_data; i++) {
        //     cout << (data[i] & ((1ULL << bw_append) - 1)) << " ";
        // }
        // cout << endl;
    }

    xt->z_extend(num_data, data, data, bw_tot, bw_tot + 1); // for compare, we need one more bit
    // cout << endl;
    // for (int i = 0; i < num_data; i++) {
    //     cout << (data[i] & ((1ULL << (bw_tot + 1)) - 1)) << " ";
    // }
    // cout << endl;
    // cout << endl;
    // for (int i = 0; i < num_data; i++) {
    //     cout << (data[i] & ((1ULL << bw_tot) - 1)) << " ";
    // }
    // cout << endl;
    // quick sort
    topKSelect(data, 0, num_data - 1, k, bw_tot + 1);
    // cout << endl;
    // for (int i = 0; i < num_data; i++) {
    //     cout << (data[i] & ((1ULL << bw_tot) - 1)) << " ";
    // }
    // cout << endl;
    trunc->truncate_and_reduce(k, data, res, bw_append, bw_tot);
    // cout << endl;
    // for (int i = 0; i < num_data; i++) {
    //     cout << (res[i] & ((1ULL << bw_data) - 1)) << " ";
    // }
    // cout << endl;
    delete[] perm;
}

// tree
void Frequency::count_kth(uint64_t *res, uint64_t *frequency, int k, int num_stand, int num_data,
                          int32_t bw_data, int32_t bw_res) {
    // bw_res leads to different network pack efficiency !
    // cout << endl << "K: " << k << endl;
    res[0] = 0;
    bool flag = false;
    if (num_stand % 2 == 0) {
        num_stand += 1;
        flag = true;
    }
    uint64_t *endpoints = new uint64_t[num_stand]();
    uint64_t *endpoints_temp = new uint64_t[(num_stand + 1) / 2]();
    uint64_t *endpoints_temp2 = new uint64_t[(num_stand + 1) / 2]();
    endpoints[0] = frequency[0];
    if (flag) {
        for (int i = 1; i < num_stand - 1; i++) {
            endpoints[i] = endpoints[i - 1] + frequency[i];
            // endpoints[i] &= (1ULL << bw_res) - 1;
        }
        xt->z_extend(num_stand - 1, endpoints, endpoints, bw_res, bw_res + 1);
        endpoints[num_stand - 1] = endpoints[num_stand - 2];
    } else {
        for (int i = 1; i < num_stand; i++) {
            endpoints[i] = endpoints[i - 1] + frequency[i];
            // endpoints[i] &= (1ULL << bw_res) - 1;
        }
        xt->z_extend(num_stand, endpoints, endpoints, bw_res, bw_res + 1);
    }
    bw_res += 1;

    // cout << "endpoints: " << endl;
    // for (int i = 0; i < num_stand; i++) {
    //     cout << endpoints[i] << ", ";
    // }
    // cout << endl;

    int depth = 0;
    int size = num_stand;
    int mid_idx = size / 2; // choose the middle element
    uint64_t temp = 0;
    uint64_t temptwo = 0;
    uint8_t *temp_t = new uint8_t[num_stand]();
    uint8_t temp_eq;
    uint8_t *res_tree = new uint8_t[num_stand];
    while (true) {
        if (size <= 3) {
            uint8_t *res_t = new uint8_t[size]();
            uint8_t *res_eq = new uint8_t[size]();
            uint8_t *res_t1 = new uint8_t[size]();
            uint8_t *res_tb2a = new uint8_t[size]();
            uint64_t *temp_t64 = new uint64_t[size]();
            uint64_t *temp_t264 = new uint64_t[size]();


            for (int i = 0; i < size; i++) {
                if (party == ALICE) {
                    temp_t64[i] = (endpoints[i] - k) & ((1ULL << (bw_res)) - 1);
                    temp_t264[i] = temp_t64[i] & ((1ULL << (bw_res - 1)) - 1);
                } else {
                    temp_t64[i] = endpoints[i] & ((1ULL << (bw_res)) - 1);
                    temp_t264[i] = (1ULL << (bw_res - 1)) - (temp_t64[i] & ((1ULL << (bw_res - 1)) - 1));
                }
            }

            aux->mill_and_eq->compare_with_eq(res_t, res_eq, temp_t264, size, bw_res);

            for (int i = 0; i < size; i++) {
                if (party == ALICE) {
                    res_t[i] ^= 1; // this used to reverse the cmp result to get wa <= t/2 - wb
                    res_t[i] ^= (temp_t64[i] < (1ULL << (bw_res - 1)));
                } else {
                    res_t[i] ^= (temp_t64[i] < (1ULL << (bw_res - 1)));
                }
            }

            for (int i = 0; i < size; i++) {
                if (party == ALICE) {
                    res_t[i] ^= 1;
                    res_eq[i] ^= 1;
                }
            }

            aux->AND(res_t, res_eq, res_t, size);

            for (int i = 0; i < size; i++) {
                if (party == ALICE) {
                    res_t[i] ^= 1;
                }
            }

            // cout << "---------------" << endl;
            // for (int i = 0; i < size; i++) {
            //     cout << (int) res_t[i] << ", ";
            // }
            // cout << endl;

            if (party == ALICE) {
                for (int i = 0; i < size - 1; i++) {
                    res_t1[i] = res_t[i] ^ 1;
                }
                res_tb2a[0] = res_t[0];
            } else {
                for (int i = 0; i < size - 1; i++) {
                    res_t1[i] = res_t[i];
                }
                res_tb2a[0] = res_t[0];
            }


            // res_t: k <= x0, k <= x1, k <= x2
            // res_t1: k > x0, k > x2
            aux->AND(res_t + 1, res_t1, res_tb2a + 1, size - 1);
            // res_tb2a: k <= x0 ; x0 < k <= x1; x1 < k <= x2;

            // cout << "---------------" << endl;
            // for (int i = 0; i < size; i++) {
            //     cout << (int) res_tb2a[i] << ", ";
            // }
            // cout << endl;

            aux->B2A(res_tb2a, temp_t64, size, bw_res - 1); // bw_res should equals to the log(num_stand)

            for (int i = 0; i < size; i++) {
                res[0] += (temp_t64[i] * i);
            }
            // res[0] &= (1ULL << (bw_res - 1)) - 1;
            delete[] temp_t64;
            delete[] temp_t264;
            delete[] res_t;
            delete[] res_t1;
            delete[] res_tb2a;
            break;
        }

        if (party == ALICE) {
            temp = (endpoints[mid_idx] - k) & ((1ULL << (bw_res)) - 1);
            temptwo = temp & ((1ULL << (bw_res - 1)) - 1);
            aux->mill_and_eq->compare_with_eq(&res_tree[depth], &temp_eq, &temptwo, 1, bw_res);
            res_tree[depth] ^= (temp < (1ULL << (bw_res - 1)));
        } else {
            temp = endpoints[mid_idx] & ((1ULL << (bw_res)) - 1);
            temptwo = (1ULL << (bw_res - 1)) - (temp & ((1ULL << (bw_res - 1)) - 1));
            aux->mill_and_eq->compare_with_eq(&res_tree[depth], &temp_eq, &temptwo, 1, bw_res);
            res_tree[depth] ^= 1;
            res_tree[depth] ^= (temp < (1ULL << (bw_res - 1)));
        }

        if (party == ALICE) {
            res_tree[depth] ^= 1;
            temp_eq ^= 1;
        }

        aux->AND(&res_tree[depth], &temp_eq, &res_tree[depth], 1);

        if (party == ALICE) {
            res_tree[depth] ^= 1;
        }

        // cout << endl;
        // cout << "mid_value: " << endpoints[mid_idx] << endl;
        // cout << "tree decision: " << (int) res_tree[depth] << endl;
        // cout << endl;

        // mux two array with length (size + 1) / 2
        for (int j = 0; j < (size + 1) / 2; j++) {
            endpoints_temp[j] = endpoints[size / 2 + j] - endpoints[j];
        }
        if (party == ALICE) {
            std::fill(temp_t, temp_t + (size + 1) / 2, res_tree[depth]);
        } else {
            std::fill(temp_t, temp_t + (size + 1) / 2, res_tree[depth] ^ 1);
            // we change for the convenience of array mux temp_t = 0 to choose the left half tree
        }
        // we here use k independent mux, which leads to k\lambda, but we can actually reduce it to \lambda
        aux->multiplexer(temp_t, endpoints_temp, endpoints_temp2, (size + 1) / 2, bw_res, bw_res);

        for (int j = 0; j < (size + 1) / 2; j++) {
            endpoints[j] += endpoints_temp2[j];
            // endpoints[j] &= (1ULL << bw_res) - 1;
        }
        size = (size + 1) / 2;
        if (size % 2 == 0) {
            endpoints[size] = endpoints[size - 1];
            size += 1;
        }
        mid_idx = size / 2;
        depth += 1;

        // cout << endl;
        // for (int j = 0; j < size; j++) {
        //     cout << endpoints[j] << ", ";
        // }
        // cout << endl;
    }

    // cout << endl;
    // cout << "res[0]: " << res[0] << endl;

    // temp2 used for the tree decision value append
    cout << "depth: " << depth << endl;
    uint64_t *temp2 = new uint64_t[depth]();
    uint64_t *temp3 = new uint64_t[depth]();
    int t = 0;
    size = num_stand;
    while (true) {
        if (size <= 3) {
            break;
        }
        temp2[t] = size / 2;
        size = (size + 1) / 2;
        if (size % 2 == 0) {
            size += 1;
        }
        t += 1;
    }

    if (party == ALICE) {
        for (int i1 = 0; i1 < depth; i1++) {
            res_tree[i1] ^= 1;
        }
    }


    // cout << endl << "append value: " << endl;
    // for (int i1 = 0; i1 < depth; i1++) {
    //     cout << temp2[i1] << ", ";
    // }
    // cout << endl;

    // cout << endl << "res_tree: " << endl;
    // for (int i1 = 0; i1 < depth; i1++) {
    //     cout << (int) res_tree[i1] << ", ";
    // }
    // cout << endl;

    aux->B2A(res_tree, temp3, depth, bw_res - 1);

    // cout << endl << "temp3: " << endl;
    // for (int i1 = 0; i1 < depth; i1++) {
    //     cout << temp3[i1] << ", ";
    // }
    // cout << endl;

    for (int i = 0; i < depth; i++) {
        res[0] += (temp2[i] * temp3[i]);
    }


    // res[0] &= (1ULL << (bw_res - 1)) - 1;


    // cout << endl;
    // cout << "res[0]: " << res[0] << endl;


    delete [] endpoints;
    delete [] endpoints_temp;
    delete [] endpoints_temp2;
    delete[] temp_t;
    delete [] res_tree;
    delete[] temp2;
    delete[] temp3;
}


// linear
// void Frequency::count_kth_linear(uint64_t *res, uint64_t *frequency, int k, int num_stand, int num_data,
//                           int32_t bw_data, int32_t bw_res) {
//     cout << endl << "K: " << k << endl;
//     bool flag = false;
//     if (num_stand % 2 == 0) {
//         num_stand += 1;
//         flag = true;
//     }
//     uint64_t *endpoints = new uint64_t[num_stand]();
//     uint64_t *endpoints_temp = new uint64_t[(num_stand + 1) / 2]();
//     uint64_t *endpoints_temp2 = new uint64_t[(num_stand + 1) / 2]();
//     endpoints[0] = frequency[0];
//     if (flag) {
//         for (int i = 1; i < num_stand - 1; i++) {
//             endpoints[i] = endpoints[i - 1] + frequency[i];
//             endpoints[i] &= (1ULL << bw_res) - 1;
//         }
//         xt->z_extend(num_stand - 1, endpoints, endpoints, bw_res, bw_res + 1);
//         endpoints[num_stand - 1] = endpoints[num_stand - 2];
//     } else {
//         for (int i = 1; i < num_stand; i++) {
//             endpoints[i] = endpoints[i - 1] + frequency[i];
//             endpoints[i] &= (1ULL << bw_res) - 1;
//         }
//         xt->z_extend(num_stand, endpoints, endpoints, bw_res, bw_res + 1);
//     }
//     bw_res += 1;
//
//     cout << "endpoints: " << endl;
//     for (int i = 0; i < num_stand; i++) {
//         cout << endpoints[i] << ", ";
//     }
//     cout << endl;
//
//     int depth = 0;
//     int size = num_stand;
//     int mid_idx = size / 2; // choose the middle element
//     uint64_t temp = 0;
//     uint8_t *temp_t = new uint8_t[num_stand]();
//     uint8_t *res_tree = new uint8_t[num_stand];
//     while (true) {
//         if (size <= 3) {
//             uint8_t *res_t = new uint8_t[size]();
//             uint64_t *temp_t64 = new uint64_t[size]();
//             for (int i = 0; i < size; i++) {
//                 if (party == ALICE) {
//                     temp_t64[i] = endpoints[i] - k;
//                 } else {
//                     temp_t64[i] = endpoints[i];
//                 }
//             }
//             aux->wrap_computation(temp_t64, res_t, size, bw_res - 1);
//             for (int i = 0; i < size; i++) {
//                 if (party == ALICE) {
//                     res_t[i] ^= ((endpoints[i] >> (bw_res - 1)) & 1);
//                 } else {
//                     res_t[i] ^= ((endpoints[i] >> (bw_res - 1)) & 1);
//                     res_t[i] ^= 1;
//                 }
//             }
//
//
//             cout << endl << "the last 3 elements: " << endl;
//             for (int i = 0; i < size; i++) {
//                 cout << (int)res_t[i] << ", ";
//             }
//             cout << endl;
//
//
//             aux->B2A(res_t, temp_t64, size, bw_res); // bw_res should equals to the log(num_stand)
//             for (int i = 0; i < size; i++) {
//                 res[0] += temp_t64[i];
//             }
//             res[0] &= (1ULL << bw_res) - 1;
//             delete[] temp_t64;
//             delete[] res_t;
//             break;
//         }
//         if (party == ALICE) {
//             temp = endpoints[mid_idx] - k;
//             aux->wrap_computation(&temp, &res_tree[depth], 1, bw_res - 1);
//             res_tree[depth] ^= ((endpoints[mid_idx] >> (bw_res - 1)) & 1);
//         } else {
//             temp = endpoints[mid_idx];
//             aux->wrap_computation(&temp, &res_tree[depth], 1, bw_res - 1);
//             res_tree[depth] ^= ((endpoints[mid_idx] >> (bw_res - 1)) & 1);
//             res_tree[depth] ^= 1;
//         }
//
//         cout << endl;
//         cout << "mid_value: " << endpoints[mid_idx] << endl;
//         cout << "tree decision: " << (int)res_tree[depth] << endl;
//         cout << endl;
//
//         // mux two array with length (size + 1) / 2
//         for (int j = 0; j < (size + 1) / 2; j++) {
//             endpoints_temp[j] = endpoints[size / 2 + j] - endpoints[j];
//         }
//         if (party == ALICE) {
//             std::fill(temp_t, temp_t + (size + 1) / 2, res_tree[depth]);
//         } else {
//             std::fill(temp_t, temp_t + (size + 1) / 2, res_tree[depth] ^ 1);
//         }
//         aux->multiplexer(temp_t, endpoints_temp, endpoints_temp2, (size + 1) / 2, bw_res, bw_res);
//
//         for (int j = 0; j < (size + 1) / 2; j++) {
//             endpoints[j] += endpoints_temp2[j];
//             endpoints[j] &= (1ULL << bw_res) - 1;
//         }
//         size = (size + 1) / 2;
//         if (size % 2 == 0) {
//             endpoints[size] = endpoints[size - 1];
//             size += 1;
//         }
//         mid_idx = size / 2;
//         depth += 1;
//
//         cout << endl;
//         for (int j = 0; j < size; j++) {
//             cout << endpoints[j] << ", ";
//         }
//         cout << endl;
//     }
//
//     // temp2 used for the tree decision value append
//     uint64_t *temp2 = new uint64_t[depth]();
//     int t = 0;
//     size = num_stand;
//     while (true) {
//         if (size <= 3) {
//             break;
//         }
//         temp2[t] = size / 2;
//         size = (size + 1) / 2;
//         if (size % 2 == 0) {
//             size += 1;
//         }
//         t += 1;
//     }
//
//     if (party == ALICE) {
//         for (int i1 = 0; i1 < depth; i1++) {
//             res_tree[i1] ^= 1;
//         }
//     }
//
//
//     cout << endl << "append value: " << endl;
//     for (int i1 = 0; i1 < depth; i1++) {
//         cout << temp2[i1] << ", ";
//     }
//     cout << endl;
//
//     cout << endl << "res_tree: " << endl;
//     for (int i1 = 0; i1 < depth; i1++) {
//         cout << (int) res_tree[i1] << ", ";
//     }
//     cout << endl;
//
//     aux->multiplexer_two_plain(res_tree, temp2, temp2, depth, bw_res, bw_res);
//     for (int i = 0; i < depth; i++) {
//         res[0] += temp2[i];
//     }
//
//     if (party == ALICE) {
//         res[0] -= 1;
//     }
//
//     res[0] &= (1ULL << bw_res) - 1;
//
//
//     cout << endl;
//     cout << "res[0]: " << res[0] << endl;
//
//
//     delete [] endpoints;
//     delete [] endpoints_temp;
//     delete [] endpoints_temp2;
//     delete[] temp_t;
//     delete [] res_tree;
//     delete[] temp2;
// }
