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

// to count the frequency of specific number.
// the input is in Z/{num_stand}Z and the result is in Z/{2^{bw_res}}Z.
void Frequency::count_shift(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                            int32_t bw_data, int32_t bw_res) {
    uint8_t *t_shifted = new uint8_t[num_stand];
    uint64_t mask_data = (bw_data == 64 ? -1 : ((1ULL << bw_data) - 1));
    uint64_t mask_res = (bw_res == 64 ? -1 : ((1ULL << bw_res) - 1));
    uint64_t *t = new uint64_t[num_stand];
    for (int i = 0; i < num_data; i++) {
        aux->uniShare_naive_bool(t_shifted, num_stand, data[i] % num_stand);
        aux->B2A(t_shifted, t, num_stand, bw_res);
        for (int j = 0; j < num_stand; j++) {
            res[j] += t[j];
            res[j] &= mask_res;
        }
    }
    delete[] t_shifted;
    delete[] t;
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
            count_eq(count, data, num_data, stand, num_stand, bw_data, bw_res + 1);
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
        count_shift(count, data, num_data, stand, num_stand, bw_data, bw_res + 1);
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


void Frequency::mode_CRT(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                         int32_t bw_data, int32_t bw_res, uint8_t eq = 0) {
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
            count_shift(res_t[i], data_fragments[i], num_data, stand, fragment_modulus[i], bw_data_t, bw_res + 1);
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
            count_shift(res_t[i], data_fragments[i], num_data, stand, fragment_modulus[i], bw_data_t, bw_res + 1);
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
