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

void print_block_(const block128 &b) {
    uint8_t bytes[16];
    _mm_storeu_si128((__m128i *) bytes, b); // copy to byte array
    for (int i = 0; i < 16; ++i)
        std::cout << std::hex << std::setw(2) << std::setfill('0') << (int) bytes[i];
    std::cout << std::endl;
}

void Frequency::count_sort(uint64_t *res, uint64_t *frequency, int num_stand, int num_data,
                           int32_t bw_data, int32_t bw_res) {
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

    // share conversion 2^bw_data -> num_data
    uint8_t *t = new uint8_t[num_stand]();
    uint64_t *appendix = new uint64_t[num_stand]();
    this->aux->wrap_computation(endpoints, t, num_stand, bw_data);
    this->aux->multiplexer_two_plain(t, (1ULL << bw_data) % num_data, appendix, num_stand, bw_data, bw_data);
    for (int i = 0; i < num_stand; i++) {
        endpoints[i] -= appendix[i];
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
    auto *prg = new PRG128[num_stand * num_data * 2]();
    auto **x = new uint8_t *[num_stand];
    for (int i = 0; i < num_stand; i++) {
        x[i] = new uint8_t[num_data * 2]();
    }
    auto **uniShr = new uint8_t *[num_stand];
    for (int i = 0; i < num_stand; i++) {
        uniShr[i] = new uint8_t[num_data * 2]();
    }
    auto **interval_indicate = new uint8_t *[num_stand];
    for (int i = 0; i < num_stand; i++) {
        interval_indicate[i] = new uint8_t[num_data]();
    }
    int save_memory = 8;
    if (party == sci::ALICE) {
        // set the offset directly in the index
        for (int i = 0; i < num_stand; i++) {
            std::fill(x[i], x[i] + endpoints[i], 1);
            std::fill(x[i] + num_data + endpoints[i], x[i] + 2 * num_data, 1);
        }

        // do not need to generate a 2*num_data seed, we can set the second half seed directly as offset < num_data
        this->aux->nMinus1OUTNOT_batch(seeds, num_stand, 2 * num_data, nullptr);
        // generate the shift translation shares


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
        auto **a = new uint8_t *[num_stand];
        for (int i = 0; i < num_stand; i++) {
            a[i] = new uint8_t[num_data * 2]();
        }
        auto **b = new uint8_t *[num_stand];
        for (int i = 0; i < num_stand; i++) {
            b[i] = new uint8_t[num_data * 2]();
        }

        for (int i = 0; i < num_stand; i++) {
            for (int j = 0; j < num_data * 2; j++) {
                prg[i * num_data * 2 + j].reseed(&seeds[i * num_data * 2 + j]);
            }
        }

        for (int i = 0; i < num_stand; i++) {
            for (int k = 0; k < (num_data * 2 / save_memory); k++) {
                for (int j = 0; j < num_data * 2; j++) {
                    prg[i * num_data * 2 + j].random_bool((bool *) shift_translate[j], save_memory);
                }

                cout << "shift_translate" << endl;
                for (int j = 0; j < num_data * 2; j++) {
                    for (int l = 0; l < save_memory; l++) {
                        cout << (int) shift_translate[j][l] << " ";
                    }
                    cout << endl;
                }

                for (int j = 0; j < num_data * 2; j++) {
                    for (int l = 0; l < save_memory; l++) {
                        a[i][l + k * save_memory] ^= shift_translate[j][l];
                        b[i][j] ^= shift_translate[(j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
                    }
                }
            }
            int last = num_data * 2 % save_memory;
            if (last == 0) {
                continue;
            }
            int k = (num_data * 2 / save_memory);
            for (int j = 0; j < num_data * 2; j++) {
                prg[i * num_data * 2 + j].random_bool((bool *) shift_translate[j], last);
            }
            cout << "shift_translate" << endl;
            for (int j = 0; j < num_data * 2; j++) {
                for (int l = 0; l < last; l++) {
                    cout << (int) shift_translate[j][l] << " ";
                }
                cout << endl;
            }
            for (int j = 0; j < num_data * 2; j++) {
                for (int l = 0; l < last; l++) {
                    a[i][l + k * save_memory] ^= shift_translate[j][l];
                    b[i][j] ^= shift_translate[(j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
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


        for (int i = 0; i < num_stand; i++) {
            for (int j = 0; j < num_data * 2; j++) {
                x[i][j] = x[i][j] ^ a[i][j];
            }
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
        pack_bits(x, x_packed, num_stand, 2 * num_data);
        iopack->io->send_data(x_packed, ((2 * num_data * num_stand + 7) / 8) * sizeof(uint8_t));


        // allocate the uniShr, and the
        auto **mux_choice = new uint8_t *[num_stand];
        for (int i = 0; i < num_stand; i++) {
            mux_choice[i] = new uint8_t[num_data * 2]();
            for (int j = 0; j < num_data; j++) {
                mux_choice[i][j] = b[i][j] ^ b[i][j + num_data];
            }
        }
        uint64_t *temp = new uint64_t[num_stand];
        for (int i = 0; i < num_stand; i++) {
            temp[i] = num_data - endpoints[i];
        }
        this->aux->mill->compare(t, temp, num_stand, bw_data);
        // mux_choice padding a 0-string
        this->aux->multiplexer_bShr(t, mux_choice, uniShr, num_stand, num_data * 2);
        for (int i = 0; i < num_stand; i++) {
            for (int j = 0; j < num_data; j++) {
                uniShr[i][j] ^= b[i][j + num_data];
                uniShr[i][j] ^= uniShr[i][j + num_data];
            }
        }
        for (int i = 0; i < num_data * 2; i++) {
            delete[] shift_translate[i];
        }
        delete[] shift_translate;
        for (int i = 0; i < num_stand; i++) {
            delete[] a[i];
            delete[] b[i];
            delete[] mux_choice[i];
        }
        delete[] x_packed;
        delete[] mux_choice;
        delete[] temp;
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
        auto **c = new uint8_t *[num_stand];
        for (int i = 0; i < num_stand; i++) {
            c[i] = new uint8_t[num_data * 2]();
        }

        for (int i = 0; i < num_stand; i++) {
            for (int j = 0; j < num_data * 2; j++) {
                prg[i * num_data * 2 + j].reseed(&seeds[i * num_data * 2 + j]);
            }
        }

        for (int i = 0; i < num_stand; i++) {
            for (int k = 0; k < (num_data * 2 / save_memory); k++) {
                for (int j = 0; j < num_data * 2; j++) {
                    prg[i * num_data * 2 + j].random_bool((bool *) shift_translate[j], save_memory);
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
                        c[i][j] ^= shift_translate[(j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
                        c[i][(l + k * save_memory + endpoints[i] + 2 * num_data) % (2 * num_data)] ^= shift_translate[
                            j][l];
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
                prg[i * num_data * 2 + j].random_bool((bool *) shift_translate[j], last);
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
                    c[i][j] ^= shift_translate[(j - l - k * save_memory + num_data * 2) % (num_data * 2)][l];
                    c[i][(l + k * save_memory + endpoints[i] + 2 * num_data) % (2 * num_data)] ^= shift_translate[j][l];
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

        auto **xMa = new uint8_t *[num_stand];
        for (int i = 0; i < num_stand; i++) {
            xMa[i] = new uint8_t[num_data * 2]();
        }

        auto **shift_xMa = new uint8_t *[num_stand];
        for (int i = 0; i < num_stand; i++) {
            shift_xMa[i] = new uint8_t[num_data * 2]();
        }
        auto *xMa_packed = new uint8_t[(2 * num_data * num_stand + 7) / 8];
        iopack->io->recv_data(xMa_packed, ((2 * num_data * num_stand + 7) / 8) * sizeof(uint8_t));
        unpack_bits(xMa_packed, xMa, num_stand, 2 * num_data);
        delete[] xMa_packed;

        // cout << endl;
        // cout << "------------xMa-------------" << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         cout << (int) xMa[i][j] << ", ";
        //     }
        //     cout << endl;
        // }


        for (int i = 0; i < num_stand; i++) {
            std::memcpy(shift_xMa[i], xMa[i] + (2 * num_data - endpoints[i]), endpoints[i]);
            std::memcpy(shift_xMa[i] + endpoints[i], xMa[i], 2 * num_data - endpoints[i]);
        }

        // cout << endl;
        // cout << "------------shift_xMa-------------" << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         cout << (int) shift_xMa[i][j] << ", ";
        //     }
        //     cout << endl;
        // }

        for (int i = 0; i < num_stand; i++) {
            for (int j = 0; j < num_data * 2; j++) {
                x[i][j] = c[i][j] ^ shift_xMa[i][j];
            }
        }


        // cout << endl;
        // cout << "---------c xor shift_xMa-----------" << endl;
        // for (int i = 0; i < num_stand; i++) {
        //     for (int j = 0; j < num_data * 2; j++) {
        //         cout << (int) x[i][j] << ", ";
        //     }
        //     cout << endl;
        // }


        // seems choose the wrong side, and need some modification to solve the situation
        uint8_t *first_half = new uint8_t[num_stand]();
        auto **mux_choice = new uint8_t *[num_stand];
        for (int i = 0; i < num_stand; i++) {
            mux_choice[i] = new uint8_t[num_data * 2](); // *2 due to the shifted vector is shared
            for (int j = num_data; j < 2 * num_data; j++) {
                mux_choice[i][j] = x[i][j - num_data] ^ x[i][j];
            }
        }
        this->aux->mill->compare(t, endpoints, num_stand, bw_data);
        // mux_choice padding a 0-string
        this->aux->multiplexer_bShr(first_half, mux_choice, uniShr, num_stand, num_data * 2);
        for (int i = 0; i < num_stand; i++) {
            for (int j = num_data; j < 2 * num_data; j++) {
                uniShr[i][j] ^= x[i][j];
                uniShr[i][j - num_data] ^= uniShr[i][j];
            }
        }
        for (int i = 0; i < 2 * num_data; i++) {
            delete[] shift_translate[i];
        }
        delete[] shift_translate;
        for (int i = 0; i < num_stand; i++) {
            delete[] c[i];
            delete[] mux_choice[i];
            delete[] xMa[i];
            delete[] shift_xMa[i];
        }
        delete[] c;
        delete[] xMa;
        delete[] shift_xMa;
        delete[] mux_choice;
        delete[] first_half;
    }

    // for (int i = 0; i < num_stand; i++) {
    //     cout << endl;
    //     for (int j = 0; j < num_data; j++) {
    //         cout << (int) uniShr[i][j] << ", ";
    //     }
    // }
    // cout << endl;

    // xor unishr
    if (party == ALICE) {
        for (int j = 0; j < num_data; j++) {
            interval_indicate[0][j] = 1 ^ uniShr[0][j];
        }
    }
    for (int i = 1; i < num_stand; i++) {
        for (int j = 0; j < num_data; j++) {
            interval_indicate[i][j] = uniShr[i][j] ^ uniShr[i - 1][j];
        }
    }
    uint64_t *temp = new uint64_t[num_data]();
    // MUX, can batch
    for (int i = 0; i < num_stand; i++) {
        aux->multiplexer_two_plain(interval_indicate[i], i, temp, num_data, bw_data, bw_data);
        for (int j = 0; j < num_data; j++) {
            res[j] += temp[j];
        }
    }

    delete[] endpoints;
    delete[] t;
    delete[] appendix;
    delete[] seeds;
    delete[] prg;
    for (int i = 0; i < num_stand; i++) {
        delete[] x[i];
        delete[] uniShr[i];
        delete[] interval_indicate[i];
    }
    delete[] x;
    delete[] uniShr;
    delete[] interval_indicate;
}

void random_permutation(uint32_t *perm, int n) {
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

void Frequency::shuffle_sort(uint64_t *res, uint64_t *data, int num_stand, int num_data,
                             int32_t bw_data, int32_t bw_res) {
    // oblivious shuffle
    // 1. element extend
    // 2. opv
    // 3. shuffle
    uint32_t *perm = new uint32_t[num_data]();
    random_permutation(perm, num_data);
    for (int i = 0; i < num_data; i++) {
        cout << perm[i] << " ";
    }
    cout << endl;

    // quick sort
}
