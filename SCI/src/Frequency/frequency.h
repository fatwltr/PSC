//
// Created by a1141 on 25-7-15.
//

#ifndef FREQUENCT_H
#define FREQUENCT_H

#include "BuildingBlocks/aux-protocols.h"
#include "Millionaire/millionaire.h"
#include "BuildingBlocks/truncation.h"
#include "BuildingBlocks/value-extension.h"
#include "LinearOT/linear-ot.h"
#include "OT/emp-ot.h"
#include "GC/emp-sh2pc.h"
#include "NonLinear/maxpool.h"

class Frequency {
public:
    int party;
    sci::IOPack *iopack;
    sci::OTPack *otpack;
    AuxProtocols *aux;
    XTProtocol *xt;
    Truncation *trunc;
    LinearOT *mult;
    Equality *equality;
    MaxPoolProtocol<uint64_t> *maxpool_oracle;

    Frequency(int party, sci::IOPack *iopack, sci::OTPack *otpack);

    ~Frequency();

    void count_eq(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                  int32_t bw_data, int32_t bw_res);

    void count_eq_batch(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                        int32_t bw_data, int32_t bw_res);

    void count_shift(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand, int32_t bw_data,
                     int32_t bw_res);

    void count_shift_batch(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand,
                           int32_t bw_data, int32_t bw_res);

    void count_eq_inner_self(uint64_t *res, uint64_t *data, int num_data, int num_stand, int32_t bw_data,
                             int32_t bw_res);

    void mode_naive(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand, int32_t bw_data,
                    int32_t bw_res, uint8_t eq);

    void mode_CRT_eq(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand, int32_t bw_data,
                     int32_t bw_res);

    void mode_CRT_shift(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand, int32_t bw_data,
                        int32_t bw_res);

    void count_sort(uint64_t *res, uint64_t *frequency, int num_stand, int num_data,
                           int32_t bw_data, int32_t bw_res);
    void shuffle_sort(uint64_t *res, uint64_t *data, int num_stand, int num_data,
                           int32_t bw_data, int32_t bw_res);
    int partition(uint64_t *arr, int low, int high, int bw);
    void quickSort(uint64_t *arr, int low, int high, int bw);

    void oblivious_shuffle(uint64_t *res, uint64_t *perm, uint64_t* data, int num_data, int32_t bw_data);
    void oblivious_shuffle_reverse(uint64_t *res, uint64_t *perm, uint64_t* data, int num_data, int32_t bw_data);

    void pack_bits(uint8_t **x, uint8_t* x_packed, int rows, int cols);

    void unpack_bits(uint8_t* x_packed, uint8_t **x, int rows, int cols) ;

    void pack_data(uint64_t *x, uint64_t *x_packed, int length, int bw);

    void unpack_data(uint64_t *x_packed, uint64_t *x, int length, int bw);
};

#endif
