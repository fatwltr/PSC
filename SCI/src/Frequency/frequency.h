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

    void count_shift(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand, int32_t bw_data,
                     int32_t bw_res);

    void count_eq_inner_self(uint64_t *res, uint64_t *data, int num_data, int num_stand, int32_t bw_data,
                             int32_t bw_res);

    void mode_naive(uint64_t *res, uint64_t *data, int num_data, uint64_t *stand, int num_stand, int32_t bw_data,
                    int32_t bw_res, uint8_t eq);

    void mode_CRT(uint64_t* res, uint64_t *data, int num_data, uint64_t *stand, int num_stand, int32_t bw_data,
                  int32_t bw_res, uint8_t eq);
};

#endif
