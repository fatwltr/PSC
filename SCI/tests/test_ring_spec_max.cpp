#include <iostream>
#include <string>

#include "globals.h"
#include "utils/emp-tool.h"
#include "BuildingBlocks/aux-protocols.h"
#include "Millionaire/equality.h"

using namespace std;
using namespace sci;

int party = 0, port = 8000;
string address = "127.0.0.1";
int dim = 3;
int bw = 4;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of elements");
    amap.arg("bw", bw, "Bitwidth of input");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    IOPack *io = new IOPack(party, port, address);
    sci::OTPack *otpack = new sci::OTPack(io, party);

    auto aux = new AuxProtocols(party, iopack, otpack);
    auto max = new MaxPoolProtocol<uint64_t>(party, RING, io, bw, 4, 0, otpack);

    uint64_t *data = new uint64_t[dim];
    uint64_t res;
    uint64_t resIdx;
    uint64_t bw_idx = 4;

    PRG128 prg;
    prg.random_data(data, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; ++i) {
        data[i] &= (1ULL << (bw - 1)) - 1;
        cout << data[i] << " ";
    }
    cout << endl;

    max->funcMaxMPC_spcial(dim, data, &res, &resIdx, bw_idx);

    cout << endl;
    cout << "max value: " << res << endl;
    cout << "max idx: " << resIdx << endl;

    delete[] data;
    delete max;
    delete aux;
    delete otpack;
    delete io;
    return 0;
}
