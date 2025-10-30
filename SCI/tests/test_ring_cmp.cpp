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
int dim = 1 << 16;
int bw = 32;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of elements");
    amap.arg("bw", bw, "Bitwidth of input");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    IOPack *iopack = new IOPack(party, port, address);
    sci::OTPack * otpack = new sci::OTPack(iopack, party);
    MillionaireProtocol *mill = new MillionaireProtocol(party, iopack, otpack);

    uint64_t *data = new uint64_t[dim];
    uint8_t *res_eq = new uint8_t[dim];

    PRG128 prg;
    prg.random_data(data, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; ++i) {
        data[i] &= (1ULL << bw) - 1;
    }
    auto start = clock_start();
    uint64_t comm = iopack->get_comm();
    mill->compare(res_eq, data, dim, bw);

    comm = iopack->get_comm() - comm;
    long long t = time_from(start);

    cout << "Time\t" << t / (1000.0) << " ms" << endl;
    cout << "Bytes Sent\t" << comm << " bytes" << endl;
    // for (int i = 0; i < std::min(dim, 10); ++i) {
    //     cout << "data[" << i << "]=" << data[i] << ", eq[" << i << "]=" << int(res_eq[i]) << endl;
    // }

    delete[] data;
    delete[] res_eq;
    delete mill;
    delete otpack;
    delete io;
    return 0;
}
