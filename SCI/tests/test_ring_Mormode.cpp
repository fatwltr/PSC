//
// Created by a1141 on 25-7-15.
//
#include <iostream>
#include <string>

#include "globals.h"
#include "utils/emp-tool.h"
#include "BuildingBlocks/aux-protocols.h"
#include "Millionaire/equality.h"
#include "Frequency/frequency.h"

using namespace std;
using namespace sci;

int party = 0, port = 8000;
string address = "127.0.0.1";
int num_data = 16;
int num_stand = 8;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("nd", num_data, "Number of elements");
    amap.arg("ns", num_stand, "Number of stands");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    IOPack *iopack = new IOPack(party, port, address);
    sci::OTPack *otpack = new sci::OTPack(iopack, party);
    Frequency *frequency = new Frequency(party, iopack, otpack);
    uint64_t *data = new uint64_t[num_data];
    uint64_t *stand = new uint64_t[num_stand];
    uint64_t res = 0;
    for (int i = 0; i < num_stand; ++i) {
        stand[i] = i;
    }
    int bw_data = ceil(log2(num_stand + 1));
    int bw_res = ceil(log2(num_data + 1));

    PRG128 prg;
    prg.random_data(data, num_data * sizeof(uint64_t));
    for (int i = 0; i < num_data; ++i) {
        if (i < num_data / 2 + 2) {
            data[i] = 1;
        }
        data[i] &= (1ULL << bw_data) - 1;
    }
    cout << endl;

    auto start = clock_start();


    uint64_t comm = iopack->get_comm();

    // frequency->mode_Mor(&res, data, num_data, stand, num_stand, bw_data, bw_res);
    frequency->mode_Mor_batch(&res, data, num_data, stand, num_stand, bw_data, bw_res);



    comm = iopack->get_comm() - comm;
    long long t = time_from(start);

    cout << "Time\t" << t / (1000.0) << " ms" << endl;
    cout << "Bytes Sent\t" << comm << " bytes" << endl;

    // if (party == ALICE) {
    //     iopack->io->send_data(data, num_data * sizeof(uint64_t));
    //     iopack->io->send_data(&res, 1 * sizeof(uint64_t));
    // } else {
    //     uint64_t t_res = 0;
    //     uint64_t *t_data = new uint64_t[num_data];
    //     iopack->io->recv_data(t_data, num_data * sizeof(uint64_t));
    //     iopack->io->recv_data(&t_res, 1 * sizeof(uint64_t));
    //
    //     for (int i = 0; i < num_data; ++i) {
    //         data[i] += t_data[i];
    //         data[i] &= (1ULL << bw_data) - 1;
    //         cout << data[i] << ", ";
    //     }
    //     cout << endl;
    //     res += t_res;
    //
    //     cout << "the mode is: " << res << endl;
    //     delete [] t_data;
    // }

    delete[] data;
    delete[] stand;
    delete otpack;
    delete iopack;
    delete frequency;
    return 0;
}
