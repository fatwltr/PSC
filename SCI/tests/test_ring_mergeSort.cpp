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
int num_data = 8;
int num_stand = 4;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("nd", num_data, "Number of elements");
    amap.arg("ns", num_stand, "Number of stands, or the k in shuffleTopk");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    IOPack *iopack = new IOPack(party, port, address);
    sci::OTPack *otpack = new sci::OTPack(iopack, party);
    Frequency *frequency = new Frequency(party, iopack, otpack);
    uint64_t *data = new uint64_t[num_data];
    uint64_t *stand = new uint64_t[num_stand];
    for (int i = 0; i < num_stand; ++i) {
        stand[i] = i;
    }
    uint64_t *res = new uint64_t[num_data]();
    int bw_data = ceil(log2(num_stand + 1));
    int bw_res = ceil(log2(num_data + 1));

    PRG128 prg;
    prg.random_data(data, num_data * sizeof(uint64_t));
    for (int i = 0; i < num_data; ++i) {
        data[i] &= (1ULL << bw_data) - 1;
    }

    auto start = clock_start();
    uint64_t comm = iopack->get_comm();

    frequency->batcher_network_sort(data, num_stand, num_data, bw_data, bw_res);
    // frequency->shuffle_topk(res, data, num_stand, num_data, bw_data, bw_res); // the res here should be num_data as the code inner shuffle_topk reuse the res.

    comm = iopack->get_comm() - comm;
    long long t = time_from(start);

    cout << "sort Time\t" << t / (1000.0) << " ms" << endl;
    cout << "sort Bytes Sent\t" << comm << " bytes" << endl;


    delete[] data;
    delete[] res;
    delete[] stand;
    delete otpack;
    delete iopack;
    delete frequency;
    return 0;
}
