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
int eq = 0;
string address = "127.0.0.1";
int num_data = 8;
int num_stand = 4;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("eq", eq, "eq=1 use count_eq");
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
    for (int i = 0; i < num_stand; ++i) {
        stand[i] = i;
    }
    uint64_t *res = new uint64_t[num_stand]();
    int bw_data = ceil(log2(num_stand + 1));
    int bw_res = ceil(log2(num_data + 1));

    PRG128 prg;
    prg.random_data(data, num_data * sizeof(uint64_t));
    for (int i = 0; i < num_data; ++i) {
        data[i] %= num_stand;
    }

    auto start = clock_start();
    uint64_t comm = iopack->get_comm();
    if (eq == 1) {
        frequency->count_eq_batch(res, data, num_data, stand, num_stand, bw_data, bw_res);
    } else {
        // frequency->count_shift(res, data, num_data, stand, num_stand, bw_data, bw_res);
        frequency->count_shift_batch(res, data, num_data, stand, num_stand, bw_data, bw_res);
    }

    comm = iopack->get_comm() - comm;
    long long t = time_from(start);

    cout << "frequency Time\t" << t / (1000.0) << " ms" << endl;
    cout << "frequency Bytes Sent\t" << comm << " bytes" << endl;

    //
    // if (party == ALICE) {
    //     iopack->io->send_data(data, num_data * sizeof(uint64_t));
    //     iopack->io->send_data(res, num_stand * sizeof(uint64_t));
    // } else {
    //     uint64_t *t_res = new uint64_t[num_stand];
    //     uint64_t *t_data = new uint64_t[num_data];
    //     iopack->io->recv_data(t_data, num_data * sizeof(uint64_t));
    //     iopack->io->recv_data(t_res, num_stand * sizeof(uint64_t));
    //     for (int i = 0; i < num_data; ++i) {
    //         cout << t_data[i] << " ";
    //     }
    //     cout << endl;
    //     for (int i = 0; i < num_data; ++i) {
    //         cout << data[i] << " ";
    //     }
    //     cout << endl;
    //     for (int i = 0; i < num_stand; ++i) {
    //         t_res[i] += res[i];
    //         t_res[i] &= (1ULL << bw_res) - 1;
    //     }
    //     for (int i = 0; i < num_data; ++i) {
    //         t_data[i] += data[i];
    //         t_data[i] %= num_stand;
    //     }
    //     for (int i = 0; i < num_data; ++i) {
    //         cout << t_data[i] << " ";
    //     }
    //     cout << endl;
    //     for (int i = 0; i < num_stand; ++i) {
    //         cout << t_res[i] << " ";
    //     }
    //     cout << endl;
    //
    //     delete [] t_data;
    //     delete [] t_res;
    // }

    uint64_t* sorted_array = new uint64_t[num_data]();
    auto sort_start = clock_start();
    uint64_t comm_second = iopack->get_comm();

    frequency->count_sort(sorted_array, res, num_stand, num_data, bw_data, bw_res);

    comm_second = iopack->get_comm() - comm_second;
    long long sort_time = time_from(sort_start);

    cout << endl;
    cout << "Ex-sort Time\t" << sort_time / (1000.0) << " ms" << endl;
    cout << "Ex-sort Bytes Sent\t" << comm_second << " bytes" << endl;

    delete[] data;
    delete[] res;
    delete[] stand;
    delete otpack;
    delete iopack;
    delete frequency;
    return 0;
}
