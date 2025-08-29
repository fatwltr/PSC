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

#include "BuildingBlocks/aux-protocols.h"
#include "utils/emp-tool.h"
#include <iostream>
using namespace sci;
using namespace std;

int party, port = 8000, dim = 16;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
AuxProtocols *aux;

void test_sec_shift() {
    int bw_in = 32;
    PRG128 prg;
    uint64_t x = 3;
    uint8_t *ori_Bit_vector = new uint8_t[dim]();
    uint8_t *shifted = new uint8_t[dim]();

    block128 *seed = new block128[dim]();
    if (party == ALICE) {
        aux->uniShare_naive_bool(ori_Bit_vector, dim, x);
    }
    else {
        aux->uniShare_naive_bool(ori_Bit_vector, dim, x + 1);
    }
    // aux->uniShare_CRT_bool(ori_Bit_vector, dim, x);
    // aux->nMinus1OUTNOT(seed, dim, x);

    if (party == ALICE) {
        std::cout << std::endl;
        std::cout << "Alice vector" << std::endl;
        for (int i = 0; i < dim; i++) {
            std::cout << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(ori_Bit_vector[i]) << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << std::endl;
        std::cout << "Bob vector" << std::endl;
        for (int i = 0; i < dim; i++) {
            std::cout << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(ori_Bit_vector[i]) << " ";
        }
        std::cout << std::endl;
    }
    delete[] ori_Bit_vector;
    delete[] shifted;
    delete[] seed;
}


int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("d", dim, "Size of vector");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    iopack = new IOPack(party, port, "127.0.0.1");
    otpack = new OTPack(iopack, party);
    uint64_t num_rounds;

    aux = new AuxProtocols(party, iopack, otpack);

    test_sec_shift();

    return 0;
}
