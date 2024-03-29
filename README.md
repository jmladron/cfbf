# CFBF
Cuckoo Filter with an Integrated Bloom Filter (CFBF) 

## Description

The goal of CFBF is to integrate a Bloom filter in the Cuckoo Filter to improve its performance during burst of hundreds of insertions. This code is used in the paper "CFBF: Reducing the Insertion Time of Cuckoo Filters with an Integrated Bloom Filter" by P. Reviriego, J. Martínez and S. Pontarelli

## Compilation

This Visual Studio C++ project consists of 3 files: CF.hpp (declaration file), CF.cpp (source file) and main.cpp (test bench)

## Command line arguments

Command line arguments for CFBF are:

- m: test mode, 1:CFBF filter test with insertion burst, 2:CFBF filter FPR having p keys in the BF, default value is 1
- s: filter size, default value is 8192
- c: filter cells in bucket, default value is 4
- o: filter occupancy, default value is 95%
- t: bf threshold, default value is 1000
- b: insertion burst, default value is 1% of filter size x filter cells in bucket
- p: keys in BF, default value is p = {40, 81, 122, 163, 204, 245, 286, 327, 368, 409}
- r: runs (trials), default value is 1000
- f: fingerprint_bits, default value is f = {12}

#### example: cfbf.exe o=94 t=10 b=491 p=40 r=5000 f=12

## Authors

The Cuckoo Filter was developed by S. Pontarelli and the CFBF was developed by J. Martinez

## License

MIT
