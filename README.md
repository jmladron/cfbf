# cfbf
Cuckoo Filter with an Integrated Bloom Filter (CFBF)

## Description

The goal of CFBF is to integrate a Bloom filter in the Cuckoo Filter to improve its performance during burst of hundreds of insertions

## Compilation

This Visual Studio C++ project consists of 3 files: CF.hpp (declaration file), CF.cpp (source file) and main.cpp (testbench)

## Usage

Command line arguments for CFBF are:

s: filter size, default value is 8192

c: filter cells in bucket, default value is 4

o: filter occupancy, default value is 95%

t: bf threshold, default value is 1000

b: insertion burst, default value is 1% of filter size x filter cells in bucket

r: runs (trials), default value is 1000

f: fingerprint_bits, default value is f = {12}

example: cfbf.exe o=94 t=10 b=491 r=5000 f=12

## Authors

The Cuckoo Filter was developed by S. Pontarelli (salvatore.pontarelli@uniroma2.it) and the enhanced CFBF was developed by J. Martinez (jorge.martinez@ufv.es)

## License

?
