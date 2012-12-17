#!/bin/bash

gcc -L/home/dcampora/temp/tbb41_20121003oss/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21 -ltbb -lrt -pthread -O2 cpu-parallel.cpp HelperSeq.cpp SeqAlg.cpp float_SeqAlg.cpp
