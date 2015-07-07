#!/bin/bash
# Intel Compiler Version

# Linux Intel Compiler
#icpc --std=c++11 -qopenmp -O3 -ipo -no-prec-div -fp-model fast=2 -D USE_MPI -I /opt/intel/impi/5.0.2.044/include64 -L /opt/intel/impi/5.0.2.044/lib64 -lpthread -lmpi -lrt r2ks.cpp -o r2ks

# MacOSX Clang
clang-omp++ -stdlib=libc++ --std=c++11 -fopenmp -O3 -D USE_MPI -lmpi -lpthread r2ks.cpp -o r2ks
