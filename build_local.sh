#!/bin/bash
icpc -g --std=c++11 -qopenmp -O3 -ipo -no-prec-div -fp-model fast=2 -I /opt/intel/impi/5.0.2.044/include64 -L /opt/intel/impi/5.0.2.044/lib64 -lpthread -lmpi -lrt r2ks.cpp -o r2ks
