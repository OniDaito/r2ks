# R2kS - A  Novel Measure for Comparing Gene Expression Based on Ranked Gene Lists

This code is an MPI based version of the algorithm described in [this paper by Ni and Vingron](http://pubman.mpdl.mpg.de/pubman/item/escidoc:1702988/component/escidoc:1711821/Ni.pdf).

Using dynamic programming to compute the Kolmogorov-Smirnov statistic, we can quickly see how similar one list is to another, in both forward and reverse modes. This is useful for measuring how similar diseases or medicines are in terms of their gene expression.


## Building

To build, you need either OpenMPI, MPICH or Intel MPI installed and running. At some point I'll release a the fork for local OpenMP versions. 

This code has been tested on GCC and ICC on Ubuntu linux.

Execute the script '''build_local.sh''' or run  something like the following:

    icpc --std=c++11 -qopenmp -O3 -ipo -no-prec-div -fp-model fast=2 -I /opt/intel/impi/5.0.2.044/include64 \
    -L /opt/intel/impi/5.0.2.044/lib64 -lpthread -lmpi -lrt r2ks.cpp -o r2ks

## Running

You should have an MPI setup installed and running

    mpirun -n 4 ./r2ks -f <filename>

Options include '''-t''' for a two tailed test and  '''-w <weight>''' to use the weight.

## File Format

Files should be of the format

    <num genes> <num lists>
    <list 0 of gene expressions>
    <list 1 of gene expressions>
    ...
    ...
    <list n of gene expressions>

## TODO

* More options and better UI support
* re-introduce the OpenMP option for local machines and different types of clusters
* Experiment with a CUDA version

## Author / Contact

Benjamin Blundell - b.blundell@qmul.ac.uk
