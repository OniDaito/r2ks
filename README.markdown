# R2kS - A  Novel Measure for Comparing Gene Expression Based on Ranked Gene Lists

This code is an MPI based version of the algorithm described in [this paper by Ni and Vingron](http://pubman.mpdl.mpg.de/pubman/item/escidoc:1702988/component/escidoc:1711821/Ni.pdf).

Using dynamic programming to compute the Kolmogorov-Smirnov statistic, we can quickly see how similar one list is to another, in both forward and reverse modes. From the paper:

Bioinformatics analyses frequently yield results in the form of lists of genes sorted by, for example, sequence similarity to a query sequence or degree of differential expression of a gene upon a change of cellular condition. Comparison of such results may depend strongly on the particular scoring system throughout the entire list, although the crucial information resides in which genes are ranked at the top of the list. Here, we propose to reduce the lists to the mere ranking of the genes and to compare only the ranked lists. To this end, we introduce a measure of similarity between ranked lists. Our measure puts particular emphasis on finding the same items near the top of the list, while the genes further down should not have a strong influence. Our approach can be understood as a special version of a two-dimensional Kolmogorov-Smirnov statistic. We present a dynamic pro- gramming algorithm for its computation and study the distribution of the similarity values. The performance on simulated and on real biological data is studied in comparison to other available measures.


## Building

To build, you need either OpenMPI, MPICH or Intel MPI installed and running. At some point I'll release a the fork for local OpenMP versions. 

This code has been tested on GCC and ICC on Ubuntu linux.

Execute the script **build_local.sh** or run  something like the following:

    icpc --std=c++11 -qopenmp -O3 -ipo -no-prec-div -fp-model fast=2 -D USE_MPI -I /opt/intel/impi/5.0.2.044/include64 \
    -L /opt/intel/impi/5.0.2.044/lib64 -lpthread -lmpi -lrt r2ks.cpp -o r2ks

There is an option -D USE_MPI which turns on MPI. Remove this switch to use OpenMP locally.

Look inside the build_local.sh for build options on Linux and OSX. You'll need clang-omp++ for building on OSX.

## Running

You should have an MPI setup installed and running

    mpirun -n 4 ./r2ks -f <filename>

Or running locally with OpenMP

    OMP_PLACES=cores OMP_PROC_BIND=spread OMP_NUM_THREADS=4 ./r2ks -f test100.t


Options include **-t** for a two tailed test and  **-w <pivot>** to use the weight with the given pivot point (index into the list)

## File Format

Files should be of the format

    <length N of all lists> <num lists>
    <list 0 of length N>
    <list 1 of length N>
    ...
    ...
    <list n of length N>

## TODO

* More options and better UI support
* re-introduce the OpenMP option for local machines and different types of clusters
* Experiment with a CUDA version

## Author / Contact

Benjamin Blundell - b.blundell@qmul.ac.uk
