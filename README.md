# Adaptive Mesh Refinement Solvers


# Dependencies

* p4est and sc:

$ ./bootstrap


Building debug version:

$ ./configure CFLAGS="-O0 -g -Wall -Wuninitialized" --enable-debug --enable-mpi --disable-shared

$ make && make install


# Install

First, need to build p4est and sc libraries and link with project.

$ mkdir build && cd build
 
$ cmake .. -DP4EST_LIBRARY=../../p4est/lib -DP4EST_INCLUDE_DIR=../../p4est/include -DSC_LIBRARY=../../p4est/lib -DSC_INCLUDE_DIR=../../p4est/include

$ make && make install
