# Adaptive Mesh Refinement Solvers


## Dependencies

p4est and sc:  
$ ./bootstrap

Build debug:  
$ ./configure CFLAGS="-O0 -g -Wall -Wuninitialized" --enable-debug --enable-mpi --enable-shared  
$ make && make install


## Install

First, need to build p4est and sc libraries and link with project.

$ mkdir build && cd build  
$ cmake .. -DP4EST_LIBRARY=PATH_TO_LIB_DIR -DP4EST_INCLUDE_DIR=PATH_TO_INCLUDE_DIR -DSC_LIBRARY=PATH_TO_LIB_DIR -DSC_INCLUDE_DIR=PATH_TO_INCLUDE_DIR  
$ make && make install
