# Adaptive Mesh Refinement Solvers


## Dependencies  
cmake >= 2.6  
lua  
autoconf >= 1.11  
mpich >= 1.15  
  

Need to build p4est and sc libraries and link with project.  

$ export DEBUG="-O0 -g -Wall -Wuninitialized"  
$ export FAST="-O2 -w"  

$ ./bootstrap  

Build debug:  
$ ./configure CFLAGS=$DEBUG CC=mpicc CXX=mpicxx --enable-mpi --enable-shared --without-blas  
$ make && make install  


## Build project

Determine p4est and sc libraries path  

$ export P4EST_LIBRARY=/usr/local/lib  
$ export P4EST_INCLUDE=/usr/local/include  
  
$ export SC_LIBRARY=/usr/local/lib  
$ export SC_INCLUDE=/usr/local/include  


Build target  

$ mkdir build && cd build  
$ cmake .. -DP4EST_LIBRARY=$P4EST_LIBRARY -DP4EST_INCLUDE_DIR=$P4EST_INCLUDE -DSC_LIBRARY=$P4EST_LIBRARY -DSC_INCLUDE_DIR=$P4EST_INCLUDE  
$ make amr  
