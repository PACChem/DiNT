#!/bin/bash

#TINKLIB="/home/ajasper/KTP/dint/tinker/tinker/source"
cd build
cmake -DCMAKE_Fortran_FLAGS=" -O3 -g -C" -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../exe ../
make
#make install
