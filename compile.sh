#!/bin/bash

#rm -r build/*
cd build/
#cmake -DCMAKE_Fortran_FLAGS=" -O3 -g -L/home/ajasper/KTP/dint/sprng/lib -L/usr/lib64 -L/home/moberg/lib -lmlfg" -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../ ../
cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../ ../
make
make install

#cd src/
#make MakefileOLD
