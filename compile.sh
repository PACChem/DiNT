#!/bin/bash

#rm -r build/*
cd build/
cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../ ../
make VERBOSE=1
make install

#cd src/
#make MakefileOLD
