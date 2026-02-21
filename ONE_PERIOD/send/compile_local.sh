#!/bin/bash

rm -f ANNT.exe *.o *.mod

#gfortran-13 -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow \
#  -fopenmp -ffree-line-length-512 \
#  toolbox_CE.f90 mod_parameters.f90 Toolbox.f90 ModuleSAVE.f90 FUNCTIONS.f90 ModuleSIMULATION.f90 MAIN.f90 \
#  -llapack -lblas \
#  -o ANNT.exe
#./ANNT.exe

#gfortran-13 -O0 -g \
#  -fcheck=all \
#  -fbacktrace \
#  -ffpe-trap=invalid,zero,overflow \
#  -Wall -Wextra -Wno-tabs \
#  -fopenmp \
#  -ffree-line-length-512 \
#  toolbox_CE.f90 mod_parameters.f90 Toolbox.f90 ModuleSAVE.f90 FUNCTIONS.f90 ModuleSIMULATION.f90 MAIN.f90 \
#  -llapack -lblas \
#  -o ANNT.exe

gfortran -o ANNT.exe toolbox_CE.f90 mod_parameters.f90 Toolbox.f90 ModuleSAVE.f90 FUNCTIONS.f90  ModuleSIMULATION.f90 MAIN.f90 \
-O3 -march=native -funroll-loops -fopenmp -ffree-line-length-512  \
  -llapack -lblas

./ANNT.exe  

