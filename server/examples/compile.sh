#!/bin/bash
#SBATCH --job-name=CpH
# normal-32g or normal-64g
#SBATCH --partition=normal-64g
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=28gb

module swap gnu12 intel/2023.2.1

rm ANNT.exe 

ifort $F90FLAGS -o ANNT.exe toolbox_CE.f90 mod_parameters.f90 sub_discretize.f90 Toolbox.f90 ModuleSAVE.f90 FUNCTIONS.f90 ModuleSIMULATION.f90 MAIN.f90 -mkl -mcmodel=large -shared-intel -qopenmp 