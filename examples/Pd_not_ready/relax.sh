#!/bin/bash
#SBATCH --ntasks-per-node=32
##SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --partition=long

ulimit -s unlimited
ulimit -v unlimited
#ulimit -a

# fix this according to your installation path
QE_PATH=~/soft/qe-7.1/build/bin

mpirun ${QE_PATH}/pw.x < relax.in > relax.out
