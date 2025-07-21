#!/bin/bash
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00

# fix this according to your installation path
QE_PATH=~/soft/qe-7.1/build/bin

mpirun ${QE_PATH}/pw.x < scf.in > scf.out
mpirun ${QE_PATH}/ph.x -pd .true. < ph.in > ph.out
python pp.py
mpirun ${QE_PATH}/pw.x < nscf.in > nscf.out
mpirun ${QE_PATH}/epw.x -input epw.in > epw.out


