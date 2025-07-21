#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --time=7-00:00:00

# add path to your installation folder to use 'skies' program
#
export PATH=/path/to/skies-binary-folder:${PATH}

# do not forget to add paths to dynamic libs
# (an example is given)
#
export LD_LIBRARY_PATH=~/soft/tbb/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=~/local/gcc-11.2/lib64:${LD_LIBRARY_PATH}

####################################################################
#              Evaluates the Thermal Conductivity                  #
####################################################################
#
# NOTE: the low T approxumation requires SpecFunc_pp_xx.dat and SpecFunc_mm_xx.dat files
# evaluated for one single electronic energy value (the Fermi energy)
#
skies thermal-cond --range=[10,1400] < epw.skies.in
