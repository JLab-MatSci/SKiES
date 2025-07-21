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

###################################################################################################
#          Evaluates the Transport Spectral Function in the Low Temperature Approximation         #
###################################################################################################
#
# NOTE: --bands=[5,5] is used to only consider the 5th wide band in the calculation.
# It significantly speeds up the calculations.
#
skies a2f --kgrid=[16,16,16] --qgrid=[8,8,8] --signs=[1,1] --bands=[5,5] --eF=12.9148 --tetra < epw.skies.in
