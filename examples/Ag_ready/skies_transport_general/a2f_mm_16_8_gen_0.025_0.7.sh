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

# NOTE: --epsilons option is used for general Allen's formulas calculations.
# The number of bins in the interval is set equal to 56 to obtain integration step of 0.025 eV.
# Additionally the range of phonon energies is supplied (in meV)
#
skies a2f --kgrid=[16,16,16] --qgrid=[8,8,8] --signs=[-1,-1] --bands=[5,5] --eF=12.9148 --tetra --epsilons=[-0.7,0.7] --eps_bins=56 --omegas=[0.1,25] < epw.skies.in
