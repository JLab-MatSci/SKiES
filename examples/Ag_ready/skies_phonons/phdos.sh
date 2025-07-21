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

######################################################################
#          Evaluates the Phonon Density of States (PhDOS)            #
######################################################################
#
#
skies phdos --grid=[30,30,30] --tetra < epw.skies.in
