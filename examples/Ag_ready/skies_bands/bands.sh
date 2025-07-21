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
#                    Evaluates the Band Structure                    #
######################################################################
#
# NOTE: a file called Kpath must be supplied in this folder
# An example is given for a fcc lattice
#
skies bands --eF=12.9148 < epw.skies.in
