/**
 @file
 @brief A basic C++ style wrapper over some MPI-routines
 @author Galtsov Ilya
 */
#include <mpi.h>

namespace skies { namespace mpi {

void init(int *argc, char ***argv);

void finalize();

int rank();

int size();

} // mpi
} // skies
