/**
 @file
 @brief A basic C++ style wrapper over some MPI-routines
 @author Galtsov Ilya
 */
#include <mpi.h>
#include <vector>

namespace skies { namespace mpi {

void init(int *argc, char ***argv);

void finalize();

int rank();

int size();

std::pair<std::vector<int>, std::vector<int>>
prepare_rcounts_displs(size_t size);

} // mpi
} // skies
