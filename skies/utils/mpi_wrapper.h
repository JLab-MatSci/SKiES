/**
 @file
 @brief A basic C++ style wrapper over some MPI-routines
 @author Galtsov Ilya
 */
#include <mpi.h>
#include <skies/common/ndimarrays.h>

namespace skies { namespace mpi {

void init(int *argc, char ***argv);

void finalize();

int rank();

int size();

std::pair<std::vector<int>, std::vector<int>>
prepare_rcounts_displs(size_t size);

arrays::array2D sum_one_from_many(const arrays::array2D& part);
arrays::array3D sum_one_from_many(const arrays::array3D& part);

} // mpi
} // skies
