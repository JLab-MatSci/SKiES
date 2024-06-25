#include <mpi.h>
#include <numeric>

#include <skies/common/alg.h>
#include <skies/utils/mpi_wrapper.h>

#include <iostream>

namespace skies { namespace mpi {

void init(int *argc, char ***argv)
{
    MPI_Init(argc, argv);
}

void finalize()
{
    MPI_Finalize();
}

int rank()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int size()
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

std::pair<std::vector<int>, std::vector<int>>
prepare_rcounts_displs(size_t size)
{
    int nproc = mpi::size();
    if (nproc <= size)
        return get_rcounts_displs(size, nproc);

    std::vector<int> rcounts(size, 1);
    std::vector<int> displs(size, 0);
    std::iota(displs.begin(), displs.end(), 0);
    return std::make_pair(rcounts, displs);
}

} // mpi
} // skies
