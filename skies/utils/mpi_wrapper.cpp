#include <mpi.h>
#include <numeric>

#include <skies/common/alg.h>
#include <skies/utils/mpi_wrapper.h>

#include <iostream>

namespace skies { namespace mpi {

using namespace arrays;

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
    if (nproc <= static_cast<int>(size))
        return get_rcounts_displs(size, nproc);

    std::vector<int> rcounts(size, 1);
    std::vector<int> displs(size, 0);
    std::iota(displs.begin(), displs.end(), 0);
    return std::make_pair(rcounts, displs);
}

array2D sum_one_from_many(const array2D& part)
{
    auto n1 = part.size();
    auto n2 = part[0].size();
    auto size = n1 * n2;
    auto flatten_part = flatten(part);
    auto flatten = array1D(size, 0.0);
    MPI_Reduce(flatten_part.data(), flatten.data(), size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(flatten.data(), size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return reshape(flatten, n1, n2);
}

array3D sum_one_from_many(const array3D& part)
{
    auto n1 = part.size();
    auto n2 = part[0].size();
    auto n3 = part[0][0].size();
    auto size = n1 * n2 * n3;
    auto flatten_part = flatten(part);
    auto flatten = array1D(size, 0.0);
    MPI_Reduce(flatten_part.data(), flatten.data(), size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(flatten.data(), size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return reshape(flatten, n1, n2, n3);
}

} // mpi
} // skies
