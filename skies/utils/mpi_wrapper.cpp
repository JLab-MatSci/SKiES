#include <mpi.h>
#include "mpi_wrapper.h"

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

} // mpi
} // skies
