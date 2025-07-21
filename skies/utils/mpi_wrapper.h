#pragma once

#include <mpi.h>

namespace skies {
namespace utils {
namespace mpi   {


struct mpi_handler {
	mpi_handler() {
		MPI_Init(NULL, NULL);
	}

	~mpi_handler() {
		MPI_Finalize();
	}
};

inline int rank(MPI_Comm comm = MPI_COMM_WORLD) {
	int rank;
	MPI_Comm_rank(comm, &rank);
	return rank;
}

inline int size(MPI_Comm comm = MPI_COMM_WORLD) {
	int size;
	MPI_Comm_size(comm, &size);
	return size;
}

inline bool is_root(MPI_Comm comm = MPI_COMM_WORLD) {
	return rank() == 0;
}

} // mpi
} // utils
} // skies
