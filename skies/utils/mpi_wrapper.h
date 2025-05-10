// #pragma once
#ifndef H_SKIES_MPI
#define H_SKIES_MPI

#include <mpi.h>

namespace skies {
namespace utils {
namespace mpi   {

#ifndef SKIES_MPI

inline bool is_root() {
	return true;
}

inline int rank() {
	return 0;
}

inline int size() {
	return 1;
}

#else

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

template <class T>
constexpr MPI_Datatype mpi_t();

template <>
constexpr MPI_Datatype mpi_t<bool>() {
    return MPI_C_BOOL;
}

template <>
constexpr MPI_Datatype mpi_t<int>() {
    return MPI_INT;
}

template <>
constexpr MPI_Datatype mpi_t<unsigned int>() {
    return MPI_UNSIGNED;
}

template <>
constexpr MPI_Datatype mpi_t<long long int>() {
    return MPI_LONG_LONG_INT;
}

template <>
constexpr MPI_Datatype mpi_t<unsigned long long>() {
    return MPI_UNSIGNED_LONG_LONG;
}

template <>
constexpr MPI_Datatype mpi_t<std::size_t>() {
    return MPI_UNSIGNED_LONG_LONG;
}

template <>
constexpr MPI_Datatype mpi_t<double>() {
    return MPI_DOUBLE;
}

#endif

} // mpi
} // utils
} // skies

#endif