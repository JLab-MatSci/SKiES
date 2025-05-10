#include <skies/utils/mpi_wrapper.h>
#include <skies/lattices/decompose.h>

namespace skies {

using namespace skies::arrays;
using namespace skies::utils;

MPIDecomposer::MPIDecomposer(std::size_t npt, int rank, int size)
	: rank_(rank), size_(size) {
		for (std::size_t i = rank; i < npt; i += size) {
			inds_loc_.push_back(i);
		}
}

array2D MPIDecomposer::grid_loc(const array2D& grid_glob) const {
	array2D res;
	res.reserve(inds_loc_.size());
	for (auto idx : inds_loc_)
		res.push_back(grid_glob[idx]);
	return res;
}

void MPIDecomposer::reduce(array1D& inout) const {
	array1D recv_buffer(inout.size(), 0.0);
#ifdef SKIES_MPI
	MPI_Reduce(inout.data(), recv_buffer.data(), inout.size(),
			   mpi::mpi_t<double>(), MPI_SUM, 0, MPI_COMM_WORLD);
#endif
	if (rank_ == 0)
		inout = std::move(recv_buffer);
}


} // skies
