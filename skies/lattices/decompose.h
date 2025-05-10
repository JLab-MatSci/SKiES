#pragma once

#include <numeric>

#include <skies/common/ndimarrays.h>

namespace skies {

class Decomposer {
public:
	virtual ~Decomposer() = default;

	virtual std::size_t npt_loc() const = 0;

	virtual std::size_t index_glob(std::size_t idx_loc) const = 0;

	virtual arrays::array2D grid_loc(const arrays::array2D& grid_glob) const = 0;

	virtual void reduce(arrays::array1D& inout) const = 0;

	virtual bool is_parallel() const = 0;

	virtual const std::vector<std::size_t>& inds_loc() const = 0;
};

class SerialDecomposer : public Decomposer {
public:
	explicit SerialDecomposer(std::size_t npt) : npt_(npt), inds_loc_(npt_) {
		std::iota(inds_loc_.begin(), inds_loc_.end(), 0);
	}

	std::size_t npt_loc() const override { return npt_; }
	std::size_t index_glob(std::size_t idx_loc) const override { return idx_loc; }
	arrays::array2D grid_loc(const arrays::array2D& grid_glob) const override { return grid_glob; }
	void reduce(arrays::array1D& inout) const override {}
	bool is_parallel() const override { return false; }

	const std::vector<std::size_t>& inds_loc() const { return inds_loc_; }

private:
	std::size_t npt_;
	std::vector<std::size_t> inds_loc_;
};


class MPIDecomposer : public Decomposer {
public:

	// creates a number of local subgrids in a round-robin fashion
	MPIDecomposer(std::size_t npt, int rank, int size);

	std::size_t npt_loc() const override { return 1000; }
	std::size_t index_glob(std::size_t idx_loc) const override { return idx_loc; }
	arrays::array2D grid_loc(const arrays::array2D& grid_glob) const override;
	void reduce(arrays::array1D& inout) const override;
	bool is_parallel() const override { if (size_ == 1) return false; return true; }

	const std::vector<std::size_t>& inds_loc() const { return inds_loc_; }

private:
	int rank_, size_;
	std::vector<std::size_t> inds_loc_;
};

} // skies