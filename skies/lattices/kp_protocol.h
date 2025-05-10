/*-----------------------------------------------------------------------
    * SKiES - Solver of Kinetic Equation for Solids
    * 
    * (C) 2025 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)
    *
    * SKiES may only be utilized for non-profit research.
    * Citing appropriate sources is required when using SKiES.
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */
#pragma once

#include <memory>
#include <unordered_map>

#include <skies/common/ndimarrays.h>
#include <skies/lattices/decompose.h>

namespace skies {

/**
 * \brief Class which contains basic information on reciprocal lattice
*/
class KPprotocol {
public:
    KPprotocol();
    /**
     * \brief Constructor of MP-scheme
     * @param n1 number of k-points in kx direction
     * @param n2 number of k-points in ky direction
     * @param n3 number of k-points in kz direction
    */
    //KPprotocol(size_t n1, size_t n2, size_t n3);
    //KPprotocol(const KPprotocol& oth);
    //KPprotocol& operator= (const KPprotocol& that);

    KPprotocol(std::size_t n1, std::size_t n2, std::size_t n3,
               std::shared_ptr<Decomposer> decomposer = nullptr);

    std::tuple<size_t, size_t, size_t> mesh() const;

    /**
     * @brief creates local subcell of a given k-point
     * @param ind k-point index
     * @return vector of 8 integers corresponding to subcell points around a given point.
     *         Additionally exploits periodic boundary conditions at edge points
     */
    std::vector<size_t> local_subcell(size_t ind) const;

    const std::vector<size_t>& range() const &;

    arrays::array1D find_vd_from_ind(size_t ind) const;
    std::vector<int> find_vi_from_ind(size_t ind) const;
    size_t find_ind_from_vi(const std::vector<int>&) const;

    size_t nkpt() const { return nkpt_; }
    const arrays::array2D& grid() const & { return grid_; }
    arrays::array2D&& grid() && { return std::move(grid_); }
    const std::vector<std::vector<int>>& igrid() const & { return igrid_; }

    std::size_t npt_loc() const { return decomposer_->npt_loc(); }
    std::size_t index_glob(std::size_t idx_loc) const { return decomposer_->index_glob(idx_loc); }
    arrays::array2D grid_loc() const { return decomposer_->grid_loc(grid_); }
    void reduce(arrays::array1D& inout) const { decomposer_->reduce(inout); }
    bool is_parallel() const { return decomposer_->is_parallel(); }

    const std::vector<std::size_t>& inds_loc() const { return decomposer_->inds_loc(); }

private:
    void apply_pbc(std::vector<std::vector<int>>& subcell_v) const;

private:
    size_t n1_;
    size_t n2_;
    size_t n3_;
    size_t nkpt_;          // num of k-points
    arrays::array2D grid_; // contains crystal coords of kpts
    std::vector<std::vector<int>> igrid_; // contains integer coords of kpts (MP-scheme)

    std::vector<size_t> kprange_;
    std::unordered_map<size_t, arrays::array1D> vd_from_ind_;
    std::unordered_map<size_t, std::vector<int>> vi_from_ind_;

    friend bool operator== (const KPprotocol&, const KPprotocol&);

    std::shared_ptr<Decomposer> decomposer_;
};

bool operator== (const KPprotocol& kp1, const KPprotocol& kp2);

} // skies
