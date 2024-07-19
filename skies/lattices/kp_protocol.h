/**
 @file
 @brief Description of k-point mesh class used in spectral function calculation
 @author Galtsov Ilya
 */
#pragma once

#include <unordered_map>
#include <skies/common/ndimarrays.h>

namespace skies {

/**
 * \brief Class which contains basic information on reciprocal lattice
*/
class KPprotocol {
private:
    size_t n1_;
    size_t n2_;
    size_t n3_;
public:
    size_t nkpt;          // num of k-points
    arrays::array2D grid; // contains crystal coords of kpts
    std::vector<std::vector<int>> igrid; // contains integer coords of kpts (MP-scheme)
public:
    KPprotocol() {}

    /**
     * \brief Constructor of MP-scheme
     * @param n1 number of k-points in kx direction
     * @param n2 number of k-points in ky direction
     * @param n3 number of k-points in kz direction
    */
    KPprotocol(int n1, int n2, int n3);

    std::tuple<size_t, size_t, size_t> mesh() const;

    /**
     * @brief creates local subcell of a given k-point
     * @param ind k-point index
     * @return vector of 8 integers corresponding to subcell points around a given point.
     *         Additionally exploits periodic boundary conditions at edge points
     */
    std::vector<size_t> local_subcell(size_t ind) const;

    std::vector<size_t> range() const;

    arrays::array1D find_vd_from_ind(size_t ind) const;
    std::vector<int> find_vi_from_ind(size_t ind) const;
    size_t find_ind_from_vi(const std::vector<int>&) const;

private:
    std::vector<size_t> kprange_;
    std::unordered_map<size_t, arrays::array1D> vd_from_ind_;
    std::unordered_map<size_t, std::vector<int>> vi_from_ind_;

    void apply_pbc(std::vector<std::vector<int>>& subcell_v) const;
};

} // skies
