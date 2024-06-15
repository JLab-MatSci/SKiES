/**
 @file
 @brief Description of k-point mesh class used in spectral function calculation
 @author Galtsov Ilya
 */
#pragma once

#include <skies/common/ndimarrays.h>

namespace skies {

/**
 * \brief Class which contains basic information on reciprocal lattice
*/
class KPprotocol {
public:
    size_t nkpt;          // num of k-points
    arrays::array2D grid; // contains coords of kpts in format [x, y, z]
public:
    KPprotocol() {}

    /**
     * \brief Constructor of MP-scheme
     * @param n1 number of k-points in kx direction
     * @param n2 number of k-points in ky direction
     * @param n3 number of k-points in kz direction
    */
    KPprotocol(int n1, int n2, int n3);

};

} // skies
