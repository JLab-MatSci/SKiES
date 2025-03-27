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

#include <math.h>
#include <string>

#include <skies/common/ndimarrays.h>

/**
 * \brief External function called from EPW Fortran code  to fill in lattice info
 * Information on a latttice constant, a number of atoms in a unit cell
 * and a matrix of primitive lattice vectors is read automatically during EPW initialization
 * @param lat_const lattice constant
 * @param nat number of atoms in a unit cell
 * @param coords matrix of primitive lattice vectors
*/
#ifdef __cplusplus
extern "C" {
#endif

void fillLattice(double* lat_const, double* unit_cell_vol, double coords[3][3]);

#ifdef __cplusplus
}
#endif

namespace skies {

/**
 * \brief Class which contains basic information on lattice
*/
class Lattprotocol {
public:
    static double latt_const;
    static double latt_volume;
    static arrays::array2D latt_coords; // contains three lattice vectors 
                            // written in rows in some real units
    /**
     * Lattice type
    */
    enum class Lattice
    {
        bcc = 0, 
        fcc = 1, 
        sc  = 2,
        undefined = -1
    };

    static Lattice latt_type;

public:
    Lattprotocol() {}
    /**
     * \brief Constructor from lattice type, lattice constant and number of atoms in a unit cell
     * @param lat lattice type
     * @param a lattice constant
     * @param nat number of atoms in a unit cell
    */
    Lattprotocol(Lattice lat, double a);

    /**
     * \brief Calculates the inverse-cell matrix
    */
    static arrays::array2D calc_inv_cell();

    // /**
    //  * \brief Calculates volume of a unit cell
    // */
    // static double calc_volume() {
    //     return (calc_3d_det(latt_coords) >= 0) ? calc_3d_det(latt_coords) : -calc_3d_det(latt_coords);
    // }

// some helper-functions
private:
    /**
     * \brief Calculates determinant of a 2x2 matrix
     * @param w matrix
    */
    static double calc_2d_det(const arrays::array2D& w)
    {
        if (w.size() != 2)
            throw std::runtime_error("Matrix 2x2 is needed");
        double res = w[0][0]*w[1][1] - w[0][1]*w[1][0];
        return res;
    }

    /**
     * \brief Calculates determinant of a 3x3 matrix
     * @param w matrix
    */
    static double calc_3d_det(const arrays::array2D& w)
    {
        if (w.size() != 3)
            throw std::runtime_error("Matrix 3x3 is needed");
        double res = w[0][0] * (w[1][1]*w[2][2] - w[1][2]*w[2][1])
                   - w[0][1] * (w[1][0]*w[2][2] - w[1][2]*w[2][0])
                   + w[0][2] * (w[1][0]*w[2][1] - w[1][1]*w[2][0]);
        return res;
    }

    /**
     * \brief Calculates the inverse of a 3x3 matrix
     * @param w matrix
    */
    static arrays::array2D calc_inv_3d(const arrays::array2D& w)
    {
        if (w.size() != 3)
            throw std::runtime_error("Matrix 3x3 is needed");
        double inv_det = 1.0 / calc_3d_det(w);
        arrays::array2D inv;
        inv.resize(3, arrays::array1D(3));
        inv[0][0] =  calc_2d_det({{w[1][1], w[2][1]}, 
                                       {w[1][2], w[2][2]}});
        inv[0][1] = -calc_2d_det({{w[0][1], w[2][1]}, 
                                       {w[0][2], w[2][2]}});
        inv[0][2] =  calc_2d_det({{w[0][1], w[1][1]}, 
                                       {w[0][2], w[1][2]}});

        inv[1][0] = -calc_2d_det({{w[1][0], w[2][0]}, 
                                       {w[1][2], w[2][2]}});
        inv[1][1] =  calc_2d_det({{w[0][0], w[2][0]}, 
                                       {w[0][2], w[2][2]}});
        inv[1][2] = -calc_2d_det({{w[0][0], w[1][0]}, 
                                       {w[0][2], w[1][2]}});

        inv[2][0] =  calc_2d_det({{w[1][0], w[2][0]}, 
                                       {w[1][1], w[2][1]}});
        inv[2][1] = -calc_2d_det({{w[0][0], w[2][0]},
                                       {w[0][1], w[2][1]}});
        inv[2][2] =  calc_2d_det({{w[0][0], w[1][0]}, 
                                       {w[0][1], w[1][1]}});
        for (auto&& v : inv)
            for (auto&& x : v)
                x *= inv_det;
        return inv;
    }
};

} // skies
