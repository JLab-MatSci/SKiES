/*-----------------------------------------------------------------------
    * SKiES - Solver of Kinetic Equation for Solids
    * Version 1.0.0
    * 
    * (C) 2024 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)
    *
    * SKiES may only be utilized for non-profit research.
    * Citing appropriate sources is required when using SKiES.
    * 
    * Project homepage: https://github.com/i-Galts/SKiES
    * 
    * This file: Electronic eigenvalue calculations header
    * Declares interfaces for band energy interpolation and
    * Fermi level determination
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */

#pragma once

#include <cassert>

#include <skies/quantities/basic_quantity.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Fills the number of electrons in the system.
 * 
 * This function populates the provided pointer with the number of electrons
 * calculated in the current context. It is intended for interoperability with
 * Fortran code.
 * 
 * @param nelec Pointer to an integer where the number of electrons will be stored.
 */
void fillNelec(int* nelec);

/**
 * @brief Fills the number of bands in the system.
 * 
 * This function populates the provided pointer with the number of bands
 * calculated in the current context. It is intended for interoperability with
 * Fortran code.
 * 
 * @param nbands Pointer to an integer where the number of bands will be stored.
 */
void fillNbands(int* nbands);

/**
 * @brief Interpolates eigenvalues at the specified k-point.
 * 
 * This function interpolates the eigenvalues at a specified k-point,
 * given by three pointers to its x-, y-, and z-components. The interpolated
 * eigenvalues are saved in the output array `eigenvalues[]`, which should
 * be allocated with a size equal to the number of bands.
 * 
 * @param kx Pointer to the x-component of the k-point.
 * @param ky Pointer to the y-component of the k-point.
 * @param kz Pointer to the z-component of the k-point.
 * @param eigenvalues Output array that will contain the interpolated eigenvalues.
 */
void interpEigenValueAt(const double* kx, const double* ky, const double* kz, double eigenvalues[]);

#ifdef __cplusplus
}
#endif

namespace skies { namespace quantities {

/**
 * @class EigenValueDrawable
 * @brief Class for handling eigenvalue data.
 * 
 * This class extends the AnyQuantity class and provides methods for
 * interpolating eigenvalue data at specified points.
 */
class EigenValueDrawable : public AnyQuantity {
public:
    static size_t nbands; ///< Number of bands

    /**
     * @brief Returns the name of the quantity.
     * @return A string representing the name of the quantity.
     */
    std::string name() const override { return "EigenValue"; }

    /**
     * @brief Interpolates eigenvalues at the specified k-point.
     * @param k The k-point at which to interpolate.
     * @return An array of interpolated eigenvalues.
     */
    arrays::array1D interpolate_at(const arrays::array1D& k) const override;
};

/**
 * @class EigenValue
 * @brief Class for handling eigenvalue calculations.
 * 
 * This class provides static methods for interpolating eigenvalues
 * at specified points.
 */
class EigenValue {
public:
    static size_t nbands; ///< Number of bands
    static int nelec;     ///< Number of electrons
    static double eF;     ///< Fermi level

    /**
     * @brief Finds the Fermi level using the bisection method.
     * 
     * This function calculates the Fermi level based on the number of electrons
     * and the specified temperature. It uses a numerical root-finding method
     * to solve the equation for the Fermi level.
     * 
     * @param TeV The temperature in electron volts.
     * @param crit The convergence criterion for the root-finding method.
     */
    static void find_eF(double TeV, double crit = 1e-16);

    /**
     * @brief Interpolates eigenvalues at the specified k-point.
     * 
     * This function takes a k-point represented by a 3D array and interpolates
     * the corresponding eigenvalues. The results are returned as an array of
     * eigenvalues relative to the Fermi level.
     * 
     * @param k The k-point at which to interpolate the eigenvalues.
     * @return An array of interpolated eigenvalues.
     */
    static arrays::array1D interpolate_at(const arrays::array1D& k);
};

} // quantities
} // skies
