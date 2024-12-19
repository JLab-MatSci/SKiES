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
    * This file: Eigenfrequency calculations header
    * Declares interfaces for handling eigenfrequencies and eigenvalues,
    * used primarily in DOS and transport calculations.
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */
#pragma once

#include <skies/quantities/basic_quantity.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief External function called from EPW Fortran code to obtain the number of phonon modes in the calculation.
 * 
 * This function fills the provided pointer with the number of phonon modes
 * calculated in the current context. It is intended for interoperability with Fortran.
 * 
 * @param nmodes Pointer to an integer where the number of modes will be stored.
 */
void fillNmodes(int* nmodes);

/**
 * @brief External function called from EPW Fortran code to interpolate eigen frequency at provided k-point.
 * 
 * This function interpolates the eigen frequencies at a specified k-point,
 * given by three pointers to its x-, y-, and z-components. The interpolated
 * eigen frequencies are saved in the output array `rfreqs[]`, which should
 * be allocated with a size equal to the number of phonon bands.
 * 
 * @param kx Pointer to the x-component of the k-point.
 * @param ky Pointer to the y-component of the k-point.
 * @param kz Pointer to the z-component of the k-point.
 * @param rfreqs Output array that will contain the interpolated eigen frequencies (in Ry).
 */
void interpEigenFreqAt(const double* kx, const double* ky, const double* kz, double rfreqs[]);

/**
 * @brief Interpolates eigen frequencies in one dimension at the provided k-point.
 * 
 * This function performs interpolation of eigen frequencies at a specified k-point
 * in one dimension. The output frequencies are stored in the provided array.
 * 
 * @param kx Pointer to the x-component of the k-point.
 * @param ky Pointer to the y-component of the k-point.
 * @param kz Pointer to the z-component of the k-point.
 * @param rfreqs Output array that will contain the interpolated eigen frequencies.
 */
void interpEigenFreq1DAt(const double* kx, const double* ky, const double* kz, double rfreqs[]);

#ifdef __cplusplus
}
#endif

namespace skies { namespace quantities {

/**
 * @class EigenFrequencyDrawable
 * @brief Class for handling eigenfrequency data.
 * 
 * This class extends the AnyQuantity class and provides methods for
 * interpolating eigenfrequency data at specified points.
 */
class EigenFrequencyDrawable : public AnyQuantity {
public:
    static size_t nmodes; ///< Number of modes

    /**
     * @brief Returns the name of the quantity.
     * @return A string representing the name of the quantity.
     */
    std::string name() const override { return "EigenFrequency"; }

    /**
     * @brief Interpolates eigenfrequencies at the specified q-point.
     * @param q The q-point at which to interpolate.
     * @return An array of interpolated eigenfrequencies.
     */
    arrays::array1D interpolate_at(const arrays::array1D& q) const override;
};

/**
 * @class EigenFrequency
 * @brief Class for handling eigenfrequency calculations.
 * 
 * This class provides static methods for interpolating eigenfrequencies
 * at specified points.
 */
class EigenFrequency {
public:
    static size_t nmodes; ///< Number of modes

    /**
     * @brief Interpolates eigenfrequencies at the specified q-point.
     * @param q The q-point at which to interpolate.
     * @return An array of interpolated eigenfrequencies.
     */
    static arrays::array1D interpolate_at(const arrays::array1D& q);
};

} // quantities
} // skies
