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

#include <skies/common/ndimarrays.h>

namespace skies { namespace transport {

/**
 * @brief calculates inner integral in formulas (39) and (40) in Smirnov N. A., 2022 "Ab initio caclculations for transport ..."
 *
 * This function computes the inner integral for inelastic scattering based on the given temperature,
 * frequency array, and a2f matrix. It is used in the context of electron-phonon interactions.
 *
 * @param Temp The temperature in Kelvin [K].
 * @param freqs A 1D array representing the frequency values.
 * @param a2f A 2D array representing the electron-phonon coupling matrix.
 * @return double The computed value of the inner integral.
 */
double calc_inner_integral_inelastic(double Temp, const arrays::array1D& freqs, const arrays::array2D& a2f);

/**
 * @brief Calculates the inner integral for inelastic scattering using separate a2f_plus and a2f_minus matrices.
 *
 * This function computes the inner integral for inelastic scattering based on the given temperature,
 * frequency array, and two separate electron-phonon coupling matrices (a2f_plus and a2f_minus).
 *
 * @param Temp The temperature in Kelvin [K].
 * @param freqs A 1D array representing the frequency values.
 * @param a2f_plus A 2D array representing the positive component of the electron-phonon coupling matrix.
 * @param a2f_minus A 2D array representing the negative component of the electron-phonon coupling matrix.
 * @return double The computed value of the inner integral.
 */
double calc_inner_integral_inelastic(double Temp, const arrays::array1D& freqs,
                                     const arrays::array2D& a2f_plus, const arrays::array2D& a2f_minus);

/**
 * @brief Calculates the electronic conductivity due to inelastic scattering.
 *
 * This function computes the electronic conductivity for a range of temperatures using the provided
 * a2f file and writes the results to an output file. The unit cell volume is used for normalization.
 *
 * @param Temps A 1D array representing the range of temperatures in Kelvin [K].
 * @param a2f_fnm The file name containing the a2f data (electron-phonon coupling matrix).
 * @param cond_fnm The file name where the computed conductivity will be written.
 * @param unit_cell_vol The volume of the unit cell used for normalization.
 */
void calc_elec_cond_inelastic(const arrays::array1D& Temps, const char* a2f_fnm, const char* cond_fnm, double unit_cell_vol);

/**
 * @brief Calculates the thermal conductivity due to inelastic scattering.
 *
 * This function computes the thermal conductivity for a range of temperatures using the provided
 * a2f_plus and a2f_minus files and writes the results to an output file. The unit cell volume is used for normalization.
 *
 * @param Temps A 1D array representing the range of temperatures in Kelvin [K].
 * @param a2f_plus_fnm The file name containing the positive component of the a2f data.
 * @param a2f_minus_fnm The file name containing the negative component of the a2f data.
 * @param cond_fnm The file name where the computed thermal conductivity will be written.
 * @param unit_cell_vol The volume of the unit cell used for normalization.
 */
void calc_therm_cond_inelastic(const arrays::array1D& Temps, const char* a2f_plus_fnm, const char* a2f_minus_fnm,
                              const char* cond_fnm, double unit_cell_vol);

} // transport
} // skies
