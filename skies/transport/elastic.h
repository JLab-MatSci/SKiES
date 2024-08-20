/**
 * @file
 * @brief Calculate transport properties using Allen's formulas accounting for elastic scattering
 * @author Galtsov Ilya
 */
#pragma once

#include <skies/common/ndimarrays.h>

namespace skies { namespace transport {

using innerIntegCalculator = std::function<double(int sign, double Te, double omega,
                                           const arrays::array1D& epsilons, const arrays::array1D& transDOSes,
                                           const arrays::array1D& a2f, int s, int s_prime)>;

/**
 * \brief Calculates integral over \epsilon in (29) 2022 Smirnov N.A. Phys. Rev. B. V.106.
 * @param sign sign of The Times...
 * @param Te current temperature in eV
 * @param omega current phonon frequency in eV
 * @param epsilons array of electronic energies to integrate over
 * @param transDOSes array of smeared values of transport DOS evaluated for these electronic energies
 * @param a2f array of transport spectral functions evaluated for these electronic energies and given phonon enenrgy
*/
double calc_inner_integral_elec_elastic(int sign, double Te, double omega, const arrays::array1D& epsilons,
                                   const arrays::array1D& transDOSes, const arrays::array1D& a2f, int s = 1, int s_prime = 1);
void calc_elec_cond_elastic(const arrays::array1D& Temps, const arrays::array1D& ion_Temps, const char* a2f_fnm, const char* cond_fnm, double unit_cell_vol);
void calc_therm_cond_elastic(const arrays::array1D& Temps, const arrays::array1D& ion_Temps,
                             const char* a2f_plus_fnm, const char* a2f_minus_fnm, const char* cond_fnm,
                             double unit_cell_vol, const std::string& a2f_pm_fnm = "", const std::string& a2f_mp_fnm = "");
void calc_seebeck_elastic(const arrays::array1D& Temps, const arrays::array1D& ion_Temps,
                          const char* a2f_plus_fnm, const char* a2f_minus_fnm, const char* cond_fnm,
                          double unit_cell_vol, const std::string& a2f_pm_fnm = "", const std::string& a2f_mp_fnm = "");

// small helper function for inverting Q-matrix
inline arrays::array2D calc_inv_2d(const arrays::array2D& w)
{
    if (w.size() != 2)
        throw std::runtime_error("Matrix 2x2 is needed");
    double det = w[0][0] * w[1][1] - w[0][1] * w[1][0];
    arrays::array2D inv;
    inv.resize(2, arrays::array1D(2));
    inv[0][0] =  w[1][1];
    inv[0][1] = -w[0][1];
    inv[1][0] = -w[1][0];
    inv[1][1] =  w[0][0];
    for (auto&& v : inv)
        for (auto&& x : v)
            x /= det;
    return inv;
}

} // transport
} // skies
