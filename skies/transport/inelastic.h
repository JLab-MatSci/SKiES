#pragma once

#include <skies/common/ndimarrays.h>

namespace skies { namespace transport {

// calculates inner integral in formulas (39) and (40) in Smirnov N. A., 2022 "Ab initio caclculations for transport ..."
double calc_inner_integral_inelastic(double Temp, const arrays::array1D& freqs, const arrays::array2D& a2f);
double calc_inner_integral_inelastic(double Temp, const arrays::array1D& freqs,
                                     const arrays::array2D& a2f_plus, const arrays::array2D& a2f_minus);

// Temps are given in [K]
void calc_elec_cond_inelastic(const arrays::array1D& Temps, const char* a2f_fnm, const char* cond_fnm, double unit_cell_vol);
void calc_therm_cond_inelastic(const arrays::array1D& Temps, const char* a2f_plus_fnm, const char* a2f_minus_fnm,
                              const char* cond_fnm, double unit_cell_vol);

} // transport
} // skies
