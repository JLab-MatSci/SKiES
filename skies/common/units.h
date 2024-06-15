/**
 @file
 @brief Physical constants and conversion from Rydberg atomic units (r.a.u.) to others needed in the project
 @author Galtsov Ilya
 */
# pragma once

#include <cmath>

namespace skies { namespace units {
    // numerical constants
    const double pi = 4.0 * std::atan(1.0);

    // in SGSE system
    const double c = 299792458.;
    const double mu0 = 4.0e-7 * pi;
    const double hplanck = 6.62607015e-34;
    const double eV = 1.602176634e-19;
    const double me = 9.1093837015e-31;
    const double mp = 1.67262192369e-27;
    const double Nav = 6.02214076e23;
    const double kB = 1.380649e-23;
    const double amu = 1.66053906660e-27;
    const double eps0 = (1 / mu0 / c*c);
    const double hbar = hplanck / (2 * pi);

    // conversion
    const double Ha_in_eV = me * eV*eV*eV / 16 / pi / pi / eps0 / eps0 / hbar / hbar;
    const double eV_in_Ry = 0.073498810939358;
    const double Ry_in_eV = 1.0 / eV_in_Ry;
    const double Ry_in_J = 2.1798741 * 1e-18;
    const double A_in_bohr = 1.8897259886;
    const double bohr_in_A = 1.0 / A_in_bohr;

    const double kB_as_Ry_over_K = kB / Ry_in_J;
    const double rau_in_cm_over_W = 1.1743967859e-7;
    const double cm_over_W_in_rau = 1.0 / rau_in_cm_over_W;
    const double rau_in_m_over_s = 1.094 * 1e6; 

    const double hbar_in_eV_s = hbar * 6.241509e18;

    const double rau_in_muOm_cm = 100.0 / 2.2999241;

    const double rau_in_muV = 19.241363 * 1e6;

} // units
} // skies
