#include <cmath>
#include <fstream>

#include <skies/common/units.h>
#include <skies/transport/iohandler.h>
#include <skies/transport/inelastic.h>

#include <iostream>

namespace skies { namespace transport {

using namespace skies::arrays;

using skies::units::pi;
using skies::units::hbar;
using skies::units::A_in_bohr;
using skies::units::eV_in_Ry;
using skies::units::rau_in_muOm_cm;

using skies::units::kB;
using skies::units::kB_as_Ry_over_K;
using skies::units::rau_in_cm_over_W;
using skies::units::rau_in_m_over_s;
using skies::units::Ry_in_J;

//////////////////////  Resistivity  //////////////////////

// calculates inner integral in formulas (39) and (40) in Smirnov N. A., 2022 "Ab initio caclculations for transport ..."
double calc_inner_integral_inelastic(double TeV, const array1D& freqs, const array2D& a2f)
{
    double dw = freqs[1] - freqs[0];

    double integ = 0.0;
    for (size_t i = 1; i < freqs.size(); ++i)
    {
        double x = 0.5 * freqs[i] / TeV;
        integ += (1.0 / freqs[i]) * (x * x / sinh(x) / sinh(x)) * a2f[i][0];
    }
    integ *= dw;
    return integ;
}

void calc_elec_cond_inelastic(const array1D& Temps, const char* a2f_fnm, const char* cond_fnm, double unit_cell_vol)
{
    IHandler ihandler(a2f_fnm);
    if (ihandler.epsilons().size() > 1)
        throw std::runtime_error("There must be just Fermi level in low T formulas.");
    array1D resist;
    // unit cell volumes must be given in [bohr^3]

    for (size_t i = 0; i < Temps.size(); ++i)
    {
        double Temp = Temps[i];
        double TeV = Temp / 11606.0;
        double integ = calc_inner_integral_inelastic(TeV, ihandler.omegas(), ihandler.a2f());
        double e2 = 2.0;
        double prefactor = 2.0 * pi * unit_cell_vol * TeV * eV_in_Ry * rau_in_muOm_cm / (e2 * ihandler.transDOSes()[0] * units::Ry_in_eV);
        resist.push_back(integ * prefactor);
    }
    OHandler ohandler(a2f_fnm, cond_fnm, ResistType::Electrical);
    ohandler.dump(Temps, resist);
}

//////////////////////  Thermal Electron Resistivity  //////////////////////

double calc_inner_integral_inelastic(double TeV, const array1D& freqs,
                                     const array2D& a2f_plus, const array2D& a2f_minus)
{
    double dw = freqs[1] - freqs[0];
    double integ = 0.0;
    for (size_t i = 1; i < freqs.size(); ++i)
    {
        double x = 0.5 * freqs[i] / TeV;
        integ += (1.0 / freqs[i]) * (x * x / sinh(x) / sinh(x))
               * (a2f_plus[i][0] * (1.0 + x * x / pi / pi) + a2f_minus[i][0] * (3.0 * x * x / pi / pi));
    }
    integ *= dw;
    return integ;
}

void calc_therm_cond_inelastic(const array1D& Temps, const char* a2f_plus_fnm, const char* a2f_minus_fnm,
                        const char* cond_fnm, double unit_cell_vol)
{
    IHandler iplusHandler(a2f_plus_fnm), iminusHandler(a2f_minus_fnm);
    if (iplusHandler.epsilons().size() > 1 || iminusHandler.epsilons().size() > 1)
        throw std::runtime_error("There must be just Fermi level in low T formulas.");
    array1D kappas;

    std::cout << 6.0 * unit_cell_vol * 100.0 / (pi * iplusHandler.transDOSes()[0] * units::Ry_in_eV * kB * hbar) << std::endl;
    for (size_t i = 0; i < Temps.size(); ++i)
    {
        double Temp = Temps[i];
        double TeV = Temp / 11606.0;
        double integ = calc_inner_integral_inelastic(TeV, iplusHandler.omegas(),
                                                     iplusHandler.a2f(), iminusHandler.a2f());
        kappas.push_back(pi * iplusHandler.transDOSes()[0] * units::Ry_in_eV * kB_as_Ry_over_K / (6.0 * unit_cell_vol * integ * rau_in_cm_over_W));
    }
    
    OHandler ohandler(a2f_plus_fnm, cond_fnm, ResistType::Thermal);
    ohandler.dump(Temps, kappas);
}

} // transport
} // skies
