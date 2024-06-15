#include <cassert>
#include <numeric>
#include <functional>

#include <skies/common/units.h>
#include <skies/sampling/sampling.h>
#include <skies/transport/elastic.h>
#include <skies/transport/iohandler.h>

#include <skies/quantities/eigenvals.h>

#include <iostream>

namespace skies { namespace transport {

using namespace skies::arrays;
using namespace skies::bzsampling;

using skies::units::pi;
using skies::units::hbar;
using skies::units::A_in_bohr;
using skies::units::eV_in_Ry;
using skies::units::rau_in_muOm_cm;

using skies::units::kB;
using skies::units::rau_in_m_over_s;
using skies::units::Ry_in_J;

using skies::units::kB_as_Ry_over_K;
using skies::units::rau_in_cm_over_W;

using skies::units::rau_in_muV;

array1D fermi_dirac_factor(int sign, double Te, double omega, const array1D& epsilons)
{
    assert(epsilons.size() > 1);
    array1D factor(epsilons.size(), 0.0);
    std::transform(epsilons.begin(), epsilons.end(), factor.begin(),
                   [sign, Te, omega] (double e) { return fermi_dirac(e, Te) * (1 - fermi_dirac(e + sign * omega, Te)); });
    return factor;
}

array1D bose_einstein_factor(double Ti, const array1D& omegas, int add)
{
    assert(omegas.size() > 1);
    array1D factor(omegas.size(), 0.0);
    std::transform(omegas.begin(), omegas.end(), factor.begin(),
                   [Ti, add] (double om) { return bose_einstein(om, Ti) + add; });
    return factor;
}

array1D calc_a2f_over_transDOS(const array1D& epsilons, const array1D& transDOSes, const array1D& a2f)
{
    array1D a2f_over_transDOS(epsilons.size(), 0.0);
    std::vector<int> inds(epsilons.size(), 0); std::iota(inds.begin(), inds.end(), 0);
    std::transform(inds.begin(), inds.end(), a2f_over_transDOS.begin(),
                   [a2f, transDOSes] (size_t i) { return a2f[i] / transDOSes[i]; });
    return a2f_over_transDOS;
}

array1D calc_external_integral(const array1D& Temps, int ord, innerIntegCalculator calc_inner_integ, const char* a2f_fnm)
{
    IHandler ihandler(a2f_fnm);
    if (ihandler.epsilons().size() < 2)
        throw std::runtime_error("There must be at least two values of electron energies in inelastic formulas.");
    if (ihandler.omegas().size() < 2)
        throw std::runtime_error("There must be at least two values of phonon frequencies.");
    assert(ord == 0 || ord == 1 || ord == 2);

    auto omegas = ihandler.omegas();
    double dom  = omegas[1] - omegas[0];

    auto epsilons = ihandler.epsilons();
    auto transDOSes = ihandler.transDOSes();
    auto a2f = ihandler.a2f();

    array1D ext_integs;
    for (size_t i = 0; i < Temps.size(); ++i)
    {
        double Temp = Temps[i];
        double TeV = Temp / 11606.0;

        array1D epsilon_integs_plus, epsilon_integs_minus;
        for(size_t iom = 0; iom < omegas.size(); ++iom)
        {
            epsilon_integs_plus.push_back(calc_inner_integ(  1, TeV, omegas[iom], epsilons, transDOSes, a2f[iom]));
            epsilon_integs_minus.push_back(calc_inner_integ(-1, TeV, omegas[iom], epsilons, transDOSes, a2f[iom]));
        }

        auto bose_factor_plus  = bose_einstein_factor(TeV, omegas, 0);
        auto bose_factor_minus = bose_einstein_factor(TeV, omegas, 1);
        double omega_integ_plus = dom * std::inner_product(bose_factor_plus.begin(), bose_factor_plus.end(), epsilon_integs_plus.begin(), 0.0);
        double omega_integ_minus = dom * std::inner_product(bose_factor_minus.begin(), bose_factor_minus.end(), epsilon_integs_minus.begin(), 0.0);

        ext_integs.push_back((omega_integ_plus + omega_integ_minus) / std::pow(TeV, ord + 1));
    }
    return ext_integs;
}

// Q_{x0,x0}
double calc_inner_integral_elec_elastic(int sign, double Te, double omega, const array1D& epsilons,
                                   const array1D& transDOSes, const array1D& a2f)
{
    assert(epsilons.size() > 1);
    double inner_integ{ 0.0 };
    double de = epsilons[1] - epsilons[0];
    auto neps = epsilons.size();

    auto a2f_over_transDOS = calc_a2f_over_transDOS(epsilons, transDOSes, a2f);
    auto factor = fermi_dirac_factor(sign, Te, omega, epsilons);
    inner_integ = de * std::inner_product(factor.begin(), factor.end(), a2f_over_transDOS.begin(), 0.0);
    return inner_integ;
}

void calc_elec_cond_elastic(const array1D& Temps, const char* a2f_fnm, const char* cond_fnm, double unit_cell_vol)
{
    // unit cell volumes must be given in [bohr^3]
    double prefactor = 2.0 * pi * unit_cell_vol * rau_in_muOm_cm / eV_in_Ry;
    // multiply by ev_in_Ry^2 due to integration over de and dOm
    prefactor *= eV_in_Ry * eV_in_Ry;

    auto Qx0x0  = prefactor * calc_external_integral(Temps, 0, calc_inner_integral_elec_elastic, a2f_fnm);
    double e2 = 2.0;
    auto resist = Qx0x0 * (0.5 / e2); // e2 = 2.0 in r.a.u.
    
    OHandler ohandler(a2f_fnm, cond_fnm, ResistType::Electrical);
    ohandler.dump(Temps, resist);
}

// Q_{x0,x1}
double calc_inner_integral_therm_elastic_ord1(int sign, double Te, double omega, 
                                              const array1D& epsilons, const array1D& transDOSes, const array1D& a2f)
{
    assert(epsilons.size() > 1);
    double inner_integ{ 0.0 };
    double de = epsilons[1] - epsilons[0];
    auto neps = epsilons.size();

    auto a2f_over_transDOS = calc_a2f_over_transDOS(epsilons, transDOSes, a2f);

    array1D epsilon_first_moment(neps, 0.0);
    std::vector<int> inds(neps, 0); std::iota(inds.begin(), inds.end(), 0);
    std::transform(inds.begin(), inds.end(), epsilon_first_moment.begin(),
                   [sign, omega, epsilons, a2f_over_transDOS] (size_t i) { return (2.0 * epsilons[i] + sign * omega) * a2f_over_transDOS[i]; });

    auto factor = fermi_dirac_factor(sign, Te, omega, epsilons);
    inner_integ = de * std::inner_product(factor.begin(), factor.end(), epsilon_first_moment.begin(), 0.0);
    return inner_integ;
}

// Q_{x1,x1}, s = s' = +1
double calc_inner_integral_therm_elastic_ord2_p(int sign, double Te, double omega, 
                                                const array1D& epsilons, const array1D& transDOSes, const array1D& a2f)
{
    // TODO assert that sign == +1 in a2f file
    assert(epsilons.size() > 1);
    double inner_integ{ 0.0 };
    double de = epsilons[1] - epsilons[0];
    auto neps = epsilons.size();

    auto a2f_over_transDOS = calc_a2f_over_transDOS(epsilons, transDOSes, a2f);

    array1D epsilon_second_moment(neps, 0.0);
    std::vector<int> inds(neps, 0); std::iota(inds.begin(), inds.end(), 0);
    std::transform(inds.begin(), inds.end(), epsilon_second_moment.begin(),
                   [sign, omega, epsilons, a2f_over_transDOS] (size_t i) { 
                    return (2.0 * epsilons[i] + sign * omega) * (2.0 * epsilons[i] + sign * omega) * a2f_over_transDOS[i];
    });

    auto factor = fermi_dirac_factor(sign, Te, omega, epsilons);
    inner_integ = de * std::inner_product(factor.begin(), factor.end(), epsilon_second_moment.begin(), 0.0);
    return inner_integ;
}

// Q_{x1,x1}, s = s' = -1 // actually is the same as calc_inner_integral_elec_elastic
double calc_inner_integral_therm_elastic_ord2_m(int sign, double Te, double omega, 
                                                const array1D& epsilons, const array1D& transDOSes, const array1D& a2f)
{
    // TODO assert that sign == -1 in a2f file
    assert(epsilons.size() > 1);
    double inner_integ{ 0.0 };
    double de = epsilons[1] - epsilons[0];
    auto neps = epsilons.size();

    auto a2f_over_transDOS = calc_a2f_over_transDOS(epsilons, transDOSes, a2f);

    array1D epsilon_second_moment(neps, 0.0);
    std::vector<int> inds(neps, 0); std::iota(inds.begin(), inds.end(), 0);
    std::transform(inds.begin(), inds.end(), epsilon_second_moment.begin(),
                   [omega, a2f_over_transDOS] (size_t i) { 
                    return omega * omega * a2f_over_transDOS[i];
    });

    auto factor = fermi_dirac_factor(sign, Te, omega, epsilons);
    inner_integ = de * std::inner_product(factor.begin(), factor.end(), epsilon_second_moment.begin(), 0.0);
    return inner_integ;
}

void calc_therm_cond_elastic(const array1D& Temps, const char* a2f_plus_fnm,
                             const char* a2f_minus_fnm, const char* cond_fnm, double unit_cell_vol)
{
    // unit cell volumes must be given in [bohr^3]
    // each prefactor must be multiplied by eV_in_Ry which comes from integration over epsilons and omegas and from the order ord in Q matrix
    double prefactor = 2.0 * pi * unit_cell_vol * eV_in_Ry;
    auto Qx0x0 = prefactor * calc_external_integral(Temps, 0, calc_inner_integral_elec_elastic, a2f_plus_fnm);

    prefactor = std::sqrt(3.0) * unit_cell_vol * eV_in_Ry;
    auto Qx0x1 = prefactor * calc_external_integral(Temps, 1, calc_inner_integral_therm_elastic_ord1, a2f_plus_fnm);

    prefactor = 6.0 / pi * unit_cell_vol * eV_in_Ry;
    prefactor *= 0.25; // from summation over s, s' in (32)
    auto Qx1x1 = prefactor * (calc_external_integral(Temps, 2, calc_inner_integral_therm_elastic_ord2_p, a2f_plus_fnm)
                            + calc_external_integral(Temps, 2, calc_inner_integral_therm_elastic_ord2_m, a2f_minus_fnm));

    array1D kappas;
    for (size_t i = 0; i < Temps.size(); ++i)
    {
        auto TeV = Temps[i] / 11606.0;
        array2D Q = { { Qx0x0[i], Qx0x1[i] }, { Qx0x1[i], Qx1x1[i] } };
        auto invQ = calc_inv_2d( Q );

        prefactor = (2.0 * pi * pi / 3.0) * TeV * eV_in_Ry * kB_as_Ry_over_K; // in r.a.u.
        auto kappa = prefactor * (invQ[1][1] - (invQ[0][1] * invQ[0][1]) / invQ[0][0]); // still in r.a.u.
        kappas.push_back(kappa / rau_in_cm_over_W); // go to W / (cm K)
    }

    OHandler ohandler(a2f_plus_fnm, cond_fnm, ResistType::Thermal);
    ohandler.dump(Temps, kappas);
}

void calc_seebeck_elastic(const arrays::array1D& Temps, const char* a2f_plus_fnm,
                          const char* a2f_minus_fnm, const char* cond_fnm, double unit_cell_vol)
{
    // unit cell volumes must be given in [bohr^3]
    // each prefactor must be multiplied by eV_in_Ry which comes from integration over epsilons and omegas and from the order ord in Q matrix
    double prefactor = 2.0 * pi * unit_cell_vol * eV_in_Ry;
    auto Qx0x0 = prefactor * calc_external_integral(Temps, 0, calc_inner_integral_elec_elastic, a2f_plus_fnm);

    prefactor = std::sqrt(3.0) * unit_cell_vol * eV_in_Ry;
    auto Qx0x1 = prefactor * calc_external_integral(Temps, 1, calc_inner_integral_therm_elastic_ord1, a2f_plus_fnm);

    prefactor = 6.0 / pi * unit_cell_vol * eV_in_Ry;
    prefactor *= 0.25; // from summation over s, s' in (32)
    auto Qx1x1 = prefactor * (calc_external_integral(Temps, 2, calc_inner_integral_therm_elastic_ord2_p, a2f_plus_fnm)
                            + calc_external_integral(Temps, 2, calc_inner_integral_therm_elastic_ord2_m, a2f_minus_fnm));

    array1D seebeck;
    double e2 = 2.0; // in [r.a.u.]
    prefactor = (pi / std::sqrt(3.0 * e2)) * kB_as_Ry_over_K; // in r.a.u.
    for (size_t i = 0; i < Temps.size(); ++i)
    {
        auto TeV = Temps[i] / 11606.0;
        array2D Q = { { Qx0x0[i], Qx0x1[i] }, { Qx0x1[i], Qx1x1[i] } };
        auto invQ = calc_inv_2d( Q );

        auto s = prefactor * invQ[0][1] / invQ[0][0]; // still in r.a.u.
        seebeck.push_back(s * rau_in_muV); // go to (muV / K)
    }

    OHandler ohandler(a2f_plus_fnm, cond_fnm, ResistType::Seebeck);
    ohandler.dump(Temps, seebeck);
}

} // transport
} // skies
