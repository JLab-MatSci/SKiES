#include <cassert>
#include <numeric>
#include <functional>

#include <skies/common/units.h>
#include <skies/sampling/sampling.h>
#include <skies/transport/elastic.h>
#include <skies/transport/iohandler.h>

#include <skies/quantities/eigenvals.h>

#include <skies/utils/tbb_wrapper.h>

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
        [sign, Te, omega] (double e) {
            return fermi_dirac(e, Te) * (1 - fermi_dirac(e + sign * omega, Te));
    });
    return factor;
}

array1D bose_einstein_factor(double Ti, const array1D& omegas, int add)
{
    assert(omegas.size() > 1);
    array1D factor(omegas.size(), 0.0);
    std::transform(omegas.begin(), omegas.end(), factor.begin(),
        [Ti, add] (double om) {
            return bose_einstein(om, Ti) + add;
    });
    return factor;
}

array1D calc_a2f_over_transDOS(const array1D& epsilons, const array1D& transDOSes, const array1D& a2f)
{
    array1D a2f_over_transDOS(epsilons.size(), 0.0);
    std::vector<int> inds(epsilons.size(), 0); std::iota(inds.begin(), inds.end(), 0);
    std::transform(inds.begin(), inds.end(), a2f_over_transDOS.begin(),
        [a2f, transDOSes] (size_t i) {
            return a2f[i] / (transDOSes[i] * units::Ry_in_eV);
    });
    return a2f_over_transDOS;
}

array2D calc_external_integral(const array1D& Temps, const array1D& ion_Temps, int ord, innerIntegCalculator calc_inner_integ, const std::string& a2f_fnm)
{
    IHandler ihandler(a2f_fnm.c_str());
    if (ihandler.epsilons().size() < 2)
        throw std::runtime_error("There must be at least two values of electron energies in low T formulas.");
    if (ihandler.omegas().size() < 2)
        throw std::runtime_error("There must be at least two values of phonon frequencies.");
    assert(ord == 0 || ord == 1 || ord == 2);

    auto omegas = ihandler.omegas();
    double dom  = omegas[1] - omegas[0];

    auto epsilons = ihandler.epsilons();
    auto transDOSes = ihandler.transDOSes();
    auto a2f = ihandler.a2f();

    int s       = ihandler.sign();
    int s_prime = ihandler.sign_pr();

    array2D ext_integs;
    for (size_t i = 0; i < Temps.size(); ++i)
    {
        double Temp = Temps[i];
        double TeV = Temp / 11606.0;
        array1D epsilon_integs_plus, epsilon_integs_minus;

        for(size_t iom = 0; iom < omegas.size(); ++iom)
        {
            epsilon_integs_plus.push_back(calc_inner_integ(  1, TeV, omegas[iom], epsilons, transDOSes, a2f[iom], s, s_prime));
            epsilon_integs_minus.push_back(calc_inner_integ(-1, TeV, omegas[iom], epsilons, transDOSes, a2f[iom], s, s_prime));
        }

        array1D elec_temp_line;
        if (ion_Temps.empty())
        {
            auto bose_factor_plus  = bose_einstein_factor(TeV, omegas, 0);
            auto bose_factor_minus = bose_einstein_factor(TeV, omegas, 1);
            double omega_integ_plus  = dom * std::inner_product(bose_factor_plus.begin(), bose_factor_plus.end(), epsilon_integs_plus.begin(), 0.0);
            double omega_integ_minus = dom * std::inner_product(bose_factor_minus.begin(), bose_factor_minus.end(), epsilon_integs_minus.begin(), 0.0);

            elec_temp_line.push_back((omega_integ_plus + omega_integ_minus) / std::pow(TeV, ord + 1));
        }
        else
        {
            for (size_t j = 0; j < ion_Temps.size(); ++j)
            {
                double ion_Temp = ion_Temps[j];
                double ion_TeV = ion_Temp / 11606.0;

                auto bose_factor_plus  = bose_einstein_factor(ion_TeV, omegas, 0);
                auto bose_factor_minus = bose_einstein_factor(ion_TeV, omegas, 1);
                double omega_integ_plus  = dom * std::inner_product(bose_factor_plus.begin(), bose_factor_plus.end(), epsilon_integs_plus.begin(), 0.0);
                double omega_integ_minus = dom * std::inner_product(bose_factor_minus.begin(), bose_factor_minus.end(), epsilon_integs_minus.begin(), 0.0);

                elec_temp_line.push_back((omega_integ_plus + omega_integ_minus) / std::pow(TeV, ord + 1));
            }
        }
        ext_integs.push_back(elec_temp_line);
    }
    return ext_integs;
}

// Q_{x0,x0}
double calc_inner_integral_elec_elastic(int sign, double Te, double omega, const array1D& epsilons,
                                   const array1D& transDOSes, const array1D& a2f,
                                   int s, int s_prime) // needed by signature of innerIntegCalculator
{
    assert((s == 1) && (s_prime == 1));
    assert(epsilons.size() > 1);
    double inner_integ{ 0.0 };
    double de = epsilons[1] - epsilons[0];

    auto a2f_over_transDOS = calc_a2f_over_transDOS(epsilons, transDOSes, a2f);
    auto factor = fermi_dirac_factor(sign, Te, omega, epsilons);
    inner_integ = de * std::inner_product(factor.begin(), factor.end(), a2f_over_transDOS.begin(), 0.0);
    return inner_integ;
}

void calc_elec_cond_elastic(const array1D& Temps, const array1D& ion_Temps, const char* a2f_fnm, const char* cond_fnm, double unit_cell_vol)
{
    // unit cell volumes must be given in [bohr^3]
    double prefactor = 2.0 * pi * unit_cell_vol * rau_in_muOm_cm / eV_in_Ry;
    // multiply by ev_in_Ry^2 due to integration over de and dOm
    prefactor *= eV_in_Ry * eV_in_Ry;

    auto Qx0x0  = prefactor * calc_external_integral(Temps, ion_Temps, 0, calc_inner_integral_elec_elastic, a2f_fnm);
    double e2 = 2.0;
    auto resist = Qx0x0 * (0.5 / e2); // e2 = 2.0 in r.a.u.

    OHandler ohandler(a2f_fnm, cond_fnm, ResistType::Electrical, ion_Temps);
    ohandler.dump(Temps, resist);
}

// Q_{x0,x1}
double calc_inner_integral_therm_elastic_ord1(int sign, double Te, double omega, 
                                              const array1D& epsilons, const array1D& transDOSes,
                                              const array1D& a2f, int s, int s_prime)
{
    assert(s == 1); // the odd correction for Q_01 is only from s_prime = -1, but s == 1 any way
    assert(epsilons.size() > 1);
    double inner_integ{ 0.0 };
    double de = epsilons[1] - epsilons[0];
    auto neps = epsilons.size();

    auto a2f_over_transDOS = calc_a2f_over_transDOS(epsilons, transDOSes, a2f);

    array1D epsilon_first_moment(neps, 0.0);
    std::vector<int> inds(neps, 0); std::iota(inds.begin(), inds.end(), 0);
    std::transform(inds.begin(), inds.end(), epsilon_first_moment.begin(), [&] (size_t i) {
        return ((epsilons[i] + s_prime * epsilons[i]) + sign * omega) * a2f_over_transDOS[i];
    });

    auto factor = fermi_dirac_factor(sign, Te, omega, epsilons);
    inner_integ = de * std::inner_product(factor.begin(), factor.end(), epsilon_first_moment.begin(), 0.0);
    return inner_integ;
}

// Q_{x1,x1}
double calc_inner_integral_therm_elastic_ord2(int sign, double Te, double omega, 
                                              const array1D& epsilons, const array1D& transDOSes, const array1D& a2f,
                                              int s, int s_prime)
{
    // TODO assert that sign == +1 in a2f file
    assert(epsilons.size() > 1);
    double inner_integ{ 0.0 };
    double de = epsilons[1] - epsilons[0];
    auto neps = epsilons.size();

    auto a2f_over_transDOS = calc_a2f_over_transDOS(epsilons, transDOSes, a2f);

    array1D epsilon_second_moment(neps, 0.0);
    std::vector<int> inds(neps, 0); std::iota(inds.begin(), inds.end(), 0);
    std::transform(inds.begin(), inds.end(), epsilon_second_moment.begin(), [&] (size_t i) {
        return ((epsilons[i] + s       * epsilons[i]) + sign * omega)
             * ((epsilons[i] + s_prime * epsilons[i]) + sign * omega) * a2f_over_transDOS[i];
    });

    auto factor = fermi_dirac_factor(sign, Te, omega, epsilons);
    inner_integ = de * std::inner_product(factor.begin(), factor.end(), epsilon_second_moment.begin(), 0.0);
    return inner_integ;
}

void calc_therm_cond_elastic(const array1D& Temps, const array1D& ion_Temps, const char* a2f_plus_fnm,
                             const char* a2f_minus_fnm, const char* cond_fnm, double unit_cell_vol,
                             const std::string& a2f_pm_fnm, const std::string& a2f_mp_fnm)
{
    // unit cell volumes must be given in [bohr^3]
    // each prefactor must be multiplied by eV_in_Ry which comes from integration over epsilons and omegas and from the order ord in Q matrix
    double prefactor = 2.0 * pi * unit_cell_vol * eV_in_Ry;
    auto Qx0x0 = prefactor * calc_external_integral(Temps, ion_Temps, 0, calc_inner_integral_elec_elastic, a2f_plus_fnm);

    prefactor = std::sqrt(3.0) * unit_cell_vol * eV_in_Ry;
    auto Qx0x1 = prefactor * calc_external_integral(Temps, ion_Temps, 1, calc_inner_integral_therm_elastic_ord1, a2f_plus_fnm);
    if (a2f_pm_fnm.size())
        Qx0x1 = Qx0x1
               + prefactor * calc_external_integral(Temps, ion_Temps, 1, calc_inner_integral_therm_elastic_ord1, a2f_pm_fnm);

    prefactor = 6.0 / pi * unit_cell_vol * eV_in_Ry;
    prefactor *= 0.25; // from summation over s, s' in (32)
    auto Qx1x1 = prefactor * (calc_external_integral(Temps, ion_Temps, 2, calc_inner_integral_therm_elastic_ord2, a2f_plus_fnm)
                            + calc_external_integral(Temps, ion_Temps, 2, calc_inner_integral_therm_elastic_ord2, a2f_minus_fnm));
    if (a2f_pm_fnm.size())
        Qx1x1 = Qx1x1
               + prefactor * (calc_external_integral(Temps, ion_Temps, 2, calc_inner_integral_therm_elastic_ord2, a2f_pm_fnm)
                            + calc_external_integral(Temps, ion_Temps, 2, calc_inner_integral_therm_elastic_ord2, a2f_mp_fnm));
    array2D kappas;
    for (size_t i = 0; i < Temps.size(); ++i)
    {
        auto TeV = Temps[i] / 11606.0;
        array1D elec_temp_line;
        if (ion_Temps.empty())
        {
            array2D Q = { { Qx0x0[i][0], Qx0x1[i][0] }, { Qx0x1[i][0], Qx1x1[i][0] } };
            auto invQ = calc_inv_2d( Q );

            prefactor = (2.0 * pi * pi / 3.0) * TeV * eV_in_Ry * kB_as_Ry_over_K; // in r.a.u.
            auto kappa = prefactor * (invQ[1][1] - (invQ[0][1] * invQ[0][1]) / invQ[0][0]); // still in r.a.u.

            elec_temp_line.push_back(kappa / rau_in_cm_over_W);
        }
        else
        {
            for (size_t j = 0; j < ion_Temps.size(); ++j)
            {
                array2D Q = { { Qx0x0[i][j], Qx0x1[i][j] }, { Qx0x1[i][j], Qx1x1[i][j] } };
                auto invQ = calc_inv_2d( Q );

                prefactor = (2.0 * pi * pi / 3.0) * TeV * eV_in_Ry * kB_as_Ry_over_K; // in r.a.u.
                auto kappa = prefactor * (invQ[1][1] - (invQ[0][1] * invQ[0][1]) / invQ[0][0]); // still in r.a.u.
                elec_temp_line.push_back(kappa / rau_in_cm_over_W); // go to W / (cm K)
            }
        }
        kappas.push_back(elec_temp_line);
    }
    OHandler ohandler(a2f_plus_fnm, cond_fnm, ResistType::Thermal, ion_Temps);
    ohandler.dump(Temps, kappas);
}

array2D evaluate_Qa0b0(char alpha, char beta, double unit_cell_vol, const array1D& Temps, const array1D& ion_Temps);
array2D evaluate_Qa0b1(char alpha, char beta, double unit_cell_vol, const array1D& Temps, const array1D& ion_Temps);
array2D evaluate_Qa1b1(char alpha, char beta, double unit_cell_vol, const array1D& Temps, const array1D& ion_Temps);

void calc_seebeck_elastic(char alpha,
 			  const arrays::array1D& Temps, const arrays::array1D& ion_Temps,
			  const char* cond_fnm, double unit_cell_vol)
{
    // unit cell volumes must be given in [bohr^3]
    // each prefactor must be multiplied by eV_in_Ry which comes from integration over epsilons and omegas and from the order ord in Q matrix
    auto Qa0x0 = evaluate_Qa0b0(alpha, 'x', unit_cell_vol, Temps, ion_Temps);
    auto Qa0x1 = evaluate_Qa0b1(alpha, 'x', unit_cell_vol, Temps, ion_Temps);
    auto Qa1x1 = evaluate_Qa1b1(alpha, 'x', unit_cell_vol, Temps, ion_Temps);

    auto Qa0y0 = evaluate_Qa0b0(alpha, 'y', unit_cell_vol, Temps, ion_Temps);
    auto Qa0y1 = evaluate_Qa0b1(alpha, 'y', unit_cell_vol, Temps, ion_Temps);
    auto Qa1y1 = evaluate_Qa1b1(alpha, 'y', unit_cell_vol, Temps, ion_Temps);

    auto Qa0z0 = evaluate_Qa0b0(alpha, 'z', unit_cell_vol, Temps, ion_Temps);
    auto Qa0z1 = evaluate_Qa0b1(alpha, 'z', unit_cell_vol, Temps, ion_Temps);
    auto Qa1z1 = evaluate_Qa1b1(alpha, 'z', unit_cell_vol, Temps, ion_Temps);


    array2D seebeck;
    double e2 = 2.0; // in [r.a.u.]
    double prefactor = (pi / std::sqrt(3.0 * e2)) * kB_as_Ry_over_K; // in r.a.u.
    for (size_t i = 0; i < Temps.size(); ++i)
    {
        array1D elec_temp_line;
        if (ion_Temps.empty())
        {
            array2D Qax = { { Qa0x0[i][0], Qa0x1[i][0] }, { Qa0x1[i][0], Qa1x1[i][0] } };
            auto invQax = calc_inv_2d( Qax );
            auto Sax = prefactor * Qax[0][1] * invQax[1][1]; // still in r.a.u.

	    array2D Qay = { { Qa0y0[i][0], Qa0y1[i][0] }, { Qa0y1[i][0], Qa1y1[i][0] } };
            auto invQay = calc_inv_2d( Qay );
            auto Say = prefactor * Qay[0][1] * invQay[1][1]; // still in r.a.u.

	    array2D Qaz = { { Qa0z0[i][0], Qa0z1[i][0] }, { Qa0z1[i][0], Qa1z1[i][0] } };
            auto invQaz = calc_inv_2d( Qaz );
            auto Saz = prefactor * Qaz[0][1] * invQaz[1][1]; // still in r.a.u.

            elec_temp_line.push_back((Sax + Say + Saz) * rau_in_muV); // go to (muV / K)
        }
        else
        {
            for (size_t j = 0; j < ion_Temps.size(); ++j)
            {
		array2D Qax = { { Qa0x0[i][j], Qa0x1[i][j] }, { Qa0x1[i][j], Qa1x1[i][j] } };
            	auto invQax = calc_inv_2d( Qax );
            	auto Sax = prefactor * Qax[0][1] * invQax[1][1]; // still in r.a.u.

            	array2D Qay = { { Qa0y0[i][j], Qa0y1[i][j] }, { Qa0y1[i][j], Qa1y1[i][j] } };
            	auto invQay = calc_inv_2d( Qay );
            	auto Say = prefactor * Qay[0][1] * invQay[1][1]; // still in r.a.u.

            	array2D Qaz = { { Qa0z0[i][j], Qa0z1[i][j] }, { Qa0z1[i][j], Qa1z1[i][j] } };
            	auto invQaz = calc_inv_2d( Qaz );
            	auto Saz = prefactor * Qaz[0][1] * invQaz[1][1]; // still in r.a.u.

            	elec_temp_line.push_back((Sax + Say + Saz) * rau_in_muV); // go to (muV / K)
            }
        }
        seebeck.push_back(elec_temp_line);
    }

    std::string a2f_pp_aa_fnm = "SpecFunc_pp_";
    a2f_pp_aa_fnm += std::string{alpha};
    a2f_pp_aa_fnm += std::string{alpha};
    a2f_pp_aa_fnm += std::string(".dat");

    OHandler ohandler(a2f_pp_aa_fnm.c_str(), cond_fnm, ResistType::Seebeck, ion_Temps);
    ohandler.dump(Temps, seebeck);
}

array2D evaluate_Qa0b0(char alpha, char beta, double unit_cell_vol, const array1D& Temps, const array1D& ion_Temps)
{
    std::string a2f_pp_ab_fnm = "SpecFunc_pp_";
    a2f_pp_ab_fnm += std::string{alpha};
    a2f_pp_ab_fnm += std::string{beta};
    a2f_pp_ab_fnm += std::string(".dat");

    double prefactor = 2.0 * pi * unit_cell_vol * eV_in_Ry;
    auto Qa0b0 = prefactor * calc_external_integral(Temps, ion_Temps, 0, calc_inner_integral_elec_elastic, a2f_pp_ab_fnm);
    return Qa0b0;
}

array2D evaluate_Qa0b1(char alpha, char beta, double unit_cell_vol, const array1D& Temps, const array1D& ion_Temps)
{
    std::string a2f_pp_ab_fnm = "SpecFunc_pp_";
    a2f_pp_ab_fnm += std::string{alpha};
    a2f_pp_ab_fnm += std::string{beta};
    a2f_pp_ab_fnm += std::string(".dat");

    double prefactor = std::sqrt(3.0) * unit_cell_vol * eV_in_Ry;
    auto Qa0b1 = prefactor * calc_external_integral(Temps, ion_Temps, 1, calc_inner_integral_therm_elastic_ord1, a2f_pp_ab_fnm);
    /*if (a2f_pm_fnm.size())
    {
        std::cout << a2f_pm_fnm << std::endl;
        Qx0x1 = Qx0x1
               + prefactor * calc_external_integral(Temps, ion_Temps, 1, calc_inner_integral_therm_elastic_ord1, a2f_pm_fnm);
    }*/

    return Qa0b1;

}

array2D evaluate_Qa1b1(char alpha, char beta, double unit_cell_vol, const array1D& Temps, const array1D& ion_Temps)
{
    std::string a2f_pp_ab_fnm = "SpecFunc_pp_";
    a2f_pp_ab_fnm += std::string{alpha};
    a2f_pp_ab_fnm += std::string{beta};
    a2f_pp_ab_fnm += std::string(".dat");

    std::string a2f_mm_ab_fnm = "SpecFunc_mm_";
    a2f_mm_ab_fnm += std::string{alpha};
    a2f_mm_ab_fnm += std::string{beta};
    a2f_mm_ab_fnm += std::string(".dat");

    double prefactor = 6.0 / pi * unit_cell_vol * eV_in_Ry;
    prefactor *= 0.25; // from summation over s, s' in (32)
    auto Qa1b1 = prefactor * (calc_external_integral(Temps, ion_Temps, 2, calc_inner_integral_therm_elastic_ord2, a2f_pp_ab_fnm)
                            + calc_external_integral(Temps, ion_Temps, 2, calc_inner_integral_therm_elastic_ord2, a2f_mm_ab_fnm));
    /*if (a2f_pm_fnm.size())
        Qx1x1 = Qx1x1
               + prefactor * (calc_external_integral(Temps, ion_Temps, 2, calc_inner_integral_therm_elastic_ord2, a2f_pm_fnm)
                            + calc_external_integral(Temps, ion_Temps, 2, calc_inner_integral_therm_elastic_ord2, a2f_mp_fnm));
    */
    return Qa1b1;
}

} // transport
} // skies
