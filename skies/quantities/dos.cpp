/*-----------------------------------------------------------------------
    * SKiES - Solver of Kinetic Equation for Solids
    * (C) 2024 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)
    *
    * SKiES may only be utilized for non-profit research.
    * Citing appropriate sources is required when using SKiES.
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */
#include <skies/quantities/dos.h>
#include <skies/utils/tbb_wrapper.h>

namespace skies { namespace quantities {

using namespace arrays;
using namespace bzsampling;

double evaluate_trdos_at_value(double value,
                               double smearing,
                               SamplingFunc sampling,
                               const array2D& energies,
                               const array2D& velocs_squared)
{

    return evaluate_dos_at_value<EigenValue>(value, smearing, sampling, energies, velocs_squared);
}

double evaluate_smeared_trdos_at_value(double value,
                                       double smearing,
                                       SamplingFunc sampling,
                                       const array2D& energies,
                                       const array2D& velocs_squared,
                                       double elec_temp)
{
    array1D trDOSes_tmp(32, 0.0);
    double de = 2.0; // TODO: TRY TO FIX IT
    auto epsilons = create_range(value - de, value + de, 30);
    std::transform(epsilons.begin(), epsilons.end(), trDOSes_tmp.begin(),
                    [&] (double e) { return evaluate_trdos_at_value(e, smearing, sampling, energies, velocs_squared); });
    return smear_with_fd(trDOSes_tmp, epsilons, value, elec_temp);
}

void evaluate_trdos(const arrays::array2D& grid,
                    const array1D& range,
                    char alpha,
                    double smearing,
                    bzsampling::SamplingFunc sampl_type)
{
    std::ofstream os("VelocitiesDOS.dat");
    os << std::right;
    os << std::setw(12) << "# Energy, eV";
    os << std::setw(25) << "Transport DOS";
    os << std::setw(24) << " [r.a.u.]";
    os << std::endl;

    auto nkpt = grid.size();
    array2D eigenens(nkpt, array1D(EigenValue::nbands));
    std::transform(PAR grid.begin(), grid.end(), eigenens.begin(), [] (auto&& k) {
        return EigenValue::interpolate_at(k);
    });

    auto prepare_velocs = [&] (char cart, array2D& elvelocs, array2D& elvelocs_sq) {
        elvelocs.resize(nkpt, array1D(EigenValue::nbands, 0.0));
        elvelocs_sq.resize(nkpt, array1D(EigenValue::nbands, 0.0));
        std::transform(PAR grid.begin(), grid.end(), elvelocs.begin(),
                    [&] (auto&& k) { return Velocities(cart).interpolate_at(k); });
        std::transform(PAR elvelocs.begin(), elvelocs.end(), elvelocs_sq.begin(),
                [] (auto&& v) {
                    auto squared_v = array1D(EigenValue::nbands, 0.0);
                    std::transform(v.begin(), v.end(), squared_v.begin(), [] (auto&& x) { return x * x; });
                    return squared_v;
        });
    };

    array2D elvelocs_alpha, elvelocs_alpha_sq;
    prepare_velocs(alpha, elvelocs_alpha, elvelocs_alpha_sq);

    array1D trDOSes(range.size());
    std::transform(range.begin(), range.end(), trDOSes.begin(), [&] (auto&& eps) {
        return evaluate_trdos_at_value(eps, smearing, sampl_type, eigenens, elvelocs_alpha_sq);
    });
    trDOSes = trDOSes * units::Ry_in_eV; // go to [r.a.u.]

    for (size_t i = 0; i < range.size(); ++i)
        os << std::setprecision(6) << std::setw(12) << range[i] << std::setw(49) << trDOSes[i] << std::endl;
    os.close();
}

} // quantities
} // skies