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
                    double smearing,
                    bzsampling::SamplingFunc sampl_type)
{
    std::ofstream os("VelocitiesDOS.dat");
    os << std::right;
    os << std::setw(12) << "# Energy, eV";
    os << std::setw(25) << "Transport DOS";
    os << std::setw(24) << " [13.605685 * Ry bohr^2]";
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

    array2D elvelocs_x, elvelocs_y, elvelocs_z;
    array2D elvelocs_x_sq, elvelocs_y_sq, elvelocs_z_sq;
    prepare_velocs('x', elvelocs_x, elvelocs_x_sq);
    prepare_velocs('y', elvelocs_y, elvelocs_y_sq);
    prepare_velocs('z', elvelocs_z, elvelocs_z_sq);

    array1D trDOSes(range.size());
    array1D trDOSes_x(range.size()), trDOSes_y(range.size()), trDOSes_z(range.size());

    std::transform(range.begin(), range.end(), trDOSes_x.begin(), [&] (auto&& eps) {
        return evaluate_trdos_at_value(eps, smearing, sampl_type, eigenens, elvelocs_x_sq);
    });
    std::transform(range.begin(), range.end(), trDOSes_y.begin(), [&] (auto&& eps) {
        return evaluate_trdos_at_value(eps, smearing, sampl_type, eigenens, elvelocs_y_sq);
    });
    std::transform(range.begin(), range.end(), trDOSes_z.begin(), [&] (auto&& eps) {
        return evaluate_trdos_at_value(eps, smearing, sampl_type, eigenens, elvelocs_z_sq);
    });
    trDOSes = (trDOSes_x + trDOSes_y + trDOSes_z) * (1.0 / 3.0);
    trDOSes = trDOSes * units::Ry_in_eV; // go to [r.a.u.]

    int rank{ 0 };
#ifdef SKIES_MPI
    rank = skies::mpi::rank();
#endif
    if (!rank)
    {
        for (size_t i = 0; i < range.size(); ++i)
            os << std::setprecision(6) << std::setw(12) << range[i] << std::setw(49) << trDOSes[i] << std::endl;
        os.close();
    }
}

} // quantities
} // skies