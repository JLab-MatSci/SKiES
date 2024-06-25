#include <skies/quantities/dos.h>

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
                    bzsampling::SamplingFunc sampl_type,
                    char cart)
{
    std::ofstream os("VelocitiesDOS.dat");
    os << std::right;
    os << std::setw(12) << "# Energy, eV";
    os << std::setw(25) << "Transport DOS";
    os << std::setw(24) << " [13.605685 * Ry bohr^2]";
    os << std::endl;

    auto nkpt = grid.size();
    array2D eigenvals(nkpt, array1D(EigenValue::nbands));
    array2D velocs(nkpt, array1D(EigenValue::nbands));
    array2D velocs_squared(nkpt, array1D(EigenValue::nbands));

    std::transform(grid.begin(), grid.end(), eigenvals.begin(),
                [] (auto&& k) { return EigenValue::interpolate_at(k); });
    std::transform(grid.begin(), grid.end(), velocs.begin(),
                [cart] (auto&& k) { return VelocitiesDrawable(cart).interpolate_at(k); });
    std::transform(velocs.begin(), velocs.end(), velocs_squared.begin(),
                    [] (const array1D& v) {
                        auto squared_v = array1D(EigenValue::nbands, 0.0);
                        std::transform(v.begin(), v.end(), squared_v.begin(), [] (auto x) { return x * x; });
                        return squared_v;
    });

    int rank{ 0 };
    array1D trDOSes(range.size(), 0.0);
    array1D trDOSes_tmp(range.size(), 0.0);
#ifdef SKIES_MPI
    auto rcounts = mpi::prepare_rcounts_displs(range.size()).first;
    auto displs  = mpi::prepare_rcounts_displs(range.size()).second;

    rank = skies::mpi::rank();
    for (int i = displs[rank]; i < displs[rank] + rcounts[rank]; ++i)
        trDOSes_tmp[i] = evaluate_trdos_at_value(range[i], smearing, sampl_type, eigenvals, velocs_squared);

    MPI_Reduce(trDOSes_tmp.data(), trDOSes.data(), range.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(trDOSes.data(), range.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    std::transform(range.begin(), range.end(), trDOSes.begin(), [] (double v) {
        return evaluate_trdos_at_value(v, smearing, sampl_type, eigenvals, velocs_squared);
    });
#endif
    if (!rank)
    {
        for (size_t i = 0; i < range.size(); ++i)
            os << std::setprecision(6) << std::setw(12) << range[i] << std::setw(49) << trDOSes[i] << std::endl;
        os.close();
    }
    return;
}

} // quantities
} // skies