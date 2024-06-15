#include <skies/quantities/dos.h>

namespace skies { namespace quantities {

using namespace arrays;
using namespace bzsampling;

double evaluate_trdos_at_value(const array2D& grid,
                               double value,
                               double smearing,
                               SamplingFunc sampling,
                               char cart)
{
    auto nkpt = grid.size();
    array2D weights(nkpt, array1D(EigenValue::nbands));
    std::transform(grid.begin(), grid.end(), weights.begin(),
                    [cart] (const array1D& k)
                    {
                        array1D v = Velocities(cart).interpolate_at(k);
                        std::transform(v.begin(), v.end(), v.begin(), [] (double d) { return d*d; });
                        return v;
                    });
    return evaluate_dos_at_value<EigenValue>(grid, value, smearing, sampling, weights);
}

double evaluate_smeared_trdos_at_value(const array2D& grid,
                                       double value,
                                       double smearing,
                                       SamplingFunc sampling,
                                       double elec_temp,
                                       char cart)
{
    array1D trDOSes_tmp(32, 0.0);
    auto epsilons = create_range(value - 1, value + 1, 30);
    std::transform(epsilons.begin(), epsilons.end(), trDOSes_tmp.begin(),
                    [&] (double e) { return evaluate_trdos_at_value(grid, e, smearing, sampling, cart); });
    return smear_with_fd(trDOSes_tmp, epsilons, value, elec_temp);
}

} // quantities
} // skies
