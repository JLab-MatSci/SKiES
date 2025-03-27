/*-----------------------------------------------------------------------
    * SKiES - Solver of Kinetic Equation for Solids
    * 
    * (C) 2025 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)
    *
    * SKiES may only be utilized for non-profit research.
    * Citing appropriate sources is required when using SKiES.
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */
#include <cassert>

#include <skies/quantities/eigenfreqs.h>

using namespace skies::arrays;

size_t skies::quantities::EigenFrequency::nmodes;
size_t skies::quantities::EigenFrequencyDrawable::nmodes;

#ifdef __cplusplus
extern "C" {
#endif

void fillNmodes(int* nmodes) {
    skies::quantities::EigenFrequencyDrawable::nmodes = static_cast<size_t>(*nmodes);
    skies::quantities::EigenFrequency::nmodes = static_cast<size_t>(*nmodes);
}

#ifdef __cplusplus
}
#endif

namespace skies { namespace quantities {

array1D EigenFrequencyDrawable::interpolate_at(const array1D& q) const
{
    assert(q.size() == 3);
    double rfreqs[nmodes];
    interpEigenFreq1DAt(&q[0], &q[1], &q[2], rfreqs);
    auto freqs = array1D(rfreqs, rfreqs + nmodes) * units::Ry_in_eV;
    return freqs;
}

array1D EigenFrequency::interpolate_at(const array1D& q)
{
    assert(q.size() == 3);
    double rfreqs[nmodes];
    interpEigenFreqAt(&q[0], &q[1], &q[2], rfreqs);
    auto freqs = array1D(rfreqs, rfreqs + nmodes) * units::Ry_in_eV;
    return freqs;
}

} // quantities
} // skies
