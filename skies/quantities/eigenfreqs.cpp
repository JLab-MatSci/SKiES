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
