#include <numeric>

#include <skies/common/alg.h>
#include <skies/quantities/eigenvals.h>
#include <skies/lattices/kp_protocol.h>
#include <skies/utils/tbb_wrapper.h>

#include <iostream>

using namespace skies::arrays;
using namespace skies::bzsampling;

size_t skies::quantities::EigenValueDrawable::nbands;
size_t skies::quantities::EigenValue::nbands;
int skies::quantities::EigenValue::nelec;
double skies::quantities::EigenValue::eF;

#ifdef __cplusplus
extern "C" {
#endif

void fillNbands(int* nbands) {
    skies::quantities::EigenValueDrawable::nbands = static_cast<size_t>(*nbands);
    skies::quantities::EigenValue::nbands = static_cast<size_t>(*nbands);
}

void fillNelec(double* nelec) {
    skies::quantities::EigenValue::nelec = *nelec;
}

#ifdef __cplusplus
}
#endif

namespace skies { namespace quantities {

array1D EigenValueDrawable::interpolate_at(const array1D& k) const
{
    assert(k.size() == 3);
    double reigs[nbands];
    interpEigenValueAt(&k[0], &k[1], &k[2], reigs);
    auto eigs = array1D(reigs, reigs + nbands) * units::Ry_in_eV;
    eigs = eigs - array1D(nbands, EigenValue::eF); // return energies relative to Fermi level
    return eigs;
}

array1D EigenValue::interpolate_at(const arrays::array1D& k)
{
    assert(k.size() == 3);
    double reigs[nbands];
    interpEigenValueAt(&k[0], &k[1], &k[2], reigs);
    auto eigs = array1D(reigs, reigs + nbands) * units::Ry_in_eV;
    eigs = eigs - array1D(nbands, EigenValue::eF); // return energies relative to Fermi level
    return eigs;
}

void EigenValue::find_eF(double TeV, double crit)
{
    // 1. fill eigenvalues on some large grid
    array2D eigenens2D;
    size_t nk1 = 50;
    size_t nk2 = 50;
    size_t nk3 = 50;
    skies::KPprotocol kprot{nk1, nk2, nk3};
    auto kpts = kprot.grid();
    eigenens2D.resize(kpts.size(), array1D(EigenValue::nbands));
    std::transform(PAR kpts.begin(), kpts.end(), eigenens2D.begin(), [] (auto&& k) {
        return EigenValue::interpolate_at(k);
    });
    array1D eigenens;
    for (auto row : eigenens2D)
        for (auto e : row)
            eigenens.push_back(e);

    auto sum_fd = [TeV, eigenens, kpts] (double mu) -> double
    {
        array1D fd_terms(kpts.size() * EigenValue::nbands, 0.0);
        std::transform(eigenens.begin(), eigenens.end(), fd_terms.begin(),
                       [TeV, mu] (double e) { return fermi_dirac(e - mu, TeV); });
        return std::reduce(fd_terms.begin(), fd_terms.end(), 0.0) / kpts.size() - EigenValue::nelec / 2.0; // spin
    };

    // 2. we know nelec, we have to solve equation nelec - \sum_{kn} f(e_{nk}, TeV) = 0
    double a{ -100.0 }; // initial left boundary for root, in eV
    double b{  100.0 }; // initial right boundary for root, in eV
    EigenValue::eF = find_root_bisect(sum_fd, a, b, crit);
}

} // quantities
} // skies
