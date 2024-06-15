#include <cassert>

#include <skies/quantities/gmatrix.h>


using namespace skies::arrays;

namespace skies { namespace quantities {

array1D EPHMatrixDrawable::interpolate_at(const array1D& q) const {
    assert(q.size() == 3);
    double rgmat;
    interpGmat1DAtq(&k[0], &k[1], &k[2], &q[0], &q[1], &q[2], &nu, &nini, &nfin, &rgmat);
    array1D gmat(1); gmat[0] = rgmat * units::Ry_in_eV;
    return gmat;
}

double EPHMatrixSquared::interpolate_at(const arrays::array1D& k,
                                         const arrays::array1D& q,
                                         int nu,
                                         int nini,
                                         int nfin)
{
    assert(k.size() == 3);
    assert(q.size() == 3);
    assert(nini > -1);
    assert(nfin > -1);
    assert(nu > -1);
    double gmat2;
    interpGmat2Atq(&k[0], &k[1], &k[2], &q[0], &q[1], &q[2], &nu, &nini, &nfin, &gmat2);
    return gmat2 * units::Ry_in_eV * units::Ry_in_eV;
}

} // quantities
} // skies
