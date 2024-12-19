/*-----------------------------------------------------------------------
    * SKiES - Solver of Kinetic Equation for Solids
    * Version 1.0.0
    * 
    * (C) 2024 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)
    *
    * SKiES may only be utilized for non-profit research.
    * Citing appropriate sources is required when using SKiES.
    * 
    * Project homepage: https://github.com/i-Galts/SKiES
    * 
    * This file: Electron-phonon matrix elements implementation
    * Implements interpolation of electron-phonon coupling matrix elements
    * and provides squared matrix element calculations.
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */

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
