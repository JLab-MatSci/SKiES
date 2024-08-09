#pragma once

#include <string>
#include <complex>
#include <iostream>

#include <skies/common/ndimarrays.h>
#include <skies/quantities/basic_quantity.h>

#ifdef __cplusplus
extern "C" {
#endif

void interpGmat1DAtq(const double* kx, const double* ky, const double* kz,
                   const double* qx, const double* qy, const double* qz,
                   const int* nu, const int* n, const int* m, double*);

void interpGmat2Atq(const double* kx, const double* ky, const double* kz,
                   const double* qx, const double* qy, const double* qz,
                   const int* nu, const int* n, const int* m, double*);

#ifdef __cplusplus
}
#endif

namespace skies { namespace quantities {

using complex   = std::complex<double>;
using array1D_c = std::vector<complex>;
using array2D_c = std::vector<array1D_c>;
using array3D_c = std::vector<array2D_c>;
using array4D_c = std::vector<array3D_c>;
using array5D_c = std::vector<array4D_c>;


class EPHMatrixDrawable : public AnyQuantity {
public:
    arrays::array1D k = {0.0, 0.0, 0.0}; // usually everyone uses Gamma here

    int nu; // phonon branch
    int nini, nfin; // initial and final bands

public:
    EPHMatrixDrawable(const arrays::array1D& k, int nu, int n, int m)
        : k(k), nu(nu), nini(n), nfin(m) {}

    std::string name() const override { return "EPHMatrix"; }

    arrays::array1D interpolate_at(const arrays::array1D& q) const override;
};

class EPHMatrixSquared {
public:
    static double interpolate_at(const arrays::array1D& k,
                                 const arrays::array1D& q,
                                 int nu,
                                 int nini,
                                 int nfin);
};

} // quantities
} // skies
