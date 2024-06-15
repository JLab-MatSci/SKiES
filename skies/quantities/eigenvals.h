/**
 @file
 @brief Description of eigen energies class (hamiltonian matrix eigen values)
 @author Galtsov Ilya
 */
#pragma once

#include <cassert>

#include <skies/quantities/basic_quantity.h>

#ifdef __cplusplus
extern "C" {
#endif

void fillNelec(double* nelec);
void fillNbands(int* nbands);
void interpEigenValueAt(const double* kx, const double* ky, const double* kz, double[]);

#ifdef __cplusplus
}
#endif

namespace skies { namespace quantities {

class EigenValueDrawable : public AnyQuantity {
public:
    static size_t nbands;
    arrays::array2D values;
public:
    std::string name() const override { return "EigenValue"; }
    arrays::array1D interpolate_at(const arrays::array1D& k) const override;
};

class EigenValue {
public:
    static size_t nbands;
    static int nelec;
    static double eF;

    static void find_eF(double TeV, double crit = 1e-16);
    static arrays::array1D interpolate_at(const arrays::array1D& k);
};

} // quantities
} // skies
