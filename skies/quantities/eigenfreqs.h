/**
 @file
 @brief Description of eigen frequencies class (dynamical matrix eigen values)
 @author Galtsov Ilya
 */
#pragma once

#include <skies/quantities/basic_quantity.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief External function called from EPW Fortran code to obtain the number of phonon modes in the calculation
 * @param nmodes number of modes
*/
void fillNmodes(int* nmodes);

/**
 * \brief External function called from EPW Fortran code to interpolate eigen frequency at provided k-point
 * k-point is given by three pointers to its x-, y-, z-components for interoperability with Fortran.
 * The interpolated eigen frequencies are saved as output array rfreqs[] (in Ry), the size of this array equals the
 * number of phonon bands in the calculation. 
 * @param kx the first k-point component
 * @param ky the second k-point component
 * @param kz the third k-point component
 * @param rfreqs output array which contains interpolated eigen frequencies
*/
void interpEigenFreqAt(const double* kx, const double* ky, const double* kz, double rfreqs[]);

void interpEigenFreq1DAt(const double* kx, const double* ky, const double* kz, double rfreqs[]);

#ifdef __cplusplus
}
#endif

namespace skies { namespace quantities {


class EigenFrequencyDrawable : public AnyQuantity {
private:
    arrays::array1D values_;
public:
    static size_t nmodes;
public:
    std::string name() const override { return "EigenFrequency"; }
    arrays::array1D interpolate_at(const arrays::array1D& q) const override;
};

class EigenFrequency {
public:
    static size_t nmodes;
    static arrays::array1D interpolate_at(const arrays::array1D& q);
};

} // quantities
} // skies
