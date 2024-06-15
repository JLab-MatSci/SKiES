#include <vector>
#include <math.h>
#include <stdexcept>
#include <numeric>

#include <skies/sampling/sampling.h>

namespace skies { namespace bzsampling {

using namespace arrays;

double fermi_dirac(double x, double sigma)
{
    return 1 / (exp(x / sigma) + 1);
}

double bose_einstein(double x, double sigma)
{
    return 1 / (exp(x / sigma) - 1);
}

double fd_derivative(double x, double sigma)
{
    if (x / sigma > 36 || x / sigma < -36)
        return 0.0;
    return exp(x / sigma) / (exp(x / sigma) + 1) 
                               / (exp(x / sigma) + 1) / sigma;
}

double gauss(double x, double sigma)
{
    if (x / sigma > 36 || x / sigma < -36)
        return 0.0;
    return 1.0 / sqrt(2 * 4 * atan(1.0)) / sigma
                        * exp(-(x / sigma) * (x / sigma) / 2);
}

// quan must be evaluated on a given range (in [eV]), sigma is in [eV]
// returns smeared value at given e
double smear_with_fd(const array1D& quan, const array1D& range, double e, double sigma)
{
    if (range.size() < 2)
        throw std::runtime_error("Too few points in teh given range.");
    if (quan.size() != range.size())
        throw std::runtime_error("Quantity array and range array must have the same size");
    array1D smearing(quan.size(), 0.0);
    std::transform(range.begin(), range.end(), smearing.begin(), [sigma, e] (double val) { return fd_derivative(val - e, sigma); });
    double d = range[1] - range[0];
    return d * std::inner_product(quan.begin(), quan.end(), smearing.begin(), 0.0);
}

} // bzsampling
} // skies
