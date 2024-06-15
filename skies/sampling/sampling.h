#pragma once

#include <map>
#include <string>
#include <functional>
#include <stdexcept>

#include <skies/common/ndimarrays.h>

namespace skies { namespace bzsampling {

using SamplingFunc   = std::function<double(double, double)>;
using SamplingParams = std::map<const std::string, const std::string>;

enum class SamplType { gs, fd };

inline SamplType hash_type (const std::string& type)
{
    if (type == "gs") return SamplType::gs;
    if (type == "fd") return SamplType::fd;
    throw std::runtime_error("Unknown type of sampling. Please enter one of 'gs' or 'fd'");
}

double fermi_dirac(double x, double sigma);
double bose_einstein(double x, double sigma);
double fd_derivative(double x, double sigma);
double gauss(double x, double sigma);

inline SamplingFunc switch_sampling(std::string type)
{
    switch (hash_type(type))
    {
        case SamplType::gs :
            return gauss;
        case SamplType::fd :
            return fd_derivative;
        default:
            break;
    }
    return SamplingFunc{};
}

double smear_with_fd(const arrays::array1D& quan, const arrays::array1D& range, double e, double sigma);

} // bzsampling
} // skies
