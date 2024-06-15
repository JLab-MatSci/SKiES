#include <fstream>
#include <numeric>
#include <algorithm>
#include <cassert>

#include <skies/common/alg.h>
#include <skies/common/ndimarrays.h>
#include <skies/sampling/sampling.h>

using namespace skies::arrays;
using namespace skies::bzsampling;

namespace skies { namespace spectral {

double calc_lambda_tr_one(double sigma, double omega)
{
    std::ifstream ifs("lambda_tr.dat");
    std::string line;

    array1D freqs;
    array1D lambdas;

    if (ifs.good()) getline(ifs, line);
    if (ifs.good()) getline(ifs, line);
    auto splitted_line = custom_split(line, ' ');
    int nqpt = std::stoi(splitted_line[3].data(), 0);
    if (ifs.good()) getline(ifs, line);
    while (ifs.good()) {
        getline(ifs, line);
        if (!line.empty()) {
            auto splitted_line = custom_split(line, ' ');
            freqs.push_back(std::strtod(splitted_line[4].data(), 0));
            lambdas.push_back(std::strtod(splitted_line[5].data(), 0));
        }
    }

    assert(freqs.size() == lambdas.size());
    array1D deltas(freqs.size());
    std::transform(freqs.begin(), freqs.end(), deltas.begin(), 
                            [sigma, omega] (double omqnu) { return fd_derivative(omega - omqnu, sigma); });

    return std::inner_product(deltas.begin(), deltas.end(), lambdas.begin(), 0) / nqpt;
}

void calc_lambda_tr(double sigma, double begin, double end, int bins)
{
    std::ofstream os("lambda_tr_integ.dat");
    os << "Integrated mode-resolved transport coupling strength" << std::endl;
    os << "Eigenfrequency" << "            " << "integated lambda" << std::endl;
    arrays::array1D range(bins + 2);
    double d = (end - begin) / (bins + 1);
    std::generate(range.begin(), range.end(), [d, pos = begin - d] () mutable { pos += d; return pos; });
    for (auto om : range)
        os << om << "            " << calc_lambda_tr_one(sigma, om) << std::endl;
    os.close();
}

} // spectral
} // skies
