#pragma once

namespace skies { namespace spectral {

double calc_lambda_tr_one(double sigma, double omega);

void calc_lambda_tr(double sigma, double begin, double end, int bins);

} // spectral
} // skies
