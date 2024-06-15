#pragma once

#include <list>
#include <fstream>

#include <skies/quantities/eigenfreqs.h>
#include <skies/lattices/latt_protocol.h>
#include <skies/quantities/basic_quantity.h>

namespace skies { namespace interpol {

using arrays::operator+;
using arrays::operator-;
using arrays::operator*;
using arrays::matmul;
using arrays::find_norm;

// array1D must be given in cryst coords
using kpLabelCoords = std::vector<std::pair<const std::string, const arrays::array1D>>;

template <typename Quan, typename... Args>
double evaluate_between_seq_points(std::ofstream& quan_file,
                                 const arrays::array1D& begin,
                                 const arrays::array1D& end,
                                 const arrays::array2D& Bmat,
                                 double start,
                                 int bins,
                                 Args&&... args)
{
    double delta = std::sqrt(find_norm(matmul(Bmat,  end - begin)));
    double step = delta / (bins + 1);
    for (int i = 0; i < bins; ++i)
    {
        auto medium = begin + (end - begin) * ((i + 1.0) / (bins + 1.0));
        auto values = Quan(std::forward<Args>(args)...).interpolate_at(medium);
        start += step;
        quan_file << std::right;
        quan_file << std::setw(23) << std::setprecision(4) << start;
        for (auto&& v : values)
            quan_file << std::setw(12) << std::setprecision(4) << v;
        quan_file << std::endl;
    }
    return delta;
}

template <typename Quan, typename... Args>
void evaluate_along_kp_path(const kpLabelCoords& kp_path,
                            int bins = 30,
                            Args&&... args)
{
    auto Bmat = Lattprotocol::calc_inv_cell();

    double mileage{ 0.0 };
    std::ofstream klabels("KLABELS");
    klabels << std::right;
    klabels << std::setw(7) << "# Label    ";
    klabels << std::setw(19) << "Distance [Angstrom]";
    klabels << std::endl;

    std::ofstream quanFile(Quan(std::forward<Args>(args)...).name() + ".dat");
    quanFile << std::right;
    quanFile << std::setw(26) << "# Kpath [1 / Angstrom]   ";
    quanFile << std::setw(14) << Quan(std::forward<Args>(args)...).name() <<  " dispersion [eV]";
    quanFile << std::endl;

    auto kp_end = std::prev(kp_path.end());
    for (auto it = kp_path.begin(); it != kp_end; ++it) // loop over all segments
    {
        // write to KLABELS file
        auto label = it->first;
        klabels << std::setw(7) << label;
        klabels << std::setw(23) << std::setprecision(6) << mileage;
        klabels << std::endl;

        auto begin = it->second;
        auto   end = std::next(it)->second;
        double delta = evaluate_between_seq_points<Quan, Args...>(quanFile, begin, end, Bmat, mileage, bins, args...);
        mileage += delta;
    }

    auto label = kp_end->first;
    klabels << std::setw(7) << label;
    klabels << std::setw(23) << std::setprecision(6) << mileage;
    klabels << std::endl;
    klabels.close();
    quanFile.close();
}

} // interpol
} // skies
