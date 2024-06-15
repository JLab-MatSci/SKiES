/**
 @file
 @brief Routines for density of states calculations
 @author Galtsov Ilya
 */
#pragma once

#include <fstream>

#include <skies/quantities/eigenfreqs.h>
#include <skies/quantities/elvelocs.h>
#include <skies/lattices/kp_protocol.h>

#include <launch/timer.h>

namespace skies { namespace quantities {

/**
 * \brief Evaluates transport DOS at given value (in eV)
 * @param grid 2D grid of k-points
 * @param value energy to calculate DOS at
 * @param smearing energy smearing used in given type of sampling (in eV)
 * @param sampling specified type of smearing function
 * @param cart one of 'x', 'y' or 'z'
*/
double evaluate_trdos_at_value(const arrays::array2D&  grid,
                               double value,
                               double smearing,
                               bzsampling::SamplingFunc sampling,
                               char cart);

/**
 * \brief Evaluates smeared transport DOS at given value (in eV)
 * @param grid 2D grid of k-points
 * @param value energy to calculate DOS at
 * @param smearing energy smearing used in given type of sampling (in eV)
 * @param sampling specified type of smearing function
 * @param elec_temp smearing parameter (in eV)
 * @param one of 'x', 'y' or 'z'
*/
double evaluate_smeared_trdos_at_value(const arrays::array2D& grid,
                                       double value,
                                       double smearing,
                                       bzsampling::SamplingFunc sampling,
                                       double elec_temp,
                                       char cart);

/**
 * \brief Evaluates DOS at given value (in eV). The quantity of interest is provided as a template parameter
 * @param grid 2D grid of k-points
 * @param value energy to calculate DOS at
 * @param smearing energy smearing used in given type of sampling (in eV)
 * @param sampling specified type of smearing function
 * @param weights 2D array calculated over given grid to serve as weights multipliers
*/
template <typename Quan>
double evaluate_dos_at_value(const arrays::array2D& grid,
                             double value,
                             double smearing,
                             bzsampling::SamplingFunc sampling,
                             const arrays::array2D& weights)
{
    size_t nkpt = grid.size();
    size_t nbnd = std::is_same_v<Quan, EigenFrequencyDrawable> ? EigenFrequencyDrawable::nmodes
	    						    : EigenValue::nbands;
    arrays::array2D values;
    values.resize(nkpt, arrays::array1D(nbnd));
    std::transform(grid.begin(), grid.end(), values.begin(),
                        [] (const arrays::array1D& k) { return Quan().interpolate_at(k); });
    double res{ 0.0 };
    for (size_t ikpt = 0; ikpt < nkpt; ++ikpt) {
        for (size_t ibnd = 0; ibnd < nbnd; ++ibnd) {
            res += weights[ikpt][ibnd] * sampling(value - values[ikpt][ibnd], smearing);
        }
    }
    res /= nkpt; // needed by definition of DOS
    return res;
}

/**
 * \brief Evaluates DOS in a given range of energies. The quantity of interest is provided as a template parameter.
 * @param grid 2D grid of k-points at which the quantity is evaluated
 * @param begin lower energy in the range (in eV)
 * @param end   upper energy in the range (in eV)
 * @param bins  number of bins in the range
 * @param smearing energy smearing used in given type of sampling (in eV)
 * @param sampl_type specified type of smearing function
 * @param weights 2D array calculated over given grid to serve as weights multipliers
*/
template <typename Quan>
void evaluate_dos(const arrays::array2D& grid,
                  double begin,
                  double end,
                  int bins,
                  double smearing,
                  bzsampling::SamplingFunc sampl_type,
                  const arrays::array2D& weights)
{
    std::ofstream os(Quan().name() + "DOS.dat");
    auto range = arrays::create_range(begin, end, bins);
    os << std::right;
    os << std::setw(12) << "# Energy, eV";
    os << std::setw(25) << Quan().name() + " DOS";
    if (Quan().name() == "Velocities")
        os << std::setw(24) << " [13.605685 * Ry bohr^2]";
    else
        os << std::setw(9) << " [1 / eV]";
    os << std::endl;

    for (auto&& v : range)
        if (!std::is_same_v<Quan, VelocitiesDrawable>)
            os << std::setprecision(6) << std::setw(12) << v << std::setw(34) << evaluate_dos_at_value<Quan>(grid, v, smearing, sampl_type, weights) << '\n';
        else
            os << std::setprecision(6) << std::setw(12) << v << std::setw(49) << evaluate_dos_at_value<EigenValue>(grid, v, smearing, sampl_type, weights) << '\n';
    os.close();
}

} // quantities
} // skies
