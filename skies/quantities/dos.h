/**
 @file
 @brief Routines for density of states calculations
 @author Galtsov Ilya
 */
#pragma once

#include <fstream>
#include <numeric>

#ifdef SKIES_MPI
#include <skies/utils/mpi_wrapper.h>
#endif
#include <skies/utils/tbb_wrapper.h>
#include <skies/quantities/eigenfreqs.h>
#include <skies/quantities/elvelocs.h>
#include <skies/lattices/kp_protocol.h>

#include <launch/timer.h>

namespace skies { namespace quantities {

/**
 * \brief Evaluates transport DOS at given value (in eV)
 * @param value energy to calculate DOS at
 * @param smearing energy smearing used in given type of sampling (in eV)
 * @param sampling specified type of smearing function
 * @param velocs_squared for perfomance issues precalculated squared velocs must be provided
*/
double evaluate_trdos_at_value(double value,
                               double smearing,
                               bzsampling::SamplingFunc sampling,
                               const arrays::array2D& energies,
                               const arrays::array2D& velocs_squared);

/**
 * \brief Evaluates smeared transport DOS at given value (in eV)
 * @param value energy to calculate DOS at
 * @param smearing energy smearing used in given type of sampling (in eV)
 * @param sampling specified type of smearing function
 * @param elec_temp smearing parameter (in eV)
*/
double evaluate_smeared_trdos_at_value(double value,
                                       double smearing,
                                       bzsampling::SamplingFunc sampling,
                                       const arrays::array2D& energies,
                                       const arrays::array2D& squared_velocs,
                                       double elec_temp);

/**
 * \brief Evaluates DOS at given value (in eV). The quantity of interest is provided as a template parameter
 * @param value energy to calculate DOS at
 * @param smearing energy smearing used in given type of sampling (in eV)
 * @param sampling specified type of smearing function
 * @param values 2D array of precalculated values on the given grid
 * @param weights 2D array calculated over given grid to serve as weights multipliers
*/
template <typename Quan>
double evaluate_dos_at_value(double value,
                             double smearing,
                             bzsampling::SamplingFunc sampling,
                             const arrays::array2D& values, 
                             const arrays::array2D& weights)
{
    auto nkpt = values.size();
    auto nbnd = values[0].size();
    std::vector<size_t> ikpts(nkpt);
    std::iota(ikpts.begin(), ikpts.end(), 0);
    double dos = std::transform_reduce(PAR ikpts.begin(), ikpts.end(), 0.0, std::plus<double>(),
        [&] (auto&& ik) -> double {
            std::vector<size_t> ibands(nbnd);
            std::iota(ibands.begin(), ibands.end(), 0);
            return std::transform_reduce(ibands.begin(), ibands.end(), 0.0, std::plus<double>(),
                [&] (auto&& n) -> double {
                    return weights[ik][n] * sampling(value - values[ik][n], smearing);
            });
        }
    );
    dos /= nkpt;
    return dos;
}

/**
 * \brief Evaluates DOS at given value (in eV). The quantity of interest is provided as a template parameter
 * @param value energy to calculate DOS at
 * @param smearing energy smearing used in given type of sampling (in eV)
 * @param sampling specified type of smearing function
 * @param values 2D array of precalculated values on the given grid
*/
template <typename Quan>
double evaluate_dos_at_value(double value,
                             double smearing,
                             bzsampling::SamplingFunc sampling,
                             const arrays::array2D& values)
{
    auto nkpt = values.size();
    auto nbnd = values[0].size();
    std::vector<size_t> ikpts(nkpt);
    std::iota(ikpts.begin(), ikpts.end(), 0);
    double dos = std::transform_reduce(PAR ikpts.begin(), ikpts.end(), 0.0, std::plus<double>(),
        [&] (auto&& ik) -> double {
            std::vector<size_t> ibands(nbnd);
            std::iota(ibands.begin(), ibands.end(), 0);
            return std::transform_reduce(ibands.begin(), ibands.end(), 0.0, std::plus<double>(),
                [&] (auto&& n) -> double {
                    return sampling(value - values[ik][n], smearing);
            });
        }
    );
    dos /= nkpt;
    return dos;
}

void evaluate_trdos(const arrays::array2D& grid,
                    const arrays::array1D& range,
                    double smearing,
                    bzsampling::SamplingFunc sampl_type);

/**
 * \brief Evaluates DOS in a given energy range. The quantity of interest is provided as a template parameter.
 * @param grid 2D grid of k-points at which the quantity is evaluated
 * @param range energy range
 * @param smearing energy smearing used in given type of sampling (in eV)
 * @param sampl_type specified type of smearing function
 * @param weights 2D array calculated over given grid to serve as weights multipliers
*/
template <typename Quan>
void evaluate_dos(const arrays::array2D& grid,
                  const arrays::array1D& range,
                  double smearing,
                  bzsampling::SamplingFunc sampl_type,
                  const arrays::array2D& weights)
{
    if (std::is_same_v<Quan, VelocitiesDrawable>)
        throw std::runtime_error("Use evaluate_trdos instead\n");
    std::ofstream os(Quan().name() + "DOS.dat");
    os << std::right;
    os << std::setw(12) << "# Energy, eV";
    os << std::setw(25) << Quan().name() + " DOS";
        os << std::setw(9) << " [1 / eV]";
    os << std::endl;

    auto nkpt = grid.size();
    size_t nbnd = std::is_same_v<Quan, EigenFrequencyDrawable>
                                ? EigenFrequencyDrawable::nmodes
	    					    : EigenValue::nbands;
    arrays::array2D values;
    values.resize(nkpt, arrays::array1D(nbnd));
    std::transform(PAR grid.begin(), grid.end(), values.begin(), [] (auto&& k) {
        return Quan().interpolate_at(k);
    });

    arrays::array1D DOSes(range.size(), 0.0);
    arrays::array1D DOSes_tmp(range.size(), 0.0);
    std::transform(range.begin(), range.end(), DOSes.begin(), [&] (double v) {
        return evaluate_dos_at_value<Quan>(v, smearing, sampl_type, values, weights);
    });

    int rank{ 0 };
#ifdef SKIES_MPI
    rank = skies::mpi::rank();
#endif
    if (!rank)
    {
        for (size_t i = 0; i < range.size(); ++i)
            os << std::setprecision(6) << std::setw(12) << range[i] << std::setw(34) << DOSes[i] << std::endl;
        os.close();
    }
    return;
}

} // quantities
} // skies
