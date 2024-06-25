/**
 @file
 @brief Routines for density of states calculations
 @author Galtsov Ilya
 */
#pragma once

#include <fstream>

#include <skies/utils/mpi_wrapper.h>
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
    double res{ 0.0 };
    for (size_t ikpt = 0; ikpt < nkpt; ++ikpt)
        for (size_t ibnd = 0; ibnd < nbnd; ++ibnd)
            res += weights[ikpt][ibnd] * sampling(value - values[ikpt][ibnd], smearing);
    res /= nkpt; // needed by definition of DOS
    return res;
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
    double res{ 0.0 };
    for (size_t ikpt = 0; ikpt < nkpt; ++ikpt)
        for (size_t ibnd = 0; ibnd < nbnd; ++ibnd)
            res += sampling(value - values[ikpt][ibnd], smearing);
    res /= nkpt; // needed by definition of DOS
    return res;
}

void evaluate_trdos(const arrays::array2D& grid,
                    const arrays::array1D& range,
                    double smearing,
                    bzsampling::SamplingFunc sampl_type,
                    char cart);

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
    std::transform(grid.begin(), grid.end(), values.begin(),
                        [] (const arrays::array1D& k) { return Quan().interpolate_at(k); });

    int rank{ 0 };
    arrays::array1D DOSes(range.size(), 0.0);
    arrays::array1D DOSes_tmp(range.size(), 0.0);
#ifdef SKIES_MPI
    auto rcounts = mpi::prepare_rcounts_displs(range.size()).first;
    auto displs  = mpi::prepare_rcounts_displs(range.size()).second;

    rank = skies::mpi::rank();
    for (int i = displs[rank]; i < displs[rank] + rcounts[rank]; ++i)
        DOSes_tmp[i] = evaluate_dos_at_value<Quan>(range[i], smearing, sampl_type, values);

    MPI_Reduce(DOSes_tmp.data(), DOSes.data(), range.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(DOSes.data(), range.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    std::transform(range.begin(), range.end(), DOSes.begin(), [] (double v) {
        return evaluate_dos_at_value<Quan>(v, smearing, sampl_type, values);
    });
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
