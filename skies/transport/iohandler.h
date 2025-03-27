/*-----------------------------------------------------------------------
    * SKiES - Solver of Kinetic Equation for Solids
    * 
    * (C) 2025 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)
    *
    * SKiES may only be utilized for non-profit research.
    * Citing appropriate sources is required when using SKiES.
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */
#pragma once

#include <fstream>
#include <skies/common/ndimarrays.h>

namespace skies { namespace transport {

/**
 * @brief Enumeration to specify the type of resistivity being handled.
 */
enum class ResistType 
{
    Electrical, ///< Represents electrical resistivity.
    Thermal     ///< Represents thermal resistivity.
};

/**
 * @brief Input handler class for reading and managing transport-related data.
 *
 * This class is responsible for reading and storing data from an a2f file, including
 * electron-phonon coupling matrices, frequency arrays, and other parameters required
 * for transport calculations.
 */
class IHandler {
public:
    /**
     * @brief Constructs an IHandler object by reading data from the specified a2f file.
     *
     * @param a2f_fnm The file name of the a2f file containing electron-phonon coupling data.
     */
    IHandler(const char* a2f_fnm);

    /**
     * @brief Destructor for the IHandler class.
     */
    ~IHandler();

    /**
     * @brief Returns the electron-phonon coupling matrix (a2f).
     *
     * @return const arrays::array2D& Reference to the 2D array representing the a2f matrix.
     */
    const arrays::array2D& a2f() { return a2f_; }

    /**
     * @brief Returns the array of phonon frequencies (omegas).
     *
     * @return const arrays::array1D& Reference to the 1D array representing phonon frequencies.
     */
    const arrays::array1D& omegas() { return omegas_; }

    /**
     * @brief Returns the array of electronic energies (epsilons).
     *
     * @return const arrays::array1D& Reference to the 1D array representing electronic energies.
     */
    const arrays::array1D& epsilons() { return epsilons_; }

    /**
     * @brief Returns the array of transport density of states (transDOSes).
     *
     * @return arrays::array1D& Reference to the 1D array representing transport DOS values.
     */
    arrays::array1D& transDOSes() { return transDOSes_; }

    /**
     * @brief Returns the smearing parameter for electrons.
     *
     * @return double The electron smearing value.
     */
    double elec_smearing() { return elec_smearing_; }

    /**
     * @brief Returns the smearing parameter for phonons.
     *
     * @return double The phonon smearing value.
     */
    double phon_smearing() { return phon_smearing_; }

    /**
     * @brief Returns the alpha parameter.
     *
     * @return int The alpha value.
     */
    int alpha() { return alpha_; }

    /**
     * @brief Returns the beta parameter.
     *
     * @return int The beta value.
     */
    int beta() { return beta_; }

    /**
     * @brief Returns the sign parameter.
     *
     * @return int The sign value.
     */
    int sign() { return sign_; }

    /**
     * @brief Returns the sign_pr parameter.
     *
     * @return int The sign_pr value.
     */
    int sign_pr() { return sign_pr_; }

private:
    std::ifstream ifs_; ///< Input file stream for reading the a2f file.

    double elec_smearing_; ///< Smearing parameter for electrons.
    double phon_smearing_; ///< Smearing parameter for phonons.
    char alpha_;           ///< Alpha parameter.
    char beta_;            ///< Beta parameter.
    int sign_;             ///< Sign parameter.
    int sign_pr_;          ///< Sign_pr parameter.

    arrays::array1D epsilons_;   ///< Array of electronic energies.
    arrays::array1D transDOSes_; ///< Array of transport density of states.
    arrays::array1D omegas_;     ///< Array of phonon frequencies.
    arrays::array2D a2f_;        ///< Electron-phonon coupling matrix.
};

/**
 * @brief Output handler class for writing transport-related results to files.
 *
 * This class is responsible for writing computed resistivity data to output files.
 * It supports both scalar and tensorial resistivity outputs.
 */
class OHandler {
public:
    /**
     * @brief Constructs an OHandler object for writing resistivity data.
     *
     * @param a2f_fnm The file name of the a2f file (optional, for reference).
     * @param cond_fnm The file name where the resistivity data will be written.
     * @param type The type of resistivity (electrical or thermal).
     * @param ion_Temps Optional array of ion temperatures (default: empty array).
     */
    OHandler(const char* a2f_fnm, const char* cond_fnm, ResistType type, const arrays::array1D& ion_Temps = arrays::array1D());

    /**
     * @brief Destructor for the OHandler class.
     */
    ~OHandler();

    /**
     * @brief Writes scalar resistivity data to the output file.
     *
     * @param Temps Array of temperatures corresponding to the resistivity values.
     * @param resist Array of resistivity values.
     */
    void dump(const arrays::array1D& Temps, const arrays::array1D& resist);

    /**
     * @brief Writes tensorial resistivity data to the output file.
     *
     * @param Temps Array of temperatures corresponding to the resistivity tensors.
     * @param resist 2D array of resistivity tensors.
     */
    void dump(const arrays::array1D& Temps, const arrays::array2D& resist);

private:
    std::ofstream ofs_;         ///< Output file stream for writing resistivity data.
    arrays::array1D ion_Temps_; ///< Array of ion temperatures (optional).
};

} // namespace transport
} // namespace skies
