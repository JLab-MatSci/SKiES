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

#include <skies/common/ndimarrays.h>
#include <skies/sampling/sampling.h>
#include <skies/sampling/tetrahedra.h>

namespace skies { namespace spectral {

/**
 * @brief Class for calculating spectral functions.
 *
 * The `SpecFunc` class provides functionality to compute spectral functions for electronic
 * and phononic systems. It supports multiple initialization methods, including grid-based
 * sampling and file-based input. The class also includes helper functions for preparing
 * eigenvalues, velocities, and Fermi surface harmonics, as well as calculating inner sums
 * and spectral functions.
 */
class SpecFunc {
public:
    /**
     * @brief Constructs a `SpecFunc` object using grid specifications and sampling functions.
     *
     * @param kpgrid Grid dimensions for the electron k-points.
     * @param qpgrid Grid dimensions for the phonon q-points.
     * @param elec_sampling Sampling function for electrons.
     * @param phon_sampling Sampling function for phonons.
     * @param elec_smearing Smearing parameter for electrons.
     * @param phon_smearing Smearing parameter for phonons.
     * @param sign Sign parameter for velocities (+1 or -1).
     * @param sign_pr Primed counterpart of the sign parameter.
     * @param Te Electronic temperature in eV.
     * @param alpha Cartesian index of the first velocity component.
     * @param beta Cartesian index of the second velocity component.
     * @param is_tetra Flag indicating whether tetrahedron integration is used (default: false).
     */
    SpecFunc(const std::vector<size_t>& kpgrid,
            const std::vector<size_t>& qpgrid,
            bzsampling::SamplingFunc elec_sampling,
            bzsampling::SamplingFunc phon_sampling,
            double elec_smearing,
            double phon_smearing,
            int sign,
            int sign_pr,
            double Te,
            char alpha,
            char beta,
            bool is_tetra = false);

    /**
     * @brief Constructs a `SpecFunc` object with precomputed epsilon values.
     *
     * This constructor allows specifying precomputed epsilon values for the calculation.
     *
     * @param kpgrid Grid dimensions for the electron k-points.
     * @param qpgrid Grid dimensions for the phonon q-points.
     * @param epsilons Precomputed epsilon values.
     * @param elec_sampling Sampling function for electrons.
     * @param phon_sampling Sampling function for phonons.
     * @param elec_smearing Smearing parameter for electrons.
     * @param phon_smearing Smearing parameter for phonons.
     * @param sign Sign parameter for velocities (+1 or -1).
     * @param sign_pr Primed counterpart of the sign parameter.
     * @param Te Electronic temperature in eV.
     * @param alpha Cartesian index of the first velocity component.
     * @param beta Cartesian index of the second velocity component.
     * @param is_tetra Flag indicating whether tetrahedron integration is used (default: false).
     */
    SpecFunc(const std::vector<size_t>& kpgrid,
            const std::vector<size_t>& qpgrid,
            const arrays::array1D& epsilons,
            bzsampling::SamplingFunc elec_sampling,
            bzsampling::SamplingFunc phon_sampling,
            double elec_smearing,
            double phon_smearing,
            int sign,
            int sign_pr,
            double Te,
            char alpha,
            char beta,
            bool is_tetra = false);

    /**
     * @brief Constructs a `SpecFunc` object from a lambda_tr.dat file.
     *
     * This constructor initializes the object using data from a file containing precomputed
     * lambda_tr values.
     *
     * @param fname File name containing the lambda_tr data.
     */
    SpecFunc(const std::string& fname);

    /**
     * @brief Destructor for the `SpecFunc` class.
     */
    ~SpecFunc();

    /**
     * @brief Calculates the spectral function for a given frequency.
     *
     * @param Omega Frequency at which the spectral function is calculated.
     * @return arrays::array1D The computed spectral function values.
     */
    arrays::array1D calc_spec_func(double Omega);

    /**
     * @brief Returns the electronic smearing parameter.
     *
     * @return double The electronic smearing value.
     */
    double elec_smearing() const { return elec_smearing_; }

    /**
     * @brief Returns the phononic smearing parameter.
     *
     * @return double The phononic smearing value.
     */
    double phon_smearing() const { return phon_smearing_; }

    /**
     * @brief Returns the number of k-points along the x-axis.
     *
     * @return size_t Number of k-points along the x-axis.
     */
    size_t nkx() const { return std::get<0>(kprot_.mesh()); }

    /**
     * @brief Returns the number of k-points along the y-axis.
     *
     * @return size_t Number of k-points along the y-axis.
     */
    size_t nky() const { return std::get<1>(kprot_.mesh()); }

    /**
     * @brief Returns the number of k-points along the z-axis.
     *
     * @return size_t Number of k-points along the z-axis.
     */
    size_t nkz() const { return std::get<2>(kprot_.mesh()); }

    /**
     * @brief Returns the number of q-points along the x-axis.
     *
     * @return size_t Number of q-points along the x-axis.
     */
    size_t nqx() const { return std::get<0>(qprot_.mesh()); }

    /**
     * @brief Returns the number of q-points along the y-axis.
     *
     * @return size_t Number of q-points along the y-axis.
     */
    size_t nqy() const { return std::get<1>(qprot_.mesh()); }

    /**
     * @brief Returns the number of q-points along the z-axis.
     *
     * @return size_t Number of q-points along the z-axis.
     */
    size_t nqz() const { return std::get<2>(qprot_.mesh()); }

    /**
     * @brief Returns the number of modes.
     *
     * @return size_t Number of phonon modes.
     */
    size_t nmds() const { return nmds_; }

    /**
     * @brief Returns the sign parameter for velocities.
     *
     * @return int The sign parameter (+1 or -1).
     */
    int sign() const { return sign_; }

    /**
     * @brief Returns the primed counterpart of the sign parameter.
     *
     * @return int The primed sign parameter.
     */
    int sign_pr() const { return sign_pr_; }

    /**
     * @brief Returns the Cartesian index of the first velocity component.
     *
     * @return char Cartesian index ('x', 'y', or 'z').
     */
    char alpha() const { return alpha_; }

    /**
     * @brief Returns the Cartesian index of the second velocity component.
     *
     * @return char Cartesian index ('x', 'y', or 'z').
     */
    char beta() const { return beta_; }

    /**
     * @brief Returns the array of epsilon values.
     *
     * @return arrays::array1D& Reference to the array of epsilon values.
     */
    arrays::array1D& epsilons() { return epsilons_; }

    /**
     * @brief Returns the array of epsilon values (const version).
     *
     * @return const arrays::array1D& Const reference to the array of epsilon values.
     */
    const arrays::array1D& epsilons() const { return epsilons_; }

    /**
     * @brief Returns the array of density of states (DOS) values.
     *
     * @return const arrays::array1D& Const reference to the array of DOS values.
     */
    const arrays::array1D& doses() const { return DOSes_; }

    /**
     * @brief Returns the array of transport density of states (trDOS) values.
     *
     * @return const arrays::array1D& Const reference to the array of trDOS values.
     */
    const arrays::array1D& trans_doses() const { return trDOSes_; }

    /**
     * @brief Sets the type of electronic smearing.
     *
     * @param type Type of electronic smearing (e.g., "Gaussian", "Lorentzian").
     */
    void set_type_of_el_smear(const std::string& type) { type_of_el_smear_ = type; }

    /**
     * @brief Sets the type of phononic smearing.
     *
     * @param type Type of phononic smearing (e.g., "Gaussian", "Lorentzian").
     */
    void set_type_of_ph_smear(const std::string& type) { type_of_ph_smear_ = type; }

    /**
     * @brief Returns the type of electronic smearing.
     *
     * @return std::string Type of electronic smearing.
     */
    std::string get_type_of_el_smear() const { return type_of_el_smear_; }

    /**
     * @brief Returns the type of phononic smearing.
     *
     * @return std::string Type of phononic smearing.
     */
    std::string get_type_of_ph_smear() const { return type_of_ph_smear_; }

    /**
     * @brief Sets the lower band index for calculations.
     *
     * @param low_band Lower band index.
     */
    void set_low_band(int low_band) { low_band_ = low_band; }

    /**
     * @brief Sets the upper band index for calculations.
     *
     * @param high_band Upper band index.
     */
    void set_high_band(int high_band) { high_band_ = high_band; }

    /**
     * @brief Returns the lower band index.
     *
     * @return int Lower band index.
     */
    int get_low_band() const { return low_band_; }

    /**
     * @brief Returns the upper band index.
     *
     * @return int Upper band index.
     */
    int get_high_band() const { return high_band_; }

    /**
     * @brief Sets the electronic temperature.
     *
     * @param Te Electronic temperature in eV.
     */
    void set_elec_temp(double Te) { Te_ = Te; }

    /**
     * @brief Returns the electronic temperature.
     *
     * @return double Electronic temperature in eV.
     */
    double get_elec_temp() const { return Te_; }

    /**
     * @brief Checks if the calculation is set to continue from a previous state.
     *
     * @return bool True if continuing calculation, false otherwise.
     */
    bool is_continue_calc() const { return is_continue_calc_; }

    /**
     * @brief Checks if tetrahedron integration is enabled.
     *
     * @return bool True if tetrahedron integration is enabled, false otherwise.
     */
    bool is_tetra() const { return is_tetra_; }

private:
    /**
     * @brief Initializes internal variables and prepares data structures.
     */
    void init();

    /**
     * @brief Prepares eigenfrequencies for calculations.
     *
     * @param eigenfreqs Array to store eigenfrequencies.
     */
    void prepare_eigenfreqs(arrays::array2D& eigenfreqs);

    /**
     * @brief Prepares eigenenergies for calculations.
     *
     * @param eigenens Array to store eigenenergies.
     * @param q Optional q-point vector (default: zero vector).
     */
    void prepare_eigenens(arrays::array2D& eigenens, const arrays::array1D& q = {0.0, 0.0, 0.0});

    /**
     * @brief Prepares electronic velocities for calculations.
     *
     * @param cart Cartesian index ('x', 'y', or 'z').
     * @param elvelocs Array to store electronic velocities.
     * @param q Optional q-point vector (default: zero vector).
     */
    void prepare_velocs(char cart, arrays::array2D& elvelocs, const arrays::array1D& q = {0.0, 0.0, 0.0});

    /**
     * @brief Prepares squared electronic velocities for calculations.
     *
     * @param elvelocs Array of electronic velocities.
     * @param elvelocs_sq Array to store squared velocities.
     */
    void prepare_squared_velocs(const arrays::array2D& elvelocs, arrays::array2D& elvelocs_sq);

    /**
     * @brief Prepares Fermi surface harmonics for calculations.
     *
     * @param elvelocs Array of electronic velocities.
     * @param elvelocs_sq Array of squared velocities.
     * @param fsh Array to store Fermi surface harmonics.
     */
    void prepare_fsh(const arrays::array2D& elvelocs, const arrays::array2D& elvelocs_sq, arrays::array2D& fsh);

    /**
     * @brief Prepares Fermi surface harmonics using tetrahedron handlers.
     *
     * @param th_dos Tetrahedron handler for DOS.
     * @param th_trdos Tetrahedron handler for transport DOS.
     * @param elvelocs Array of electronic velocities.
     * @param fsh Array to store Fermi surface harmonics.
     */
    void prepare_fsh(const tetrahedra::TetraHandler& th_dos, const tetrahedra::TetraHandler& th_trdos,
                    const arrays::array2D& elvelocs, arrays::array2D& fsh);

    /**
     * @brief Calculates the inner sum for a subarray at low temperatures.
     *
     * @param ik Index of the k-point.
     * @param iq Index of the q-point.
     * @param imd Index of the mode.
     * @param eigenens Array of eigenenergies.
     * @param eigenens_qk Array of eigenenergies at q+k.
     * @param matels Array of matrix elements.
     */
    void calc_inner_sum_in_subarray_lowt(size_t ik, size_t iq, size_t imd,
                                        arrays::array2D& eigenens,
                                        arrays::array2D& eigenens_qk,
                                        arrays::array3D& matels);

    /**
     * @brief Calculates the inner sum for a subarray using epsilon values.
     *
     * @param ik Index of the k-point.
     * @param iq Index of the q-point.
     * @param imd Index of the mode.
     * @param eigenens Array of eigenenergies.
     * @param eigenens_qk Array of eigenenergies at q+k.
     * @param matels Array of matrix elements.
     */
    void calc_inner_sum_in_subarray_epsilons(size_t ik, size_t iq, size_t imd,
                                            arrays::array2D& eigenens,
                                            arrays::array2D& eigenens_qk,
                                            arrays::array3D& matels);

    /**
     * @brief Calculates the inner sum for a specific q-point, mode, and epsilon index.
     *
     * @param iq Index of the q-point.
     * @param imd Index of the mode.
     * @param ieps Index of the epsilon value.
     * @return double The computed inner sum.
     */
    double calc_inner_sum(size_t iq, size_t imd, size_t ieps);

    /**
     * @brief Calculates the inner sum for a subarray.
     *
     * @param iq Index of the q-point.
     * @param imd Index of the mode.
     * @param ieps Index of the epsilon value.
     * @return double The computed inner sum.
     */
    double calc_inner_sum_in_subarray(size_t iq, size_t imd, size_t ieps);

    /**
     * @brief Calculates the inner sum for a subarray using tetrahedron integration.
     *
     * @param iq Index of the q-point.
     * @param imd Index of the mode.
     * @param ieps Index of the epsilon value.
     * @return double The computed inner sum.
     */
    double calc_inner_sum_in_subarray_tetra(size_t iq, size_t imd, size_t ieps);

    /**
     * @brief Calculates the external sum for a given frequency.
     *
     * @param Omega Frequency at which the external sum is calculated.
     * @return arrays::array1D The computed external sum values.
     */
    arrays::array1D calc_exter_sum(double Omega);

    // Private member variables
    KPprotocol kprot_; ///< Protocol for k-point sampling.
    KPprotocol qprot_; ///< Protocol for q-point sampling.

    int sign_;         ///< Sign parameter for velocities (+1 or -1).
    int sign_pr_;      ///< Primed counterpart of the sign parameter.

    size_t nbnd_;      ///< Number of bands.
    arrays::array2D eigenens_; ///< Eigenenergies.

    size_t nmds_;      ///< Number of phonon modes.
    arrays::array2D eigenfreqs_; ///< Eigenfrequencies.

    arrays::array2D elvelocs_alpha_; ///< Electronic velocities for alpha direction.
    arrays::array2D elvelocs_beta_;  ///< Electronic velocities for beta direction.
    arrays::array2D elvelocs_alpha_sq_; ///< Squared electronic velocities for alpha direction.
    arrays::array2D elvelocs_beta_sq_;  ///< Squared electronic velocities for beta direction.

    tetrahedra::TetraHandler th_dos_; ///< Tetrahedron handler for DOS.
    tetrahedra::TetraHandler th_trdos_alpha_; ///< Tetrahedron handler for transport DOS (alpha).
    tetrahedra::TetraHandler th_trdos_beta_;  ///< Tetrahedron handler for transport DOS (beta).

    arrays::array2D fsh_alpha_; ///< Fermi surface harmonics for alpha direction.
    arrays::array2D fsh_beta_;  ///< Fermi surface harmonics for beta direction.

    size_t low_band_{ 0 };      ///< Lower band index.
    size_t high_band_{ quantities::EigenValue::nbands - 1 }; ///< Upper band index.

    arrays::array1D epsilons_; ///< Array of epsilon values.
    arrays::array1D DOSes_;    ///< Array of density of states (DOS) values.
    arrays::array1D trDOSes_;  ///< Array of transport density of states (trDOS) values.

    arrays::array3D inner_sum_; ///< Inner sum array.
    bool is_full_{ false };     ///< Flag indicating if the calculation is complete.
    size_t iq_cont_{ 0 };       ///< Continuation index for q-points.
    size_t imd_cont_{ 0 };      ///< Continuation index for modes.

    double Te_{ 0.258 };        ///< Electronic temperature in eV (default corresponds to 3000 K).
    char alpha_{ 'x' };         ///< Cartesian index of the first velocity component.
    char beta_{ 'x' };          ///< Cartesian index of the second velocity component.

    bool is_tetra_{ false };    ///< Flag indicating if tetrahedron integration is enabled.
    bool is_continue_calc_{ false }; ///< Flag indicating if the calculation is continuing from a previous state.

    bzsampling::SamplingFunc elec_sampling_; ///< Sampling function for electrons.
    bzsampling::SamplingFunc phon_sampling_; ///< Sampling function for phonons.
    double elec_smearing_; ///< Smearing parameter for electrons.
    double phon_smearing_; ///< Smearing parameter for phonons.

    std::string type_of_el_smear_; ///< Type of electronic smearing.
    std::string type_of_ph_smear_; ///< Type of phononic smearing.

    arrays::array2D eigenens_qk_; ///< Eigenenergies at q+k.
    arrays::array2D elvelocs_qk_alpha_; ///< Electronic velocities at q+k for alpha direction.
    arrays::array2D elvelocs_qk_beta_;  ///< Electronic velocities at q+k for beta direction.
    arrays::array2D fsh_qk_alpha_; ///< Fermi surface harmonics at q+k for alpha direction.
    arrays::array2D fsh_qk_beta_;  ///< Fermi surface harmonics at q+k for beta direction.
};

/**
 * @brief Main driver function to calculate spectral functions.
 *
 * This function drives the calculation of spectral functions for a given range of frequencies.
 *
 * @param a2f `SpecFunc` object used for the calculation.
 * @param omegas Array of frequencies at which the spectral function is calculated.
 * @param fname File name where the results will be written.
 */
void calc_spec_func(SpecFunc& a2f, const arrays::array1D& omegas, const std::string& fname);

} // namespace spectral
} // namespace skies