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
#include <skies/lattices/kp_protocol.h>

namespace skies { namespace tetrahedra {

/**
 * @brief Class for handling single tetrahedron-based Brillouin zone integration.
 *
 * The `TetraHandler` class provides functionality for calculating density of states (DOS) using tetrahedron methods.
 * It supports both k-point and q-point protocols and allows for customization of integration parameters.
 */
class TetraHandler {
public:
    /**
     * @brief Default constructor for the `TetraHandler` class.
     */
    TetraHandler();

    /**
     * @brief Constructs a `TetraHandler` object with matrix elements and eigenenergies.
     *
     * @tparam Matels Type of the matrix elements array.
     * @tparam EigenEns Type of the eigenenergies array.
     * @param A Matrix elements array (assumed to have dimensions: nbnd x nkpt).
     * @param eigenens Eigenenergies array (assumed to have dimensions: nbnd x nkpt).
     */
    template <typename Matels, typename EigenEns>
    TetraHandler(Matels&& A, EigenEns&& eigenens)
        : A_(std::forward<Matels>(A))
        , energies_(std::forward<EigenEns>(eigenens))
    {
        // Validation logic can be added here if needed.
    }

    /**
     * @brief Evaluates the density of states (DOS) at a specific value.
     *
     * @param value Energy value at which the DOS is evaluated.
     * @param use_qprot Flag indicating whether to use the q-point protocol (default: false).
     * @return double The computed DOS value.
     */
    double evaluate_dos_at_value(double value, bool use_qprot = false) const;
    double evaluate_nos_at_value(double value, bool use_qprot = false) const;

private:
    /**
     * @brief Evaluates the density of states (DOS) at a specific k-point and band index.
     *
     * @param ik Index of the k-point.
     * @param n Band index.
     * @param value Energy value at which the DOS is evaluated.
     * @param use_qprot Flag indicating whether to use the q-point protocol (default: false).
     * @return double The computed DOS value.
     */
    double evaluate_dos_at(size_t ik, size_t n, double value, bool use_qprot = false) const;

private:
    static KPprotocol kprot_; ///< Static k-point protocol.
    static KPprotocol qprot_; ///< Static q-point protocol.
    static bool phon_tag_;    ///< Static flag indicating phonon-related calculations.

public:
    /**
     * @brief Sets the k-point protocol for all `TetraHandler` instances.
     *
     * @param kprot The k-point protocol to set.
     */
    static void set_kprot(const KPprotocol& kprot);

    /**
     * @brief Sets the q-point protocol for all `TetraHandler` instances.
     *
     * @param qprot The q-point protocol to set.
     */
    static void set_qprot(const KPprotocol& qprot);

    /**
     * @brief Sets the phonon tag for all `TetraHandler` instances.
     *
     * @param phon_tag The phonon tag to set.
     */
    static void set_phon_tag(bool phon_tag);

    /**
     * @brief Returns the current k-point protocol.
     *
     * @return const KPprotocol& Reference to the k-point protocol.
     */
    static const KPprotocol& kprot();

    /**
     * @brief Returns the current q-point protocol.
     *
     * @return const KPprotocol& Reference to the q-point protocol.
     */
    static const KPprotocol& qprot();

    /**
     * @brief Checks if the `TetraHandler` class has been initialized.
     *
     * @return bool True if initialized, false otherwise.
     */
    static bool is_initialized();

    /**
     * @brief Checks the phonon tag status.
     *
     * @return bool True if phonon-related calculations are enabled, false otherwise.
     */
    static bool phon_tag();

private:
    arrays::array2D A_;       ///< Matrix elements array (dimensions: nbnd x nkpt).
    arrays::array2D energies_; ///< Eigenenergies array (dimensions: nbnd x nkpt).

    /**
     * @brief Creates tetrahedra for integration within a subcell.
     *
     * @param subcell Subcell indices defining the integration region.
     * @return std::vector<std::vector<size_t>> List of tetrahedra indices.
     */
    std::vector<std::vector<size_t>>
    create_tetrahedra(const std::vector<size_t>& subcell) const;
};

/**
 * @brief Class for handling doubly constrained tetrahedron-based Brillouin zone integration.
 *
 * The `DoubleTetraHandler` class implements the doubly constrained BZ integration method proposed by P. B. Allen
 * (physica status solidi (b), 1983, V.120, N.2, P.529-538). It is used for evaluating joint density of states (JDOS)
 * and related quantities.
 */
class DoubleTetraHandler {
public:
    /**
     * @brief Constructs a `DoubleTetraHandler` object with global matrix elements, eigenenergies, and frequencies.
     *
     * @param A_glob Global matrix elements array (dimensions: nbnd x nkpt x nmodes).
     * @param epsilons_glob Global eigenenergies array (dimensions: nbnd x nkpt).
     * @param omegas_glob Global frequencies array (dimensions: nmodes x nkpt).
     */
    DoubleTetraHandler(arrays::array3D&& A_glob,
                    arrays::array2D&& epsilons_glob,
                    arrays::array2D&& omegas_glob);

    /**
     * @brief Evaluates the density of states (DOS) at specific indices and values.
     *
     * @param ik Index of the k-point.
     * @param n Band index.
     * @param m Mode index.
     * @param E Energy value.
     * @param O Frequency value.
     * @return double The computed DOS value.
     */
    double evaluate_dos_at(size_t ik, size_t n, size_t m, double E, double O) const;

    /**
     * @brief Evaluates the density of states (DOS) at specific energy and frequency values.
     *
     * @param E Energy value.
     * @param O Frequency value.
     * @return double The computed DOS value.
     */
    double evaluate_dos_at_values(double E, double O) const;

private:
    static KPprotocol kprot_; ///< Static k-point protocol.

public:
    /**
     * @brief Sets the k-point protocol for all `DoubleTetraHandler` instances.
     *
     * @param kprot The k-point protocol to set.
     */
    static void set_kprot(const KPprotocol& kprot);

    /**
     * @brief Returns the current k-point protocol.
     *
     * @return const KPprotocol& Reference to the k-point protocol.
     */
    static const KPprotocol& kprot();

    /**
     * @brief Checks if the `DoubleTetraHandler` class has been initialized.
     *
     * @return bool True if initialized, false otherwise.
     */
    static bool is_initialized();

private:
    arrays::array3D A_glob_;      ///< Global matrix elements array (dimensions: nbnd x nkpt x nmodes).
    arrays::array2D epsilons_glob_; ///< Global eigenenergies array (dimensions: nbnd x nkpt).
    arrays::array2D omegas_glob_;   ///< Global frequencies array (dimensions: nmodes x nkpt).

    /**
     * @brief Creates tetrahedra for integration within a subcell.
     *
     * @param subcell Subcell indices defining the integration region.
     * @return std::vector<std::vector<size_t>> List of tetrahedra indices.
     */
    std::vector<std::vector<size_t>>
    create_tetrahedra(const std::vector<size_t>& subcell) const;

    /**
     * @brief Evaluates the integral for case 0.
     *
     * @param EE Energy value.
     * @param E Array of energy values.
     * @param OO Frequency value.
     * @param O Array of frequency values.
     * @param A Array of matrix elements.
     * @return double The computed integral value.
     */
    double evaluate_integral_0(double EE, const arrays::array1D& E,
                            double OO, const arrays::array1D& O,
                            const arrays::array1D& A) const;

    /**
     * @brief Evaluates the integral for case 1.
     *
     * @param EE Energy value.
     * @param E Array of energy values.
     * @param OO Frequency value.
     * @param O Array of frequency values.
     * @param A Array of matrix elements.
     * @return double The computed integral value.
     */
    double evaluate_integral_1(double EE, const arrays::array1D& E,
                            double OO, const arrays::array1D& O,
                            const arrays::array1D& A) const;

    /**
     * @brief Evaluates the integral for case 3.
     *
     * @param EE Energy value.
     * @param E Array of energy values.
     * @param OO Frequency value.
     * @param O Array of frequency values.
     * @param A Array of matrix elements.
     * @return double The computed integral value.
     */
    double evaluate_integral_3(double EE, const arrays::array1D& E,
                            double OO, const arrays::array1D& O,
                            const arrays::array1D& A) const;
};

/**
 * @brief Evaluates the density of states (DOS) over a range of values.
 *
 * @param range Array of energy or frequency values.
 */
void evaluate_dos(const arrays::array1D& range);

/**
 * @brief Evaluates the density of states (DOS) at a specific value.
 *
 * @param eigenens Eigenenergies array.
 * @param value Energy or frequency value.
 * @return double The computed DOS value.
 */
double evaluate_dos(const arrays::array2D& eigenens, double value);

/**
 * @brief Evaluates the density of states (DOS) at a specific value with optional q-point protocol.
 *
 * @param eigenens Eigenenergies array.
 * @param range Array of energy or frequency values.
 * @param value Energy or frequency value.
 * @param use_qprot Flag indicating whether to use the q-point protocol (default: false).
 * @return double The computed DOS value.
 */
double evaluate_dos(const arrays::array2D& eigenens, const arrays::array2D& range, double value, bool use_qprot = false);

/**
 * @brief Evaluates the phonon density of states (phDOS) over a range of values.
 *
 * @param range Array of frequency values.
 */
void evaluate_phdos(const arrays::array1D& range);

/**
 * @brief Evaluates the transport density of states (trDOS) over a range of values.
 *
 * @param range Array of energy values.
 * @param alpha Cartesian index of the velocity component.
 */
void evaluate_trdos(const arrays::array1D& range, char alpha);

} // namespace tetrahedra
} // namespace skies