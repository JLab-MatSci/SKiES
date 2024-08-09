#pragma once

#include <skies/common/ndimarrays.h>
#include <skies/lattices/kp_protocol.h>

namespace skies { namespace tetrahedra {

class TetraHandler;
class DoubleTetraHandler;

void evaluate_dos(const arrays::array1D& range);
double evaluate_dos(const arrays::array2D&, double value);
double evaluate_dos(const arrays::array2D&, const arrays::array2D&, double value, bool use_qprot = false);

void evaluate_phdos(const arrays::array1D& range);

void evaluate_trdos(const arrays::array1D& range);

class TetraHandler {
public:
    TetraHandler();

    // A and eigenens are assumed to have a form (nbnd x nkpt)
    template <typename Matels, typename EigenEns>
    TetraHandler(Matels&& A, EigenEns&& eigenens)
        : A_(std::forward<Matels>(A))
        , energies_(std::forward<EigenEns>(eigenens))
    {
        // if ((kprot_.nkpt() != A_[0].size()) || (kprot_.nkpt() != energies_[0].size()))
        //     throw std::runtime_error("Arrays of matrix elements A_k, energies e_k must contain kprot.nkpt elements");
    }

    double evaluate_dos_at_value(double value, bool use_qprot = false) const;

private:
    double evaluate_dos_at(size_t ik, size_t n, double value, bool use_qprot = false) const;

private:
    static KPprotocol kprot_;
    static KPprotocol qprot_;
    static bool phon_tag_;
public:
    static void set_kprot(const KPprotocol& kprot);
    static void set_qprot(const KPprotocol& qprot);
    static void set_phon_tag(bool phon_tag);
    static const KPprotocol& kprot();
    static const KPprotocol& qprot();
    static bool phon_tag();
    static bool is_initialized();

private:
    arrays::array2D A_;
    arrays::array2D energies_;

    std::vector<std::vector<size_t>>
    create_tetrahedra(const std::vector<size_t>& subcell) const;
};

/**
 * Class for handling doubly constrained BZ integration as proposed by P. B. Allen:
 *               physica status solidi (b). 1983. V.120. N.2. P.529-538.
 */
class DoubleTetraHandler {
public:
    DoubleTetraHandler(arrays::array3D&& A_glob,
                       arrays::array2D&& epsilons_glob,
                       arrays::array2D&& omegas_glob);
    double evaluate_dos_at(size_t ik, size_t n, size_t m, double E, double O) const;
    double evaluate_dos_at_values(double E, double O) const;

private:
    static KPprotocol kprot_;
public:
    static void set_kprot(const KPprotocol& kprot);
    static const KPprotocol& kprot();
    static bool is_initialized();

private:
    arrays::array3D A_glob_;
    arrays::array2D epsilons_glob_;
    arrays::array2D omegas_glob_;

    std::vector<std::vector<size_t>>
    create_tetrahedra(const std::vector<size_t>& subcell) const;

    double evaluate_integral_0(double EE, const arrays::array1D& E,
                               double OO, const arrays::array1D& O,
                               const arrays::array1D& A) const;
    double evaluate_integral_1(double EE, const arrays::array1D& E,
                               double OO, const arrays::array1D& O,
                               const arrays::array1D& A) const;
    double evaluate_integral_3(double EE, const arrays::array1D& E,
                               double OO, const arrays::array1D& O,
                               const arrays::array1D& A) const;
};

} // tetrahedra
} // skies