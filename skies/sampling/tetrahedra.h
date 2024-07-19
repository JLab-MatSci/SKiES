#pragma once

#include <skies/common/ndimarrays.h>
#include <skies/lattices/kp_protocol.h>

namespace skies { namespace tetrahedra {

void evaluate_dos(const KPprotocol& kprot,
                  const arrays::array1D& range);

void evaluate_phdos(const KPprotocol& kprot,
                  const arrays::array1D& range);

void evaluate_trdos(const KPprotocol& kprot,
                    const arrays::array1D& range,
                    char cart);

double evaluate_dos_at_value(const arrays::array1D& A,
                    const arrays::array1D& energies,
                    const KPprotocol& kprot,
                    double value);

double evaluate_dos_at_value(const arrays::array2D& A,
                    const arrays::array2D& energies,
                    const KPprotocol& kprot,
                    double value);

double evaluate_dos_at_values(const arrays::array1D& A,
                    const arrays::array1D& epsilons,
                    const arrays::array1D& omegas,
                    const KPprotocol& kprot,
                    size_t start,
                    size_t finish,
                    double E,
                    double O);

double evaluate_dos_at_values(const arrays::array3D& A,
                              const arrays::array2D& epsilons,
                              const arrays::array2D& omegas,
                              const KPprotocol& kprot,
                              size_t start,
                              size_t finish,
                              double E,
                              double O);

class TetraHandler {
public:
    TetraHandler(const arrays::array1D& A, const arrays::array1D& energies, const KPprotocol& kprot);
    double evaluate_dos_at(size_t ik, double value) const;

private:
    arrays::array1D A_;
    arrays::array1D energies_;
    KPprotocol kprot_;

    std::vector<std::vector<size_t>>
    create_tetrahedra(const std::vector<size_t>& subcell) const;
};

/**
 * Class for handling doubly constrained BZ integration as proposed by P. B. Allen:
 *               physica status solidi (b). 1983. V.120. N.2. P.529-538.
 */
class DoubleTetraHandler {
public:
    DoubleTetraHandler(const arrays::array1D& A_glob,
                       const arrays::array1D& epsilons_glob,
                       const arrays::array1D& omegas_glob,
                       const KPprotocol& kprot);
    double evaluate_dos_at(size_t ik, double E, double O) const;

private:
    arrays::array1D A_glob_;
    arrays::array1D epsilons_glob_;
    arrays::array1D omegas_glob_;
    KPprotocol kprot_;

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