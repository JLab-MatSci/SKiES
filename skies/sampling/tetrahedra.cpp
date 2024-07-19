#include <cassert>
#include <numeric>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <skies/sampling/tetrahedra.h>
#include <skies/quantities/elvelocs.h>
#include <skies/quantities/eigenfreqs.h>

#include <skies/utils/mpi_wrapper.h>

namespace skies { namespace tetrahedra {

using namespace arrays;
using namespace quantities;

double evaluate_dos_at_value(const array1D& A,
                    const array1D& energies,
                    const KPprotocol& kprot,
                    double value)
{
    if ((kprot.nkpt != A.size()) || (kprot.nkpt != energies.size()))
        throw std::runtime_error("Arrays of matrix elements A_k, energies e_k must contain kprot.nkpt elements");
    
    double dos{ 0.0 };
    TetraHandler th(A, energies, kprot);
#ifdef SKIES_MPI
    int  rank = mpi::rank();
    auto rcounts_displs = mpi::prepare_rcounts_displs(kprot.nkpt);
    auto rcounts = rcounts_displs.first;
    auto displs  = rcounts_displs.second;
    auto count{ rcounts[rank] };
    auto displ{ displs[rank]  };

    double dos_part{ 0.0 };
    for (auto ik = displ; ik < displ + count; ++ik)
        dos_part += th.evaluate_dos_at(ik, value);
    MPI_Reduce(&dos_part, &dos, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dos, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    for (auto&& ik : kprot.range())
        dos += th.evaluate_dos_at(ik, value);
#endif
    return dos;
}

// assumed that A and energies have shape (nkpt x nbnd)
double evaluate_dos_at_value(const array2D& A,
                             const array2D& energies,
                             const KPprotocol& kprot,
                             double value)
{
    if ((kprot.nkpt != A.size()) || (kprot.nkpt != energies.size()))
        throw std::runtime_error("Arrays of matrix elements A_k, energies e_k must contain kprot.nkpt elements");
    double dos{ 0.0 };
    auto AT = transpose(A);
    auto eT = transpose(energies);
    for (size_t i = 0; i < AT.size(); ++i)
        dos += evaluate_dos_at_value(AT[i], eT[i], kprot, value);
    dos /= kprot.nkpt;
    return dos;
}

void evaluate_dos(const KPprotocol& kprot,
                  const arrays::array1D& range)
{
    std::ofstream os("EigenValueDOS.dat");
    os << std::right;
    os << std::setw(12) << "# Energy, eV";
    os << std::setw(25) << "EigenValue DOS";
    os << std::setw(9) << " [1 / eV]";
    os << std::endl;

    auto grid = kprot.grid;
    auto nkpt = grid.size();
    array2D weights(nkpt, array1D(EigenValueDrawable::nbands, 1));
    array2D energies(nkpt, array1D(EigenValueDrawable::nbands, 0.0));
    std::transform(grid.begin(), grid.end(), energies.begin(),
                    [] (auto&& k) { return EigenValueDrawable().interpolate_at(k); });

    array1D dos(range.size(), 0.0);
    std::transform(range.begin(), range.end(), dos.begin(), [&] (double v) {
        return evaluate_dos_at_value(weights, energies, kprot, v);
    });

    for (size_t i = 0; i < range.size(); ++i)
        os << std::setprecision(6) << std::setw(12) << range[i] << std::setw(34) << dos[i] << std::endl;
    os.close();
}

void evaluate_phdos(const KPprotocol& kprot,
                    const arrays::array1D& range)
{
    std::ofstream os("EigenFrequencyDOS.dat");
    os << std::right;
    os << std::setw(12) << "# Energy, eV";
    os << std::setw(25) << "EigenFrequency DOS";
    os << std::setw(9) << " [1 / eV]";
    os << std::endl;

    auto grid = kprot.grid;
    auto nkpt = grid.size();
    array2D weights(nkpt, array1D(EigenFrequencyDrawable::nmodes, 1));
    array2D energies(nkpt, array1D(EigenFrequencyDrawable::nmodes, 0.0));
    std::transform(grid.begin(), grid.end(), energies.begin(),
                    [] (auto&& k) { return EigenFrequencyDrawable().interpolate_at(k); });

    array1D dos(range.size(), 0.0);
    std::transform(range.begin(), range.end(), dos.begin(), [&] (double v) {
        return evaluate_dos_at_value(weights, energies, kprot, v);
    });

    for (size_t i = 0; i < range.size(); ++i)
        os << std::setprecision(6) << std::setw(12) << range[i] << std::setw(34) << dos[i] << std::endl;
    os.close();
}

void evaluate_trdos(const KPprotocol& kprot,
                    const arrays::array1D& range,
                    char cart)
{
    std::ofstream os("VelocitiesDOS.dat");
    os << std::right;
    os << std::setw(12) << "# Energy, eV";
    os << std::setw(25) << "Transport DOS";
    os << std::setw(24) << " [13.605685 * Ry bohr^2]";
    os << std::endl;

    auto grid = kprot.grid;
    auto nkpt = grid.size();
    array2D eigenvals(nkpt, array1D(EigenValue::nbands, 0.0));
    array2D velocs(nkpt, array1D(EigenValue::nbands, 0.0));
    array2D velocs_squared(nkpt, array1D(EigenValue::nbands, 0.0));

    std::transform(grid.begin(), grid.end(), eigenvals.begin(),
                [] (auto&& k) { return EigenValue::interpolate_at(k); });
    std::transform(grid.begin(), grid.end(), velocs.begin(),
                [cart] (auto&& k) { return VelocitiesDrawable(cart).interpolate_at(k); });
    std::transform(velocs.begin(), velocs.end(), velocs_squared.begin(),
                    [] (const array1D& v) {
                        auto squared_v = array1D(EigenValue::nbands, 0.0);
                        std::transform(v.begin(), v.end(), squared_v.begin(), [] (auto x) { return x * x; });
                        return squared_v;
    });

    array1D trDOSes(range.size(), 0.0);
    std::transform(range.begin(), range.end(), trDOSes.begin(), [&] (double v) {
        return evaluate_dos_at_value(velocs_squared, eigenvals, kprot, v);
    });

    for (size_t i = 0; i < range.size(); ++i)
        os << std::setprecision(6) << std::setw(12) << range[i] << std::setw(34) << trDOSes[i] << std::endl;
    os.close();
}

// it is assumed that energies and matrix elements are precalculated
// on a grid in kprot and A[0] corresponds to 0th k-point
TetraHandler::TetraHandler(const array1D& A,
                           const array1D& energies,
                           const KPprotocol& kprot)
    : A_(A)
    , energies_(energies)
    , kprot_(kprot)
    {}

double TetraHandler::evaluate_dos_at(size_t ik, double value) const
{
    auto subcell = kprot_.local_subcell(ik);
    auto tetrahedra = create_tetrahedra(subcell);
    assert(tetrahedra.size() == 6);
    double dos{ 0.0 };
    for (auto&& t : tetrahedra)
    {
        // c - array of 4 elements according to (B1) - (B4) from
        // Lambin P., Vigneron J. P. // Physical Review B. 1984. V.29. N.6. P.3430.
        array1D c(4, 0.0);
        array1D local_energies, local_matels;
        for (auto&& i : t)
        {
            local_matels.push_back(A_[i]);
            local_energies.push_back(energies_[i]);
        }
        assert(local_matels.size() == 4 && local_energies.size() == 4);
        std::sort(local_energies.begin(), local_energies.end());

        auto E1 = local_energies[0];
        auto E2 = local_energies[1];
        auto E3 = local_energies[2];
        auto E4 = local_energies[3];

        auto E = value;

        if ((E1 <= E) && (E <= E2))
        {
            if ((E1 == E2) || (E1 == E3) || (E1 == E4))
            {
                c = {0.0, 0.0, 0.0, 0.0};
            }
            else
            {
                c[0]  = (E2 - E)/(E2 - E1) + (E3 - E)/(E3 - E1) + (E4 - E)/(E4 - E1);
                c[0] *= (E - E1)*(E - E1)/(E4 - E1)/(E3 - E1)/(E2 - E1);
                c[1]  = (E - E1)*(E - E1)*(E - E1)/(E2 - E1)/(E2 - E1)/(E3 - E1)/(E4 - E1);
                c[2]  = (E - E1)*(E - E1)*(E - E1)/(E2 - E1)/(E3 - E1)/(E3 - E1)/(E4 - E1);
                c[3]  = (E - E1)*(E - E1)*(E - E1)/(E2 - E1)/(E3 - E1)/(E4 - E1)/(E4 - E1);
            }
        } else
        if ((E2 <= E) && (E <= E3))
        {
            if ((E1 == E3) || (E2 == E4) || (E2 == E3) || (E1 == E4))
            {
                c = {0.0, 0.0, 0.0, 0.0};
            }
            else
            {
                auto cc1 = (E3 - E)/(E3 - E1)/(E3 - E1);
                cc1 *= ((E3 - E)*(E - E2)/(E4 - E2)/(E3 - E2) + (E4 - E)*(E - E1)/(E4 - E1)/(E4 - E2) + (E3 - E)*(E - E1)/(E3 - E2)/(E4 - E1));
                auto cc2 = (E4 - E)/(E4 - E1)/(E4 - E1);
                cc2 *= ((E4 - E)*(E - E1)/(E4 - E2)/(E3 - E1) + (E4 - E)*(E - E2)/(E4 - E2)/(E3 - E2) + (E3 - E)*(E - E1)/(E3 - E1)/(E3 - E2));
                c[0] = 0.5 * (cc1 + cc2);

                auto cc3 = (E3 - E)/(E3 - E2)/(E3 - E2);
                cc3 *= ((E3 - E)*(E - E2)/(E4 - E2)/(E3 - E1) + (E4 - E)*(E - E2)/(E4 - E2)/(E4 - E1) + (E3 - E)*(E - E1)/(E3 - E1)/(E4 - E1));
                auto cc4 = (E4 - E)/(E4 - E2)/(E4 - E2);
                cc4 *= ((E3 - E)*(E - E2)/(E3 - E2)/(E3 - E1) + (E4 - E)*(E - E1)/(E4 - E1)/(E3 - E1) + (E4 - E)*(E - E2)/(E3 - E2)/(E4 - E1));
                c[1] = 0.5 * (cc3 + cc4);

                auto cc5 = (E - E2)/(E3 - E2)/(E3 - E2);
                cc5 *= ((E3 - E)*(E - E2)/(E4 - E2)/(E3 - E1) + (E4 - E)*(E - E2)/(E4 - E2)/(E4 - E1) + (E3 - E)*(E - E1)/(E3 - E1)/(E4 - E1));
                auto cc6 = (E - E1)/(E3 - E1)/(E3 - E1);
                cc6 *= ((E3 - E)*(E - E2)/(E4 - E2)/(E3 - E2) + (E4 - E)*(E - E1)/(E4 - E1)/(E4 - E2) + (E3 - E)*(E - E1)/(E3 - E2)/(E4 - E1));
                c[2] = 0.5 * (cc5 + cc6);

                auto cc7 = (E - E2)/(E4 - E2)/(E4 - E2);
                cc7 *= ((E3 - E)*(E - E2)/(E3 - E2)/(E3 - E1) + (E4 - E)*(E - E1)/(E4 - E1)/(E3 - E1) + (E4 - E)*(E - E2)/(E3 - E2)/(E4 - E1));
                auto cc8 = (E - E1)/(E4 - E1)/(E4 - E1);
                cc8 *= ((E4 - E)*(E - E1)/(E4 - E2)/(E3 - E1) + (E4 - E)*(E - E2)/(E4 - E2)/(E3 - E2) + (E3 - E)*(E - E1)/(E3 - E1)/(E3 - E2));
                c[3] = 0.5 * (cc7 + cc8);
            }
        } else
        if ((E3 <= E) && (E <= E4))
        {
            if ((E1 == E4) || (E2 == E4) || (E3 == E4))
            {
                c = {0.0, 0.0, 0.0, 0.0};
            }
            else
            {
                c[0]  = (E4 - E)*(E4 - E)*(E4 - E)/(E4 - E1)/(E4 - E1)/(E4 - E2)/(E4 - E3);
                c[1]  = (E4 - E)*(E4 - E)*(E4 - E)/(E4 - E1)/(E4 - E2)/(E4 - E2)/(E4 - E3);
                c[2]  = (E4 - E)*(E4 - E)*(E4 - E)/(E4 - E1)/(E4 - E2)/(E4 - E3)/(E4 - E3);
                c[3]  = (E - E3)/(E4 - E3) + (E - E2)/(E4 - E2) + (E - E1)/(E4 - E1);
                c[3] *= (E4 - E)*(E4 - E)/(E4 - E1)/(E4 - E2)/(E4 - E3);
            }
        }
        else
        {
            c = {0.0, 0.0, 0.0, 0.0};
        }

        double dos_part{ 0.0 };
        if ((E1 == E2) && (E1 == E3) && (E1 == E4) && (E == E1))
            dos_part = 0.25;
        else
            dos_part = std::inner_product(local_matels.begin(), local_matels.end(), c.begin(), 0.0);

        dos += dos_part;
    }
    
    if (dos < 1.0e-12) dos = 0.0;

    dos /= 6;
    return dos;
}

double evaluate_dos_at_values(const array1D& A_glob,
                    const array1D& epsilons_glob,
                    const array1D& omegas_glob,
                    const KPprotocol& kprot,
                    size_t start,
                    size_t finish,
                    double E,
                    double O)
{
    double dos{ 0.0 };
    DoubleTetraHandler dth(A_glob, epsilons_glob, omegas_glob, kprot);
    for (size_t ik = start; ik < finish; ++ik)
        dos += dth.evaluate_dos_at(ik, E, O);
    return dos;
}

// here it is assumed that shape of A is (nbnd x nbnd x nkpt)
// but espilons/omegas have shape (nkpt x nbnd) 
double evaluate_dos_at_values(const array3D& A_glob,
                              const array2D& epsilons_glob,
                              const array2D& omegas_glob,
                              const KPprotocol& kprot,
                              size_t start,
                              size_t finish,
                              double E,
                              double O)
{
    if ((kprot.nkpt != A_glob[0][0].size()) || (kprot.nkpt != epsilons_glob.size()) || (kprot.nkpt != omegas_glob.size()))
        throw std::runtime_error("Arrays of matrix elements A_k, epsilons e_k and omegas om_k must contain kprot.nkpt elements");
    double dos{ 0.0 };
    auto nmax = A_glob.size();
    auto mmax = A_glob[0].size();
    auto eTg = transpose(epsilons_glob);
    auto oTg = transpose(omegas_glob);
    for (size_t n = 0; n < nmax ; ++n)
        for (size_t m = 0; m < mmax; ++m)
            dos += evaluate_dos_at_values(A_glob[n][m], eTg[n], oTg[m], kprot, start, finish, E, O);
    return dos;
}

std::vector<std::vector<size_t>>
TetraHandler::create_tetrahedra(const std::vector<size_t>& sc) const
{
    std::vector<std::vector<size_t>> tetra(6);
    tetra[0] = { sc[0], sc[1], sc[2], sc[5] };
    tetra[1] = { sc[1], sc[2], sc[3], sc[5] };
    tetra[2] = { sc[2], sc[3], sc[5], sc[7] };
    tetra[3] = { sc[2], sc[5], sc[6], sc[7] };
    tetra[4] = { sc[2], sc[4], sc[5], sc[6] };
    tetra[5] = { sc[0], sc[2], sc[4], sc[5] };
    return tetra;
}

std::vector<std::vector<size_t>>
DoubleTetraHandler::create_tetrahedra(const std::vector<size_t>& sc) const
{
    std::vector<std::vector<size_t>> tetra(6);
    tetra[0] = { sc[0], sc[1], sc[2], sc[5] };
    tetra[1] = { sc[1], sc[2], sc[3], sc[5] };
    tetra[2] = { sc[2], sc[3], sc[5], sc[7] };
    tetra[3] = { sc[2], sc[5], sc[6], sc[7] };
    tetra[4] = { sc[2], sc[4], sc[5], sc[6] };
    tetra[5] = { sc[0], sc[2], sc[4], sc[5] };
    return tetra;
}

DoubleTetraHandler::DoubleTetraHandler(const array1D& A_glob,
                                       const array1D& epsilons_glob,
                                       const array1D& omegas_glob,
                                       const KPprotocol& kprot)
    : A_glob_(A_glob)
    , epsilons_glob_(epsilons_glob)
    , omegas_glob_(omegas_glob)
    , kprot_(kprot)
{}

double DoubleTetraHandler::evaluate_integral_0(double EE, const array1D& E,
                           double OO, const array1D& O,
                           const array1D& A) const
{
    double h{ 0.0 };
    double midpoint{ 0.0 };
    double f{ 0.0 }, g{ 0.0 };

    double a  = (EE - E[0]) / (OO - O[0]);
    array1D aa(4, 0.0), pp(4, 0.0);
    aa[1] = (E[1] - E[0]) / (O[1] - O[0]);
    aa[2] = (E[2] - E[0]) / (O[2] - O[0]);
    aa[3] = (E[3] - E[0]) / (O[3] - O[0]);
    pp[1] = (A[1] - A[0]) / (O[1] - O[0]);
    pp[2] = (A[2] - A[0]) / (O[2] - O[0]);
    pp[3] = (A[3] - A[0]) / (O[3] - O[0]);
    double as{ 0.0 }; double ps{ 0.0 };
    double am{ 0.0 }; double pm{ 0.0 };
    double al{ 0.0 }; double pl{ 0.0 };
    std::vector<size_t> ii{1, 2, 3};
    do {
        if ( (aa[ii[0]] <= aa[ii[1]]) && (aa[ii[1]] <= aa[ii[2]]) )
        {
            as = aa[ii[0]]; ps = pp[ii[0]];
            am = aa[ii[1]]; pm = pp[ii[1]];
            al = aa[ii[2]]; pl = pp[ii[2]];
            break;
        }
    } while (std::next_permutation(ii.begin(), ii.end()));
    assert((as <= am) && (am <= al));

    if ((as <= a) && (a <= am)) // case AI
    {
        if (!((as == am) || (as == al)))
        {
            h = 2*ps + (pm - ps)*(a - as)/(am - as) + (pl - ps)*(a - as)/(al - as);
            g = (a - as)/(am - as)/(al - as);
        }
    }
    else
    if ((am <= a) && (a <= al)) // case AII
    {
        if (!((am == al) || (as == al)))
        {
            h = 2*pl + (pm - pl)*(al - a)/(al - am) + (ps - pl)*(al - a)/(al - as);
            g = (al - a)/(al - am)/(al - as);
        }
    }
    midpoint = A[0] + 0.5 * (OO - O[0]) * h;
    f = (OO - O[0])/(O[1] - O[0])/(O[2] - O[0])/(O[3] - O[0]);
    return 6.0 * f * g * midpoint;
}

double DoubleTetraHandler::evaluate_integral_1(double EE, const array1D& E,
                           double OO, const array1D& O,
                           const array1D& A) const
{
    double h{ 0.0 };
    double midpoint{ 0.0 };
    double f{ 0.0 }, g{ 0.0 };

    double b  = (EE - E[1]) / (OO - O[1]);
    array1D bb(4, 0.0), qq(4, 0.0);
    bb[0] = (E[0] - E[1]) / (O[0] - O[1]);
    bb[2] = (E[2] - E[1]) / (O[2] - O[1]);
    bb[3] = (E[3] - E[1]) / (O[3] - O[1]);
    qq[0] = (A[0] - A[1]) / (O[0] - O[1]);
    qq[2] = (A[2] - A[1]) / (O[2] - O[1]);
    qq[3] = (A[3] - A[1]) / (O[3] - O[1]);
    double bs{ 0.0 }; double qs{ 0.0 };
    double bm{ 0.0 }; double qm{ 0.0 };
    double bl{ 0.0 }; double ql{ 0.0 };
    std::vector<size_t> ii{0, 2, 3};
    do {
        if ( (bb[ii[0]] <= bb[ii[1]]) && (bb[ii[1]] <= bb[ii[2]]) )
        {
            bs = bb[ii[0]]; qs = qq[ii[0]];
            bm = bb[ii[1]]; qm = qq[ii[1]];
            bl = bb[ii[2]]; ql = qq[ii[2]];
            break;
        }
    } while (std::next_permutation(ii.begin(), ii.end()));
    assert((bs <= bm) && (bm <= bl));

    if ((bs <= b) && (b <= bm)) // case BI
    {
        if (!((bm == bs) || (bl == bs)))
        {
            h = 2*qs + (qm - qs)*(b - bs)/(bm - bs) + (ql - qs)*(b - bs)/(bl - bs);
            g = (b - bs)/(bm - bs)/(bl - bs);
        }
    }
    else
    if ((bm <= b) && (b <= bl)) // case BII
    {
        if (!((bl == bm) || (bl == bs)))
        {
            h = 2*ql + (qm - ql)*(bl - b)/(bl - bm) + (qs - ql)*(bl - b)/(bl - bs);
            g = (bl - b)/(bl - bm)/(bl - bs);
        }
    }

    midpoint = A[1] + 0.5 * (OO - O[1]) * h;
    f = (OO - O[1])/(O[1] - O[0])/(O[2] - O[1])/(O[3] - O[1]);
    return 6.0 * f * g * midpoint;
}

double DoubleTetraHandler::evaluate_integral_3(double EE, const array1D& E,
                           double OO, const array1D& O,
                           const array1D& A) const
{
    double h{ 0.0 };
    double midpoint{ 0.0 };
    double f{ 0.0 }, g{ 0.0 };

    double c  = (EE - E[3]) / (OO - O[3]);
    array1D cc(4, 0.0), rr(4, 0.0);
    cc[0] = (E[0] - E[3]) / (O[0] - O[3]);
    cc[1] = (E[1] - E[3]) / (O[1] - O[3]);
    cc[2] = (E[2] - E[3]) / (O[2] - O[3]);
    rr[0] = (A[0] - A[3]) / (O[0] - O[3]);
    rr[1] = (A[1] - A[3]) / (O[1] - O[3]);
    rr[2] = (A[2] - A[3]) / (O[2] - O[3]);
    double cs{ 0.0 }; double rs{ 0.0 };
    double cm{ 0.0 }; double rm{ 0.0 };
    double cl{ 0.0 }; double rl{ 0.0 };
    std::vector<size_t> ii{0, 1, 2};
    do {
        if ( (cc[ii[0]] <= cc[ii[1]]) && (cc[ii[1]] <= cc[ii[2]]) )
        {
            cs = cc[ii[0]]; rs = rr[ii[0]];
            cm = cc[ii[1]]; rm = rr[ii[1]];
            cl = cc[ii[2]]; rl = rr[ii[2]];
            break;
        }
    } while (std::next_permutation(ii.begin(), ii.end()));
    assert((cs <= cm) && (cm <= cl));

    if ((cs <= c) && (c <= cm)) // case CI
    {
        if (!((cm == cs) || (cl == cs)))
        {
            h = 2*rs + (rm - rs)*(c - cs)/(cm - cs) + (rl - rs)*(c - cs)/(cl - cs);
            g = (c - cs)/(cm - cs)/(cl - cs);
        }
    }
    else
    if ((cm <= c) && (c <= cl)) // case CII
    {
        if (!((cl == cm) || (cl == cs)))
        {
            h = 2*rl + (rm - rl)*(cl - c)/(cl - cm) + (rs - rl)*(cl - c)/(cl - cs);
            g = (cl - c)/(cl - cm)/(cl - cs);
        }
    }

    midpoint = A[3] + 0.5 * (OO - O[3]) * h;
    f = (O[3] - OO)/(O[3] - O[0])/(O[3] - O[1])/(O[3] - O[2]);
    return 6.0 * f * g * midpoint;
}

double DoubleTetraHandler::evaluate_dos_at(size_t ik, double E, double O) const
{
    auto subcell = kprot_.local_subcell(ik); // global index ik!
    auto tetrahedra = create_tetrahedra(subcell);
    assert(tetrahedra.size() == 6);
    double dos{ 0.0 };
    for (auto&& t : tetrahedra)
    {
        array1D local_epsilons, local_omegas, local_matels;
        for (auto&& i : t) // global index i!
        {
            local_matels.push_back(A_glob_[i]);
            local_omegas.push_back(omegas_glob_[i]);
            local_epsilons.push_back(epsilons_glob_[i]);
        }
        assert((local_matels.size() == 4) && (local_epsilons.size() == 4) && (local_omegas.size() == 4));
        std::sort(local_omegas.begin(), local_omegas.end());

        auto O0 = local_omegas[0];
        auto O1 = local_omegas[1];
        auto O2 = local_omegas[2];
        auto O3 = local_omegas[3];

        if ((O0 == O) && (O1 == O) && (O2 == O) && (O3 == O))
        {
            dos += 0.25;
            continue;
        }

        if ((O0 < O) && (O <= O1)) // case A
        {
            if (!((O0 == O1) || (O0 == O2) || (O0 == O3) || (O0 == O)))
                dos += evaluate_integral_0(E, local_epsilons, O, local_omegas, local_matels);
        } else
        if ((O1 <= O) && (O <= O2)) // case B
        {
            if (!((O1 == O0) || (O1 == O2) || (O1 == O3) || (O1 == O) || (O0 == O1) || (O0 == O2) || (O0 == O3) || (O0 == O)))
                dos += evaluate_integral_0(E, local_epsilons, O, local_omegas, local_matels)
                    -  evaluate_integral_1(E, local_epsilons, O, local_omegas, local_matels);
        } else
        if ((O2 <= O) && (O <= O3)) // case C
        {
            if (!((O3 == O0) || (O3 == O1) || (O3 == O2) || (O3 == O)))
                dos += evaluate_integral_3(E, local_epsilons, O, local_omegas, local_matels);
        }
    }
    
    if (dos < 1.0e-12) dos = 0.0;

    dos /= 6;
    return dos;
}

} // tetrahedra
} // skies