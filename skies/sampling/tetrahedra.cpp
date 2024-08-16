#include <cassert>
#include <numeric>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <skies/sampling/tetrahedra.h>
#include <skies/quantities/elvelocs.h>
#include <skies/quantities/eigenfreqs.h>

#ifdef SKIES_MPI
#include <skies/utils/mpi_wrapper.h>
#endif

#include <skies/utils/tbb_wrapper.h>

namespace skies { namespace tetrahedra {

using namespace arrays;
using namespace quantities;

KPprotocol TetraHandler::kprot_;
KPprotocol TetraHandler::qprot_;
bool TetraHandler::phon_tag_ = false;
KPprotocol DoubleTetraHandler::kprot_;

TetraHandler::TetraHandler() {}

void TetraHandler::set_kprot(const KPprotocol& kprot)
{
    kprot_ = kprot;
}

void TetraHandler::set_qprot(const KPprotocol& qprot)
{
    qprot_ = qprot;
}

void TetraHandler::set_phon_tag(bool phon_tag)
{
    phon_tag_ = phon_tag;
}

const KPprotocol& TetraHandler::kprot()
{
    return kprot_;
}

const KPprotocol& TetraHandler::qprot()
{
    return qprot_;
}

bool TetraHandler::phon_tag()
{
    return phon_tag_;
}

bool TetraHandler::is_initialized()
{
    bool kprot_ok = kprot_.nkpt() > 0;
    if (!phon_tag_)
        return kprot_ok;
    bool qprot_ok = qprot_.nkpt() > 0;
    return kprot_ok && qprot_ok;
}

void DoubleTetraHandler::set_kprot(const KPprotocol& kprot)
{
    kprot_ = kprot;
}

const KPprotocol& DoubleTetraHandler::kprot()
{
    return kprot_;
}

bool DoubleTetraHandler::is_initialized()
{
    return kprot_.nkpt() > 0;
}

double evaluate_dos(const array2D& eigenens, double value)
{
    assert(TetraHandler::is_initialized());
    auto grid = TetraHandler::kprot().grid();
    array2D weights(EigenValueDrawable::nbands, array1D(grid.size(), 1));
    TetraHandler th(std::move(weights), transpose(eigenens));
    return th.evaluate_dos_at_value(value);
}

double evaluate_dos(const array2D& weights,
                    const array2D& eigenens,
                    double value,
                    bool use_qprot)
{
    assert(TetraHandler::is_initialized());
    TetraHandler th(weights, transpose(eigenens));
    return th.evaluate_dos_at_value(value, use_qprot);
}

double TetraHandler::evaluate_dos_at_value(double value, bool use_qprot) const
{
    std::vector<size_t> ikpts;
    if (use_qprot)
        ikpts = qprot_.range();
    else
        ikpts = kprot_.range();
    double dos = std::transform_reduce(PAR ikpts.begin(), ikpts.end(), 0.0, std::plus<double>(),
        [&] (auto&& ik) -> double {
            std::vector<size_t> ibands(A_.size());
            std::iota(ibands.begin(), ibands.end(), 0);
            return std::transform_reduce(ibands.begin(), ibands.end(), 0.0, std::plus<double>(),
                [&] (auto&& n) -> double {
                    return evaluate_dos_at(ik, n, value, use_qprot);
            });
        }
    );
    dos /= kprot_.nkpt();
    return dos;
}

void evaluate_dos(const arrays::array1D& range)
{
    std::ofstream os("EigenValueDOS.dat");
    os << std::right;
    os << std::setw(12) << "# Energy, eV";
    os << std::setw(25) << "EigenValue DOS";
    os << std::setw(9) << " [1/eV/spin/cell]";
    os << std::endl;

    assert(TetraHandler::is_initialized());
    auto grid = TetraHandler::kprot().grid();
    auto nkpt = grid.size();
    array2D weights(EigenValueDrawable::nbands, array1D(nkpt, 1));
    array2D energies(nkpt, array1D(EigenValueDrawable::nbands, 0.0));
    std::transform(grid.begin(), grid.end(), energies.begin(),  [] (auto&& k) {
        return EigenValueDrawable().interpolate_at(k);
    });

    array1D dos(range.size(), 0.0);
    TetraHandler th(std::move(weights), transpose(energies));
    std::transform(range.begin(), range.end(), dos.begin(), [&] (double v) {
        return th.evaluate_dos_at_value(v);
    });

    int rank{ 0 };
#ifdef SKIES_MPI
    rank = skies::mpi::rank();
#endif
    if (!rank)
    {
        for (size_t i = 0; i < range.size(); ++i)
            os << std::setprecision(6) << std::setw(12) << range[i] << std::setw(34) << dos[i] << std::endl;
        os.close();
    }
}

void evaluate_phdos(const arrays::array1D& range)
{
    std::ofstream os("EigenFrequencyDOS.dat");
    os << std::right;
    os << std::setw(12) << "# Energy, eV";
    os << std::setw(25) << "EigenFrequency DOS";
    os << std::setw(9) << " [1 / eV]";
    os << std::endl;

    assert(TetraHandler::is_initialized());
    auto grid = TetraHandler::kprot().grid();
    auto nkpt = grid.size();
    array2D weights(EigenFrequencyDrawable::nmodes, array1D(nkpt, 1));
    array2D energies(nkpt, array1D(EigenFrequencyDrawable::nmodes, 0.0));
    std::transform(grid.begin(), grid.end(), energies.begin(), [] (auto&& k) {
        return EigenFrequencyDrawable().interpolate_at(k); // CAUTION! threads
    });

    array1D dos(range.size(), 0.0);
    TetraHandler th(std::move(weights), transpose(energies));
    std::transform(range.begin(), range.end(), dos.begin(), [&] (double v) {
        return th.evaluate_dos_at_value(v);
    });

    int rank{ 0 };
#ifdef SKIES_MPI
    rank = skies::mpi::rank();
#endif
    if (!rank)
    {
        for (size_t i = 0; i < range.size(); ++i)
            os << std::setprecision(6) << std::setw(12) << range[i] << std::setw(34) << dos[i] << std::endl;
        os.close();
    }
}

void evaluate_trdos(const arrays::array1D& range)
{
    std::ofstream os("VelocitiesDOS.dat");
    os << std::right;
    os << std::setw(12) << "# Energy, eV";
    os << std::setw(25) << "Transport DOS";
    os << std::setw(24) << " [r.a.u.]";
    os << std::endl;

    assert(TetraHandler::is_initialized());
    auto grid = TetraHandler::kprot().grid();
    auto nkpt = grid.size();
    array2D eigenens(nkpt, array1D(EigenValue::nbands, 0.0));
    std::transform(PAR grid.begin(), grid.end(), eigenens.begin(), [] (auto&& k) {
        return EigenValue::interpolate_at(k);
    });

    auto prepare_velocs = [&] (char cart, array2D& elvelocs, array2D& elvelocs_sq) {
        elvelocs.resize(nkpt, array1D(EigenValue::nbands, 0.0));
        elvelocs_sq.resize(nkpt, array1D(EigenValue::nbands, 0.0));
        std::transform(PAR grid.begin(), grid.end(), elvelocs.begin(),
                    [&] (auto&& k) { return Velocities(cart).interpolate_at(k); });
        std::transform(PAR elvelocs.begin(), elvelocs.end(), elvelocs_sq.begin(),
                [] (auto&& v) {
                    auto squared_v = array1D(EigenValue::nbands, 0.0);
                    std::transform(v.begin(), v.end(), squared_v.begin(), [] (auto&& x) { return x * x; });
                    return squared_v;
        });
    };

    array2D elvelocs_x, elvelocs_y, elvelocs_z;
    array2D elvelocs_x_sq, elvelocs_y_sq, elvelocs_z_sq;
    prepare_velocs('x', elvelocs_x, elvelocs_x_sq);
    prepare_velocs('y', elvelocs_y, elvelocs_y_sq);
    prepare_velocs('z', elvelocs_z, elvelocs_z_sq);

    array1D trDOSes(range.size());
    array1D trDOSes_x(range.size()), trDOSes_y(range.size()), trDOSes_z(range.size());

    tetrahedra::TetraHandler th_x(transpose(elvelocs_x_sq), transpose(eigenens));
    std::transform(range.begin(), range.end(), trDOSes_x.begin(), [th = std::move(th_x)] (auto&& eps) {
        return th.evaluate_dos_at_value(eps);
    });
    tetrahedra::TetraHandler th_y(transpose(elvelocs_y_sq), transpose(eigenens));
    std::transform(range.begin(), range.end(), trDOSes_y.begin(), [th = std::move(th_y)] (auto&& eps) {
        return th.evaluate_dos_at_value(eps);
    });
    tetrahedra::TetraHandler th_z(transpose(elvelocs_z_sq), transpose(eigenens));
    std::transform(range.begin(), range.end(), trDOSes_z.begin(), [th = std::move(th_z)] (auto&& eps) {
        return th.evaluate_dos_at_value(eps);
    });
    trDOSes = (trDOSes_x + trDOSes_y + trDOSes_z) * (1.0 / 3.0);
    trDOSes = trDOSes * units::Ry_in_eV; // go to [r.a.u.]

    int rank{ 0 };
#ifdef SKIES_MPI
    rank = skies::mpi::rank();
#endif
    if (!rank)
    {
        for (size_t i = 0; i < range.size(); ++i)
            os << std::setprecision(6) << std::setw(12) << range[i] << std::setw(34) << trDOSes[i] << std::endl;
        os.close();
    }
}

double TetraHandler::evaluate_dos_at(size_t ik, size_t n, double value, bool use_qprot) const
{
    std::vector<size_t> subcell;
    if (use_qprot)
        subcell = qprot_.local_subcell(ik); // global index ik
    else
        subcell = kprot_.local_subcell(ik); // global index ik
    auto tetrahedra = create_tetrahedra(subcell);
    assert(tetrahedra.size() == 6);
    double dos{ 0.0 };
    for (auto&& t : tetrahedra)
    {
        // c - array of 4 elements according to (B1) - (B4) from
        // Lambin P., Vigneron J. P. // Physical Review B. 1984. V.29. N.6. P.3430.
        array1D c(4, 0.0);
        array1D local_energies(4, 0.0), local_matels(4, 0.0);
        for (size_t i = 0; i < 4; ++i)
        {
            auto iik = t[i]; // global index iik!
            local_matels[i] = A_[n][iik];
            local_energies[i] = energies_[n][iik];
        }
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

double DoubleTetraHandler::evaluate_dos_at_values(double E, double O) const
{
    auto nmax = A_glob_.size();
    auto mmax = A_glob_[0].size();
    auto ikpts = kprot_.range();
    double dos = std::transform_reduce(PAR ikpts.begin(), ikpts.end(), 0.0, std::plus<double>(),
        [&] (auto&& ik) -> double {
            std::vector<size_t> ibands(nmax);
            std::iota(ibands.begin(), ibands.end(), 0);
            return std::transform_reduce(ibands.begin(), ibands.end(), 0.0, std::plus<double>(),
                [&] (auto&& n) -> double {
                    std::vector<size_t> jbands(mmax);
                    std::iota(jbands.begin(), jbands.end(), 0);
                    return std::transform_reduce(jbands.begin(), jbands.end(), 0.0, std::plus<double>(),
                        [&] (auto&& m) -> double {
                            return evaluate_dos_at(ik, n, m, E, O);
                    });
            });
        }
    );
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

DoubleTetraHandler::DoubleTetraHandler(array3D&& A_glob,
                                       array2D&& epsilons_glob,
                                       array2D&& omegas_glob)
    : A_glob_(std::move(A_glob))
    , epsilons_glob_(std::move(epsilons_glob))
    , omegas_glob_(std::move(omegas_glob))
{
    if ((kprot_.nkpt() != A_glob_[0][0].size()) || (kprot_.nkpt() != epsilons_glob_.size()) || (kprot_.nkpt() != omegas_glob_.size()))
        throw std::runtime_error("Arrays of matrix elements A_k, epsilons e_k and omegas om_k must contain kprot.nkpt elements");
}

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

double DoubleTetraHandler::evaluate_dos_at(size_t ik, size_t n, size_t m, double E, double O) const
{
    auto subcell = kprot_.local_subcell(ik); // global index ik!
    auto tetrahedra = create_tetrahedra(subcell);
    assert(tetrahedra.size() == 6);
    double dos{ 0.0 };
    for (auto&& t : tetrahedra)
    {
        array1D local_epsilons(4, 0.0), local_omegas(4, 0.0), local_matels(4, 0.0);
        for (size_t i = 0; i < 4; ++i)
        {
            auto&& iik = t[i]; // global index iik!
            local_matels[i] = A_glob_[n][m][iik];
            local_omegas[i] = omegas_glob_[iik][m];
            local_epsilons[i] = epsilons_glob_[iik][n];
        }
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