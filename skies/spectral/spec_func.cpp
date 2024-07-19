#include <map>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <iostream>
#include <cstring>
#include <cstdlib>

#include <skies/common/alg.h>
#include <skies/quantities/dos.h>
#include <skies/quantities/gmatrix.h>
#include <skies/spectral/spec_func.h>

#include <skies/sampling/tetrahedra.h>

#include <launch/timer.h>

#ifdef SKIES_MPI
#include <skies/utils/mpi_wrapper.h>
#endif

using namespace skies::arrays;
using namespace skies::bzsampling;
using namespace skies::quantities;

#ifdef SKIES_MPI
using namespace skies::mpi;
#endif

namespace skies { namespace spectral {

void SpecFunc::init(const std::vector<size_t>& kpgrid,
                    const std::vector<size_t>& qpgrid,
                    bzsampling::SamplingFunc elec_sampling,
                    bzsampling::SamplingFunc phon_sampling,
                    double elec_smearing,
                    double phon_smearing
) {
    nkx_ = kpgrid[0]; nqx_ = qpgrid[0];
    nky_ = kpgrid[1]; nqy_ = qpgrid[1];
    nkz_ = kpgrid[2]; nqz_ = qpgrid[2];

    kpts_ = KPprotocol(nkx_, nky_, nkz_).grid;
    qpts_ = KPprotocol(nqx_, nqy_, nqz_).grid;

    nqpt_ = qpts_.size();
    nkpt_ = kpts_.size();

    nbnd_ = EigenValue::nbands;
    nmds_ = EigenFrequency::nmodes;

    low_band_ = 0;
    high_band_ = nbnd_ - 1;

    eigenens_.resize(nkpt_, array1D(nbnd_));
    elvelocs_.resize(nkpt_, array1D(nbnd_));
    eigenfreqs_.resize(nqpt_, array1D(nmds_));

    std::transform(kpts_.begin(), kpts_.end(), eigenens_.begin(),
                    [] (auto&& k) { return EigenValue::interpolate_at(k); });
    std::transform(kpts_.begin(), kpts_.end(), elvelocs_.begin(),
                    [this] (auto&& k) { return Velocities(cart_).interpolate_at(k); });
    std::transform(qpts_.begin(), qpts_.end(), eigenfreqs_.begin(),
                    [] (auto&& q) { return EigenFrequencyDrawable().interpolate_at(q); });

    elec_sampling_ = elec_sampling;
    phon_sampling_ = phon_sampling;
    elec_smearing_ = elec_smearing;
    phon_smearing_ = phon_smearing;
}

SpecFunc::SpecFunc(const std::vector<size_t>& kpgrid,
                   const std::vector<size_t>& qpgrid,
                   bzsampling::SamplingFunc elec_sampling,
                   bzsampling::SamplingFunc phon_sampling,
                   double elec_smearing,
                   double phon_smearing,
                   int sign,
                   double Te,
                   char cart,
                   bool is_tetra
)
    : sign_(sign)
    , sign_pr_(sign)
    , Te_(Te)
    , cart_(cart)
    , is_tetra_(is_tetra)
{
    init(kpgrid, qpgrid, elec_sampling, phon_sampling, elec_smearing, phon_smearing);
    epsilons_.push_back(0);

    launch::Timer t;
    t.start("========= Started transport DOS calculations...");
    array2D velocs_squared(nkpt_, array1D(EigenValue::nbands));
    std::transform(elvelocs_.begin(), elvelocs_.end(), velocs_squared.begin(),
                    [] (const array1D& v) {
                        auto squared_v = array1D(EigenValue::nbands, 0.0);
                        std::transform(v.begin(), v.end(), squared_v.begin(), [] (auto x) { return x * x; });
                        return squared_v;               
    });
    double trDOS{ 0.0 };
    if (is_tetra_)
    {
        KPprotocol kprot(nkx_, nky_, nkz_);
        trDOS = tetrahedra::evaluate_dos_at_value(velocs_squared, eigenens_, kprot, 0);
    }
    else
        trDOS = evaluate_smeared_trdos_at_value(0, elec_smearing_, elec_sampling_, eigenens_, velocs_squared, Te_);
    trDOSes_.push_back(trDOS);
    dump_trdos_file();

    t.stop("========= Transport DOS is evaluated");
    t.print_elapsed("\t  Transport DOS evaluation time: ");
}

SpecFunc::SpecFunc(const std::vector<size_t>& kpgrid,
               const std::vector<size_t>& qpgrid,
               const arrays::array1D& epsilons,
               bzsampling::SamplingFunc elec_sampling,
               bzsampling::SamplingFunc phon_sampling,
               double elec_smearing,
               double phon_smearing,
               int sign,
               int sign_pr,
               double Te,
               char cart,
               bool is_tetra
)
    : sign_(sign)
    , sign_pr_(sign_pr)
    , epsilons_(epsilons)
    , Te_(Te)
    , cart_(cart)
    , is_tetra_(is_tetra)
{
    init(kpgrid, qpgrid, elec_sampling, phon_sampling, elec_smearing, phon_smearing);
    trDOSes_.resize(epsilons_.size(), 0.0);
// #ifdef SKIES_MPI
//     auto trDOSes_tmp = array1D(epsilons_.size(), 0.0);
// #endif

    launch::Timer t;
    t.start("\t  Started transport DOS calculations...");

    array2D velocs_squared(nkpt_, array1D(EigenValue::nbands));
    std::transform(elvelocs_.begin(), elvelocs_.end(), velocs_squared.begin(),
                    [] (const array1D& v) {
                        auto squared_v = array1D(EigenValue::nbands, 0.0);
                        std::transform(v.begin(), v.end(), squared_v.begin(), [] (auto x) { return x * x; });
                        return squared_v;
    });
// #ifdef SKIES_MPI
//     auto rcounts = mpi::prepare_rcounts_displs(epsilons_.size()).first;
//     auto displs  = mpi::prepare_rcounts_displs(epsilons_.size()).second;

//     int rank = skies::mpi::rank();
//     for (int i = displs[rank]; i < displs[rank] + rcounts[rank]; ++i)
//     {
//         double e = epsilons_[i];
//         trDOSes_tmp[i] = evaluate_smeared_trdos_at_value(e, elec_smearing_, elec_sampling_, eigenens_, velocs_squared, Te_);
//     }

//     MPI_Reduce(trDOSes_tmp.data(), trDOSes_.data(), epsilons_.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//     MPI_Bcast(trDOSes_.data(), epsilons_.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
// #else

    if (is_tetra_)
    {
        KPprotocol kprot(nkx_, nky_, nkz_);
        std::transform(epsilons_.begin(), epsilons_.end(), trDOSes_.begin(), [this, &velocs_squared, &kprot] (double v) {
            return tetrahedra::evaluate_dos_at_value(velocs_squared, eigenens_, kprot, v);
        });
    }
    else
    {
        std::transform(epsilons_.begin(), epsilons_.end(), trDOSes_.begin(),
            [this, &velocs_squared] (double e)
            {
                auto sm_trDOS = evaluate_smeared_trdos_at_value(e, elec_smearing_, elec_sampling_, eigenens_, velocs_squared, Te_); 
                return sm_trDOS;
            });
    }
    dump_trdos_file();
//#endif
    t.stop("\t  Transport DOS is evaluated. Results are written to VelocitiesDOS.dat");
    std::cout << "\t  Time elapsed: " << t.elapsed() << " ms" << std::endl;
}

void SpecFunc::dump_trdos_file()
{
     int rank{ 0 };
#ifdef SKIES_MPI
    rank = mpi::rank();
#endif
    std::ofstream os;
    std::string fname = "VelocitiesDOS.dat";
    if (!rank)
    {
        os.open(fname);
        std::string header;
        header += "# Generated in spectral function construnction\n";
        if (is_tetra_)
            header += "# smearing: tetrahedra\n";
        else
            header +=  "# smearing: " + std::to_string(elec_smearing_) + " ev\n";
        header += "# velocity component: " + std::to_string(0)
               +  "\n# electron energy [eV]               transport DOS [r.a.u.]\n";
        os << header;
        for (size_t i = 0; i < epsilons_.size(); ++i)
            os << std::setprecision(6) << std::setw(12) << epsilons_[i] << std::setw(34) << trDOSes_[i] << std::endl;
        os.close();
    }
}

SpecFunc::SpecFunc(const std::string& fname)
{
    std::ifstream ifs(fname);
    std::string line;

    if (fname != "LambdaTr_plus.dat" && fname != "LambdaTr_minus.dat")
        throw std::runtime_error("To continue an interrupted calculation 'LambdaTr_plus(minus).dat' file must be in the working dir");

    if (ifs.fail())
        throw std::runtime_error("The file for continuation does not exist");

    if (ifs.good()) getline(ifs, line);
    if (ifs.good()) getline(ifs, line);
    auto splitted_line = custom_split(line, ' ');

    auto xnq = splitted_line[3].data();
    nqx_ = std::stoi(custom_split(xnq, 'x')[0].data());
    nqy_ = std::stoi(custom_split(xnq, 'x')[1].data());
    nqz_ = std::stoi(custom_split(xnq, 'x')[2].data());

    auto xnk = splitted_line[7].data();
    nkx_ = std::stoi(custom_split(xnk, 'x')[0].data());
    nky_ = std::stoi(custom_split(xnk, 'x')[1].data());
    nkz_ = std::stoi(custom_split(xnk, 'x')[2].data());

    if (ifs.good()) getline(ifs, line);
    nmds_ = std::stod(custom_split(line, ' ')[3].data());

    if (ifs.good()) getline(ifs, line);

    double elec_smearing{ 1.0 };
    bzsampling::SamplingFunc elec_sampling;
    if (ifs.good()) getline(ifs, line);
    std::string word = custom_split(line, ' ')[2].data();
    if (word != "tetrahedra")
    {
        is_tetra_ = false;
        elec_sampling = switch_sampling(custom_split(line, ' ')[2].data());
        elec_smearing = std::stod(custom_split(line, ' ')[3].data());
    }
    else
    {
        is_tetra_ = true;
    }

    double phon_smearing{ 1.0 };
    bzsampling::SamplingFunc phon_sampling;
    if (ifs.good()) getline(ifs, line);
    word = custom_split(line, ' ')[2].data();
    if (word != "tetrahedra")
    {
        is_tetra_ = false;
        phon_sampling = switch_sampling(custom_split(line, ' ')[2].data());
        phon_smearing = std::stod(custom_split(line, ' ')[3].data());
    }
    else
    {
        is_tetra_ = true;
    }

    if (ifs.good()) getline(ifs, line);
    sign_    = std::stoi(custom_split(line, ' ')[1].data());

    if (ifs.good()) getline(ifs, line);
    sign_pr_ = std::stoi(custom_split(line, ' ')[1].data());

    if (ifs.good()) getline(ifs, line);
    cart_ = *custom_split(line, ' ')[2].data();

    if (ifs.good()) getline(ifs, line);
    auto epsilons_line = custom_split(line, ' ');
    if (epsilons_line.size() < 5)
        throw std::runtime_error("Electron energy list must contain at least one value.");
    for (size_t i = 4; i < epsilons_line.size(); ++i)
        epsilons_.push_back(std::stod(epsilons_line[i]));

    if (ifs.good()) getline(ifs, line);
    auto transDOSes_line = custom_split(line, ' ');
    if (transDOSes_line.size() < 5)
        throw std::runtime_error("Transport DOS list must contain at least one value.");
    for (size_t i = 4; i < transDOSes_line.size(); ++i)
        trDOSes_.push_back(std::stod(transDOSes_line[i]) * units::eV_in_Ry);

    if (ifs.good()) getline(ifs, line);

    size_t cnt{0};
    size_t iq{0}, imd{0};
    nqpt_ = nqx_ * nqy_ * nqz_;
    inner_sum_.resize(epsilons_.size(), array2D(nqpt_, array1D(nmds_, 0.0)));
    while (ifs.good())
    {
        getline(ifs, line);
        if (imd == nmds_) { imd = 0; iq++; }
        if (!line.empty()) {
            for (size_t ieps = 0, s = epsilons_.size(); ieps < s; ++ieps)
                inner_sum_[ieps][iq][imd] = std::stod(custom_split(line, ' ')[7 + ieps * 2].data());
            cnt++;
            imd++;
        }
    }
    if (cnt < nqpt_ * nmds_) is_full_ = false;
    else is_full_ = true;
    iq_cont_  = iq;
    imd_cont_ = imd;
    if (imd == nmds_) { iq_cont_++; imd_cont_ = 0; }
    for (size_t ieps = 0, s = epsilons_.size(); ieps < s; ++ieps)
        inner_sum_[ieps].resize(iq_cont_ + 1);

    ifs.close();

    std::vector<size_t> kpgrid = {nkx_, nky_, nkz_};
    std::vector<size_t> qpgrid = {nqx_, nqy_, nqz_};

    init(kpgrid, qpgrid, elec_sampling, phon_sampling, elec_smearing, phon_smearing);

    is_continue_calc_ = true;
}

SpecFunc::~SpecFunc() {}

array1D SpecFunc::calc_spec_func(double Omega)
{
    array1D a2f = 0.5 * calc_exter_sum(Omega) * ((is_tetra_) ? (1.0 / nkpt_) : (1.0 / nkpt_ / nqpt_));
    assert(a2f.size() == epsilons_.size());
    for (size_t ieps = 0, s = epsilons_.size(); ieps < s; ++ieps)
        a2f[ieps] /= ( trDOSes_[ieps] / (is_continue_calc_ ? skies::units::Ry_in_eV : 1) );
    return a2f;
}

// calculates inner sum with fixed q and \nu
double SpecFunc::calc_inner_sum(size_t iq, size_t imd, size_t ieps)
{
#if defined(SKIES_MPI)
// pure MPI version
    int rank = mpi::rank();
    int nproc = mpi::size();

    auto rcounts = std::vector<int>(nproc);
    auto  displs = std::vector<int>(nproc);

    rcounts = get_rcounts_displs(nkpt_, nproc).first;
    displs  = get_rcounts_displs(nkpt_, nproc).second;

    int count{ rcounts[rank] };
    int displ{ displs[rank]  };

    double locInnerSum{ 0.0 };
    if (is_tetra_)
        locInnerSum = calc_inner_sum_in_subarray_tetra(iq, imd, ieps, displ, displ + count, low_band_, high_band_);
    else
        locInnerSum = calc_inner_sum_in_subarray(iq, imd, ieps, displ, displ + count, low_band_, high_band_);

    // final sum count
    double sum_of_sums{0.0};
    MPI_Reduce(&locInnerSum, &sum_of_sums, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sum_of_sums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    inner_sum_[ieps][iq][imd] = sum_of_sums;
    return sum_of_sums;
#endif
// pure serial version
    double innerSum = calc_inner_sum_in_subarray(iq, imd, ieps, 0, nkpt_, low_band_, high_band_);
    inner_sum_[ieps][iq][imd] = innerSum;
    return innerSum;
}

double
SpecFunc::calc_inner_sum_in_subarray(size_t iq, size_t imd, size_t ieps,
                                     size_t start, size_t finish,
                                     size_t low_band, size_t high_band)
{
    assert(low_band <= high_band);
    double matel2{ 0.0 };
    double innerSum{ 0.0 };
    size_t ik{ start };
    for (; ik < finish; ++ik)
    {
        auto k = kpts_[ik];
        auto q = qpts_[iq];
        auto qk = k + q;
        auto eps = epsilons_[ieps]; // current electron energy level
        // first loop for quantities at initial k-point and band n
        for (size_t n = low_band; n < high_band + 1; ++n)
        {
            // k point is handled beforehand
            double delta_ekn  = elec_sampling_(eigenens_[ik][n] - eps, elec_smearing_);
            double vkn  = elvelocs_[ik][n];

            // k+q point is handled on the fly
            auto qkbands = EigenValue::interpolate_at(qk); 
            auto qkvels  = Velocities(cart_).interpolate_at(qk);

            // inner loop for quantities at final k+q-point and band m
            for (size_t m = low_band; m < high_band + 1; ++m)
            {
                double delta_eqkm = elec_sampling_(qkbands[m] - eps, elec_smearing_);
                matel2 = EPHMatrixSquared::interpolate_at(k, q, imd, n, m);
                innerSum += matel2 * (vkn - sign_*qkvels[m])*(vkn - sign_pr_*qkvels[m]) * delta_ekn * delta_eqkm;
            }
        }
    }

    return innerSum;
}

double
SpecFunc::calc_inner_sum_in_subarray_tetra(size_t iq, size_t imd, size_t ieps,
                                     size_t start, size_t finish,
                                     size_t low_band, size_t high_band)
{
    // similar to calc_inner_sum_in_subarray, but uses tetrahedron method for BZ sampling
    assert(low_band <= high_band);
    double matel2{ 0.0 };
    array3D matels;
    auto nbnd = high_band - low_band + 1;

    size_t ik{ start };
    array2D eigenens(nkpt_, array1D(nbnd, 0.0));
    array2D eigenens_qk(nkpt_, array1D(nbnd, 0.0));
    matels.resize(nbnd, array2D(nbnd, array1D(nkpt_, 0.0)));
    for (; ik < finish; ++ik)
    {
        auto k = kpts_[ik];
        auto q = qpts_[iq];
        auto qk = k + q;
        auto tmp_qk = EigenValue::interpolate_at(qk);
        auto qkvels  = Velocities(cart_).interpolate_at(qk);
        for (size_t n = low_band; n < high_band + 1; ++n)
        {
            eigenens[ik][n - low_band] = eigenens_[ik][n];
            eigenens_qk[ik][n - low_band] = tmp_qk[n];
            double vkn  = elvelocs_[ik][n];

            for (size_t m = low_band; m < high_band + 1; ++m)
            {
                double vkqm = qkvels[m];
                matel2 = EPHMatrixSquared::interpolate_at(k, q, imd, n, m);
                auto veloc_squared = (vkn - sign_*vkqm)*(vkn - sign_pr_*vkqm);
                matels[n - low_band][m - low_band][ik] = matel2 * veloc_squared;
            }
        }
    }
#ifdef SKIES_MPI
    eigenens = mpi::sum_one_from_many(eigenens);
    eigenens_qk = mpi::sum_one_from_many(eigenens_qk);
    matels = mpi::sum_one_from_many(matels);
#endif 
    auto eps = epsilons_[ieps];
    KPprotocol kprot(nkx_, nky_, nkz_);
    double innerSum = tetrahedra::evaluate_dos_at_values(matels, eigenens, eigenens_qk, kprot, start, finish, eps, eps);
    return innerSum;
}

void dump_header_lambda_file(const SpecFunc& a2f, std::ofstream& os)
{
    os << "Mode-resolved transport coupling strength" << std::endl;

    os << "num. of q-points: "   << a2f.nqx() << "x" << a2f.nqy() << "x" << a2f.nqz()
       << ", num. of k-points: " << a2f.nkx() << "x" << a2f.nky() << "x" << a2f.nkz() << std::endl; 
    os << "num. of modes: " << a2f.nmds() << std::endl;

    os << "Fermi level: " << EigenValue::eF << " eV" << std::endl;

    if (a2f.is_tetra())
    {
        os << "electron sampling: " << "tetrahedra" << std::endl;
        os << "phonon sampling: "   << "tetrahedra" << std::endl;
    }
    else
    {
        os << "electron sampling: " << a2f.get_type_of_el_smear() << " " << a2f.elec_smearing() << " eV" << std::endl;
        os << "phonon sampling: "   << a2f.get_type_of_ph_smear() << " " << a2f.phon_smearing() << " eV" << std::endl;
    }

    os << "sign: "    << a2f.sign()    << std::endl;
    os << "sign': " << a2f.sign_pr() << std::endl;

    os << "velocity component: " << a2f.cart() << std::endl;

    os << "electron energy list [eV]: ";
    for (auto&& e : a2f.epsilons()) os << e << " ";
    os << std::endl;
    os << "transport DOS list [r.a.u.]: ";
    for (auto&& e : a2f.trans_doses()) os << e * skies::units::Ry_in_eV  << " "; // go to [r.a.u.]
    os << std::endl;

    os << std::setw(5) << std::left << "iq" << std::right
       << std::setw(8)  << "qx"
       << std::setw(8)  << "qy"
       << std::setw(8)  << "qz"
       << std::setw(6)  << "nu"
       << std::setw(15) << "om_q_nu";
    for (size_t i = 0; i < a2f.epsilons().size(); ++i)
    {
        os << std::setw(15) << "lambda_q_nu"
           << std::setw(15) << "inner_sum";
    }
    os << std::endl;
}

array1D SpecFunc::calc_exter_sum(double Omega)
{
    int rank{ 0 };
#ifdef SKIES_MPI
    rank = mpi::rank();
#endif
    array1D exterSum(epsilons_.size());
    KPprotocol qprot(nqx_, nqy_, nqz_);
    if (!is_full_)
    {
        std::ofstream os;
        std::string fname = "LambdaTr";
        std::string suffix = sign() > 0 ? "_plus.dat" : "_minus.dat";
        fname += suffix; 

        if (is_continue_calc_)
        {
            if (!rank)
                os.open(fname, std::ios_base::app);
            inner_sum_.resize(epsilons_.size(), array2D(nqpt_, array1D(nmds_)));
        }
        else
        {
            if (!rank)
            {
                os.open(fname);
                dump_header_lambda_file(*this, os);
            }
            inner_sum_.resize(epsilons_.size(), array2D(nqpt_, array1D(nmds_, 0.0)));
        }

        for (size_t iq = iq_cont_; iq < nqpt_; ++iq) 
        {
            auto modes = EigenFrequency::interpolate_at(qpts_[iq]);

            size_t imd = is_continue_calc_ ? imd_cont_ : 0;
            for (; imd < nmds_; ++imd)
            {
                if (is_tetra_)
                {
                    for (size_t ieps = 0; ieps < epsilons_.size(); ++ieps)
                    {
                        // just accumulate inner matrix elements inner_sum_{q\nu} for the next use of tetrahedron method
                        calc_inner_sum(iq, imd, ieps); 
                    }
                }
                else
                {
                    double delta_omql = phon_sampling_(modes[imd] - Omega, phon_smearing_);
                    for (size_t ieps = 0; ieps < epsilons_.size(); ++ieps)
                        exterSum[ieps] += delta_omql * calc_inner_sum(iq, imd, ieps);
                }
                double eigfreq = modes[imd];
                array1D lambda(epsilons_.size());
                for (size_t ieps = 0; ieps < epsilons_.size(); ++ieps)
                    lambda[ieps] = inner_sum_[ieps][iq][imd] / nkpt_ / trDOSes_[ieps] / eigfreq;
                if (!rank)
                {
                    os << std::setw(5) << std::left << iq + 1 << std::right
                    << std::setw(8) << std::setprecision(3) << qpts_[iq][0]
                    << std::setw(8) << std::setprecision(3) << qpts_[iq][1]
                    << std::setw(8) << std::setprecision(3) << qpts_[iq][2]
                    << std::setw(6) << imd + 1
                    << std::setprecision(3) << std::setw(15) << eigfreq;
                    for (size_t ieps = 0; ieps < epsilons_.size(); ++ieps) {
                    os << std::setprecision(5) << std::setw(15) << lambda[ieps]
                        << std::setprecision(5) << std::setw(15) << inner_sum_[ieps][iq][imd];
                    }
                    os << std::endl;
                }
            } // modes
            is_continue_calc_ = false;
        } // qpts
        if (is_tetra_)
        {
            // evaluate for this Omega. At this step matrix elements inner_sum_{q\nu} have been accumulated
            for (size_t ieps = 0; ieps < epsilons_.size(); ++ieps)
                exterSum[ieps] = tetrahedra::evaluate_dos_at_value(inner_sum_[ieps], eigenfreqs_, qprot, Omega);
        }
        if (!rank)
            os.close();
        is_full_ = true;
    }
    else // !is_full_
    {
        if (is_tetra_)
        {
            // here inner matrix elements inner_sum_{q\nu} have been accumulated
            for (size_t ieps = 0; ieps < epsilons_.size(); ++ieps)
                exterSum[ieps] = tetrahedra::evaluate_dos_at_value(inner_sum_[ieps], eigenfreqs_, qprot, Omega);
            return exterSum;
        }
        for (size_t iq = 0; iq < nqpt_; ++iq)
        {
            for (size_t imd = 0; imd < nmds_; ++imd)
            {
                double delta_omql = phon_sampling_(eigenfreqs_[iq][imd] - Omega, phon_smearing_);
                for (size_t ieps = 0; ieps < epsilons_.size(); ++ieps)
                    exterSum[ieps] += delta_omql * inner_sum_[ieps][iq][imd];
            }
        }
    }
    return exterSum;
}

///////////////////////////////////////////////////////////////////////////////
////      calc_spec_func is the main driver function to calculate          ////
////                 in a range of omegas and epsilons                     ////
///////////////////////////////////////////////////////////////////////////////

void calc_spec_func(SpecFunc& a2f, const array1D& omegas, const std::string& fname)
{
    int rank{ 0 };
#ifdef SKIES_MPI
    rank = mpi::rank();
#endif
    std::ofstream os;
    if (!rank)
    {
        os.open(fname);
        auto array1D_to_string = [] (const array1D& arr) {
            std::string out;
            for (auto&& e : arr) {
                out += std::to_string(e);
                out += " ";
            }
            return out;
        };
        std::string header = 
                  "# elec_smearing: " + ((a2f.is_tetra()) ? "tetrahedra\n" : (std::to_string(a2f.elec_smearing()) + "eV\n"))
                + "# phon_smearing: " + ((a2f.is_tetra()) ? "tetrahedra\n" : (std::to_string(a2f.phon_smearing()) + "eV\n"))
                + "# sign: " + std::to_string(a2f.sign())
                + "\n# velocity component: " + std::to_string(0)
                + "\n# electron energy list [eV]: " + array1D_to_string(a2f.epsilons())
                + "\n# transport DOS for energy list [r.a.u.]: "
                + array1D_to_string(a2f.trans_doses() * ((a2f.is_continue_calc()) ? 1 : skies::units::Ry_in_eV))
                + "\n\n# frequency [eV]               spectral function\n";
        os << header;
    }

    // make calculation for specific Omega
    for (auto&& Om : omegas)
    {
        auto vals = a2f.calc_spec_func(Om);
        if (!rank)
        {
            os << std::left << std::setw(20) << Om << std::right;
            for (size_t ieps = 0; ieps < a2f.epsilons().size(); ++ieps)
                os << std::setw(12) << std::setprecision(6) << vals[ieps];
            os << std::endl;
        }
    }

    if (!rank)
        os.close();
}

} // spectral
} // skies
