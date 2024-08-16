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

#include <skies/utils/tbb_wrapper.h>

using namespace skies::arrays;
using namespace skies::bzsampling;
using namespace skies::quantities;

#ifdef SKIES_MPI
using namespace skies::mpi;
#endif

namespace skies { namespace spectral {

void SpecFunc::init()
{
    launch::Timer t;
    t.start("========= Started transport spectral function initialization...");
    if (is_tetra_)
    {
        tetrahedra::TetraHandler::set_kprot(kprot_);
        tetrahedra::TetraHandler::set_phon_tag(true);
        tetrahedra::TetraHandler::set_qprot(qprot_);
        tetrahedra::DoubleTetraHandler::set_kprot(kprot_);
    }

    // eigenens_ at (n,k)-grid are filled only once here
    eigenens_.resize(kprot_.nkpt(), array1D(nbnd_));
    std::transform(PAR kprot_.grid().begin(), kprot_.grid().end(), eigenens_.begin(), [] (auto&& k) {
        return EigenValueDrawable().interpolate_at(k);
    });

    // eigefreqs_ at (\nu,q)-grid are filled only once here
    eigenfreqs_.resize(qprot_.nkpt(), array1D(nmds_));
    std::transform(PAR qprot_.grid().begin(), qprot_.grid().end(), eigenfreqs_.begin(), [] (auto&& q) {
        return EigenFrequencyDrawable().interpolate_at(q);
    });

    // elvelocs at (n,k) are filled only once here
    fsh_alpha_.resize(kprot_.nkpt(), array1D(nbnd_, 0.0));
    prepare_squared_velocs(alpha_, elvelocs_alpha_, elvelocs_alpha_sq_);
    elvelocs_beta_ = elvelocs_alpha_;
    if (alpha_ != beta_)
    {
        fsh_beta_.resize(kprot_.nkpt(), array1D(nbnd_, 0.0));
        prepare_squared_velocs(beta_, elvelocs_beta_, elvelocs_beta_sq_);
    }

    // fsh_alpha_ and fsh_beta_ are initialized only once here
    array2D weights(EigenValueDrawable::nbands, array1D(kprot_.nkpt(), 1));
    tetrahedra::TetraHandler th_dos(std::move(weights), transpose(eigenens_));
    th_dos_ = std::move(th_dos);
    if (is_tetra_)
    {
        tetrahedra::TetraHandler th_trdos_alpha(transpose(elvelocs_alpha_sq_), transpose(eigenens_));
        th_trdos_alpha_ = std::move(th_trdos_alpha);
        if (epsilons_.size() > 1)
            prepare_fsh(th_dos_, th_trdos_alpha_, elvelocs_alpha_, fsh_alpha_);
        if (alpha_ != beta_)
        {
            tetrahedra::TetraHandler th_trdos_beta(transpose(elvelocs_beta_sq_), transpose(eigenens_));
            th_trdos_beta_ = std::move(th_trdos_beta);
            if (epsilons_.size() > 1)
                prepare_fsh(th_dos_, th_trdos_beta_, elvelocs_beta_, fsh_beta_);
        }
        else
        {
            fsh_beta_ = fsh_alpha_; // by default when alpha_ == beta_
        }
    }
    else // is_tetra
    {
        if (epsilons_.size() > 1)
            prepare_fsh(elvelocs_alpha_, elvelocs_alpha_sq_, fsh_alpha_);
        if (alpha_ != beta_)
        {
            if (epsilons_.size() > 1)
                prepare_fsh(elvelocs_beta_, elvelocs_beta_sq_, fsh_beta_);
        }
        else
        {
            fsh_beta_ = fsh_alpha_; // by default when alpha_ == beta_
        }
    }

    // dump fsh into file to continue calculations in future
    dump_array(std::string{"FSH_"} + std::string{alpha_} + std::string{".dat"}, fsh_alpha_);
    dump_array(std::string{"FSH_"} + std::string{ beta_} + std::string{".dat"},  fsh_beta_);

    // calculate DOS for range of given epsilons
    DOSes_.resize(epsilons_.size(), 0.0);
    if (is_tetra_)
    {
        std::transform(epsilons_.begin(), epsilons_.end(), DOSes_.begin(), [this] (auto&& eps) {
            return th_dos_.evaluate_dos_at_value(eps); // in [1/eV/spin/cell]
        });
    }
    else
    {
        std::transform(epsilons_.begin(), epsilons_.end(), DOSes_.begin(), [this] (auto&& eps) {
            return evaluate_dos_at_value<EigenValue>(eps, elec_smearing_, elec_sampling_, eigenens_); // in [1/eV/spin/cell]
        });
    }

    // calculate trDOS for range of given epsilons
    // explicitly evaluate transport DOS for each of x, y and z components and take average of the three
    trDOSes_.resize(epsilons_.size(), 0.0);
    array2D elvelocs_x, elvelocs_y, elvelocs_z;
    array2D elvelocs_x_sq, elvelocs_y_sq, elvelocs_z_sq;
    prepare_squared_velocs('x', elvelocs_x, elvelocs_x_sq);
    prepare_squared_velocs('y', elvelocs_y, elvelocs_y_sq);
    prepare_squared_velocs('z', elvelocs_z, elvelocs_z_sq);
    array1D trDOSes_x(epsilons_.size()), trDOSes_y(epsilons_.size()), trDOSes_z(epsilons_.size());
    if (is_tetra_)
    {
        tetrahedra::TetraHandler th_x(transpose(elvelocs_x_sq), transpose(eigenens_));
        std::transform(epsilons_.begin(), epsilons_.end(), trDOSes_x.begin(), [th = std::move(th_x)] (auto&& eps) {
            return th.evaluate_dos_at_value(eps);
        });
        tetrahedra::TetraHandler th_y(transpose(elvelocs_y_sq), transpose(eigenens_));
        std::transform(epsilons_.begin(), epsilons_.end(), trDOSes_y.begin(), [th = std::move(th_y)] (auto&& eps) {
            return th.evaluate_dos_at_value(eps);
        });
        tetrahedra::TetraHandler th_z(transpose(elvelocs_z_sq), transpose(eigenens_));
        std::transform(epsilons_.begin(), epsilons_.end(), trDOSes_z.begin(), [th = std::move(th_z)] (auto&& eps) {
            return th.evaluate_dos_at_value(eps);
        });
    }
    else // is_tetra
    {
        std::transform(epsilons_.begin(), epsilons_.end(), trDOSes_x.begin(), [&] (auto&& eps) {
            return evaluate_trdos_at_value(eps, elec_smearing_, elec_sampling_, eigenens_, elvelocs_x_sq);
        });
        std::transform(epsilons_.begin(), epsilons_.end(), trDOSes_y.begin(), [&] (auto&& eps) {
            return evaluate_trdos_at_value(eps, elec_smearing_, elec_sampling_, eigenens_, elvelocs_y_sq);
        });
        std::transform(epsilons_.begin(), epsilons_.end(), trDOSes_z.begin(), [&] (auto&& eps) {
            return evaluate_trdos_at_value(eps, elec_smearing_, elec_sampling_, eigenens_, elvelocs_z_sq);
        });
    }
    trDOSes_ = (trDOSes_x + trDOSes_y + trDOSes_z) * (1.0 / 3.0);

    t.stop("\t  Transport spectral function is initialized");
    int rank{ 0 };
#ifdef SKIES_MPI
    rank = mpi::rank();
#endif
    if (!rank)
        std::cout << "\t  Time elapsed: " << t.elapsed() << " ms" << std::endl;
}

void SpecFunc::prepare_squared_velocs(char cart, array2D& elvelocs, array2D& elvelocs_sq)
{
    elvelocs.resize(kprot_.nkpt(), array1D(nbnd_, 0.0));
    elvelocs_sq.resize(kprot_.nkpt(), array1D(nbnd_, 0.0));
    std::transform(PAR kprot_.grid().begin(), kprot_.grid().end(), elvelocs.begin(),
                [&] (auto&& k) { return Velocities(cart).interpolate_at(k); });
    std::transform(PAR elvelocs.begin(), elvelocs.end(), elvelocs_sq.begin(),
            [] (auto&& v) {
                auto squared_v = array1D(EigenValue::nbands, 0.0);
                std::transform(v.begin(), v.end(), squared_v.begin(), [] (auto&& x) { return x * x; });
                return squared_v;
    });
}

void SpecFunc::prepare_fsh(const array2D& elvelocs, const array2D& elvelocs_sq, array2D& fsh)
{
    std::transform(PAR eigenens_.begin(), eigenens_.end(), elvelocs.begin(), fsh.begin(),
        [&] (auto&& ee, auto&& vv) {
            auto fsh_line = array1D(EigenValue::nbands, 0.0);
            std::transform(ee.begin(), ee.end(), vv.begin(), fsh_line.begin(),
                [&] (auto&& e, auto&& v) {
                    auto DOS_nk = evaluate_dos_at_value<EigenValue>(e, elec_smearing_, elec_sampling_, eigenens_); // in [1/eV]
                    auto trDOS_nk = evaluate_trdos_at_value(e, elec_smearing_, elec_sampling_, eigenens_, elvelocs_sq); // in [Ry^2 * bohr^2 / eV]
                    return v / std::sqrt(trDOS_nk / DOS_nk); // unitless, v in [Ry * bohr]
            });
            return fsh_line;
    });
}

void SpecFunc::prepare_fsh(const tetrahedra::TetraHandler& th_dos,
                           const tetrahedra::TetraHandler& th_trdos,
                           const array2D& elvelocs,
                           array2D& fsh)
{
    std::transform(PAR eigenens_.begin(), eigenens_.end(), elvelocs.cbegin(), fsh.begin(),
        [&] (auto&& ee, auto&& vv) {
            auto fsh_line = array1D(EigenValue::nbands, 0.0);
            std::transform(ee.begin(), ee.end(), vv.cbegin(), fsh_line.begin(),
                [&] (auto&& e, auto&& v) {
                    auto DOS_nk = th_dos.evaluate_dos_at_value(e); // in [1/eV]
                    auto trDOS_nk = th_trdos.evaluate_dos_at_value(e); // in [Ry^2 * bohr^2 / eV]
                    return v / std::sqrt(trDOS_nk / DOS_nk); // unitless, v in [Ry * bohr]
            });
            return fsh_line;
    });
}

SpecFunc::SpecFunc(const std::vector<size_t>& kpgrid,
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
                   bool is_tetra
)
    : kprot_(KPprotocol{kpgrid[0], kpgrid[1], kpgrid[2]})
    , qprot_(KPprotocol{qpgrid[0], qpgrid[1], qpgrid[2]})
    , sign_(sign)
    , sign_pr_(sign_pr)
    , nbnd_(EigenValue::nbands)
    , nmds_(EigenFrequency::nmodes)
    , high_band_(nbnd_ - 1)
    , Te_(Te)
    , alpha_(alpha)
    , beta_(beta)
    , is_tetra_(is_tetra)
{
    epsilons_.push_back(0);
    elec_sampling_ = elec_sampling;
    phon_sampling_ = phon_sampling;
    elec_smearing_ = elec_smearing;
    phon_smearing_ = phon_smearing;
    init();
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
               char alpha,
               char beta,
               bool is_tetra
)
    : kprot_(KPprotocol{kpgrid[0], kpgrid[1], kpgrid[2]})
    , qprot_(KPprotocol{qpgrid[0], qpgrid[1], qpgrid[2]})
    , sign_(sign)
    , sign_pr_(sign_pr)
    , nbnd_(EigenValue::nbands)
    , nmds_(EigenFrequency::nmodes)
    , high_band_(nbnd_ - 1)
    , epsilons_(epsilons)
    , Te_(Te)
    , alpha_(alpha)
    , beta_(beta)
    , is_tetra_(is_tetra)
{
    elec_sampling_ = elec_sampling;
    phon_sampling_ = phon_sampling;
    elec_smearing_ = elec_smearing;
    phon_smearing_ = phon_smearing;
    init();
}


/*void SpecFunc::dump_trdos_file()
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
}*/

SpecFunc::SpecFunc(const std::string& fname)
{
    std::ifstream ifs(fname);
    std::string line;

    if (ifs.fail())
        throw std::runtime_error("The file for continuation does not exist");

    if (ifs.good()) std::getline(ifs, line);
    if (ifs.good()) std::getline(ifs, line);
    auto splitted_line = custom_split(line, ' ');

    auto xnq = splitted_line[4].data();
    size_t nqx = std::stoi(custom_split(xnq, 'x')[0].data());
    size_t nqy = std::stoi(custom_split(xnq, 'x')[1].data());
    size_t nqz = std::stoi(custom_split(xnq, 'x')[2].data());

    auto xnk = splitted_line[8].data();
    size_t nkx = std::stoi(custom_split(xnk, 'x')[0].data());
    size_t nky = std::stoi(custom_split(xnk, 'x')[1].data());
    size_t nkz = std::stoi(custom_split(xnk, 'x')[2].data());

    if (ifs.good()) std::getline(ifs, line);
    size_t nmds = std::stoi(custom_split(line, ' ')[4].data());

    if (ifs.good()) std::getline(ifs, line);
    size_t low_band = std::stoi(custom_split(line, ' ')[2].data());

    if (ifs.good()) std::getline(ifs, line);
    size_t high_band = std::stoi(custom_split(line, ' ')[2].data());

    if (ifs.good()) std::getline(ifs, line);
    double Te = std::stod(custom_split(line, ' ')[2].data());

    if (ifs.good()) std::getline(ifs, line);
    EigenValue::eF = std::stod(custom_split(line, ' ')[3].data());

    double elec_smearing{ 1.0 };
    bzsampling::SamplingFunc elec_sampling = [] (double, double) { return 0.0; };
    if (ifs.good()) std::getline(ifs, line);
    std::string word = custom_split(line, ' ')[3].data();
    bool is_tetra = false;
    if (word != "tetrahedra")
    {
        is_tetra = false;
        elec_sampling = switch_sampling(custom_split(line, ' ')[3].data());
        elec_smearing = std::stod(custom_split(line, ' ')[4].data());
    }
    else
    {
        is_tetra = true;
    }

    double phon_smearing{ 1.0 };
    bzsampling::SamplingFunc phon_sampling = [] (double, double) { return 0.0; };
    if (ifs.good()) std::getline(ifs, line);
    std::string word2 = custom_split(line, ' ')[3].data();
    if ((word2 == "tetrahedra" && word != "tetrahedra") ||
        (word2 != "tetrahedra" && word == "tetrahedra"))
        throw std::runtime_error("Error in continuation file: both samplings must be tetrahedra");
    if (word != "tetrahedra")
    {
        phon_sampling = switch_sampling(custom_split(line, ' ')[3].data());
        phon_smearing = std::stod(custom_split(line, ' ')[4].data());
    }

    if (ifs.good()) std::getline(ifs, line);
    int sign    = std::stoi(custom_split(line, ' ')[2].data());

    if (ifs.good()) std::getline(ifs, line);
    int sign_pr = std::stoi(custom_split(line, ' ')[2].data());

    if (ifs.good()) std::getline(ifs, line);
    char alpha = *custom_split(line, ' ')[4].data();

    if (ifs.good()) std::getline(ifs, line);
    char beta = *custom_split(line, ' ')[4].data();

    // initialize class fields
    kprot_ = KPprotocol{nkx, nky, nkz};
    qprot_ = KPprotocol{nqx, nqy, nqz};
    sign_ = sign;
    sign_pr_ = sign_pr;
    nbnd_ = EigenValue::nbands;
    nmds_ = nmds;
    low_band_ = low_band;
    high_band_ = high_band;

    if (ifs.good()) std::getline(ifs, line);
    auto epsilons_line = custom_split(line, ' ');
    if (epsilons_line.size() < 6)
        throw std::runtime_error("Electron energy list must contain at least one value.");
    for (size_t i = 5; i < epsilons_line.size(); ++i)
        epsilons_.push_back(std::stod(epsilons_line[i]));

    if (ifs.good()) std::getline(ifs, line);
    auto DOSes_line = custom_split(line, ' ');
    if (DOSes_line.size() < 6)
        throw std::runtime_error("DOS list must contain at least one value.");
    for (size_t i = 5; i < DOSes_line.size(); ++i)
        DOSes_.push_back(std::stod(DOSes_line[i])); // in file DOSes are given in [1/eV/spin/cell]

    if (ifs.good()) std::getline(ifs, line);
    auto trDOSes_line = custom_split(line, ' ');
    if (trDOSes_line.size() < 6)
        throw std::runtime_error("trDOS list must contain at least one value.");
    for (size_t i = 5; i < trDOSes_line.size(); ++i)
        trDOSes_.push_back(std::stod(trDOSes_line[i])); // in file trDOSes are given in [Ry^2 * bohr^2 / eV]

    if (ifs.good()) std::getline(ifs, line);

    size_t cnt{ 0 };
    size_t iq{ 0 }, imd{ 0 };
    inner_sum_.resize(epsilons_.size(), array2D(qprot_.nkpt(), array1D(nmds_, 0.0)));
    while (ifs.good())
    {
        std::getline(ifs, line);
        if (imd == nmds_)
        {
            imd = 0;
            iq++;
        }
        if (!line.empty())
        {
            for (size_t ieps = 0, s = epsilons_.size(); ieps < s; ++ieps)
                inner_sum_[ieps][iq][imd] = std::stod(custom_split(line, ' ')[7 + ieps * 2].data());
            cnt++;
            imd++;
        }
    }
    if (cnt < qprot_.nkpt() * nmds_)
        is_full_ = false;
    else
        is_full_ = true;
    iq_cont_  = iq;
    imd_cont_ = imd;
    if (imd == nmds_)
    {
        iq_cont_++;
        imd_cont_ = 0;
    }

    ifs.close();

    Te_ = Te;
    alpha_ = alpha;
    beta_  = beta;
    is_tetra_ = is_tetra;

    // load fsh if needed
    if (epsilons_.size() > 1)
    {
        load_array(std::string{"FSH_"} + std::string{alpha_} + std::string{".dat"}, fsh_alpha_);
        load_array(std::string{"FSH_"} + std::string{ beta_} + std::string{".dat"},  fsh_beta_);
    }
    elec_sampling_ = elec_sampling;
    phon_sampling_ = phon_sampling;
    elec_smearing_ = elec_smearing;
    phon_smearing_ = phon_smearing;
    init();

    is_continue_calc_ = true;
}

SpecFunc::~SpecFunc() {}

array1D SpecFunc::calc_spec_func(double Omega)
{ 
    array1D a2f = 0.5 * calc_exter_sum(Omega) * ((is_tetra_)
                                              ? (1.0 / kprot_.nkpt())
                                              : (1.0 / kprot_.nkpt() / qprot_.nkpt()));
    assert(a2f.size() == epsilons_.size());
    for (size_t ieps = 0, s = epsilons_.size(); ieps < s; ++ieps)
        a2f[ieps] /= (epsilons_.size() > 1)
                   ? DOSes_[ieps] // general: in [1/eV/spin/cell]
                   : trDOSes_[ieps]; // low-T formulas: in [Ry^2 * bohr^2 / eV]
    return a2f; // unitless
    // note about units:
    // 1) general: [1/eV] / [1/eV] -> unitless
    // 2) low-T: [eV^2] [bohr*Ry]^2 [1/eV] [1/eV] [1/eV] / [Ry^2 * bohr^2 / eV] -> unitless
}

// calculates inner sum with fixed q and \nu
double SpecFunc::calc_inner_sum(size_t iq, size_t imd, size_t ieps)
{
    double inner_sum{ 0.0 }; // dos, not divided by nkpt
    if (is_tetra_)
        inner_sum = calc_inner_sum_in_subarray_tetra(iq, imd, ieps, low_band_, high_band_);
    else
        inner_sum = calc_inner_sum_in_subarray(iq, imd, ieps, low_band_, high_band_);
    inner_sum_[ieps][iq][imd] = inner_sum;
    return inner_sum;
}


double
SpecFunc::calc_inner_sum_in_subarray(size_t iq, size_t imd, size_t ieps,
                                     size_t low_band, size_t high_band)
{
    assert(low_band <= high_band);
    auto nbnd = high_band - low_band + 1;

    const array2D& kpts = kprot_.grid(); 
    const array2D& qpts = qprot_.grid(); 

    return std::transform_reduce(PAR kprot_.range().begin(), kprot_.range().end(), 0.0, std::plus<double>(),
    [&] (auto&& ik)
    {
        double inner_sum{ 0.0 };
        auto k = kpts[ik];
        auto q = qpts[iq];
        auto qk = k + q;
        auto eps = epsilons_[ieps]; // current electron energy level

        // k+q point is handled on the fly
        auto tmp_qk = EigenValue::interpolate_at(qk);
        auto elvelocs_qk_alpha = Velocities(alpha_).interpolate_at(qk);
        auto fsh_qk_alpha = array1D(EigenValue::nbands, 0.0);
        if (epsilons_.size() > 1)
            std::transform(tmp_qk.begin(), tmp_qk.end(), elvelocs_qk_alpha.begin(), fsh_qk_alpha.begin(),
                [&] (auto&& e, auto&& v) {
                    auto DOS_mqk = evaluate_dos_at_value<EigenValue>(e, elec_smearing_, elec_sampling_, eigenens_);
                    auto trDOS_mqk = evaluate_trdos_at_value(e, elec_smearing_, elec_sampling_, eigenens_, elvelocs_alpha_sq_);
                    return v / std::sqrt(trDOS_mqk / DOS_mqk);
            });

        auto elvelocs_qk_beta = array1D(EigenValue::nbands, 0.0);
        auto fsh_qk_beta = array1D(EigenValue::nbands, 0.0);
        if (alpha_ != beta_)
        {
            elvelocs_qk_beta = Velocities(beta_).interpolate_at(qk);
            if (epsilons_.size() > 1)
                std::transform(tmp_qk.begin(), tmp_qk.end(), elvelocs_qk_beta.begin(), fsh_qk_beta.begin(),
                    [&] (auto&& e, auto&& v) {
                        auto DOS_mqk = evaluate_dos_at_value<EigenValue>(e, elec_smearing_, elec_sampling_, eigenens_);
                        auto trDOS_mqk = evaluate_trdos_at_value(e, elec_smearing_, elec_sampling_, eigenens_, elvelocs_beta_sq_);
                        return v / std::sqrt(trDOS_mqk / DOS_mqk);
                });
        }
        else
        {
            elvelocs_qk_beta = elvelocs_qk_alpha;
            fsh_qk_beta = fsh_qk_alpha;
        }

        // first loop for quantities at initial k-point and band n
        for (size_t n = 0; n < nbnd; ++n)
        {
            double delta_ekn = elec_sampling_(eigenens_[ik][n + low_band_] - eps, elec_smearing_);

            if (epsilons_.size() > 1) // general formulas
            {
                double fsh_kn_alpha = fsh_alpha_[ik][n + low_band_];
                double fsh_kn_beta  = fsh_beta_[ik][n + low_band_];
                for (size_t m = 0; m < nbnd; ++m)
                {
                    double delta_eqkm = elec_sampling_(tmp_qk[m + low_band_] - eps, elec_smearing_);
                    double fsh_qkm_alpha = fsh_qk_alpha[m + low_band_];
                    double fsh_qkm_beta  = fsh_qk_beta[m + low_band_];
                    double matel2 = EPHMatrixSquared::interpolate_at(k, q, imd, n + low_band_, m + low_band_);
                    double fsh_factor = (fsh_kn_alpha - sign_*fsh_qkm_alpha)*(fsh_kn_beta - sign_pr_*fsh_qkm_beta);
                    double delta_factor = delta_ekn * delta_eqkm;
                    inner_sum += matel2 * fsh_factor * delta_factor;
                }
            }
            else // low-T formulas
            {
                double vkn_alpha = elvelocs_alpha_[ik][n + low_band_];
                double vkn_beta  = elvelocs_beta_[ik][n + low_band_];
                for (size_t m = 0; m < nbnd; ++m)
                {
                    double delta_eqkm = elec_sampling_(tmp_qk[m + low_band_] - eps, elec_smearing_);
                    double vkqm_alpha = elvelocs_qk_alpha[m + low_band_];
                    double vkqm_beta  = elvelocs_qk_beta[m + low_band_];
                    double matel2 = EPHMatrixSquared::interpolate_at(k, q, imd, n + low_band_, m + low_band_);
                    double veloc_factor = (vkn_alpha - sign_*vkqm_alpha)*(vkn_beta - sign_pr_*vkqm_beta);
                    double delta_factor = delta_ekn * delta_eqkm;
                    inner_sum += matel2 * veloc_factor * delta_factor;
                }
            }
        }
        return inner_sum;
    });
}

double
SpecFunc::calc_inner_sum_in_subarray_tetra(size_t iq, size_t imd, size_t ieps,
                                     size_t low_band, size_t high_band)
{
    // similar to calc_inner_sum_in_subarray, but uses tetrahedron method for BZ sampling
    assert(low_band <= high_band);
    auto nbnd = high_band - low_band + 1;

    // size_t ik{ 0 };
    // size_t finish{ kprot_.nkpt() };
    array2D eigenens(kprot_.nkpt(), array1D(nbnd, 0.0));
    array2D eigenens_qk(kprot_.nkpt(), array1D(nbnd, 0.0));
    array3D matels(nbnd, array2D(nbnd, array1D(kprot_.nkpt(), 0.0)));
    const array2D& kpts = kprot_.grid();
    const array2D& qpts = qprot_.grid();

    // for (; ik < finish; ++ik)
    // {
    std::for_each(PAR kprot_.range().begin(), kprot_.range().end(), [&] (auto&& ik)
    {
        auto k = kpts[ik];
        auto q = qpts[iq];
        auto qk = k + q;
        auto tmp_qk = EigenValue::interpolate_at(qk);

        auto elvelocs_qk_alpha = Velocities(alpha_).interpolate_at(qk);
        auto fsh_qk_alpha = array1D(EigenValue::nbands, 0.0);
        if (epsilons_.size() > 1)
            std::transform(tmp_qk.begin(), tmp_qk.end(), elvelocs_qk_alpha.begin(), fsh_qk_alpha.begin(),
                [&] (auto&& e, auto&& v) {
                    auto DOS_mqk = th_dos_.evaluate_dos_at_value(e);
                    auto trDOS_mqk = th_trdos_alpha_.evaluate_dos_at_value(e);
                    return v / std::sqrt(trDOS_mqk / DOS_mqk);
            });

        auto elvelocs_qk_beta = array1D(EigenValue::nbands, 0.0);
        auto fsh_qk_beta = array1D(EigenValue::nbands, 0.0);
        if (alpha_ != beta_)
        {
            elvelocs_qk_beta = Velocities(beta_).interpolate_at(qk);
            if (epsilons_.size() > 1)
                std::transform(tmp_qk.begin(), tmp_qk.end(), elvelocs_qk_beta.begin(), fsh_qk_beta.begin(),
                    [&] (auto&& e, auto&& v) {
                        auto DOS_mqk = th_dos_.evaluate_dos_at_value(e);
                        auto trDOS_mqk = th_trdos_alpha_.evaluate_dos_at_value(e);
                        return v / std::sqrt(trDOS_mqk / DOS_mqk);
                });
        }
        else
        {
            elvelocs_qk_beta = elvelocs_qk_alpha;
            fsh_qk_beta = fsh_qk_alpha;
        }

        for (size_t n = 0; n < nbnd; ++n)
        {
            eigenens[ik][n] = eigenens_[ik][n + low_band_];
            eigenens_qk[ik][n] = tmp_qk[n + low_band_];

            if (epsilons_.size() > 1) // general formulas
            {
                double fsh_kn_alpha = fsh_alpha_[ik][n + low_band_];
                double fsh_kn_beta  = fsh_beta_[ik][n + low_band_];
                for (size_t m = 0; m < nbnd; ++m)
                {
                    double fsh_qkm_alpha = fsh_qk_alpha[m + low_band_];
                    double fsh_qkm_beta  = fsh_qk_beta[m + low_band_];
                    double matel2 = EPHMatrixSquared::interpolate_at(k, q, imd, n + low_band_, m + low_band_);
                    double fsh_factor = (fsh_kn_alpha - sign_*fsh_qkm_alpha)*(fsh_kn_beta - sign_pr_*fsh_qkm_beta);
                    matels[n][m][ik] = matel2 * fsh_factor;
                }
            }
            else // low-T formulas
            {
                double vkn_alpha = elvelocs_alpha_[ik][n + low_band_];
                double vkn_beta  = elvelocs_beta_[ik][n + low_band_];
                for (size_t m = 0; m < nbnd; ++m)
                {
                    double vkqm_alpha = elvelocs_qk_alpha[m + low_band_];
                    double vkqm_beta  = elvelocs_qk_beta[m + low_band_];
                    double matel2 = EPHMatrixSquared::interpolate_at(k, q, imd, n + low_band_, m + low_band_);
                    double veloc_factor = (vkn_alpha - sign_*vkqm_alpha)*(vkn_beta - sign_pr_*vkqm_beta);
                    matels[n][m][ik] = matel2 * veloc_factor;
                }
            }
        }
    });
    // } // kpts
    auto eps = epsilons_[ieps];
    tetrahedra::DoubleTetraHandler dth(std::move(matels), std::move(eigenens), std::move(eigenens_qk));
    return dth.evaluate_dos_at_values(eps, eps);
}

void dump_header_lambda_file(const SpecFunc& a2f, std::ofstream& os)
{
    os << "# Mode-resolved transport coupling strength" << std::endl;

    os << "# num. of q-points: "   << a2f.nqx() << "x" << a2f.nqy() << "x" << a2f.nqz()
       << ", num. of k-points: " << a2f.nkx() << "x" << a2f.nky() << "x" << a2f.nkz() << std::endl; 
    os << "# num. of modes: " << a2f.nmds() << std::endl;

    os << "# low_band: " << a2f.get_low_band() << std::endl;
    os << "# high_band: " << a2f.get_high_band() << std::endl;
    os << "# elec_temp: " << a2f.get_elec_temp() << std::endl;

    os << "# Fermi level: " << EigenValue::eF << " eV" << std::endl;

    if (a2f.is_tetra())
    {
        os << "# electron sampling: " << "tetrahedra" << std::endl;
        os << "# phonon sampling: "   << "tetrahedra" << std::endl;
    }
    else
    {
        os << "# electron sampling: " << a2f.get_type_of_el_smear() << " " << a2f.elec_smearing() << " eV" << std::endl;
        os << "# phonon sampling: "   << a2f.get_type_of_ph_smear() << " " << a2f.phon_smearing() << " eV" << std::endl;
    }

    os << "# sign: "    << a2f.sign()    << std::endl;
    os << "# sign': " << a2f.sign_pr() << std::endl;

    os << "# velocity component alpha: " << a2f.alpha() << std::endl;
    os << "# velocity component  beta: " << a2f.beta()  << std::endl;

    os << "# electron energy list [eV]: ";
    for (auto&& e : a2f.epsilons()) os << e << " ";
    os << std::endl;
    os << "# electronic DOS list [1/eV/spin/cell]: ";
    for (auto&& e : a2f.doses()) os << e << " ";
    os << std::endl;
    os << "# transport DOS list [Ry^2*bohr^2/eV]: ";
    for (auto&& e : a2f.trans_doses()) os << e << " ";
    os << std::endl;

    os << std::setw(5) << std::left << "# iq" << std::right
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
    if (!is_full_)
    {
        std::ofstream os;
        std::string filename = "LambdaTr_";
        if (sign_ > 0)   filename += std::string{'p'};
        else		     filename += std::string{'m'};
        if (sign_pr_ > 0) filename += std::string{'p'};
        else		     filename += std::string{'m'};
        filename += std::string{'_'};
        filename += std::string{alpha_};
        filename += std::string{beta_};
        filename += std::string{".dat"};

        const array2D& qpts = qprot_.grid();

        if (is_continue_calc_)
        {
            if (!rank)
                os.open(filename, std::ios_base::app);
            inner_sum_.resize(epsilons_.size(), array2D(qprot_.nkpt(), array1D(nmds_)));
        }
        else
        {
            if (!rank)
            {
                os.open(filename);
                dump_header_lambda_file(*this, os);
            }
            inner_sum_.resize(epsilons_.size(), array2D(qprot_.nkpt(), array1D(nmds_, 0.0)));
        }

        for (size_t iq = iq_cont_; iq < qprot_.nkpt(); ++iq) 
        {
            auto modes = EigenFrequency::interpolate_at(qpts[iq]); // only main thread here

            size_t imd = (is_continue_calc_) ? imd_cont_ : 0;
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
                    lambda[ieps] = inner_sum_[ieps][iq][imd] / kprot_.nkpt() / trDOSes_[ieps] / eigfreq;
                if (!rank)
                {
                    os << std::setw(5) << std::left << iq + 1 << std::right
                    << std::setw(8) << std::setprecision(3) << qpts[iq][0]
                    << std::setw(8) << std::setprecision(3) << qpts[iq][1]
                    << std::setw(8) << std::setprecision(3) << qpts[iq][2]
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
                exterSum[ieps] = tetrahedra::evaluate_dos(transpose(inner_sum_[ieps]), eigenfreqs_, Omega, true);
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
                exterSum[ieps] = tetrahedra::evaluate_dos(transpose(inner_sum_[ieps]),  eigenfreqs_, Omega, true);
            return exterSum; // already divided by nqpt!
        }
        else
        for (size_t iq = 0; iq < qprot_.nkpt(); ++iq)
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
    // note about final units:
    // 1) general: |g^2| (fx - fx)*(fy - fy) * d(eps)*d(eps)*d(Omega) -> [eV^2] [1] [1] [1/eV] [1/eV] [1/eV] = [1/eV]
    // 2) low-T: |g^2| (vx - vx)*(vy - vy) * d(eps)*d(eps)*d(Omega) -> [eV^2] [bohr*Ry]^2 [1/eV] [1/eV] [1/eV]
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
                  "# elec_smearing: " + ((a2f.is_tetra()) ? "tetrahedra\n" : (std::to_string(a2f.elec_smearing()) + " eV\n"))
                + "# phon_smearing: " + ((a2f.is_tetra()) ? "tetrahedra\n" : (std::to_string(a2f.phon_smearing()) + " eV\n"))
                + "# sign: " + std::to_string(a2f.sign())
                + "\n# sign_pr: " + std::to_string(a2f.sign_pr())
                + "\n# velocity component alpha: " + a2f.alpha()
                + "\n# velocity component beta: " + a2f.beta()
                + "\n# electron energy list [eV]: " + array1D_to_string(a2f.epsilons())
                + "\n# DOS for energy list [1/eV/spin/cell]: "
                + array1D_to_string(a2f.doses())
                + "\n# transport DOS for energy list [Ry^2*bohr^2/eV]: "
                + array1D_to_string(a2f.trans_doses())
                + "\n#\n# Frequency [eV]        Transport Spectral Function\n";
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
                os << std::setw(15) << std::setprecision(6) << vals[ieps];
            os << std::endl;
        }
    }

    if (!rank)
        os.close();
}

} // spectral
} // skies
