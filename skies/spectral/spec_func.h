#pragma once

#include <skies/common/ndimarrays.h>
#include <skies/sampling/sampling.h>
#include <skies/sampling/tetrahedra.h>

namespace skies { namespace spectral {

class SpecFunc {
public:
    SpecFunc(const std::vector<size_t>& kpgrid,
               const std::vector<size_t>& qpgrid,
               bzsampling::SamplingFunc elec_sampling,
               bzsampling::SamplingFunc phon_sampling,
               double elec_smearing,
               double phon_smearing,
               int sign,
               double Te,
               char alpha,
               char beta,
               bool is_tetra = false);

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

    // constructor from lambda_tr.dat file
    SpecFunc(const std::string& fname);

    ~SpecFunc();

    arrays::array1D calc_spec_func(double Omega);
               
private:
    void init();
    void prepare_squared_velocs(char cart, arrays::array2D& elvelocs, arrays::array2D& elvelocs_sq);
    void prepare_fsh(const arrays::array2D& elvelocs, const arrays::array2D& elvelocs_sq, arrays::array2D& fsh);
    void prepare_fsh(const tetrahedra::TetraHandler& th_dos, const tetrahedra::TetraHandler& th_trdos,
                     const arrays::array2D& elvelocs, arrays::array2D& fsh);

    KPprotocol kprot_;
    KPprotocol qprot_;

    // signs in brackets with velocities
    // sign = +1 <=> a^2F(+1, \Omega)
    // sign_pr is primed counterpart
    int sign_;
    int sign_pr_;

    size_t nbnd_;
    arrays::array2D eigenens_;

    size_t nmds_;
    arrays::array2D eigenfreqs_;

    arrays::array2D elvelocs_alpha_;
    arrays::array2D elvelocs_beta_;
    arrays::array2D elvelocs_alpha_sq_;
    arrays::array2D elvelocs_beta_sq_;

    skies::tetrahedra::TetraHandler th_dos_;
    skies::tetrahedra::TetraHandler th_trdos_alpha_;
    skies::tetrahedra::TetraHandler th_trdos_beta_;

    // Fermi Surface Harmonics at (n,k)-grid
    arrays::array2D fsh_alpha_;
    arrays::array2D fsh_beta_;

    size_t low_band_{ 0 };
    size_t high_band_{ quantities::EigenValue::nbands - 1 };

    arrays::array1D epsilons_;
    arrays::array1D DOSes_;
    arrays::array1D trDOSes_;

    arrays::array3D inner_sum_;
    bool is_full_{ false };
    size_t iq_cont_{ 0 };
    size_t imd_cont_{ 0 };

    double Te_{ 0.258 }; // electronic temperature in [eV], needed for trDOS smearing: default corr. to 3000 K
    char alpha_{ 'x' }; // cartesian index of velocities
    char beta_{ 'x' };  // cartesian index of velocities

    bool is_tetra_{ false };
    bool is_continue_calc_{ false };

    bzsampling::SamplingFunc elec_sampling_;
    bzsampling::SamplingFunc phon_sampling_;
    double elec_smearing_;
    double phon_smearing_;

    std::string type_of_el_smear_;
    std::string type_of_ph_smear_;
public:
    double elec_smearing() const { return elec_smearing_; }
    double phon_smearing() const { return phon_smearing_; }

    size_t nkx() const { return std::get<0>(kprot_.mesh()); }
    size_t nky() const { return std::get<1>(kprot_.mesh()); }
    size_t nkz() const { return std::get<2>(kprot_.mesh()); }
    size_t nqx() const { return std::get<0>(qprot_.mesh()); }
    size_t nqy() const { return std::get<1>(qprot_.mesh()); }
    size_t nqz() const { return std::get<2>(qprot_.mesh()); }

    size_t nmds() const { return nmds_; }
    int sign()    const { return sign_; }
    int sign_pr() const { return sign_pr_; }
    char alpha() const { return alpha_; }
    char beta()  const { return beta_; }
    arrays::array1D& epsilons() { return epsilons_; } 
    const arrays::array1D& epsilons() const { return epsilons_; }
    const arrays::array1D& doses() const { return DOSes_; }
    const arrays::array1D& trans_doses() const { return trDOSes_; }

    void set_type_of_el_smear(const std::string& type) { type_of_el_smear_ = type; }
    void set_type_of_ph_smear(const std::string& type) { type_of_ph_smear_ = type; }
    std::string get_type_of_el_smear() const { return type_of_el_smear_; }
    std::string get_type_of_ph_smear() const { return type_of_ph_smear_; }

    void set_low_band(int low_band) { low_band_ = low_band; }
    void set_high_band(int high_band) { high_band_ = high_band; }
    int  get_low_band() const { return low_band_; }
    int  get_high_band() const { return high_band_; }

    void   set_elec_temp(double Te) { Te_ = Te; }
    double get_elec_temp() const { return Te_; }

    bool is_continue_calc() const { return is_continue_calc_; }
    bool is_tetra() const { return is_tetra_; }

private:
    double calc_inner_sum(size_t iq, size_t imd, size_t ieps);
    double calc_inner_sum_in_subarray(size_t iq, size_t imd, size_t ieps, size_t low_band, size_t high_band);
    double calc_inner_sum_in_subarray_tetra(size_t iq, size_t imd, size_t ieps, size_t low_band, size_t high_band);
    arrays::array1D calc_exter_sum(double Omega);

    //void dump_trdos_file();
};

// main driver functions
void calc_spec_func(SpecFunc& a2f, const arrays::array1D& omegas, const std::string& fname);

} // spectral 
} // skies
