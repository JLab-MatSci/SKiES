#pragma once

#include <skies/common/ndimarrays.h>
#include <skies/sampling/sampling.h>

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
               char cart);

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
               char cart);

    // constructor from lambda_tr.dat file
    SpecFunc(const std::string& fname);

    ~SpecFunc();

    arrays::array1D calc_spec_func(double Omega);
               
private:
    void init(const std::vector<size_t>& kpgrid,
              const std::vector<size_t>& qpgrid,
              bzsampling::SamplingFunc elec_sampling,
              bzsampling::SamplingFunc phon_sampling,
              double elec_smearing,
              double phon_smearing);

    bzsampling::SamplingFunc elec_sampling_;
    bzsampling::SamplingFunc phon_sampling_;
    double elec_smearing_;
    double phon_smearing_;

    // signs in brackets with velocities
    // sign = +1 <=> a^2F(+1, \Omega)
    // sign_pr is primed counterpart
    int sign_;
    int sign_pr_;

    size_t nkpt_;
    size_t nkx_, nky_, nkz_;
    arrays::array2D kpts_;

    size_t nqpt_;
    size_t nqx_, nqy_, nqz_;
    arrays::array2D qpts_;

    size_t nbnd_;
    arrays::array2D eigenens_;
    arrays::array2D elvelocs_;

    size_t nmds_;
    arrays::array2D eigenfreqs_;

    std::map<arrays::array1D, int> cached_indices_;

    arrays::array3D inner_sum_;
    bool is_full_{ false };
    size_t iq_cont_{ 0 };
    size_t imd_cont_{ 0 };

    arrays::array1D epsilons_;
    arrays::array1D trDOSes_;

    int low_band_{ 0 };
    int high_band_;

    double Te_{ 0.258 }; // electronic temperature in [eV], needed for trDOS smearing: default corr. to 3000 K

    char cart_{ 'x' }; // cartesian index of velocities

public:
    double elec_smearing() const { return elec_smearing_; }
    double phon_smearing() const { return phon_smearing_; }
    arrays::array1D trans_doses() const { return trDOSes_; }

    int nkx() const { return nkx_; }
    int nky() const { return nky_; }
    int nkz() const { return nkz_; }
    int nqx() const { return nqx_; }
    int nqy() const { return nqy_; }
    int nqz() const { return nqz_; }
    size_t nkpt() const { return nkpt_; }
    size_t nqpt() const { return nqpt_; }

    size_t nmds() const { return nmds_; }
    int sign()    const { return sign_; }
    int sign_pr() const { return sign_pr_; }
    char cart() const { return cart_; }
    arrays::array1D& epsilons() { return epsilons_; } 
    const arrays::array1D& epsilons() const { return epsilons_; }

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

    bool is_continue_calc() { return is_continue_calc_; }

private:
    double calc_inner_sum(size_t iq, size_t imd, size_t ieps);
    double calc_inner_sum_in_subarray(size_t iq, size_t imd, size_t ieps, size_t start, size_t finish, size_t low_band, size_t high_band);
    arrays::array1D calc_exter_sum(double Omega);
    
    std::string type_of_el_smear_;
    std::string type_of_ph_smear_;

    bool is_continue_calc_{ false };
};

// main driver functions
void calc_spec_func(SpecFunc& a2f, const arrays::array1D& omegas, const std::string& fname);

} // spectral 
} // skies
