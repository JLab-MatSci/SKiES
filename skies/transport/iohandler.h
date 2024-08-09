#pragma once

#include <fstream>

#include <skies/common/ndimarrays.h>

namespace skies { namespace transport {

enum class ResistType 
{
    Electrical,
    Thermal,
    Seebeck
};

class IHandler {
public:
    IHandler(const char* a2f_fnm);
    ~IHandler();

    const arrays::array2D& a2f() { return a2f_; };
    const arrays::array1D& omegas() { return omegas_; }
    const arrays::array1D& epsilons() { return epsilons_; }
    arrays::array1D&  transDOSes() { return transDOSes_; }
    double elec_smearing() { return elec_smearing_; }
    double phon_smearing() { return phon_smearing_; }
    int    cartes_ind() { return cartes_ind_; }

private:
    std::ifstream ifs_;

    double elec_smearing_;
    double phon_smearing_;
    int    cartes_ind_;

    arrays::array1D epsilons_;
    arrays::array1D transDOSes_;
    arrays::array1D omegas_;
    arrays::array2D a2f_;
};

class OHandler {
public:
    OHandler(const char* a2f_fnm,  const char* cond_fnm, ResistType type, const arrays::array1D& ion_Temps = arrays::array1D());
    ~OHandler();

    void dump(const arrays::array1D& Temps, const arrays::array1D& resist);
    void dump(const arrays::array1D& Temps, const arrays::array2D& resist);

private:
    std::ofstream ofs_;
    arrays::array1D ion_Temps_;
};

} // transport
} // skies
