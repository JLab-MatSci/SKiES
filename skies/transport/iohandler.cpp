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
#include <fstream>

#include <skies/common/alg.h>
#include <skies/common/units.h>
#include <skies/transport/iohandler.h>

#include <iomanip>

#include <iostream>

namespace skies { namespace transport {

using namespace arrays;

using skies::units::rau_in_m_over_s;
using skies::units::Ry_in_J;

IHandler::IHandler(const char* a2f_fnm)
    : ifs_(a2f_fnm)
{
    std::string line;
    bool header_flag = false;

    if (ifs_.fail())
        throw std::runtime_error("The input file " + std::string{a2f_fnm} + " does not exist");

    while (ifs_.good())
    {
        while (!header_flag)
        {
            getline(ifs_, line);
            if (line.find("elec_smearing:") != std::string::npos)
            {
                std::string word_or_num = custom_split(line, ' ')[2].data();
                if (word_or_num == "tetrahedra")
                    elec_smearing_ = 0.0;
                else
                    elec_smearing_ = std::stod(word_or_num);
            }
            if (line.find("phon_smearing:") != std::string::npos)
            {
                std::string word_or_num = custom_split(line, ' ')[2].data();
                if (word_or_num == "tetrahedra")
                    phon_smearing_ = 0.0;
                else
                    phon_smearing_ = std::stod(word_or_num);
            }
            if (line.find("sign:") != std::string::npos)
                sign_ = std::stoi(custom_split(line, ' ').back().data(), 0);
            if (line.find("sign_pr:") != std::string::npos)
                sign_pr_ = std::stoi(custom_split(line, ' ').back().data(), 0);
            if (line.find("alpha") != std::string::npos)
                alpha_ = *custom_split(line, ' ').back().data();
            if (line.find("beta") != std::string::npos)
                beta_ = *custom_split(line, ' ').back().data();
            if (line.find("electron") != std::string::npos)
            {
                auto splitted_line = custom_split(line, ' ');
                for (size_t i = 5; i < splitted_line.size(); ++i)
                    epsilons_.push_back(std::strtod(splitted_line[i].data(), 0));
            }
            if (line.find("transport") != std::string::npos)
            {
                auto splitted_line = custom_split(line, ' ');
                for (size_t i = 7; i < splitted_line.size(); ++i)
                    transDOSes_.push_back(std::strtod(splitted_line[i].data(), 0));
            }
            if (line.find("Frequency") != std::string::npos)
                header_flag = true;
        }
        getline(ifs_, line);
        if (!line.empty())
        {
            auto splitted_line = custom_split(line, ' ');
            omegas_.push_back(std::strtod(splitted_line[0].data(), 0));
            array1D epsilons_at_omega;
            for (size_t ieps = 0; ieps < epsilons_.size(); ++ieps)
                epsilons_at_omega.push_back(std::strtod(splitted_line[ieps + 1].data(), 0));
            a2f_.push_back(epsilons_at_omega);
        }
    }
}

IHandler::~IHandler()
{
    ifs_.close();
}

OHandler::OHandler(const char* a2f_fnm, const char* cond_fnm, ResistType type, const array1D& ion_Temps)
    : ofs_(cond_fnm)
    , ion_Temps_(ion_Temps)
{
    IHandler ihandler(a2f_fnm);
    bool is_tetra = (ihandler.elec_smearing() == 0);
    switch (type)
    {
    case  ResistType::Electrical:
        if (is_tetra)
        {
            ofs_ << "#  elec_smearing: " << "tetrahedra" << std::endl;
            ofs_ << "#  phon_smearing: " << "tetrahedra" << std::endl;
        }
        else
        {
            ofs_ << "#  elec_smearing: " << ihandler.elec_smearing() << " [eV]" << std::endl;
            ofs_ << "#  phon_smearing: " << ihandler.elec_smearing() << " [eV]" << std::endl;
        }
        ofs_ << "#  Transport DOS per spin [r.a.u.] in energy list: ";
        for (size_t ieps = 0; ieps < ihandler.epsilons().size(); ++ieps)
            ofs_ << ihandler.transDOSes()[ieps] * units::Ry_in_eV << " ";
        ofs_ << "\n#" << std::endl;
        if (ion_Temps_.empty())
        {
            ofs_ << std::left << "# Te [K]" << std::setw(30) << "      Resistivity [muOm cm]" << std::endl;
        }
        else
        {
            ofs_ << "# Resistivity [muOhm cm]" << std::endl;
            ofs_ << "#  Te [K]  \\  Ti [K]:";
            for (auto&& Ti : ion_Temps_)
                ofs_ << /*std::setprecision(1) <<*/ std::setw(12) << Ti;
            ofs_ << std::endl;
        }
        break;

    case ResistType::Thermal:
        if (is_tetra)
        {
            ofs_ << "#  elec_smearing: " << "tetrahedra" << std::endl;
            ofs_ << "#  phon_smearing: " << "tetrahedra" << std::endl;
        }
        else
        {
            ofs_ << "#  elec_smearing: " << ihandler.elec_smearing() << " [eV]" << std::endl;
            ofs_ << "#  phon_smearing: " << ihandler.elec_smearing() << " [eV]" << std::endl;
        }
        ofs_ << "#  Transport DOS per spin [r.a.u.] in energy list: ";
        for (size_t ieps = 0; ieps < ihandler.epsilons().size(); ++ieps)
            ofs_ << ihandler.transDOSes()[ieps] * units::Ry_in_eV << " ";
        ofs_ << "\n#" << std::endl;
        if (ion_Temps_.empty())
        {
            ofs_ << "# Te [K]" << std::left << std::setw(35) << "      Thermal Resistivity [W/cm/K]" << std::endl;
        }
        else
        {
            ofs_ << "# Thermal Conductivity [W/cm/K]" << std::endl;
            ofs_ << "#  Te [K]  \\  Ti [K]:";
            for (auto&& Ti : ion_Temps_)
                ofs_ << std::setw(12) << Ti;
            ofs_ << std::endl;
        }
        break;

    default:
        break;
    }
}

OHandler::~OHandler()
{
    ofs_.close();
}

void OHandler::dump(const arrays::array1D& Temps, const arrays::array1D& resist)
{
    assert(resist.size() == Temps.size());
    for (size_t itemp = 0; itemp < Temps.size(); ++itemp)
            ofs_ << std::setprecision(6) << std::setw(14) << Temps[itemp]
                 << std::setprecision(6) << std::setw(14) << resist[itemp] << std::endl;
}

void OHandler::dump(const arrays::array1D& Temps, const arrays::array2D& resist)
{
    assert(resist.size() == Temps.size());
    assert((resist[0].size() == ion_Temps_.size()) || (ion_Temps_.size() == 0));
    if (ion_Temps_.empty())
        for (size_t itemp = 0; itemp < Temps.size(); ++itemp)
            ofs_ << std::setprecision(6) << std::setw(14) << Temps[itemp]
                 << std::setprecision(6) << std::setw(14) << resist[itemp][0] << std::endl;
    else
    {
        for (size_t itemp = 0; itemp < Temps.size(); ++itemp)
        {
            ofs_ << std::setw(12) << std::setprecision(6) << Temps[itemp] << "              ";
            for (size_t jtemp = 0; jtemp < ion_Temps_.size(); ++jtemp)
                ofs_ << std::setw(12) << std::setprecision(6) << resist[itemp][jtemp];
            ofs_ << std::endl;
        }
    }
}

} // transport
} // skies
