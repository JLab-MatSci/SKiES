#include <fstream>

#include <skies/common/alg.h>
#include <skies/common/units.h>
#include <skies/transport/iohandler.h>

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
        throw std::runtime_error("The input file does not exist");

    while (ifs_.good())
    {
        while (!header_flag)
        {
            getline(ifs_, line);
            if (line.find("elec_smearing:") != std::string::npos)
                elec_smearing_ = std::strtod(custom_split(line, ' ')[2].data(), 0);
            if (line.find("phon_smearing:") != std::string::npos)
                phon_smearing_ = std::strtod(custom_split(line, ' ')[2].data(), 0);
            if (line.find("velocity") != std::string::npos)
                cartes_ind_ = std::strtod(custom_split(line, ' ').back().data(), 0);
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
            if (line.find("frequency") != std::string::npos)
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

OHandler::OHandler(const char* a2f_fnm, const char* cond_fnm, ResistType type)
    : ofs_(cond_fnm)
{
    IHandler ihandler(a2f_fnm);
    switch (type)
    {
    case  ResistType::Electrical:
        ofs_ << "#  elec_smearing: " << ihandler.elec_smearing() << " [eV]" << std::endl;
        ofs_ << "#  phon_smearing: " << ihandler.phon_smearing() << " [eV]" << std::endl;
        ofs_ << "#  velocity component: " << ihandler.cartes_ind() << std::endl;
        ofs_ << "#  Transport DOS per spin [r.a.u.] in energy list: ";
        for (size_t ieps = 0; ieps < ihandler.epsilons().size(); ++ieps) ofs_ << ihandler.transDOSes()[ieps];
        ofs_ << "\n#" << std::endl;
        ofs_ << "#  Temperature [K]           Resistivity [muOm cm]" << std::endl;
        break; 

    case ResistType::Thermal:
        ofs_ << "#  elec_smearing: " << ihandler.elec_smearing() << " eV" << std::endl;
        ofs_ << "#  phon_smearing: " << ihandler.phon_smearing() << " eV" << std::endl;
        ofs_ << "#  velocity component: " << ihandler.cartes_ind() << std::endl;
        for (size_t ieps = 0; ieps < ihandler.epsilons().size(); ++ieps) ofs_ << ihandler.transDOSes()[ieps];
        ofs_ << "\n#" << std::endl;
        ofs_ << "#  Temperature [K]           Thermal Resistivity [W/cm/K]" << std::endl;
        break;

    case ResistType::Seebeck:
        ofs_ << "#  elec_smearing: " << ihandler.elec_smearing() << " eV" << std::endl;
        ofs_ << "#  phon_smearing: " << ihandler.phon_smearing() << " eV" << std::endl;
        ofs_ << "#  velocity component: " << ihandler.cartes_ind() << std::endl;
        for (size_t ieps = 0; ieps < ihandler.epsilons().size(); ++ieps) ofs_ << ihandler.transDOSes()[ieps];
        ofs_ << "\n#" << std::endl;
        ofs_ << "#  Temperature [K]           Seebeck coefficient [muV/K]" << std::endl;
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
    for (size_t itemp = 0; itemp < Temps.size(); ++itemp)
        ofs_ << Temps[itemp] << "    " << resist[itemp] << std::endl;
}

} // transport
} // skies
