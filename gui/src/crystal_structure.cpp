#include "crystal_structure.h"
#include <cmath>

namespace skies {
namespace gui {
namespace crystals {

std::vector<Atom> CrystalStructure::generateAtoms() const
{
    std::vector<Atom> atoms;

    switch(type) {
    case LatticeType::bcc:
        atoms.push_back({0, 0, 0});
        atoms.push_back({0.5f, 0.5f, 0.5f});
        break;
    case LatticeType::fcc:
        atoms.push_back({0, 0, 0});
        atoms.push_back({0.5f, 0.5f, 0});
        atoms.push_back({0.5f, 0, 0.5f});
        atoms.push_back({0, 0.5f, 0.5f});
        break;
    case LatticeType::hcp:
        atoms.push_back({0, 0, 0});
        atoms.push_back({0.5f, static_cast<float>(sqrt(3) / 2), 0});
        atoms.push_back({0.5f, static_cast<float>(sqrt(3) / 6), static_cast<float>(c / 2)});
        atoms.push_back({0, static_cast<float>(sqrt(3) * 2 / 3), static_cast<float>(c / 2)});
        break;
    }

    for(auto& atom : atoms) {
        atom.x *= a;
        atom.y *= a;
        atom.z *= (type == LatticeType::hcp) ? c : a;
    }

    return atoms;
}

std::vector<Bond> CrystalStructure::generateBonds() const {
    std::vector<Bond> bonds;
    auto atoms = generateAtoms();
    const float cutoff = (type == LatticeType::hcp) ? std::min(a,c) * 1.1f : a * 1.1f;

    for(size_t i = 0; i < atoms.size(); ++i) {
        for(size_t j = i + 1; j < atoms.size(); ++j) {
            QVector3D delta = atoms[i].position() - atoms[j].position();
            if (delta.length() < cutoff) {
                bonds.push_back({atoms[i].position(), atoms[j].position()});
            }
        }
    }

    return bonds;
}

} // crystals
} // gui
} // skies
