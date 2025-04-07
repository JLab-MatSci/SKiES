#pragma once

#include <QVector3D>
#include <vector>

namespace skies {
namespace gui {
namespace crystals {

enum class LatticeType {
    bcc,
    fcc,
    hcp
};

struct Atom {
    float x, y, z;
    QVector3D position() const { return QVector3D(x, y, z); }
};

struct Bond {
    QVector3D start, end;
};

class CrystalStructure {
public:
    LatticeType type = LatticeType::bcc;
    double a = 1.0;
    double c = 1.0;

    std::vector<Atom> generateAtoms() const;
    std::vector<Bond> generateBonds() const;
};

} // crystals
} // gui
} // skies

Q_DECLARE_METATYPE(skies::gui::crystals::LatticeType)
