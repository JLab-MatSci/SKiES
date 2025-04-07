#pragma once

#include "gl_widget.h"
#include "crystal_structure.h"

namespace skies {
namespace gui {
namespace crystals {

class CrystalGLWidget : public SkiesOpenGLWidget {
    Q_OBJECT
public:
    explicit CrystalGLWidget(QWidget *parent = nullptr);

    void updateStructure(const CrystalStructure& structure);

    void setAtomSizeScale(float scale);
    void setAtomColor(const QColor& color);

protected:
    void initializeGL() override;
    void paintGL() override;

private:
    CrystalStructure currentStructure;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    bool showBonds = true;
    float bondThickness = 2.0f;
    GLuint sphereDisplayList;
    float atomSizeScale;
    QColor atomColor;

    void drawAtom(const QVector3D& pos, float radius);
    void drawBond(const QVector3D& start, const QVector3D& end);
};

} // crystals
} // gui
} // skies
