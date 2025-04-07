#include "crystal_gl_widget.h"

#include <QOpenGLFunctions>

#include <GL/glu.h>

#include <cmath>

namespace skies {
namespace gui {
namespace crystals {

CrystalGLWidget::CrystalGLWidget(QWidget *parent) 
    : SkiesOpenGLWidget(parent)
    , sphereDisplayList(0)
    , atomSizeScale(1.0f)
    , atomColor(Qt::red)
{}

void CrystalGLWidget::setAtomSizeScale(float scale) {
    atomSizeScale = scale;
    update();
}

void CrystalGLWidget::setAtomColor(const QColor& color) {
    atomColor = color;
    update();
}

void CrystalGLWidget::initializeGL() {
    SkiesOpenGLWidget::initializeGL();

    sphereDisplayList = glGenLists(1);
    glNewList(sphereDisplayList, GL_COMPILE);
    {
        GLUquadric* quad = gluNewQuadric();
        gluSphere(quad, 1.0, 16, 16);
        gluDeleteQuadric(quad);
    }
    glEndList();

    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
}

void CrystalGLWidget::paintGL() {
    SkiesOpenGLWidget::paintGL();

    glDisable(GL_LIGHTING);

    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    
    glColor3f(atomColor.redF(), atomColor.greenF(), atomColor.blueF());

    float atomRadius = 0.15f * currentStructure.a * atomSizeScale;
    for(const auto& atom : atoms) {
        drawAtom(atom.position(), atomRadius);
    }
    
    if(showBonds) {
        glDisable(GL_LIGHTING);
        glColor3f(0.2f, 0.8f, 0.2f);
        for(const auto& bond : bonds) {
            drawBond(bond.start, bond.end);
        }
        glEnable(GL_LIGHTING);
    }
}

void CrystalGLWidget::drawAtom(const QVector3D& pos, float radius) {
    glPushMatrix();
    glTranslatef(pos.x(), pos.y(), pos.z());
    glScalef(radius, radius, radius);
    glCallList(sphereDisplayList);
    glPopMatrix();
}

void CrystalGLWidget::drawBond(const QVector3D& start, const QVector3D& end) {
    glLineWidth(bondThickness);
    glBegin(GL_LINES);
    glVertex3f(start.x(), start.y(), start.z());
    glVertex3f(end.x(), end.y(), end.z());
    glEnd();
}

void CrystalGLWidget::updateStructure(const CrystalStructure& structure) {
    currentStructure = structure;
    atoms = structure.generateAtoms();
    bonds = structure.generateBonds();
    update();
}

} // crystals
} // gui
} // skies
