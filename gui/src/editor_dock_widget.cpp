#include "editor_dock_widget.h"
#include "crystal_structure.h"

#include <QFormLayout>
#include <QLabel>
#include <QComboBox>
#include <QDoubleSpinBox>

namespace skies {
namespace gui {
namespace crystals {

EditorDockWidget::EditorDockWidget( QWidget *parent)
    : QDockWidget("Crystal Parameters", parent)
{
    setObjectName("EditorDock");
    setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    setFeatures(DockWidgetMovable | DockWidgetFloatable);

    QWidget *content = new QWidget(this);
    setWidget(content);

    setupUI();
}

void EditorDockWidget::setupUI() {
    QWidget *content = widget();
    QFormLayout *form = new QFormLayout(content);

    qRegisterMetaType<LatticeType>("LatticeType");

    latticeCombo = new QComboBox(content);
    latticeCombo->addItem("BCC", QVariant::fromValue(LatticeType::bcc));
    latticeCombo->addItem("FCC", QVariant::fromValue(LatticeType::fcc));
    latticeCombo->addItem("HCP", QVariant::fromValue(LatticeType::hcp));
    form->addRow("Lattice Type:", latticeCombo);
    
    aParam = new QDoubleSpinBox(content);
    aParam->setRange(0.1, 10.0);
    aParam->setValue(1.0);
    aParam->setSingleStep(0.1);
    form->addRow("a (Å):", aParam);
    
    aSlider = new QSlider(Qt::Horizontal, content);
    aSlider->setRange(10, 100);
    form->addRow(aSlider);
    
    cParam = new QDoubleSpinBox(content);
    cParam->setRange(0.1, 10.0);
    cParam->setValue(1.0);
    cParam->setSingleStep(0.1);
    form->addRow("c (Å):", cParam);
    
    cSlider = new QSlider(Qt::Horizontal, content);
    cSlider->setRange(10, 100);
    form->addRow(cSlider);

    connect(latticeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &EditorDockWidget::updateStructure);
    connect(aParam, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &EditorDockWidget::updateStructure);
    connect(cParam, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &EditorDockWidget::updateStructure);

    updateStructure();
}

void EditorDockWidget::updateStructure() {
    structure.type = static_cast<LatticeType>(latticeCombo->currentData().toInt());
    structure.a = aParam->value();
    structure.c = cParam->value();

    bool isHCP = (structure.type == LatticeType::hcp);
    cParam->setVisible(isHCP);
    cSlider->setVisible(isHCP);

    emit structureUpdated(structure);
}

} // crystals
} // gui
} // skies
