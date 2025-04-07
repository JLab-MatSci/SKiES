#include "crystal_viewer_dialog.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QSlider>
#include <QPushButton>
#include <QColorDialog>
#include <QLabel>

namespace skies {
namespace gui {
namespace crystals {

CrystalViewerDialog::CrystalViewerDialog(QWidget *parent)
    : QDialog(parent)
{
    setWindowTitle("Crystal Structure Viewer");
    setMinimumSize(900, 600);

    setupUI();
}

void CrystalViewerDialog::setupUI() {
    QHBoxLayout *mainLayout = new QHBoxLayout(this);

    viewport = new CrystalGLWidget(this);

    editorDock = new EditorDockWidget(this);
    editorDock->setFeatures(QDockWidget::NoDockWidgetFeatures);

    connect(editorDock, &EditorDockWidget::structureUpdated,
            viewport, &CrystalGLWidget::updateStructure);

    buttons = new QDialogButtonBox(
        QDialogButtonBox::Close, 
        Qt::Horizontal, 
        this
    );
    connect(buttons, &QDialogButtonBox::rejected, 
            this, &QDialog::reject);
    
    QWidget *dockContainer = new QWidget(this);
    QVBoxLayout *dockLayout = new QVBoxLayout(dockContainer);

    QSlider *atomSizeSlider = new QSlider(Qt::Horizontal, this);
    atomSizeSlider->setRange(1, 500);
    atomSizeSlider->setValue(100);
    connect(atomSizeSlider, &QSlider::valueChanged, [this](int value) {
        viewport->setAtomSizeScale(value / 100.0f);
    });

    QPushButton *colorPickerButton = new QPushButton("Pick Atom Color", this);
    connect(colorPickerButton, &QPushButton::clicked, [this]() {
        QColor color = QColorDialog::getColor(Qt::red, this, "Select Atom Color");
        if (color.isValid()) {
            viewport->setAtomColor(color);
        }
    });

    dockLayout->addWidget(new QLabel("Atom Size Scale:"));
    dockLayout->addWidget(atomSizeSlider);
    dockLayout->addWidget(colorPickerButton);

    dockLayout->addWidget(editorDock);
    dockLayout->addWidget(buttons);
    dockLayout->setContentsMargins(0, 0, 0, 0);
    
    mainLayout->addWidget(viewport, 1);
    mainLayout->addWidget(dockContainer, 0);
    mainLayout->setContentsMargins(5, 5, 5, 5);
}

} // crystals
} // gui
} // skies
