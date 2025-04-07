#pragma once

#include <QDockWidget>
#include "crystal_structure.h"

class QComboBox;
class QDoubleSpinBox;
class QSlider;

namespace skies {
namespace gui {
namespace crystals {

class EditorDockWidget : public QDockWidget {
    Q_OBJECT
public:
    explicit EditorDockWidget(QWidget *parent = nullptr);

signals:
    void structureUpdated(const CrystalStructure& structure);

private slots:
    void updateStructure();

private:
    void setupUI();

    CrystalStructure structure;
    QComboBox *latticeCombo;
    QDoubleSpinBox *aParam, *cParam;
    QSlider *aSlider, *cSlider;
};

} // crystals
} // gui
} // skies
