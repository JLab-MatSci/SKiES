#pragma once

#include "crystal_gl_widget.h"
#include "editor_dock_widget.h"

#include <QDialog>
#include <QDialogButtonBox>

namespace skies {
namespace gui {
namespace crystals {

class CrystalViewerDialog : public QDialog {
    Q_OBJECT
public:
    explicit CrystalViewerDialog(QWidget *parent = nullptr);
    
private:
    void setupUI();
    
    CrystalGLWidget *viewport;
    EditorDockWidget *editorDock;
    QDialogButtonBox *buttons;
};

} // crystals
} // gui
} // skies
