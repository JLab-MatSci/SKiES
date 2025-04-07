#pragma once

#include "crystal_gl_widget.h"
#include "editor_dock_widget.h"

#include <QMainWindow>
#include <QTextEdit>
#include <QGraphicsView>
#include <QStatusBar>

namespace skies { namespace gui {

class SkiesApplicationWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit SkiesApplicationWindow(QWidget *parent = nullptr);

private slots:
    void openCrystalViewer();

private:
    void setupMenu();
    void createMenuBar();
    void createCentralWidget();

    QTextEdit *logOutput;
    QGraphicsView *graphicsView;
    QStatusBar *statusBar;
    SkiesOpenGLWidget *openGLWidget;

    crystals::CrystalGLWidget *crysViewport;
    crystals::EditorDockWidget *editorDock;

    QAction *showCrystalViewerAct;
};

} // gui
} // skies
