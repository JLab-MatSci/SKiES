#include "window.h"
#include "gl_widget.h"
#include "crystal_viewer_dialog.h"

#include <QMenuBar>
#include <QVBoxLayout>
#include <QLineEdit>
#include <QLabel>
#include <QTimer>

namespace skies {
namespace gui {

SkiesApplicationWindow::SkiesApplicationWindow(QWidget *parent)
    : QMainWindow(parent)
    , logOutput(nullptr)
    , graphicsView(nullptr)
    , statusBar(new QStatusBar(this))
    , openGLWidget(nullptr)
    , crysViewport(nullptr)
    , editorDock(nullptr)
{
    setWindowTitle("SKiES");
    resize(800, 600);

    createMenuBar();

    createCentralWidget();

    statusBar->showMessage("Ready");
    setStatusBar(statusBar);
}

void SkiesApplicationWindow::setupMenu() {
    QMenu *viewMenu = menuBar()->addMenu(tr("&View"));

    showCrystalViewerAct = viewMenu->addAction(tr("&Crystal Viewer"));
    showCrystalViewerAct->setShortcut(tr("Ctrl+V"));
    showCrystalViewerAct->setStatusTip(tr("Open crystal structure viewer"));
    connect(showCrystalViewerAct, &QAction::triggered, 
            this, &SkiesApplicationWindow::openCrystalViewer);
}

void SkiesApplicationWindow::openCrystalViewer() {
    crystals::CrystalViewerDialog dialog(this);
    dialog.exec();
}

void SkiesApplicationWindow::createMenuBar() {
    QMenuBar *menuBar = new QMenuBar(this);

    QMenu *fileMenu = menuBar->addMenu(tr("&File"));
    QAction *exitAction = fileMenu->addAction(tr("Exit"));
    connect(exitAction, &QAction::triggered, this, &QMainWindow::close);

    QMenu *helpMenu = menuBar->addMenu(tr("&Help"));
    QAction *aboutAction = helpMenu->addAction(tr("About"));
    connect(aboutAction, &QAction::triggered, this, [this]() {
        statusBar->showMessage(tr("About clicked"));
    });

    setMenuBar(menuBar);

    setupMenu();
}

void SkiesApplicationWindow::createCentralWidget() {
    QWidget *centralWidget = new QWidget(this);
    QVBoxLayout *layout = new QVBoxLayout(centralWidget);
    layout->setContentsMargins(2, 2, 2, 2);
    layout->setSpacing(5);

    QLineEdit *inputField = new QLineEdit(centralWidget);
    inputField->setPlaceholderText(tr("Enter some data here"));
    layout->addWidget(new QLabel(tr("Input Data:"), centralWidget));
    layout->addWidget(inputField);

    logOutput = new QTextEdit(centralWidget);
    logOutput->setReadOnly(true);
    logOutput->append(tr("Log started..."));
    layout->addWidget(new QLabel(tr("Log Output:"), centralWidget));
    layout->addWidget(logOutput);

    connect(inputField, &QLineEdit::returnPressed, this, [this, inputField]() {
        logOutput->append(tr("Input: ") + inputField->text());
        inputField->clear();
    });

    setCentralWidget(centralWidget);
}

} // gui
} // skies
