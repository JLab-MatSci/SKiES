#include "window.h"
#include "opengl_widget.h"

#include <QMenuBar>
#include <QVBoxLayout>
#include <QLineEdit>
#include <QLabel>
#include <QTimer>

SkiesApplicationWindow::SkiesApplicationWindow(QWidget *parent)
    : QMainWindow(parent) {

    setWindowTitle("SKiES");
    resize(800, 600);

    createMenuBar();

    createCentralWidget();

    statusBar = new QStatusBar(this);
    statusBar->showMessage("Ready");
    setStatusBar(statusBar);
}

void SkiesApplicationWindow::createMenuBar() {
    QMenuBar *menuBar = new QMenuBar(this);

    QMenu *fileMenu = menuBar->addMenu("&File");
    QAction *exitAction = fileMenu->addAction("Exit");
    connect(exitAction, &QAction::triggered, this, &QMainWindow::close);

    QMenu *helpMenu = menuBar->addMenu("&Help");
    QAction *aboutAction = helpMenu->addAction("About");
    connect(aboutAction, &QAction::triggered, this, [this]() {
        statusBar->showMessage("About clicked");
    });

    setMenuBar(menuBar);
}

void SkiesApplicationWindow::createCentralWidget() {
    QWidget *centralWidget = new QWidget(this);

    QVBoxLayout *layout = new QVBoxLayout(centralWidget);
    layout->setContentsMargins(2, 2, 2, 2);
    layout->setSpacing(5);

    QLineEdit *inputField = new QLineEdit(centralWidget);
    inputField->setPlaceholderText("Enter some data here");
    layout->addWidget(new QLabel("Input Data:", centralWidget));
    layout->addWidget(inputField);

    logOutput = new QTextEdit(centralWidget);
    logOutput->setReadOnly(true);
    logOutput->append("Log started...");
    layout->addWidget(new QLabel("Log Output:", centralWidget));
    layout->addWidget(logOutput);

    SkiesOpenGLWidget *openGLWidget = new SkiesOpenGLWidget(centralWidget);
    openGLWidget->setMinimumSize(400, 400);
    
    openGLWidget->setAttribute(Qt::WA_NativeWindow);
    openGLWidget->setAttribute(Qt::WA_DontCreateNativeAncestors);
    
    layout->addWidget(openGLWidget, 1);
    
    setCentralWidget(centralWidget);

    connect(inputField, &QLineEdit::returnPressed, this, [this, inputField]() {
        logOutput->append("Input: " + inputField->text());
        inputField->clear();
    });
}