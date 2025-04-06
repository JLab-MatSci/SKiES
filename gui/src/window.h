#pragma once

#include <QMainWindow>
#include <QTextEdit>
#include <QGraphicsView>
#include <QStatusBar>

class SkiesApplicationWindow: public QMainWindow {
    Q_OBJECT
public:
    explicit SkiesApplicationWindow(QWidget *parent = nullptr);

private:
    void createMenuBar();
    void createCentralWidget();

    QTextEdit *logOutput;
    QGraphicsView *graphicsView;
    QStatusBar *statusBar;
};