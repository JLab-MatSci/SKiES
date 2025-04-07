#include <QApplication>
#include <QSurfaceFormat>
#include <QOpenGLContext>
#include <QMessageBox>

#include "src/window.h"

int main(int argc, char *argv[]) {
    QSurfaceFormat format;
    format.setVersion(2, 1);
    format.setProfile(QSurfaceFormat::CompatibilityProfile);
    QSurfaceFormat::setDefaultFormat(format);

    QApplication app(argc, argv);

    skies::gui::SkiesApplicationWindow win;
    win.show();

    return app.exec();
}