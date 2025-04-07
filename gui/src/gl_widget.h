#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector3D>

#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>

namespace skies {
namespace gui {

class SkiesOpenGLWidget : public QOpenGLWidget, protected QOpenGLFunctions {
    Q_OBJECT

public:
    explicit SkiesOpenGLWidget(QWidget *parent = nullptr);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;

    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;

private:
    void updateProjectionMatrix();

    QMatrix4x4 projection;
    QMatrix4x4 view;
    QMatrix4x4 model;

    QVector3D rotationAxis;
    qreal angularSpeed;
    QQuaternion rotation;

    QPoint lastMousePosition;
    float zoomLevel;
};

} // gui
} // skies
