#include "opengl_widget.h"

#include <math.h>

#include <QMouseEvent>
#include <QWheelEvent>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QTimer>

SkiesOpenGLWidget::SkiesOpenGLWidget(QWidget *parent)
    : QOpenGLWidget(parent), angularSpeed(0), zoomLevel(-5.0f)
{
    rotationAxis = QVector3D(0, 1, 0);

    QTimer *timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, [this](){
        angularSpeed = 0.1f;
        update();
    });
    timer->start(16);
}

void SkiesOpenGLWidget::initializeGL() {
    initializeOpenGLFunctions();
    if (!this->context()->isValid()) {
        qFatal("OpenGL context is invalid!");
    }
    
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    model.setToIdentity();
}

void SkiesOpenGLWidget::resizeGL(int w, int h) {
    updateProjectionMatrix();
    glViewport(0, 0, w, h);
}


void SkiesOpenGLWidget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (angularSpeed > 0) {
        qreal angle = angularSpeed * 16;
        rotation = QQuaternion::fromAxisAndAngle(rotationAxis, angle) * rotation;
        angularSpeed *= 0.99f;
    }
    model.setToIdentity();
    model.rotate(rotation);

    QMatrix4x4 mvp = projection * view * model;

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixf(mvp.data());

    glBegin(GL_QUADS);

    glColor3f(1.0f, 0.0f, 0.0f); glVertex3f(-1.0f, -1.0f, 1.0f);
    glColor3f(0.0f, 1.0f, 0.0f); glVertex3f( 1.0f, -1.0f, 1.0f);
    glColor3f(0.0f, 0.0f, 1.0f); glVertex3f( 1.0f,  1.0f, 1.0f);
    glColor3f(1.0f, 1.0f, 0.0f); glVertex3f(-1.0f,  1.0f, 1.0f);

    glColor3f(1.0f, 0.0f, 0.0f); glVertex3f(-1.0f, -1.0f, -1.0f);
    glColor3f(0.0f, 1.0f, 0.0f); glVertex3f(-1.0f,  1.0f, -1.0f);
    glColor3f(0.0f, 0.0f, 1.0f); glVertex3f( 1.0f,  1.0f, -1.0f);
    glColor3f(1.0f, 1.0f, 0.0f); glVertex3f( 1.0f, -1.0f, -1.0f);

    glColor3f(1.0f, 0.0f, 0.0f); glVertex3f(-1.0f, 1.0f, -1.0f);
    glColor3f(0.0f, 1.0f, 0.0f); glVertex3f(-1.0f, 1.0f,  1.0f);
    glColor3f(0.0f, 0.0f, 1.0f); glVertex3f( 1.0f, 1.0f,  1.0f);
    glColor3f(1.0f, 1.0f, 0.0f); glVertex3f( 1.0f, 1.0f, -1.0f);

    glColor3f(1.0f, 0.0f, 0.0f); glVertex3f(-1.0f, -1.0f, -1.0f);
    glColor3f(0.0f, 1.0f, 0.0f); glVertex3f( 1.0f, -1.0f, -1.0f);
    glColor3f(0.0f, 0.0f, 1.0f); glVertex3f( 1.0f, -1.0f,  1.0f);
    glColor3f(1.0f, 1.0f, 0.0f); glVertex3f(-1.0f, -1.0f,  1.0f);

    glColor3f(1.0f, 0.0f, 0.0f); glVertex3f(1.0f, -1.0f, -1.0f);
    glColor3f(0.0f, 1.0f, 0.0f); glVertex3f(1.0f,  1.0f, -1.0f);
    glColor3f(0.0f, 0.0f, 1.0f); glVertex3f(1.0f,  1.0f,  1.0f);
    glColor3f(1.0f, 1.0f, 0.0f); glVertex3f(1.0f, -1.0f,  1.0f);

    glColor3f(1.0f, 0.0f, 0.0f); glVertex3f(-1.0f, -1.0f, -1.0f);
    glColor3f(0.0f, 1.0f, 0.0f); glVertex3f(-1.0f, -1.0f,  1.0f);
    glColor3f(0.0f, 0.0f, 1.0f); glVertex3f(-1.0f,  1.0f,  1.0f);
    glColor3f(1.0f, 1.0f, 0.0f); glVertex3f(-1.0f,  1.0f, -1.0f);

    glEnd();
}

void SkiesOpenGLWidget::mousePressEvent(QMouseEvent *event) {
    lastMousePosition = event->pos();
}

void SkiesOpenGLWidget::mouseMoveEvent(QMouseEvent *event) {
    int dx = event->x() - lastMousePosition.x();
    int dy = event->y() - lastMousePosition.y();

    if (event->buttons() & Qt::LeftButton) {
        rotationAxis = QVector3D(dy, dx, 0).normalized();
        angularSpeed = std::sqrt(dx * dx + dy * dy) / 100.0f;
    }

    lastMousePosition = event->pos();
    update();
}

void SkiesOpenGLWidget::wheelEvent(QWheelEvent *event) {
    zoomLevel += event->angleDelta().y() / 120.0f * 0.5f;
    updateProjectionMatrix();
    update();
}

void SkiesOpenGLWidget::updateProjectionMatrix() {
    projection.setToIdentity();
    float aspectRatio = static_cast<float>(width()) / static_cast<float>(height());
    projection.perspective(45.0f, aspectRatio, 0.1f, 100.0f);
    view.setToIdentity();
    view.translate(0, 0, zoomLevel);
}