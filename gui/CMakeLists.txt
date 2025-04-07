set(SOURCES
    "src/window.cpp"
    "src/gl_widget.cpp"
    "src/crystal_gl_widget.cpp"
    "src/crystal_structure.cpp"
    "src/crystal_viewer_dialog.cpp"
    "src/editor_dock_widget.cpp"
)

set(HEADERS
  "src/window.h"
  "src/gl_widget.h"
  "src/crystal_gl_widget.h"
  "src/crystal_structure.h"
  "src/crystal_viewer_widget.h"
  "src/editor_dock_widget.h"
)

find_package(Qt5
  COMPONENTS
  Core
  Gui
  Widgets
#   Svg
#   PrintSupport
#   Xml
  OpenGL
  LinguistTools
  REQUIRED
)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

find_package(OpenGL COMPONENTS OpenGL REQUIRED)

add_executable(skies-gui
    "app_main.cpp"
    "${SOURCES}"
)

target_link_libraries(skies-gui
  OpenGL::GL
  OpenGL::GLU
  Qt5::Core
  Qt5::Gui
  # Qt5::PrintSupport
  Qt5::OpenGL
  # Qt5::Svg
  -lX11
  -lGL
)

install(TARGETS skies-gui
    	DESTINATION bin
    	EXPORT skies-targets
)