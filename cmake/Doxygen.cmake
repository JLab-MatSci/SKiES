function(CreateDocs in out latex)
  find_package(Doxygen)
  if (NOT DOXYGEN_FOUND)
    add_custom_target(dox COMMAND FALSE
      COMMENT "Doxygen is not found")
    return()
  endif()
  set(DOXYGEN_GENERATE_HTML YES)
  set(DOXYGEN_HTML_OUTPUT ${PROJECT_BINARY_DIR}/${out})
  set(DOXYGEN_GENERATE_LATEX YES)
  set(DOXYGEN_LATEX_OUTPUT ${PROJECT_BINARY_DIR}/${latex})
  UseDoxygenAwesomeCss()
  doxygen_add_docs(dox ${PROJECT_SOURCE_DIR}/${in})
endfunction()

macro(UseDoxygenAwesomeCss)
 include(FetchContent)
 FetchContent_Declare(doxygen-awesome-css
      GIT_REPOSITORY
      https://github.com/jothepro/doxygen-awesome-css.git
      GIT_TAG
      v2.3.3
 )
 FetchContent_MakeAvailable(doxygen-awesome-css)
 set(DOXYGEN_GENERATE_TREEVIEW YES)
 set(DOXYGEN_HAVE_DOT YES)
 set(DOXYGEN_DOT_IMAGE_FORMAT svg)
 set(DOXYGEN_DOT_TRANSPARENT YES)
 set(DOXYGEN_DISABLE_INDEX NO)
 set(DOXYGEN_SIDEBAR NO)
 set(DOXYGEN_HTML_COLORSTYLE LIGHT)
 set(DOXYGEN_HTML_EXTRA_STYLESHEET
      ${doxygen-awesome-css_SOURCE_DIR}/doxygen-awesome.css
      ${doxygen-awesome-css_SOURCE_DIR}/doxygen-awesome-sidebar-only.css)
endmacro()