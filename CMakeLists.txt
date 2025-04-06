cmake_policy(SET CMP0048 NEW)
project(SKiES VERSION 1.0.0 LANGUAGES CXX C Fortran)
cmake_minimum_required(VERSION 3.15)

if(PROJECT_BINARY_DIR STREQUAL PROJECT_SOURCE_DIR)
    message(FATAL_ERROR
    "Please never do in-source builds.\n"
    "The correct way to build a Cmake project is\n"
    "\t cmake -B <path-to-build-dir>\n")
endif()

set(CMAKE_CXX_FLAGS_DEBUG      "-g")
set(CMAKE_CXX_FLAGS_RELEASE   "-O3")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Default build type at CMake configure time")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SKIES_VERSION_MAJOR      1)
set(SKIES_VERSION_MINOR      0)
set(SKIES_VERSION_PATCH      0)
set(SKIES_VERSION "${SKIES_VERSION_MAJOR},${SKIES_VERSION_MINOR},${SKIES_VERSION_PATCH}")

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    SET(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -ggdb -fno-omit-frame-pointer -Wno-deprecated-copy -Wno-ignored-attributes")
    add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
    add_compile_options(-Wno-error=parentheses)
endif()

####====================================================####
                   # Define basic options
####====================================================####

option(SKIES_ENABLE_DOC                   "Generate docs" OFF)
option(SKIES_ENABLE_TBB            "Parallel TBB version"  ON)
option(SKIES_ENABLE_TEST                  "Tests support" OFF)
option(SKIES_ENABLE_GUI                     "GUI support" OFF)

find_package(Mbd REQUIRED)
find_package(qe REQUIRED)

# add TBB library to use C++17 parallelization facilities
if (SKIES_ENABLE_TBB)
    find_package(TBB REQUIRED)
    add_definitions(-DSKIES_TBB)
endif()

if (SKIES_ENABLE_DOC)
    include(Doxygen)
    CreateDocs(skies docs docs/latex)
    install(DIRECTORY
        ${PROJECT_BINARY_DIR}/docs
        DESTINATION ${CMAKE_INSTALL_PREFIX}
    )
endif()

# check of C/Fortran compilers compatibility is important
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
include_directories(${BLAS_INCLUDE_DIR})
include_directories(${LAPACK_INCLUDE_DIR})

add_subdirectory(skies)
add_subdirectory(launch)

if (SKIES_ENABLE_GUI)
add_subdirectory(gui)
endif()

if (SKIES_ENABLE_TEST)
    add_subdirectory(tests)
endif()

include(CMakePackageConfigHelpers)
configure_package_config_file(
          ${CMAKE_CURRENT_SOURCE_DIR}/cmake/skiesConfig.cmake.in
          "${CMAKE_CURRENT_BINARY_DIR}/skiesConfig.cmake"
          INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/share/cmake
          NO_SET_AND_CHECK_MACRO
          NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/skiesConfig.cmake"
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/cmake
)
