cmake_minimum_required(VERSION 3.13)
project(miniproject CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
add_compile_options(-Wall)

# Point at emdw's source tree (for headers) and prebuilt shared library.
# Override with -DEMDW_ROOT=... or -DEMDW_BUILD=... if either moves.
set(EMDW_ROOT  "$ENV{HOME}/devel/emdw"        CACHE PATH "emdw source tree")
set(EMDW_BUILD "${EMDW_ROOT}/build"           CACHE PATH "emdw build tree")

if(NOT EXISTS "${EMDW_BUILD}/src/libemdw.so")
  message(FATAL_ERROR
    "libemdw.so not found at ${EMDW_BUILD}/src/libemdw.so — build emdw first "
    "(cd ${EMDW_BUILD} && cmake --build .)")
endif()

# Reuse emdw's Find modules and the same third-party deps it uses.
list(APPEND CMAKE_MODULE_PATH "${EMDW_ROOT}/cmake/Modules")
find_package(Boost    REQUIRED)
find_package(gLinear  REQUIRED)
find_package(prlite   REQUIRED)

# emdw's public headers live throughout its src/ tree.
set(EMDW_SRC "${EMDW_ROOT}/src")
include_directories(
  ${EMDW_SRC}
  ${EMDW_SRC}/emdw-base
  ${EMDW_SRC}/emdw-beliefprop
  ${EMDW_SRC}/emdw-clustering
  ${EMDW_SRC}/emdw-factors
  ${EMDW_SRC}/emdw-graphing
  ${EMDW_SRC}/utils-dependancyglue
  ${EMDW_SRC}/utils-informationtheory
  ${EMDW_SRC}/utils-matrix
  ${EMDW_SRC}/utils-misc
  ${EMDW_SRC}/utils-sigmapoints
  ${EMDW_SRC}/utils-vector
  ${EMDW_SRC}/vmp
  ${BOOST_INCLUDE_DIRS}
  ${GLINEAR_INCLUDE_DIRS}
  ${PRLITE_INCLUDE_DIRS}
)

# Import the prebuilt emdw shared library as a target.
add_library(emdw SHARED IMPORTED)
set_target_properties(emdw PROPERTIES
  IMPORTED_LOCATION "${EMDW_BUILD}/src/libemdw.so")
target_link_libraries(emdw INTERFACE
  ${GLINEAR_LIBRARIES} ${PRLITE_LIBRARIES} ${BOOST_LIBRARIES}
  lapack blas ssl crypto z)

add_executable(miniproject miniproject.cc miniproject_helpers.cc)
target_link_libraries(miniproject PRIVATE emdw)
target_compile_definitions(miniproject PRIVATE
  MINIPROJECT_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(miniproject_q3 miniproject_q3.cc miniproject_helpers.cc)
target_link_libraries(miniproject_q3 PRIVATE emdw)
target_compile_definitions(miniproject_q3 PRIVATE
  MINIPROJECT_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
