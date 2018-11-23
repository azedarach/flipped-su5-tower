cmake_minimum_required(VERSION 3.1)

project(flipped_su5_fits LANGUAGES CXX Fortran)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(FlexibleSUSYColors)
include(FlexibleSUSYBuildOptions)
include(FlexibleSUSYOptionalPackages)

if(ENABLE_TEST)
  find_package(Boost 1.37.0 COMPONENTS unit_test_framework REQUIRED)
else()
  find_package(Boost 1.37.0 REQUIRED)
endif()
find_package(Eigen3 3.1 REQUIRED)
find_package(GSL REQUIRED)

include(CheckIncludeFileCXX)
check_include_file_cxx(random ENABLE_RANDOM)

find_himalaya()
find_lapack()
find_looptools()
find_threads()
find_odeint()
find_debug()
find_sqlite()

if(ENABLE_LIBRARYLINK)
  find_package(Mathematica 7.0 COMPONENTS MathLink WolframLibrary REQUIRED)
endif(ENABLE_LIBRARYLINK)

if(ENABLE_COMPILER_WARNINGS)
  include(add_cxx_flag_if_supported)
  include(CompilerWarnings)
endif(ENABLE_COMPILER_WARNINGS)

add_subdirectory(config)
add_subdirectory(fflite)
add_subdirectory(slhaea)
add_subdirectory(src)
add_subdirectory(model_specific)
add_subdirectory(models)

if(ENABLE_TEST AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test)
  enable_testing()
  add_subdirectory(test)
endif()
