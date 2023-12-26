cmake_minimum_required(VERSION 3.21)

set(major_version 1)
set(minor_version 2)
set(patch_version 0)
set(version ${major_version}.${minor_version}.${patch_version})

project(collo VERSION ${version} DESCRIPTION "Header-only library collo: collocation integrator for ODEs" HOMEPAGE_URL "https://github.com/shvak/collo" LANGUAGES CXX)

set(COLLO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(Eigen3 3.4 REQUIRED)

add_library(collo INTERFACE)

target_compile_features(collo INTERFACE cxx_std_20)
target_include_directories(collo INTERFACE
  $<BUILD_INTERFACE:${COLLO_INCLUDE_DIR}> ${EIGEN3_INCLUDE_DIRS}
)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  add_subdirectory(examples)
  include(CTest)
  if(BUILD_TESTING)
    add_subdirectory(tests)
  endif()
endif()
