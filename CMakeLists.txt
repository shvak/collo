cmake_minimum_required(VERSION 3.24)

set(MAJOR_VERSION 1)
set(MINOR_VERSION 3)
set(PATCH_VERSION 0)
set(VERSION ${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION})
project(collo
  VERSION ${VERSION}
  DESCRIPTION "Header-only library collo: collocation integrator for ODEs"
  HOMEPAGE_URL "https://github.com/shvak/collo"
  LANGUAGES CXX
)

find_package(Eigen3 3.4 REQUIRED)

include(GNUInstallDirs)

set(COLLO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(collo INTERFACE)
include(cmake/ListHeaders.cmake) # set(COLLO_HEADERS ...)
target_sources(collo INTERFACE FILE_SET HEADERS BASE_DIRS ${COLLO_INCLUDE_DIR}
  FILES ${COLLO_HEADERS})
target_compile_features(collo INTERFACE cxx_std_20)
target_link_libraries(collo INTERFACE Eigen3::Eigen)
target_include_directories(collo INTERFACE
  $<BUILD_INTERFACE:${COLLO_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
set_target_properties(collo PROPERTIES VERSION ${VERSION})

include(CMakePackageConfigHelpers)

set(COLLO_INSTALL_CONFIG ${CMAKE_INSTALL_LIBDIR}/cmake/collo)

install(TARGETS collo
  EXPORT colloTargets
  FILE_SET HEADERS COMPONENT Development
)
install(EXPORT colloTargets
  FILE colloTargets.cmake
  NAMESPACE collo::
  DESTINATION ${COLLO_INSTALL_CONFIG}
)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/colloConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/colloConfig.cmake
  INSTALL_DESTINATION ${COLLO_INSTALL_CONFIG}
)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/colloConfigVersion.cmake
  VERSION ${VERSION}
  COMPATIBILITY SameMajorVersion
)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/colloConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/colloConfigVersion.cmake
  DESTINATION ${COLLO_INSTALL_CONFIG}
)

# if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
if(PROJECT_IS_TOP_LEVEL)
  add_subdirectory(examples)
  include(CTest)
  if(BUILD_TESTING)
    add_subdirectory(tests)
  endif()
endif()
