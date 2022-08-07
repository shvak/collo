cmake_minimum_required(VERSION 3.21)

find_package(fmt REQUIRED)

add_executable(isochrone isochrone.cpp)
target_compile_features(isochrone PRIVATE cxx_std_20)
target_include_directories(isochrone PRIVATE ${CMAKE_SOURCE_DIR}/include ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(isochrone PRIVATE fmt::fmt)

add_custom_command(TARGET isochrone POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/start.dat
    ${CMAKE_CURRENT_BINARY_DIR}/start.dat)
