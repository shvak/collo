cmake_minimum_required(VERSION 3.21)

find_package(fmt REQUIRED)

add_executable(kepler kepler.cpp)
target_compile_features(kepler PRIVATE cxx_std_20)
target_include_directories(kepler PRIVATE ${CMAKE_SOURCE_DIR}/include ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(kepler PRIVATE fmt::fmt)

add_custom_command(TARGET kepler POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/start.dat
    ${CMAKE_CURRENT_BINARY_DIR}/start.dat)
