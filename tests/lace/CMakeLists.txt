cmake_minimum_required(VERSION 3.21)

find_package(doctest REQUIRED)

add_executable(vector vector.cpp)

target_compile_features(vector PRIVATE cxx_std_20)
target_include_directories(vector PRIVATE ${COLLO_INCLUDE_DIR})
target_link_libraries(vector PRIVATE doctest::doctest)

add_executable(matrix matrix.cpp)

target_compile_features(matrix PRIVATE cxx_std_20)
target_include_directories(matrix PRIVATE ${COLLO_INCLUDE_DIR})
target_link_libraries(matrix PRIVATE doctest::doctest)

include(doctest)
doctest_discover_tests(vector)
doctest_discover_tests(matrix)
