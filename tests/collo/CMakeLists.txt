cmake_minimum_required(VERSION 3.21)

find_package(doctest REQUIRED)

add_executable(lobatto lobatto.cpp)

target_compile_features(lobatto PRIVATE cxx_std_20)
target_include_directories(lobatto PRIVATE ${COLLO_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(lobatto PRIVATE doctest::doctest)

add_executable(gauss gauss.cpp)

target_compile_features(gauss PRIVATE cxx_std_20)
target_include_directories(gauss PRIVATE ${COLLO_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(gauss PRIVATE doctest::doctest)

include(doctest)
doctest_discover_tests(lobatto)
doctest_discover_tests(gauss)
