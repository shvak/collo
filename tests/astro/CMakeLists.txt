cmake_minimum_required(VERSION 3.21)

find_package(doctest REQUIRED)

add_executable(kepler-test kepler.cpp)

target_compile_features(kepler-test PRIVATE cxx_std_20)
target_include_directories(kepler-test PRIVATE ${COLLO_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(kepler-test PRIVATE doctest::doctest)

include(doctest)
doctest_discover_tests(kepler-test)
