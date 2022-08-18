cmake_minimum_required(VERSION 3.21)

find_package(doctest REQUIRED)

add_executable(base base.cpp)

target_compile_features(base PRIVATE cxx_std_20)
target_include_directories(base PRIVATE ${COLLO_INCLUDE_DIR})
target_link_libraries(base PRIVATE doctest::doctest)

add_executable(roots roots.cpp)

target_compile_features(roots PRIVATE cxx_std_20)
target_include_directories(roots PRIVATE ${COLLO_INCLUDE_DIR})
target_link_libraries(roots PRIVATE doctest::doctest)

add_executable(legendre legendre.cpp)

target_compile_features(legendre PRIVATE cxx_std_20)
target_include_directories(legendre PRIVATE ${COLLO_INCLUDE_DIR})
target_link_libraries(legendre PRIVATE doctest::doctest)

include(doctest)
doctest_discover_tests(base)
doctest_discover_tests(roots)
doctest_discover_tests(legendre)
