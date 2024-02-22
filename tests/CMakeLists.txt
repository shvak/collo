cmake_minimum_required(VERSION 3.21)

find_package(doctest REQUIRED)
include(doctest)

function(testexe testname)
  add_executable(${testname} ${testname}.cpp)
  target_compile_features(${testname} PRIVATE cxx_std_20)
  target_link_libraries(${testname} PRIVATE doctest::doctest collo)
  doctest_discover_tests(${testname})
endfunction()

add_subdirectory(numm)
add_subdirectory(lace)
add_subdirectory(collo)
add_subdirectory(astro)
