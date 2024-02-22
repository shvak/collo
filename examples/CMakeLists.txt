cmake_minimum_required(VERSION 3.21)

function(example_exe example_name)
  add_executable(${example_name} ${example_name}.cpp)
  target_compile_features(${example_name} PRIVATE cxx_std_20)
  target_link_libraries(${example_name} PRIVATE collo)
endfunction()

add_subdirectory(kepler)
add_subdirectory(isochrone)
add_subdirectory(lotka-volterra)
