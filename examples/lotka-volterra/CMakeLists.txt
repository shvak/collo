cmake_minimum_required(VERSION 3.21)

add_executable(lotka-volterra lotka-volterra.cpp)
target_compile_features(lotka-volterra PRIVATE cxx_std_20)
target_include_directories(lotka-volterra PRIVATE ${COLLO_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIRS})
