cmake_minimum_required(VERSION 3.23)
project(Pathtracing)

file(GLOB sources src/*.cpp)

add_executable(main ${sources})
target_compile_options(main PUBLIC -std=c++17 -Wall -Wconversion -O3)
include_directories(include)
