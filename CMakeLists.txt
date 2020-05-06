cmake_minimum_required(VERSION 3.10)
project(RayTracing)

include_directories(lib/stb)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -fopenmp")

add_executable(rt main.cpp vector.h image.h objects.h mesh.h tracer.h)