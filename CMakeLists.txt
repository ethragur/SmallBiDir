cmake_minimum_required(VERSION 2.8)
project(PathTracing)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -march=native -fopenmp")

set(SOURCE_FILES
	src/Main.cpp
	src/PathTracer.cpp
	src/BiDirectionalPT.cpp
    src/scene.cpp
    )

add_executable(PathTracing ${SOURCE_FILES})
