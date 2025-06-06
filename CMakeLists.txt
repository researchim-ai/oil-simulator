cmake_minimum_required(VERSION 3.10)
project(OilSimulator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(simulator src/main.cpp src/Grid.cpp src/Solver.cpp src/Parameters.cpp src/FluidProperties.cpp) 