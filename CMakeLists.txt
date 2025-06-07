cmake_minimum_required(VERSION 3.11)
project(OilSimulator)

# --- Fetch and configure Eigen ---
include(FetchContent)
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)
include_directories(${eigen_SOURCE_DIR})

# --- OpenMP Support ---
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create a library with all simulator sources except main
add_library(oillib
    src/Grid.cpp
    src/Solver.cpp
    src/Parameters.cpp
    src/FluidProperties.cpp
    src/Well.cpp
    src/Scenario.cpp
)
# Make headers from src available
target_include_directories(oillib PUBLIC src)
# Link the library to OpenMP and Eigen
target_link_libraries(oillib PUBLIC OpenMP::OpenMP_CXX eigen)

# Create the main executable and link it to the library
add_executable(simulator src/main.cpp)
target_link_libraries(simulator PRIVATE oillib)

# Target for running tests
add_executable(test_runner tests/TestSaturationUpdate.cpp)
target_link_libraries(test_runner PRIVATE oillib) 