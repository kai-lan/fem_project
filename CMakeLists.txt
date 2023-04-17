cmake_minimum_required (VERSION 3.0)
project (fem_project)

set(CMAKE_CXX_STANDARD 14)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package (Eigen3 REQUIRED)
if(NOT Eigen3_FOUND)
    error("Eigen not found!")
endif()

# SHARED or STATIC(default)
# If A.h and A.cc co-exist, only need to add A.cc
add_library(libraries Simulation.hh Mesh.cc Energy.cc FiniteDifference.hh)
target_link_libraries (libraries PUBLIC Eigen3::Eigen)

add_executable (main main.cc)
target_link_libraries(main libraries)

add_executable (unit_tests unit_tests.cc)
target_link_libraries(unit_tests libraries)

# Pybind 11
find_package(pybind11 REQUIRED)
pybind11_add_module(pyfem pyfem.cc)
target_link_libraries(pyfem PUBLIC libraries pybind11::pybind11)
set_target_properties(pyfem PROPERTIES LIBRARY_OUTPUT_DIRECTORY $<1:${PROJECT_SOURCE_DIR}/python>)