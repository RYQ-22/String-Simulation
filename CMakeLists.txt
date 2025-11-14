cmake_minimum_required(VERSION 3.10)
project(discrete_elastic_rods)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
endif()

# Find Eigen
find_package(Eigen3 REQUIRED)

# Find Polyscope
find_package(Polyscope REQUIRED)
  
# Find OpenMP
find_package(OpenMP REQUIRED)

file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp)

add_executable(${PROJECT_NAME} ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR})

# Use C++14
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(${PROJECT_NAME} polyscope)
target_compile_options(${PROJECT_NAME} PUBLIC ${OpenMP_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME} ${OpenMP_CXX_FLAGS})
