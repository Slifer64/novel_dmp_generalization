cmake_minimum_required(VERSION 3.0)

project(osqp_lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_MODULE_PATH}
)

find_package(catkin REQUIRED)

find_package(OSQPLibrary REQUIRED)

catkin_package(
  INCLUDE_DIRS
    include
    ${OSQP_INCLUDE_DIR}
    /usr/include/eigen3
)

include_directories(
  include
  /usr/include/eigen3
  ${catkin_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
  src/csc_mat.cpp
  src/quadprog.cpp
  src/qp_solver.cpp
)

target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
  ${OSQP_LIBRARIES}
)
