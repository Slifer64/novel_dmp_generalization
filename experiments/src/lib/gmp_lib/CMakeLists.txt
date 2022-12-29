cmake_minimum_required(VERSION 2.8.3)

project(gmp_lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Set the path containing additional modules like armadillo
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../install/deps/ ${CMAKE_MODULE_PATH} )

find_package(catkin REQUIRED COMPONENTS
  osqp_lib
)

catkin_package(
  INCLUDE_DIRS include ${catkin_INCLUDE_DIRS}
  DEPENDS
)

include_directories(
  include
  /usr/include/eigen3
  ${catkin_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}

  src/math/quaternions.cpp
  src/math/math.cpp
  src/io/gmp_io.cpp
  src/io/file_io.cpp
  src/io/xml_parser.cpp

  src/TrajScale/TrajScale.cpp

  src/CanonicalSystem/CanonicalSystem.cpp

  #------------------------
  src/GMP/GMP_regressor.cpp
  #------------------------
  src/GMP/GMP.cpp
  src/GMP/GMP_Update.cpp
  src/GMP/GMP_Opt.cpp
  src/GMP/GMP_MPC.cpp
  #------------------------
  src/GMP/GMPo.cpp
  src/GMP/GMPo_Update.cpp
)

add_dependencies(${PROJECT_NAME}
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
)

target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
  armadillo
  osqp_lib
)
