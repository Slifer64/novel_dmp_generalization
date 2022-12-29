cmake_minimum_required(VERSION 2.8.3)

project(io_lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
# Set the path containing additional modules like armadillo
#set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../install/deps/ ${CMAKE_MODULE_PATH} )

find_package(catkin REQUIRED)

catkin_package(
        INCLUDE_DIRS include
)

include_directories(
  include
  /usr/include/eigen3
  ${catkin_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
  src/io_utils.cpp
  src/xml_parser.cpp
  src/file_io.cpp
  src/print_utils.cpp
)

add_dependencies(${PROJECT_NAME}
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
)

target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
  armadillo
)
