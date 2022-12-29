#!/bin/bash


CATKIN_MAKE_OPTIONS="-DCMAKE_BUILD_TYPE=Release"

if [ $# -eq 0 ]; then
  CATKIN_MAKE_OPTIONS=$@
fi

cd src/

file_name="CMakeLists.txt"
if [ ! -f $file_name ]; then
  rm -rf ../build/ ../devel/ >/dev/null 
  catkin_init_workspace
fi

cd ../
catkin_make -DCMAKE_BUILD_TYPE=Release
source devel/setup.bash




