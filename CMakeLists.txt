cmake_minimum_required(VERSION 3.0)
project(quad_control)

## Default to C++17
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

## Find catkin and any catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  genmsg
  nav_msgs
  tf2_geometry_msgs
  tf2
  mavros_msgs
  geometry_msgs)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

## Declare a catkin package
catkin_package()

## Include directories and add executable
add_executable(offb_node src/offb_node.cpp ${HEADER_FILES})

target_link_libraries(offb_node Eigen3::Eigen ${catkin_LIBRARIES})
include_directories(include ${catkin_INCLUDE_DIRS})
