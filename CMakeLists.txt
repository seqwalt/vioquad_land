cmake_minimum_required(VERSION 2.8.3)
project(quad_control)

## Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
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

## Declare a catkin package
catkin_package()

## Include directories and add executable
add_executable(offb_node src/offb_node.cpp ${HEADER_FILES})
target_link_libraries(offb_node ${catkin_LIBRARIES})

target_include_directories(offb_node PRIVATE)
include_directories(include ${catkin_INCLUDE_DIRS})
