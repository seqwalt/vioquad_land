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
  mavros
  mavros_msgs
  mavlink
  geometry_msgs
  message_generation
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

## Setup custom message
add_message_files(FILES FlatOutputs.msg)

generate_messages(
    DEPENDENCIES
    std_msgs
    geometry_msgs
)

## Declare a catkin package
catkin_package(
    INCLUDE_DIRS include
    CATKIN_DEPENDS roscpp rospy std_msgs geometry_msgs message_runtime
)

## Include directories and add executable
add_executable(mavros_cmd_node src/mavros_cmd_node.cpp)

add_dependencies(mavros_cmd_node ${catkin_EXPORTED_TARGETS})
target_link_libraries(mavros_cmd_node Eigen3::Eigen ${catkin_LIBRARIES})
include_directories(include ${catkin_INCLUDE_DIRS})
