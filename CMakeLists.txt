cmake_minimum_required(VERSION 3.0)
project(quad_control)

## Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

# flags for speed (should already be enabled by default)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fopenmp -O3")

## Find catkin and any catkin packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  nav_msgs
  tf2_geometry_msgs
  tf2
  mavros
  mavros_msgs
  mavlink
  geometry_msgs
  genmsg
  gencpp
  message_generation
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

## Setup custom message and services
add_message_files(FILES FlatOutputs.msg)
add_service_files(FILES InitSetpoint.srv)

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

add_library(acado STATIC
    acado/generated_code/acado_qpoases_interface.cpp
    acado/generated_code/acado_integrator.c
    acado/generated_code/acado_solver.c
    acado/generated_code/acado_auxiliary_functions.c
    acado/qpoases/SRC/Bounds.cpp
    acado/qpoases/SRC/Constraints.cpp
    acado/qpoases/SRC/CyclingManager.cpp
    acado/qpoases/SRC/Indexlist.cpp
    acado/qpoases/SRC/MessageHandling.cpp
    acado/qpoases/SRC/QProblem.cpp
    acado/qpoases/SRC/QProblemB.cpp
    acado/qpoases/SRC/SubjectTo.cpp
    acado/qpoases/SRC/Utils.cpp
    acado/qpoases/SRC/EXTRAS/SolutionAnalysis.cpp
)

target_include_directories(acado PUBLIC
    acado/generated_code
    acado/qpoases
    acado/qpoases/INCLUDE
    acado/qpoases/SRC
)

target_link_libraries(acado ${catkin_LIBRARIES})

## Include directories and add executable
add_executable(mavros_cmd_node src/mavros_cmd_node.cpp)
add_dependencies(mavros_cmd_node ${catkin_EXPORTED_TARGETS})
target_link_libraries(mavros_cmd_node Eigen3::Eigen ${catkin_LIBRARIES})

add_executable(trajectory_gen_node src/trajectory_gen_node.cpp)
add_dependencies(trajectory_gen_node ${catkin_EXPORTED_TARGETS})
target_link_libraries(trajectory_gen_node Eigen3::Eigen ${catkin_LIBRARIES})

add_executable(mpc_node src/mpc_node.cpp)
add_dependencies(mpc_node ${catkin_EXPORTED_TARGETS})
target_link_libraries(mpc_node acado Eigen3::Eigen ${catkin_LIBRARIES})

include_directories(include ${catkin_INCLUDE_DIRS})
