cmake_minimum_required(VERSION 3.0)
project(quad_control)

set(CMAKE_BUILD_TYPE Debug)

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
  apriltag_ros
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

## acado source files and embedded version of qpOASES
# set(codegen_dir "acado/examples/quad_mpc") # for just mpc, not PAMPC
set(codegen_dir "acado/PAMPC")

add_library(acado STATIC
    ${codegen_dir}/autogen_code/acado_qpoases_interface.cpp
    ${codegen_dir}/autogen_code/acado_integrator.c
    ${codegen_dir}/autogen_code/acado_solver.c
    ${codegen_dir}/autogen_code/acado_auxiliary_functions.c
    ${codegen_dir}/qpoases/SRC/Bounds.cpp
    ${codegen_dir}/qpoases/SRC/Constraints.cpp
    ${codegen_dir}/qpoases/SRC/CyclingManager.cpp
    ${codegen_dir}/qpoases/SRC/Indexlist.cpp
    ${codegen_dir}/qpoases/SRC/MessageHandling.cpp
    ${codegen_dir}/qpoases/SRC/QProblem.cpp
    ${codegen_dir}/qpoases/SRC/QProblemB.cpp
    ${codegen_dir}/qpoases/SRC/SubjectTo.cpp
    ${codegen_dir}/qpoases/SRC/Utils.cpp
    ${codegen_dir}/qpoases/SRC/EXTRAS/SolutionAnalysis.cpp
)
target_include_directories(acado PUBLIC
    acado/
    ${codegen_dir}/autogen_code
    ${codegen_dir}/qpoases
    ${codegen_dir}/qpoases/INCLUDE
    ${codegen_dir}/qpoases/SRC
)
target_link_libraries(acado ${catkin_LIBRARIES})

## Full and updated qpOASES solver for min-snap traj generation
set(full_qpoases_dir "extra/qpOASES")
add_library(full_qpoases STATIC
    ${full_qpoases_dir}/src/BLASReplacement.cpp
    ${full_qpoases_dir}/src/Bounds.cpp
    ${full_qpoases_dir}/src/Constraints.cpp
    ${full_qpoases_dir}/src/Flipper.cpp
    ${full_qpoases_dir}/src/Indexlist.cpp
    ${full_qpoases_dir}/src/LAPACKReplacement.cpp
    ${full_qpoases_dir}/src/Makefile
    ${full_qpoases_dir}/src/Matrices.cpp
    ${full_qpoases_dir}/src/MessageHandling.cpp
    ${full_qpoases_dir}/src/Options.cpp
    ${full_qpoases_dir}/src/OQPinterface.cpp
    ${full_qpoases_dir}/src/QProblemB.cpp
    ${full_qpoases_dir}/src/QProblem.cpp
    ${full_qpoases_dir}/src/SolutionAnalysis.cpp
    ${full_qpoases_dir}/src/SparseSolver.cpp
    ${full_qpoases_dir}/src/SQProblem.cpp
    ${full_qpoases_dir}/src/SQProblemSchur.cpp
    ${full_qpoases_dir}/src/SubjectTo.cpp
    ${full_qpoases_dir}/src/Utils.cpp
)
target_include_directories(full_qpoases PUBLIC
    ${full_qpoases_dir}
    ${full_qpoases_dir}/include
)

## Include directories and add executable
add_executable(mavros_cmd_node src/mavros_cmd_node.cpp)
add_dependencies(mavros_cmd_node ${catkin_EXPORTED_TARGETS})
target_link_libraries(mavros_cmd_node Eigen3::Eigen ${catkin_LIBRARIES})

add_executable(trajectory_gen_node src/trajectory_gen_node.cpp)
add_dependencies(trajectory_gen_node ${catkin_EXPORTED_TARGETS})
target_link_libraries(trajectory_gen_node Eigen3::Eigen ${catkin_LIBRARIES})

add_executable(mpc_node src/mpc_node.cpp)
add_dependencies(mpc_node ${catkin_EXPORTED_TARGETS})
target_link_libraries(mpc_node acado full_qpoases Eigen3::Eigen ${catkin_LIBRARIES})

set(DO_TESTS TRUE)
if(DO_TESTS)
    add_executable(min_snap_test extra/tests/min_snap_traj_test.cpp)
    target_link_libraries(min_snap_test full_qpoases Eigen3::Eigen)
    
    add_executable(min_snap_test2 extra/tests/min_snap_traj_test2.cpp)
    target_link_libraries(min_snap_test2 full_qpoases Eigen3::Eigen)
endif()

include_directories(include ${catkin_INCLUDE_DIRS})
