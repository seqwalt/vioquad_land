#!/bin/bash

#########################################
############### Functions ###############
#########################################

# Check if the input directory name matches the desired directory name
check_dir() {
  input_dir=$1    # first arg
  desired_dir=$2  # second arg
  if ! [ "$(basename "$input_dir")" == "$desired_dir" ]; then
      echo "$input_dir is not the $desired_dir directory"
      exit 1
  fi
}

# Define a function to check command success
check_command_success() {
  if ! [ $? -eq 0 ]; then
    exit 1  # Exit the script with an error code
  else
    echo "Success: $1" # print argument
    echo # new line
  fi
}

# Check if the directories in an environment variable path exist
check_path_exists() {
  path_name=$1 # argument should be path name, not path contents
  # Split the PATH into individual directories using ":" as the delimiter
  IFS=":" read -ra dirs <<< "${!path_name}"
  # Loop through each directory and check for existence
  for dir in "${dirs[@]}"; do
    # Check if the directory exists if its not an empty string
    if ! [ -d "$dir" ] && ! [ "" == "$dir" ] ; then
      echo "In ${path_name}, directory $dir does not exist"
      exit 1
    fi
  done
}

#########################################
################ Script #################
#########################################

# Check if an argument is provided to script
if [ $# -ne 1 ]; then
  echo "Usage: $0 <path/to/PX4-Autopilot>"
  exit 1
fi

# Store px4 path
px4path="$1"
check_dir $px4path "PX4-Autopilot"

# Store current path
sim_dir=$(pwd)
check_dir $sim_dir "simulation"

#### ---- Setup and start gazebo ---- ####
# See https://docs.px4.io/main/en/simulation/ros_interface.html#launching-gazebo-classic-with-ros-wrappers
cd $px4path
check_command_success "cd to PX4-Autopilot"

DONT_RUN=1 make px4_sitl_default gazebo-classic
check_command_success "DONT_RUN=1 make px4_sitl_default gazebo-classic"

source $sim_dir/../../../../devel/setup.bash
check_command_success "source workspace of quad_control"

source Tools/simulation/gazebo-classic/setup_gazebo.bash $(pwd) $(pwd)/build/px4_sitl_default
check_command_success "update Gazebo and LD paths"

export ROS_PACKAGE_PATH=$ROS_PACKAGE_PATH:$(pwd):$(pwd)/Tools/simulation/gazebo-classic/sitl_gazebo-classic
check_path_exists "ROS_PACKAGE_PATH"
check_command_success "update ROS_PACKAGE_PATH"

export GAZEBO_MODEL_PATH=$GAZEBO_MODEL_PATH:$sim_dir/models:$sim_dir/worlds
check_path_exists "GAZEBO_MODEL_PATH"
check_command_success "update GAZEBO_MODEL_PATH"

roslaunch px4 posix_sitl.launch \
  gui:=false \
  sdf:=$sim_dir/models/iris_with_cameras/iris_with_cameras.sdf \
  world:=$sim_dir/worlds/VIO_landing_pad.world
