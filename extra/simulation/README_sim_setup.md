# Setting up the simulator

## A) Generate custom RotorS IMU gazebo plugin
```
sudo apt-get install libgazebo11-dev  # if using noetic
sudo apt-get install libgazebo9-dev   # if using melodic
sudo apt install protobuf-compiler

cd .../quad_control/extra/simulation/rotors_imu_plugin
mkdir build
cd build
cmake ..
make
```

## B) Launch gazebo-classic with ROS wrappers:
1. Make sure gazebo-classic is installed
2. Add executable permissions (do once):
```
chmod +x run_gazebo-classic.sh
```
2. Run gazebo-classic:
```
./run_gazebo-classic.sh my/path/to/PX4-Autopilot/
```
3. This script follows from here: https://docs.px4.io/main/en/simulation/ros_interface.html#launching-gazebo-classic-with-ros-wrappers

## (Optional) regenerate file vi_sensor.sdf from vi_sensor.urdf.xacro

```
cd .../quad_control/extra/simulation/models/vi_camera/xacro
xacro vi_camera.urdf.xacro > vi_camera.urdf   # generate file vi_camera.urdf
gz sdf -p vi_camera.urdf > vi_camera.sdf      # generate file vi_camera.sdf
mv vi_camera.sdf ../                          # move file into vi_camera directory
```
