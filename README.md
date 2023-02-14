## Simple mavros node for a PX4 quadcopter
### Usage
1. Start mavros in either simulation or hardware
2. Run the controller
```
rosrun quad_control mavros_cmd_node
```
3. Run the trajectory publisher
```
rosrun quad_control trajectory_gen_node
```
### Behaviour
- Currently, this package allows a quadcopter to follow a trajectory given by a ```.csv``` fille
    - The ```.csv``` file should have 18 columns:
    ```
    time|x|y|z|yaw|dx|dy|dz|dyaw|ddx|ddy|ddz|dddx|dddy|dddz|ddddx|ddddy|ddddz
    ```
    - Currently only position, yaw, velocity and acceleration are used for the geometric control.
    - The trajectory ```.csv``` file is assumed to 0.01 sec time steps (i.e. 100 Hz).
- The geometric control input is computed and published to mavros everytime a new trajectory set-point is published (at 100 Hz).
- Upon running the ```rosrun``` commands the following should occur
    1. ```mav_cmd_node``` waits to receive the home position from mavros
    2. ```mav_cmd_node``` waits to receive the first trajectory stepoint from ```trajectory_gen_node```
    3. ```mav_cmd_node``` takes-off and hovers at the first trajectory setpoint (position and heading)
    4. ```trajectory_gen_node``` starts publishing the trajectory at 100 Hz
    5. ```mav_cmd_node``` computes a control input and publishes a mavros setpoint command for:
        - attitude
        - angular velocity
        - thrust
    6. The quadcopter tracks the trajectory, the lands once it is done.
