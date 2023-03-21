# Message for flat outputs (and some derivatives) of a quadcopter
# Flat outputs: (position.x, position.y, position.z, yaw)
#
# reached_goal is false if the quad has not yet reached the goal

std_msgs/Header header
bool reached_goal

geometry_msgs/Point position
geometry_msgs/Vector3 velocity
geometry_msgs/Vector3 acceleration
float64 yaw
