# Message for flat outputs (and some derivatives) of a quadcopter
# Flat outputs: (position.x, position.y, position.z, yaw)
#
# trajectoryDone is false if there is a trajectory occuring,
# and true if the trajectory is over.

std_msgs/Header header
bool trajectoryDone

geometry_msgs/Point position
geometry_msgs/Vector3 velocity
geometry_msgs/Vector3 acceleration
float64 yaw
