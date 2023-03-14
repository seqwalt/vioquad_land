import csv
import pyulog
import numpy as np
from matplotlib import pyplot as plt

def load_csv_data(csv_file, operation):
    with open(csv_file) as csvfile:
        # Create a reader object
        reader = csv.reader(csvfile)
        # Skip the header row
        next(reader)
        # Initialize lists to hold the data
        Time = []
        Data = []
        # Read each row of data
        for row in reader:
            # Extract the data
            time = float(row[0])*1e-6 # convert microseconds to seconds
            data = operation(row);
            # Add the values to the lists
            Time.append(time)
            Data.append(data)
    return [Time, Data]

def get_avg_motor_PWM(row):
    # input: row from actuator_motors csv file
    # output: average motor pwm
    data = []
    for i in range(4):
        data.append(float(row[i+2]))
    return (data[0] + data[1] + data[2] + data[3])/4

def get_normalized_thrust(row):
    # input: row from actuator_controls csv file
    # output: normalized thrust (controls[3], which is column 6)
    return float(row[5])

# Open the CSV file
# actuator_controls_file = "/home/sequoyah/Desktop/test_ulog_csv/03_44_56_actuator_controls_0_0.csv";
# data = load_csv_data(actuator_controls_file, get_normalized_thrust)
# times = data[0]
# thrust = data[1]

path = '/home/sequoyah/Software/PX4-Autopilot/build/px4_sitl_default/rootfs/log/'
log = pyulog.ULog(path+'2023-03-13/03_44_56.ulg')

# show available datasets
LOG = log.data_list
for dataset in LOG:
    print(dataset.name)

# PX4 Normalized thrust [0,1]
actuator_ctrl_dataset = log.get_dataset('actuator_controls_0')
#print(dir(actuator_ctrl_dataset)) # show attributes of object
ctrl_dict = actuator_ctrl_dataset.data # dictionary
times1 = ctrl_dict['timestamp']*1e-6 # convert to seconds
thrust = ctrl_dict['control[3]']

# Mass-normalized thrust (units of m/s^2)
accel_dataset = log.get_dataset('vehicle_acceleration')
accel_dict = accel_dataset.data # dictionary
times2 = accel_dict['timestamp']*1e-6 # convert to seconds
x_acc = accel_dict['xyz[0]'] # units of m/s^2
y_acc = accel_dict['xyz[1]']
z_acc = accel_dict['xyz[2]']

att_dataset = log.get_dataset('vehicle_attitude_groundtruth')
att_dict = att_dataset.data
times3 = att_dict['timestamp']*1e-6 # convert to seconds
print(att_dict)
q1 = att_dict['q[0]'] # units of m/s^2
q2 = att_dict['q[1]']
q3 = att_dict['q[2]']
q4 = att_dict['q[3]']

# Create a plot
'''
plt.figure(1)
plt.plot(times1, thrust)
# Add labels and title
plt.xlabel('Time (sec)')
plt.ylabel('Thrust (PX4 normalized)')
plt.title('PX4 Normalized Thrust')

plt.figure(2)
plt.plot(times2, x_acc)
plt.plot(times2, y_acc)
plt.plot(times2, z_acc)
plt.xlabel('Time (sec)')
plt.ylabel('Acceleration (m/s^2)')
plt.title('Acceleration')

# Display the plots
plt.show()
'''
