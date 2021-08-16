# Processes conflog
import numpy as np
import os
import matplotlib.pyplot as plt

files = os.listdir('Attempt_3')
filtered_files = []

#Get needed files
for i, file in enumerate(files):
    if file[0:4] == 'CONF':
        filtered_files.append(file)

MVP_data = []
ORCA_data = []
VO_data = []
NONE_data = []
AIR_data = []

for file in filtered_files:
    # Position 22 (and maybe 23) is always the aimed-for flight density
    if file[23] !='_':
        flight_den = int(file[22:24])
    else:
        flight_den = int(file[22])
    # Extract data
    data = np.genfromtxt('Attempt_3/' + file, delimiter = ',', dtype = int)
    if data.size == 0:
        #Skip empty files
        continue
    # We kind of only want the last row
    needed_data = np.array(data[-1])
    # Insert flight density as first number
    needed_data = np.append(flight_den, needed_data)
    # Append data to appropriate place
    if '_MVP_' in file:
        MVP_data.append(needed_data)
    elif '_ORCA_' in file:
        ORCA_data.append(needed_data)
    elif '_VO_' in file:
        VO_data.append(needed_data)
    elif '_air_' in file:
        AIR_data.append(needed_data)
    elif '_NONE_' in file and '_air_' not in file:
        NONE_data.append(needed_data)
    else:
        print('Something is wrong.')

# Convert to numpy arrays
MVP_data = np.array(MVP_data)
ORCA_data = np.array(ORCA_data)
VO_data = np.array(VO_data)
NONE_data = np.array(NONE_data)
AIR_data = np.array(AIR_data)

# Divide the data by density
# Low
MVP_data_l = MVP_data[np.where(MVP_data[:,0] == 5)]
ORCA_data_l = ORCA_data[np.where(ORCA_data[:,0] == 5)]
VO_data_l = VO_data[np.where(VO_data[:,0] == 5)]
NONE_data_l = NONE_data[np.where(NONE_data[:,0] == 5)]
AIR_data_l = AIR_data[np.where(AIR_data[:,0] == 5)]

# Medium
MVP_data_m = MVP_data[np.where(MVP_data[:,0] == 8)]
ORCA_data_m = ORCA_data[np.where(ORCA_data[:,0] == 8)]
VO_data_m = VO_data[np.where(VO_data[:,0] == 8)]
NONE_data_m = NONE_data[np.where(NONE_data[:,0] == 8)]
AIR_data_m = AIR_data[np.where(AIR_data[:,0] == 8)]

# High
MVP_data_h = MVP_data[np.where(MVP_data[:,0] == 11)]
ORCA_data_h = ORCA_data[np.where(ORCA_data[:,0] == 11)]
VO_data_h = VO_data[np.where(VO_data[:,0] == 11)]
NONE_data_h = NONE_data[np.where(NONE_data[:,0] == 11)]
AIR_data_h = AIR_data[np.where(AIR_data[:,0] == 11)]

# We have data, now let's plot some stuff
plt.figure('conf_num')
plt.subplot(131)
plt.title('Low density')
plt.scatter([2]*len(MVP_data_l), MVP_data_l[:,2])
plt.scatter([4]*len(ORCA_data_l), ORCA_data_l[:,2])
plt.scatter([3]*len(VO_data_l), VO_data_l[:,2])
plt.scatter([0]*len(NONE_data_l), NONE_data_l[:,2])
plt.scatter([1]*len(AIR_data_l), AIR_data_l[:,2])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of conflicts [-]')
plt.show()

plt.subplot(132)
plt.title('Medium density')
plt.scatter([2]*len(MVP_data_m), MVP_data_m[:,2])
plt.scatter([4]*len(ORCA_data_m), ORCA_data_m[:,2])
plt.scatter([3]*len(VO_data_m), VO_data_m[:,2])
plt.scatter([0]*len(NONE_data_m), NONE_data_m[:,2])
plt.scatter([1]*len(AIR_data_m), AIR_data_m[:,2])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of conflicts [-]')
plt.show()

plt.subplot(133)
plt.title('High density')
plt.scatter([2]*len(MVP_data_h), MVP_data_h[:,2])
plt.scatter([4]*len(ORCA_data_h), ORCA_data_h[:,2])
plt.scatter([3]*len(VO_data_h), VO_data_h[:,2])
plt.scatter([0]*len(NONE_data_h), NONE_data_h[:,2])
plt.scatter([1]*len(AIR_data_h), AIR_data_h[:,2])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of conflicts [-]')
plt.show()

plt.figure('los_num')
plt.subplot(131)
plt.title('Low density')
plt.scatter([2]*len(MVP_data_l), MVP_data_l[:,3])
plt.scatter([4]*len(ORCA_data_l), ORCA_data_l[:,3])
plt.scatter([3]*len(VO_data_l), VO_data_l[:,3])
plt.scatter([0]*len(NONE_data_l), NONE_data_l[:,3])
plt.scatter([1]*len(AIR_data_l), AIR_data_l[:,3])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Losses of separation [-]')
plt.show()

plt.subplot(132)
plt.title('Medium density')
plt.scatter([2]*len(MVP_data_m), MVP_data_m[:,3])
plt.scatter([4]*len(ORCA_data_m), ORCA_data_m[:,3])
plt.scatter([3]*len(VO_data_m), VO_data_m[:,3])
plt.scatter([0]*len(NONE_data_m), NONE_data_m[:,3])
plt.scatter([1]*len(AIR_data_m), AIR_data_m[:,3])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Losses of separation [-]')
plt.show()

plt.subplot(133)
plt.title('High density')
plt.scatter([2]*len(MVP_data_h), MVP_data_h[:,3])
plt.scatter([4]*len(ORCA_data_h), ORCA_data_h[:,3])
plt.scatter([3]*len(VO_data_h), VO_data_h[:,3])
plt.scatter([0]*len(NONE_data_h), NONE_data_h[:,3])
plt.scatter([1]*len(AIR_data_h), AIR_data_h[:,3])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Losses of separation [-]')
plt.show()

plt.figure('geobreach_num')
plt.subplot(131)
plt.title('Low density')
plt.scatter([2]*len(MVP_data_l), MVP_data_l[:,4])
plt.scatter([4]*len(ORCA_data_l), ORCA_data_l[:,4])
plt.scatter([3]*len(VO_data_l), VO_data_l[:,4])
plt.scatter([0]*len(NONE_data_l), NONE_data_l[:,4])
plt.scatter([1]*len(AIR_data_l), AIR_data_l[:,4])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of geofence breaches [-]')
plt.show()

plt.subplot(132)
plt.title('Medium density')
plt.scatter([2]*len(MVP_data_m), MVP_data_m[:,4])
plt.scatter([4]*len(ORCA_data_m), ORCA_data_m[:,4])
plt.scatter([3]*len(VO_data_m), VO_data_m[:,4])
plt.scatter([0]*len(NONE_data_m), NONE_data_m[:,4])
plt.scatter([1]*len(AIR_data_m), AIR_data_m[:,4])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of geofence breaches [-]')
plt.show()

plt.subplot(133)
plt.title('High density')
plt.scatter([2]*len(MVP_data_h), MVP_data_h[:,4])
plt.scatter([4]*len(ORCA_data_h), ORCA_data_h[:,4])
plt.scatter([3]*len(VO_data_h), VO_data_h[:,4])
plt.scatter([0]*len(NONE_data_h), NONE_data_h[:,4])
plt.scatter([1]*len(AIR_data_h), AIR_data_h[:,4])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of geofence breaches [-]')
plt.show()

