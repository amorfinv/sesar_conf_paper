#%% Processes CONFLOG
import numpy as np
import os
import matplotlib.pyplot as plt

nm2m = 1852

results_folder = 'Attempt_7'

files = os.listdir(results_folder)
filtered_files = []

#Get needed files
for i, file in enumerate(files):
    if file[0:4] == 'CONF':
        filtered_files.append(file)

MVP_conf_data = []
ORCA_conf_data = []
VO_conf_data = []
NONE_conf_data = []
AIR_conf_data = []

for file in filtered_files:
    # Position 22 (and maybe 23) is always the aimed-for flight density
    if file[23] !='_':
        flight_den = int(file[22:24])
    else:
        flight_den = int(file[22])
    # Extract data
    data = np.genfromtxt(results_folder + '/' + file, delimiter = ',', dtype = int)
    if data.size == 0:
        #Skip empty files
        continue
    # We kind of only want the last row
    needed_data = np.array(data[-1])
    # Insert flight density as first number
    needed_data = np.append(flight_den, needed_data)
    # Append data to appropriate place
    if '_MVP_' in file:
        MVP_conf_data.append(needed_data)
    elif '_ORCA_' in file:
        ORCA_conf_data.append(needed_data)
    elif '_VO_' in file:
        VO_conf_data.append(needed_data)
    elif '_air_' in file:
        AIR_conf_data.append(needed_data)
    elif '_NONE_' in file and '_air_' not in file:
        NONE_conf_data.append(needed_data)
    else:
        print('Something is wrong.')

# Convert to numpy arrays
MVP_conf_data = np.array(MVP_conf_data)
ORCA_conf_data = np.array(ORCA_conf_data)
VO_conf_data = np.array(VO_conf_data)
NONE_conf_data = np.array(NONE_conf_data)
AIR_conf_data = np.array(AIR_conf_data)

# Divide the data by density
# Low
MVP_conf_data_l = MVP_conf_data[np.where(MVP_conf_data[:,0] == 5)]
ORCA_conf_data_l = ORCA_conf_data[np.where(ORCA_conf_data[:,0] == 5)]
VO_conf_data_l = VO_conf_data[np.where(VO_conf_data[:,0] == 5)]
NONE_conf_data_l = NONE_conf_data[np.where(NONE_conf_data[:,0] == 5)]
AIR_conf_data_l = AIR_conf_data[np.where(AIR_conf_data[:,0] == 5)]

# Medium
MVP_conf_data_m = MVP_conf_data[np.where(MVP_conf_data[:,0] == 8)]
ORCA_conf_data_m = ORCA_conf_data[np.where(ORCA_conf_data[:,0] == 8)]
VO_conf_data_m = VO_conf_data[np.where(VO_conf_data[:,0] == 8)]
NONE_conf_data_m = NONE_conf_data[np.where(NONE_conf_data[:,0] == 8)]
AIR_conf_data_m = AIR_conf_data[np.where(AIR_conf_data[:,0] == 8)]

# High
MVP_conf_data_h = MVP_conf_data[np.where(MVP_conf_data[:,0] == 11)]
ORCA_conf_data_h = ORCA_conf_data[np.where(ORCA_conf_data[:,0] == 11)]
VO_conf_data_h = VO_conf_data[np.where(VO_conf_data[:,0] == 11)]
NONE_conf_data_h = NONE_conf_data[np.where(NONE_conf_data[:,0] == 11)]
AIR_conf_data_h = AIR_conf_data[np.where(AIR_conf_data[:,0] == 11)]

'''
# We have data, now let's plot some stuff
plt.figure('conf_num')
plt.subplot(131)
plt.title('Low density')
plt.scatter([2]*len(MVP_conf_data_l), MVP_conf_data_l[:,2])
plt.scatter([4]*len(ORCA_conf_data_l), ORCA_conf_data_l[:,2])
plt.scatter([3]*len(VO_conf_data_l), VO_conf_data_l[:,2])
plt.scatter([0]*len(NONE_conf_data_l), NONE_conf_data_l[:,2])
plt.scatter([1]*len(AIR_conf_data_l), AIR_conf_data_l[:,2])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of conflicts [-]')
plt.show()

plt.subplot(132)
plt.title('Medium density')
plt.scatter([2]*len(MVP_conf_data_m), MVP_conf_data_m[:,2])
plt.scatter([4]*len(ORCA_conf_data_m), ORCA_conf_data_m[:,2])
plt.scatter([3]*len(VO_conf_data_m), VO_conf_data_m[:,2])
plt.scatter([0]*len(NONE_conf_data_m), NONE_conf_data_m[:,2])
plt.scatter([1]*len(AIR_conf_data_m), AIR_conf_data_m[:,2])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of conflicts [-]')
plt.show()

plt.subplot(133)
plt.title('High density')
plt.scatter([2]*len(MVP_conf_data_h), MVP_conf_data_h[:,2])
plt.scatter([4]*len(ORCA_conf_data_h), ORCA_conf_data_h[:,2])
plt.scatter([3]*len(VO_conf_data_h), VO_conf_data_h[:,2])
plt.scatter([0]*len(NONE_conf_data_h), NONE_conf_data_h[:,2])
plt.scatter([1]*len(AIR_conf_data_h), AIR_conf_data_h[:,2])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of conflicts [-]')
plt.show()

plt.figure('los_num')
plt.subplot(131)
plt.title('Low density')
plt.scatter([2]*len(MVP_conf_data_l), MVP_conf_data_l[:,3])
plt.scatter([4]*len(ORCA_conf_data_l), ORCA_conf_data_l[:,3])
plt.scatter([3]*len(VO_conf_data_l), VO_conf_data_l[:,3])
plt.scatter([0]*len(NONE_conf_data_l), NONE_conf_data_l[:,3])
plt.scatter([1]*len(AIR_conf_data_l), AIR_conf_data_l[:,3])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Losses of separation [-]')
plt.show()

plt.subplot(132)
plt.title('Medium density')
plt.scatter([2]*len(MVP_conf_data_m), MVP_conf_data_m[:,3])
plt.scatter([4]*len(ORCA_conf_data_m), ORCA_conf_data_m[:,3])
plt.scatter([3]*len(VO_conf_data_m), VO_conf_data_m[:,3])
plt.scatter([0]*len(NONE_conf_data_m), NONE_conf_data_m[:,3])
plt.scatter([1]*len(AIR_conf_data_m), AIR_conf_data_m[:,3])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Losses of separation [-]')
plt.show()

plt.subplot(133)
plt.title('High density')
plt.scatter([2]*len(MVP_conf_data_h), MVP_conf_data_h[:,3])
plt.scatter([4]*len(ORCA_conf_data_h), ORCA_conf_data_h[:,3])
plt.scatter([3]*len(VO_conf_data_h), VO_conf_data_h[:,3])
plt.scatter([0]*len(NONE_conf_data_h), NONE_conf_data_h[:,3])
plt.scatter([1]*len(AIR_conf_data_h), AIR_conf_data_h[:,3])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Losses of separation [-]')
plt.show()

plt.figure('geobreach_num')
plt.subplot(131)
plt.title('Low density')
plt.scatter([2]*len(MVP_conf_data_l), MVP_conf_data_l[:,4])
plt.scatter([4]*len(ORCA_conf_data_l), ORCA_conf_data_l[:,4])
plt.scatter([3]*len(VO_conf_data_l), VO_conf_data_l[:,4])
plt.scatter([0]*len(NONE_conf_data_l), NONE_conf_data_l[:,4])
plt.scatter([1]*len(AIR_conf_data_l), AIR_conf_data_l[:,4])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of geofence breaches [-]')
plt.show()

plt.subplot(132)
plt.title('Medium density')
plt.scatter([2]*len(MVP_conf_data_m), MVP_conf_data_m[:,4])
plt.scatter([4]*len(ORCA_conf_data_m), ORCA_conf_data_m[:,4])
plt.scatter([3]*len(VO_conf_data_m), VO_conf_data_m[:,4])
plt.scatter([0]*len(NONE_conf_data_m), NONE_conf_data_m[:,4])
plt.scatter([1]*len(AIR_conf_data_m), AIR_conf_data_m[:,4])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of geofence breaches [-]')
plt.show()

plt.subplot(133)
plt.title('High density')
plt.scatter([2]*len(MVP_conf_data_h), MVP_conf_data_h[:,4])
plt.scatter([4]*len(ORCA_conf_data_h), ORCA_conf_data_h[:,4])
plt.scatter([3]*len(VO_conf_data_h), VO_conf_data_h[:,4])
plt.scatter([0]*len(NONE_conf_data_h), NONE_conf_data_h[:,4])
plt.scatter([1]*len(AIR_conf_data_h), AIR_conf_data_h[:,4])
plt.xticks([0,1,2,3,4], ['None', 'Air', 'MVP', 'VO', 'ORCA'])
plt.ylabel('Number of geofence breaches [-]')
plt.show()
'''

#%% Process FLSTLOG
#Get needed files
filtered_files = []
for i, file in enumerate(files):
    if file[0:4] == 'FLST':
        filtered_files.append(file)

MVP_flst_data = []
ORCA_flst_data = []
VO_flst_data = []
NONE_flst_data = []
AIR_flst_data = []

for file in filtered_files:
    # Position 22 (and maybe 23) is always the aimed-for flight density
    if file[23] !='_':
        flight_den = int(file[22:24])
    else:
        flight_den = int(file[22])
    # Extract data
    data = np.genfromtxt(results_folder + '/' + file, delimiter = ',', dtype = 
                         ('float', 'S4', 'float', 'float', 'float', 'float', 
                          'float', 'float', 'float', 'float', 'float', 'float',
                          'float', 'bool', 'float', 'float', 'float', 'float', 
                          'int'))
    if data.size == 0:
        #Skip empty files
        continue
    # We kind of only want the last row
    needed_data = list(data)
    # Insert flight density as first number
    needed_data = list((flight_den, needed_data))
    # Append data to appropriate place
    if '_MVP_' in file:
        MVP_flst_data.append(needed_data)
    elif '_ORCA_' in file:
        ORCA_flst_data.append(needed_data)
    elif '_VO_' in file:
        VO_flst_data.append(needed_data)
    elif '_air_' in file:
        AIR_flst_data.append(needed_data)
    elif '_NONE_' in file and '_air_' not in file:
        NONE_flst_data.append(needed_data)
    else:
        print('Something is wrong.')
        
# Convert to numpy arrays
MVP_flst_data = np.array(MVP_flst_data)
ORCA_flst_data = np.array(ORCA_flst_data)
VO_flst_data = np.array(VO_flst_data)
NONE_flst_data = np.array(NONE_flst_data)
AIR_flst_data = np.array(AIR_flst_data)

# Divide the data by density
# Low
MVP_flst_data_l = MVP_flst_data[np.where(MVP_flst_data[:,0] == 5)]
ORCA_flst_data_l = ORCA_flst_data[np.where(ORCA_flst_data[:,0] == 5)]
VO_flst_data_l = VO_flst_data[np.where(VO_flst_data[:,0] == 5)]
NONE_flst_data_l = NONE_flst_data[np.where(NONE_flst_data[:,0] == 5)]
AIR_flst_data_l = AIR_flst_data[np.where(AIR_flst_data[:,0] == 5)]

# Medium
MVP_flst_data_m = MVP_flst_data[np.where(MVP_flst_data[:,0] == 8)]
ORCA_flst_data_m = ORCA_flst_data[np.where(ORCA_flst_data[:,0] == 8)]
VO_flst_data_m = VO_flst_data[np.where(VO_flst_data[:,0] == 8)]
NONE_flst_data_m = NONE_flst_data[np.where(NONE_flst_data[:,0] == 8)]
AIR_flst_data_m = AIR_flst_data[np.where(AIR_flst_data[:,0] == 8)]

# High
MVP_flst_data_h = MVP_flst_data[np.where(MVP_flst_data[:,0] == 11)]
ORCA_flst_data_h = ORCA_flst_data[np.where(ORCA_flst_data[:,0] == 11)]
VO_flst_data_h = VO_flst_data[np.where(VO_flst_data[:,0] == 11)]
NONE_flst_data_h = NONE_flst_data[np.where(NONE_flst_data[:,0] == 11)]
AIR_flst_data_h = AIR_flst_data[np.where(AIR_flst_data[:,0] == 11)]

# In this one we need more data processing
# First, compute the flight time for each
MVP_flst_flighttime_l = np.array([])
for scenario_data in MVP_flst_data_l:
    scenario_data = np.array([np.array(x) for x in scenario_data[1]])
    MVP_flst_flighttime_l = np.append(MVP_flst_flighttime_l, 
                            scenario_data[:,0] - scenario_data[:,2])

#%% Process REGLOG
#Get needed files
filtered_files = []
for i, file in enumerate(files):
    if file[0:3] == 'REG':
        filtered_files.append(file)

MVP_reg_data = []
ORCA_reg_data = []
VO_reg_data = []
NONE_reg_data = []
AIR_reg_data = []

for file in filtered_files:
    # Position 22 (and maybe 23) is always the aimed-for flight density
    if file[22] !='_':
        flight_den = int(file[21:23])
    else:
        flight_den = int(file[21])
    # Extract data
    data = np.genfromtxt(results_folder + '/' + file, delimiter = ',', dtype = int)
    if data.size == 0:
        #Skip empty files
        continue
    # We kind of only want the last row
    needed_data = np.array(data)
    # Insert flight density as first number
    needed_data = np.append(flight_den, needed_data)
    # Append data to appropriate place
    if '_MVP_' in file:
        MVP_reg_data.append(needed_data)
    elif '_ORCA_' in file:
        ORCA_reg_data.append(needed_data)
    elif '_VO_' in file:
        VO_reg_data.append(needed_data)
    elif '_air_' in file:
        AIR_reg_data.append(needed_data)
    elif '_NONE_' in file and '_air_' not in file:
        NONE_reg_data.append(needed_data)
    else:
        print('Something is wrong.')

