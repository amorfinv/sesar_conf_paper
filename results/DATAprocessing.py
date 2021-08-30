# %% Processes CONFLOG
import numpy as np
import os
import matplotlib.pyplot as plt
from cycler import cycler

nm2m = 1852

results_folder = 'Attempt_7'

files = os.listdir(results_folder)
filtered_files = []

# Set default font sizes

plt.rc('font', size=8)          # controls default text sizes
plt.rc('axes', titlesize=9)     # fontsize of the axes title
plt.rc('axes', labelsize=9)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
plt.rc('ytick', labelsize=8)    # fontsize of the tick labels
plt.rc('legend', fontsize=8)    # legend fontsize
plt.rc('figure', titlesize=9)  # fontsize of the figure title
plt.rcParams["scatter.marker"] = 'x'
plt.rcParams['lines.markersize'] = 4

# %%% Get needed files
for i, file in enumerate(files):
    if file[0:4] == 'CONF':
        filtered_files.append(file)

NONE_conf_data = []
NONE_AIR_conf_data = []
MVP_conf_data = []
MVP_AIR_conf_data = []
ORCA_conf_data = []
ORCA_AIR_conf_data = []
MVPC_conf_data = []
MVPC_AIR_conf_data = []
ORCAC_conf_data = []
ORCAC_AIR_conf_data = []

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
        if '_air_' in file:
            MVP_AIR_conf_data.append(needed_data)
        else:
            MVP_conf_data.append(needed_data)
    elif '_MVPC_' in file:
        if '_air_' in file:
            MVPC_AIR_conf_data.append(needed_data)
        else:
            MVPC_conf_data.append(needed_data)
    elif '_ORCA_' in file:
        if '_air_' in file:
            ORCA_AIR_conf_data.append(needed_data)
        else:
            ORCA_conf_data.append(needed_data)
    elif '_ORCAC_' in file:
        if '_air_' in file:
            ORCAC_AIR_conf_data.append(needed_data)
        else:
            ORCAC_conf_data.append(needed_data)
    elif '_NONE_' in file:
        if '_air_' in file:
            NONE_AIR_conf_data.append(needed_data)
        else:
            NONE_conf_data.append(needed_data)
    else:
        print('Something is wrong.')

# Convert to numpy arrays
NONE_conf_data = np.array(NONE_conf_data)
NONE_AIR_conf_data = np.array(NONE_AIR_conf_data)
MVP_conf_data = np.array(MVP_conf_data)
MVP_AIR_conf_data = np.array(MVP_AIR_conf_data)
ORCA_conf_data = np.array(ORCA_conf_data)
ORCA_AIR_conf_data = np.array(ORCA_AIR_conf_data)
MVPC_conf_data = np.array(MVPC_conf_data)
MVPC_AIR_conf_data = np.array(MVPC_AIR_conf_data)
ORCAC_conf_data = np.array(ORCAC_conf_data)
ORCAC_AIR_conf_data = np.array(ORCAC_AIR_conf_data)


globaldata_init = [NONE_conf_data, NONE_AIR_conf_data, 
              MVP_conf_data, MVP_AIR_conf_data,
              ORCA_conf_data, ORCA_AIR_conf_data,
              MVPC_conf_data, MVPC_AIR_conf_data, 
              ORCAC_conf_data, ORCAC_AIR_conf_data]

# %%% Divide the data by density
# Low
globaldata_l = []
for data in globaldata_init:
    low_data = data[np.where(data[:,0] == 5)]
    globaldata_l.append(low_data)

# Medium
globaldata_m = []
for data in globaldata_init:
    low_data = data[np.where(data[:,0] == 5)]
    globaldata_l.append(low_data)

# High
globaldata_h = []
for data in globaldata_init:
    low_data = data[np.where(data[:,0] == 5)]
    globaldata_l.append(low_data)
    
globaldata = [globaldata_l, globaldata_m, globaldata_h]

# %%% Graphs

# %%%% Number of conflicts
plt.figure('conf_num', figsize = (5,3))
ax1 = plt.subplot(131)
ax1.set_title('Low density')
ax1.scatter([2]*len(MVP_conf_data_l), MVP_conf_data_l[:,2], color = '#808080')
ax1.scatter([3]*len(ORCA_conf_data_l), ORCA_conf_data_l[:,2], color = '#808080')
ax1.scatter([0]*len(NONE_conf_data_l), NONE_conf_data_l[:,2], color = '#808080')
ax1.scatter([1]*len(AIR_conf_data_l), AIR_conf_data_l[:,2], color = '#808080')
ax1.boxplot([MVP_conf_data_l[:,2], ORCA_conf_data_l[:,2], 
             NONE_conf_data_l[:,2], AIR_conf_data_l[:,2]] , positions = [2,3,0,1])
ax1.set_xticks([0,1,2,3])
ax1.set_xticklabels(['N', 'A', 'M', 'O'])
ax1.set_ylabel('Number of conflicts [-]')

ax2 = plt.subplot(132)
ax2.set_title('Medium density')
ax2.scatter([2]*len(MVP_conf_data_m), MVP_conf_data_m[:,2], color = '#808080')
ax2.scatter([3]*len(ORCA_conf_data_m), ORCA_conf_data_m[:,2], color = '#808080')
ax2.scatter([0]*len(NONE_conf_data_m), NONE_conf_data_m[:,2], color = '#808080')
ax2.scatter([1]*len(AIR_conf_data_m), AIR_conf_data_m[:,2], color = '#808080')
ax2.boxplot([MVP_conf_data_m[:,2], ORCA_conf_data_m[:,2], 
             NONE_conf_data_m[:,2], AIR_conf_data_m[:,2]] , positions = [2,3,0,1])
ax2.set_xticks([0,1,2,3])
ax2.set_yticks([])
ax2.set_xticklabels(['N', 'A', 'M', 'O'])

ax3 = plt.subplot(133)
ax3.set_title('High density')
ax3.scatter([2]*len(MVP_conf_data_h), MVP_conf_data_h[:,2], color = '#808080')
ax3.scatter([3]*len(ORCA_conf_data_h), ORCA_conf_data_h[:,2], color = '#808080')
ax3.scatter([0]*len(NONE_conf_data_h), NONE_conf_data_h[:,2], color = '#808080')
ax3.scatter([1]*len(AIR_conf_data_h), AIR_conf_data_h[:,2], color = '#808080')
ax3.boxplot([MVP_conf_data_h[:,2], ORCA_conf_data_h[:,2], 
             NONE_conf_data_h[:,2], AIR_conf_data_h[:,2]] , positions = [2,3,0,1])
ax3.set_xticks([0,1,2,3])
ax3.set_yticks([])
ax3.set_xticklabels(['N', 'A', 'M', 'O'])

ax1.set_ylim([-2, max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
ax2.set_ylim([-2, max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
ax3.set_ylim([-2, max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
plt.savefig('Plots/conf_num.png')
plt.show()

# %%%% Number of LOS
plt.figure('los_num', figsize = (5,3))
ax1 = plt.subplot(131)
ax1.set_title('Low density')
ax1.scatter([2]*len(MVP_conf_data_l), MVP_conf_data_l[:,3], color = '#808080')
ax1.scatter([3]*len(ORCA_conf_data_l), ORCA_conf_data_l[:,3], color = '#808080')
ax1.scatter([0]*len(NONE_conf_data_l), NONE_conf_data_l[:,3], color = '#808080')
ax1.scatter([1]*len(AIR_conf_data_l), AIR_conf_data_l[:,3], color = '#808080')
ax1.boxplot([MVP_conf_data_l[:,3], ORCA_conf_data_l[:,3], 
             NONE_conf_data_l[:,3], AIR_conf_data_l[:,3]] , positions = [2,3,0,1])
ax1.set_xticks([0,1,2,3])
ax1.set_xticklabels(['N', 'A', 'M', 'O'])
ax1.set_ylabel('Number of losses of separation [-]')

ax2 = plt.subplot(132)
ax2.set_title('Medium density')
ax2.scatter([2]*len(MVP_conf_data_m), MVP_conf_data_m[:,3], color = '#808080')
ax2.scatter([3]*len(ORCA_conf_data_m), ORCA_conf_data_m[:,3], color = '#808080')
ax2.scatter([0]*len(NONE_conf_data_m), NONE_conf_data_m[:,3], color = '#808080')
ax2.scatter([1]*len(AIR_conf_data_m), AIR_conf_data_m[:,3], color = '#808080')
ax2.boxplot([MVP_conf_data_m[:,3], ORCA_conf_data_m[:,3], 
             NONE_conf_data_m[:,3], AIR_conf_data_m[:,3]] , positions = [2,3,0,1])
ax2.set_xticks([0,1,2,3])
ax2.set_yticks([])
ax2.set_xticklabels(['N', 'A', 'M', 'O'])

ax3 = plt.subplot(133)
ax3.set_title('High density')
ax3.scatter([2]*len(MVP_conf_data_h), MVP_conf_data_h[:,3], color = '#808080')
ax3.scatter([3]*len(ORCA_conf_data_h), ORCA_conf_data_h[:,3], color = '#808080')
ax3.scatter([0]*len(NONE_conf_data_h), NONE_conf_data_h[:,3], color = '#808080')
ax3.scatter([1]*len(AIR_conf_data_h), AIR_conf_data_h[:,3], color = '#808080')
ax3.boxplot([MVP_conf_data_h[:,3], ORCA_conf_data_h[:,3], 
             NONE_conf_data_h[:,3], AIR_conf_data_h[:,3]] , positions = [2,3,0,1])
ax3.set_xticks([0,1,2,3])
ax3.set_yticks([])
ax3.set_xticklabels(['N', 'A', 'M', 'O'])

ax1.set_ylim([-2, max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
ax2.set_ylim([-2, max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
ax3.set_ylim([-2, max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
plt.savefig('Plots/los_num.png')
plt.show()

# %%%% Number of geobreaches
plt.figure('geobreach_num', figsize = (5,3))
ax1 = plt.subplot(131)
ax1.set_title('Low density')
ax1.scatter([2]*len(MVP_conf_data_l), MVP_conf_data_l[:,4], color = '#808080')
ax1.scatter([3]*len(ORCA_conf_data_l), ORCA_conf_data_l[:,4], color = '#808080')
ax1.scatter([0]*len(NONE_conf_data_l), NONE_conf_data_l[:,4], color = '#808080')
ax1.scatter([1]*len(AIR_conf_data_l), AIR_conf_data_l[:,4], color = '#808080')
ax1.boxplot([MVP_conf_data_l[:,4], ORCA_conf_data_l[:,4], 
             NONE_conf_data_l[:,4], AIR_conf_data_l[:,4]] , positions = [2,3,0,1])
ax1.set_xticks([0,1,2,3])
ax1.set_xticklabels(['N', 'A', 'M', 'O'])
ax1.set_ylabel('Number of geofence breaches [-]')

ax2 = plt.subplot(132)
ax2.set_title('Medium density')
ax2.scatter([2]*len(MVP_conf_data_m), MVP_conf_data_m[:,4], color = '#808080')
ax2.scatter([3]*len(ORCA_conf_data_m), ORCA_conf_data_m[:,4], color = '#808080')
ax2.scatter([0]*len(NONE_conf_data_m), NONE_conf_data_m[:,4], color = '#808080')
ax2.scatter([1]*len(AIR_conf_data_m), AIR_conf_data_m[:,4], color = '#808080')
ax2.boxplot([MVP_conf_data_m[:,4], ORCA_conf_data_m[:,4], 
             NONE_conf_data_m[:,4], AIR_conf_data_m[:,4]] , positions = [2,3,0,1])
ax2.set_xticks([0,1,2,3])
ax2.set_yticks([])
ax2.set_xticklabels(['N', 'A', 'M', 'O'])

ax3 = plt.subplot(133)
ax3.set_title('High density')
ax3.scatter([2]*len(MVP_conf_data_h), MVP_conf_data_h[:,4], color = '#808080')
ax3.scatter([3]*len(ORCA_conf_data_h), ORCA_conf_data_h[:,4], color = '#808080')
ax3.scatter([0]*len(NONE_conf_data_h), NONE_conf_data_h[:,4], color = '#808080')
ax3.scatter([1]*len(AIR_conf_data_h), AIR_conf_data_h[:,4], color = '#808080')
ax3.boxplot([MVP_conf_data_h[:,4], ORCA_conf_data_h[:,4], 
             NONE_conf_data_h[:,4], AIR_conf_data_h[:,4]] , positions = [2,3,0,1])
ax3.set_xticks([0,1,2,3])
ax3.set_yticks([])
ax3.set_xticklabels(['N', 'A', 'M', 'O'])

ax1.set_ylim([-2, max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
ax2.set_ylim([-2, max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
ax3.set_ylim([-2, max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
plt.savefig('Plots/geobreach_num.png')
plt.show()


# %% Process FLSTLOG
# %%% Get needed files
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
MVP_flst_data = np.array(MVP_flst_data, dtype = object)
ORCA_flst_data = np.array(ORCA_flst_data, dtype = object)
VO_flst_data = np.array(VO_flst_data, dtype = object)
NONE_flst_data = np.array(NONE_flst_data, dtype = object)
AIR_flst_data = np.array(AIR_flst_data, dtype = object)

# %%% Divide the data by density
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
# %%% Flight time
# %%%% Data gathering
MVP_flst_acdiff_l = []
MVP_flst_distance_l = []
MVP_flst_avgflighttime_l = []
for scenario_data in MVP_flst_data_l:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    MVP_flst_acdiff_l.append(max_num_ac - del_num_ac)
    MVP_flst_distance_l.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    MVP_flst_avgflighttime_l.append(np.average(flight_times))

ORCA_flst_acdiff_l = []
ORCA_flst_distance_l = []
ORCA_flst_avgflighttime_l = []
for scenario_data in ORCA_flst_data_l:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    ORCA_flst_acdiff_l.append(max_num_ac - del_num_ac)
    ORCA_flst_distance_l.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    ORCA_flst_avgflighttime_l.append(np.average(flight_times))

VO_flst_acdiff_l = []
VO_flst_distance_l = []
VO_flst_avgflighttime_l = []
for scenario_data in VO_flst_data_l:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    VO_flst_acdiff_l.append(max_num_ac - del_num_ac)
    VO_flst_distance_l.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    VO_flst_avgflighttime_l.append(np.average(flight_times))

NONE_flst_acdiff_l = []
NONE_flst_distance_l = []
NONE_flst_avgflighttime_l = []
for scenario_data in NONE_flst_data_l:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    NONE_flst_acdiff_l.append(max_num_ac - del_num_ac)
    NONE_flst_distance_l.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    NONE_flst_avgflighttime_l.append(np.average(flight_times))

AIR_flst_acdiff_l = []
AIR_flst_distance_l = []
AIR_flst_avgflighttime_l = []
for scenario_data in AIR_flst_data_l:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    AIR_flst_acdiff_l.append(max_num_ac - del_num_ac)
    AIR_flst_distance_l.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    AIR_flst_avgflighttime_l.append(np.average(flight_times))
    
MVP_flst_acdiff_m = []
MVP_flst_distance_m = []
MVP_flst_avgflighttime_m = []
for scenario_data in MVP_flst_data_m:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    MVP_flst_acdiff_m.append(max_num_ac - del_num_ac)
    MVP_flst_distance_m.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    MVP_flst_avgflighttime_m.append(np.average(flight_times))

ORCA_flst_acdiff_m = []
ORCA_flst_distance_m = []
ORCA_flst_avgflighttime_m = []
for scenario_data in ORCA_flst_data_m:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    ORCA_flst_acdiff_m.append(max_num_ac - del_num_ac)
    ORCA_flst_distance_m.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    ORCA_flst_avgflighttime_m.append(np.average(flight_times))

VO_flst_acdiff_m = []
VO_flst_distance_m = []
VO_flst_avgflighttime_m = []
for scenario_data in VO_flst_data_m:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    VO_flst_acdiff_m.append(max_num_ac - del_num_ac)
    VO_flst_distance_m.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    VO_flst_avgflighttime_m.append(np.average(flight_times))

NONE_flst_acdiff_m = []
NONE_flst_distance_m = []
NONE_flst_avgflighttime_m = []
for scenario_data in NONE_flst_data_m:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    NONE_flst_acdiff_m.append(max_num_ac - del_num_ac)
    NONE_flst_distance_m.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    NONE_flst_avgflighttime_m.append(np.average(flight_times))

AIR_flst_acdiff_m = []
AIR_flst_distance_m = []
AIR_flst_avgflighttime_m = []
for scenario_data in AIR_flst_data_m:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    AIR_flst_acdiff_m.append(max_num_ac - del_num_ac)
    AIR_flst_distance_m.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    AIR_flst_avgflighttime_m.append(np.average(flight_times))
    
MVP_flst_acdiff_h = []
MVP_flst_distance_h = []
MVP_flst_avgflighttime_h = []
for scenario_data in MVP_flst_data_h:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    MVP_flst_acdiff_h.append(max_num_ac - del_num_ac)
    MVP_flst_distance_h.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    MVP_flst_avgflighttime_h.append(np.average(flight_times))

ORCA_flst_acdiff_h = []
ORCA_flst_distance_h = []
ORCA_flst_avgflighttime_h = []
for scenario_data in ORCA_flst_data_h:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    ORCA_flst_acdiff_h.append(max_num_ac - del_num_ac)
    ORCA_flst_distance_h.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    ORCA_flst_avgflighttime_h.append(np.average(flight_times))

VO_flst_acdiff_h = []
VO_flst_distance_h = []
VO_flst_avgflighttime_h = []
for scenario_data in VO_flst_data_h:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    VO_flst_acdiff_h.append(max_num_ac - del_num_ac)
    VO_flst_distance_h.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    VO_flst_avgflighttime_h.append(np.average(flight_times))

NONE_flst_acdiff_h = []
NONE_flst_distance_h = []
NONE_flst_avgflighttime_h = []
for scenario_data in NONE_flst_data_h:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    NONE_flst_acdiff_h.append(max_num_ac - del_num_ac)
    NONE_flst_distance_h.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    NONE_flst_avgflighttime_h.append(np.average(flight_times))

AIR_flst_acdiff_h = []
AIR_flst_distance_h = []
AIR_flst_avgflighttime_h = []
for scenario_data in AIR_flst_data_h:
    max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
    del_num_ac = len(scenario_data[1])
    AIR_flst_acdiff_h.append(max_num_ac - del_num_ac)
    AIR_flst_distance_h.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
    flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
    AIR_flst_avgflighttime_h.append(np.average(flight_times))

# %%%% Graphs
plt.figure('avg_flight_time', figsize = (5,3))
ax1 = plt.subplot(131)
ax1.set_title('Low density')
ax1.scatter([2]*len(MVP_flst_avgflighttime_l), MVP_flst_avgflighttime_l, color = '#808080')
ax1.scatter([3]*len(ORCA_flst_avgflighttime_l), ORCA_flst_avgflighttime_l, color = '#808080')
ax1.scatter([0]*len(NONE_flst_avgflighttime_l), NONE_flst_avgflighttime_l, color = '#808080')
ax1.scatter([1]*len(AIR_flst_avgflighttime_l), AIR_flst_avgflighttime_l, color = '#808080')
ax1.boxplot([MVP_flst_avgflighttime_l, ORCA_flst_avgflighttime_l,
             NONE_flst_avgflighttime_l, AIR_flst_avgflighttime_l], 
            positions = [2,3,0,1])
ax1.set_xticks([0,1,2,3])
ax1.set_xticklabels(['N', 'A', 'M', 'O'])
ax1.set_ylabel('Seconds [s]')

ax2 = plt.subplot(132)
ax2.set_title('Medium density')
ax2.scatter([2]*len(MVP_flst_avgflighttime_m), MVP_flst_avgflighttime_m, color = '#808080')
ax2.scatter([3]*len(ORCA_flst_avgflighttime_m), ORCA_flst_avgflighttime_m, color = '#808080')
ax2.scatter([0]*len(NONE_flst_avgflighttime_m), NONE_flst_avgflighttime_m, color = '#808080')
ax2.scatter([1]*len(AIR_flst_avgflighttime_m), AIR_flst_avgflighttime_m, color = '#808080')
ax2.boxplot([MVP_flst_avgflighttime_m, ORCA_flst_avgflighttime_m,
             NONE_flst_avgflighttime_m, AIR_flst_avgflighttime_m], 
            positions = [2,3,0,1])
ax2.set_xticks([0,1,2,3])
ax2.set_yticks([])
ax2.set_xticklabels(['N', 'A', 'M', 'O'])

ax3 = plt.subplot(133)
ax3.set_title('High density')
ax3.scatter([2]*len(MVP_flst_avgflighttime_h), MVP_flst_avgflighttime_h, color = '#808080')
ax3.scatter([3]*len(ORCA_flst_avgflighttime_h), ORCA_flst_avgflighttime_h, color = '#808080')
ax3.scatter([0]*len(NONE_flst_avgflighttime_h), NONE_flst_avgflighttime_h, color = '#808080')
ax3.scatter([1]*len(AIR_flst_avgflighttime_h), AIR_flst_avgflighttime_h, color = '#808080')
ax3.boxplot([MVP_flst_avgflighttime_h, ORCA_flst_avgflighttime_h,
             NONE_flst_avgflighttime_h, AIR_flst_avgflighttime_h], 
            positions = [2,3,0,1])
ax3.set_xticks([0,1,2,3])
ax3.set_yticks([])
ax3.set_xticklabels(['N', 'A', 'M', 'O'])

ax1.set_ylim([min(ax1.get_ylim()[0], ax2.get_ylim()[0], ax3.get_ylim()[0]), max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
ax2.set_ylim([min(ax1.get_ylim()[0], ax2.get_ylim()[0], ax3.get_ylim()[0]), max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
ax3.set_ylim([min(ax1.get_ylim()[0], ax2.get_ylim()[0], ax3.get_ylim()[0]), max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
plt.savefig('plots/avg_flight_time.png')
plt.show()

# %%% Aircraft that didn't make it
# %%%% Data gathering
# We compare number of rows in flstlog with the total number of aircraft 
# in the scenario
# Actually nothing to compare, only VO has aircraft that don't reach and
# we're dropping VO

# %%% Distance travelled
# %%%% Graphs
plt.figure('avg_distance', figsize = (5,3))
ax1 = plt.subplot(131)
ax1.set_title('Low density')
ax1.scatter([2]*len(MVP_flst_distance_l), MVP_flst_distance_l, color = '#808080')
ax1.scatter([3]*len(ORCA_flst_distance_l), ORCA_flst_distance_l, color = '#808080')
ax1.scatter([0]*len(NONE_flst_distance_l), NONE_flst_distance_l, color = '#808080')
ax1.scatter([1]*len(AIR_flst_distance_l), AIR_flst_distance_l, color = '#808080')
ax1.boxplot([MVP_flst_distance_l, ORCA_flst_distance_l,
             NONE_flst_distance_l, AIR_flst_distance_l], 
            positions = [2,3,0,1])
ax1.set_xticks([0,1,2,3])
ax1.set_xticklabels(['N', 'A', 'M', 'O'])
ax1.set_ylabel('Meters [m]')

ax2 = plt.subplot(132)
ax2.set_title('Medium density')
ax2.scatter([2]*len(MVP_flst_distance_m), MVP_flst_distance_m, color = '#808080')
ax2.scatter([3]*len(ORCA_flst_distance_m), ORCA_flst_distance_m, color = '#808080')
ax2.scatter([0]*len(NONE_flst_distance_m), NONE_flst_distance_m, color = '#808080')
ax2.scatter([1]*len(AIR_flst_distance_m), AIR_flst_distance_m, color = '#808080')
ax2.boxplot([MVP_flst_distance_m, ORCA_flst_distance_m,
             NONE_flst_distance_m, AIR_flst_distance_m], 
            positions = [2,3,0,1])
ax2.set_xticks([0,1,2,3])
ax2.set_yticks([])
ax2.set_xticklabels(['N', 'A', 'M', 'O'])

ax3 = plt.subplot(133)
ax3.set_title('High density')
ax3.scatter([2]*len(MVP_flst_distance_h), MVP_flst_distance_h, color = '#808080')
ax3.scatter([3]*len(ORCA_flst_distance_h), ORCA_flst_distance_h, color = '#808080')
ax3.scatter([0]*len(NONE_flst_distance_h), NONE_flst_distance_h, color = '#808080')
ax3.scatter([1]*len(AIR_flst_distance_h), AIR_flst_distance_h, color = '#808080')
ax3.boxplot([MVP_flst_distance_h, ORCA_flst_distance_h,
             NONE_flst_distance_h, AIR_flst_distance_h], 
            positions = [2,3,0,1])
ax3.set_xticks([0,1,2,3])
ax3.set_yticks([])
ax3.set_xticklabels(['N', 'A', 'M', 'O'])

ax1.set_ylim([min(ax1.get_ylim()[0], ax2.get_ylim()[0], ax3.get_ylim()[0]), max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
ax2.set_ylim([min(ax1.get_ylim()[0], ax2.get_ylim()[0], ax3.get_ylim()[0]), max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
ax3.set_ylim([min(ax1.get_ylim()[0], ax2.get_ylim()[0], ax3.get_ylim()[0]), max(ax1.get_ylim()[1], ax2.get_ylim()[1], ax3.get_ylim()[1])])
plt.savefig('plots/avg_distance.png')
plt.show()

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

