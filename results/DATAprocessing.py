# %% Processes CONFLOG
import numpy as np
import os
import matplotlib.pyplot as plt

nm2m = 1852

results_folder = 'Attempt_8'

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

densities = ['Low', 'Medium', 'High']
methods = ['N', 'A', 'M', 'MA', 'O', 'OA','MC','MCA', 'CR', 'CRA']
to_plot = [0,1,8,9]

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


global_conf_init = [NONE_conf_data, NONE_AIR_conf_data, 
              MVP_conf_data, MVP_AIR_conf_data,
              ORCA_conf_data, ORCA_AIR_conf_data,
              MVPC_conf_data, MVPC_AIR_conf_data, 
              ORCAC_conf_data, ORCAC_AIR_conf_data]

# %%% Divide the data by density
# Low
global_conf_l = []
for data in global_conf_init:
    if len(data) == 0:
        continue
    low_data = data[np.where(data[:,0] == 5)]
    global_conf_l.append(low_data)

# Medium
global_conf_m = []
for data in global_conf_init:
    if len(data) == 0:
        continue
    medium_data = data[np.where(data[:,0] == 8)]
    global_conf_m.append(medium_data)

# High
global_conf_h = []
for data in global_conf_init:
    if len(data) == 0:
        continue
    high_data = data[np.where(data[:,0] == 11)]
    global_conf_h.append(high_data)
    
global_conf = [global_conf_l, global_conf_m, global_conf_h]


# %%% Graphs
ylabels = ['Number of conflicts [-]', 'Number of losses of separation [-]', 'Number of geofence breaches [-]']
for k, figure in enumerate(['conf_num', 'los_num', 'geobreach_num']):
    axes = dict()
    plt.figure(figure, figsize = (5,3))
    for i, density in enumerate(densities): #Representing [5,8,11]
        # Make subplot
        axes[density] = plt.subplot(131 + i)
        axes[density].set_title(density + ' density')
        if i == 0:
            axes[density].set_ylabel(ylabels[k])
        else:
            axes[density].set_yticks([])
        counter = 0
        for j in to_plot: # Representing N, NA, M, MA, O, OA, MC, MCA, OC, OCA
            method = methods[j]
            axes[density].scatter([counter] * len(global_conf[i][j]), global_conf[i][j][:,k+2], color = '#808080')
            axes[density].boxplot(global_conf[i][j][:,k+2], positions = [counter])
            counter += 1
        axes[density].set_xticks(range(len(to_plot)))
        axes[density].set_xticklabels([methods[x] for x in to_plot])
    
        
    big_y = max([axes[x].get_ylim()[1] for x in axes])
    for ax in axes:
        axes[ax].set_ylim([-2, big_y])
    plt.savefig('Plots/'+figure+'.png')
    plt.show()


# %% Process FLSTLOG
# %%% Get needed files
filtered_files = []
for i, file in enumerate(files):
    if file[0:4] == 'FLST':
        filtered_files.append(file)

NONE_flst_data = []
NONE_AIR_flst_data = []
MVP_flst_data = []
MVP_AIR_flst_data = []
ORCA_flst_data = []
ORCA_AIR_flst_data = []
MVPC_flst_data = []
MVPC_AIR_flst_data = []
ORCAC_flst_data = []
ORCAC_AIR_flst_data = []

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
                          'int', 'float'))
    if data.size == 0:
        #Skip empty files
        continue
    # We kind of only want the last row
    needed_data = list(data)
    # Insert flight density as first number
    needed_data = list((flight_den, needed_data))
    # Append data to appropriate place
    if '_MVP_' in file:
        if '_air_' in file:
            MVP_AIR_flst_data.append(needed_data)
        else:
            MVP_flst_data.append(needed_data)
    elif '_MVPC_' in file:
        if '_air_' in file:
            MVPC_AIR_flst_data.append(needed_data)
        else:
            MVPC_flst_data.append(needed_data)
    elif '_ORCA_' in file:
        if '_air_' in file:
            ORCA_AIR_flst_data.append(needed_data)
        else:
            ORCA_flst_data.append(needed_data)
    elif '_ORCAC_' in file:
        if '_air_' in file:
            ORCAC_AIR_flst_data.append(needed_data)
        else:
            ORCAC_flst_data.append(needed_data)
    elif '_NONE_' in file:
        if '_air_' in file:
            NONE_AIR_flst_data.append(needed_data)
        else:
            NONE_flst_data.append(needed_data)
    else:
        print('Something is wrong.')

# Convert to numpy arrays
NONE_flst_data = np.array(NONE_flst_data, dtype = object)
NONE_AIR_flst_data = np.array(NONE_AIR_flst_data, dtype = object)
MVP_flst_data = np.array(MVP_flst_data, dtype = object)
MVP_AIR_flst_data = np.array(MVP_AIR_flst_data, dtype = object)
ORCA_flst_data = np.array(ORCA_flst_data, dtype = object)
ORCA_AIR_flst_data = np.array(ORCA_AIR_flst_data, dtype = object)
MVPC_flst_data = np.array(MVPC_flst_data, dtype = object)
MVPC_AIR_flst_data = np.array(MVPC_AIR_flst_data, dtype = object)
ORCAC_flst_data = np.array(ORCAC_flst_data, dtype = object)
ORCAC_AIR_flst_data = np.array(ORCAC_AIR_flst_data, dtype = object)


global_flst_init = [NONE_flst_data, NONE_AIR_flst_data, 
          MVP_flst_data, MVP_AIR_flst_data,
          ORCA_flst_data, ORCA_AIR_flst_data,
          MVPC_flst_data, MVPC_AIR_flst_data, 
          ORCAC_flst_data, ORCAC_AIR_flst_data]

# %%% Divide the data by density
# Low
global_flst_l = []
for data in global_flst_init:
    if len(data) == 0:
        continue
    low_data = data[np.where(data[:,0] == 5)]
    global_flst_l.append(low_data)

# Medium
global_flst_m = []
for data in global_flst_init:
    if len(data) == 0:
        continue
    medium_data = data[np.where(data[:,0] == 8)]
    global_flst_m.append(medium_data)

# High
global_flst_h = []
for data in global_flst_init:
    if len(data) == 0:
        continue
    high_data = data[np.where(data[:,0] == 11)]
    global_flst_h.append(high_data)
    
global_flst = [global_flst_l, global_flst_m, global_flst_h]

# In this one we need more data processing
# %%% Flight time
# %%%% Data gathering
ylabels = ['Average flight time [s]', 'Average travelled distance [m]']
for k, figure in enumerate(['avg_flight_time', 'avg_distance']):
    axes = dict()
    plt.figure(figure, figsize = (5,3))
    for i, density in enumerate(densities): #Representing [5,8,11]
        # Make subplot
        axes[density] = plt.subplot(131 + i)
        axes[density].set_title(density + ' density')
        if i == 0:
            axes[density].set_ylabel(ylabels[k])
        else:
            axes[density].set_yticks([])
        counter = 0
        for j in to_plot: # Representing N, NA, M, MA, O, OA, MC, MCA, OC, OCA
            method = methods[j]
            flst_acdiff = []
            flst_distance = []
            flst_avgflighttime = []
            for scenario_data in global_flst[i][j]:
                max_num_ac = max([int(x[1].replace(b'D', b'')) for x in scenario_data[1]])
                del_num_ac = len(scenario_data[1])
                flst_acdiff.append(max_num_ac - del_num_ac)
                flst_distance.append(np.average([x[5]*nm2m for x in scenario_data[1]]))
                flight_times = np.array([x[0] - x[2] for x in scenario_data[1]])
                flst_avgflighttime.append(np.average(flight_times))
            if k == 0:
                axes[density].scatter([counter] * len(flst_avgflighttime), flst_avgflighttime, color = '#808080')
                axes[density].boxplot(flst_avgflighttime, positions = [counter])
            else:
                axes[density].scatter([counter] * len(flst_distance), flst_distance, color = '#808080')
                axes[density].boxplot(flst_distance, positions = [counter])
            counter += 1
        axes[density].set_xticks(range(len(to_plot)))
        axes[density].set_xticklabels([methods[x] for x in to_plot])
    
        
    big_y = max([axes[x].get_ylim()[1] for x in axes])
    small_y = min([axes[x].get_ylim()[0] for x in axes])
    for ax in axes:
        axes[ax].set_ylim([small_y, big_y])
    plt.savefig('Plots/'+figure+'.png')
    plt.show()


#%% Process REGLOG
#Get needed files
filtered_files = []
for i, file in enumerate(files):
    if file[0:3] == 'REG':
        filtered_files.append(file)

NONE_reg_data = []
NONE_AIR_reg_data = []
MVP_reg_data = []
MVP_AIR_reg_data = []
ORCA_reg_data = []
ORCA_AIR_reg_data = []
MVPC_reg_data = []
MVPC_AIR_reg_data = []
ORCAC_reg_data = []
ORCAC_AIR_reg_data = []

avg_dict = dict()
avg_dict[5] = []
avg_dict[8] = []
avg_dict[11] = []
max_dict = dict()
max_dict[5] = []
max_dict[8] = []
max_dict[11] = []

for file in filtered_files:
    # Position 22 (and maybe 23) is always the aimed-for flight density
    if file[22] !='_':
        flight_den = int(file[21:23])
    else:
        flight_den = int(file[21])
    # Extract data
    data = np.genfromtxt(results_folder + '/' + file, delimiter = ',', dtype = int)
    avg_dict[flight_den].append(np.average(data[:,1]))
    max_dict[flight_den].append(max(data[:,1]))
    if data.size == 0:
        #Skip empty files
        continue
    # We kind of only want the last row
    needed_data = np.array(data)
    # Insert flight density as first number
    needed_data = np.append(flight_den, needed_data)
    # Append data to appropriate place
    if '_MVP_' in file:
        if '_air_' in file:
            MVP_AIR_reg_data.append(needed_data)
        else:
            MVP_reg_data.append(needed_data)
    elif '_MVPC_' in file:
        if '_air_' in file:
            MVPC_AIR_reg_data.append(needed_data)
        else:
            MVPC_reg_data.append(needed_data)
    elif '_ORCA_' in file:
        if '_air_' in file:
            ORCA_AIR_reg_data.append(needed_data)
        else:
            ORCA_reg_data.append(needed_data)
    elif '_ORCAC_' in file:
        if '_air_' in file:
            ORCAC_AIR_reg_data.append(needed_data)
        else:
            ORCAC_reg_data.append(needed_data)
    elif '_NONE_' in file:
        if '_air_' in file:
            NONE_AIR_reg_data.append(needed_data)
        else:
            NONE_reg_data.append(needed_data)
    else:
        print('Something is wrong.')
        
# Convert to numpy arrays
NONE_reg_data = np.array(NONE_reg_data)
NONE_AIR_reg_data = np.array(NONE_AIR_reg_data)
MVP_reg_data = np.array(MVP_reg_data)
MVP_AIR_reg_data = np.array(MVP_AIR_reg_data)
ORCA_reg_data = np.array(ORCA_reg_data)
ORCA_AIR_reg_data = np.array(ORCA_AIR_reg_data)
MVPC_reg_data = np.array(MVPC_reg_data)
MVPC_AIR_reg_data = np.array(MVPC_AIR_reg_data)
ORCAC_reg_data = np.array(ORCAC_reg_data)
ORCAC_AIR_reg_data = np.array(ORCAC_AIR_reg_data)


global_reg_init = [NONE_reg_data, NONE_AIR_reg_data, 
          MVP_reg_data, MVP_AIR_reg_data,
          ORCA_reg_data, ORCA_AIR_reg_data,
          MVPC_reg_data, MVPC_AIR_reg_data, 
          ORCAC_reg_data, ORCAC_AIR_reg_data]

# Density
# Low
global_reg_l = []
for data in global_reg_init:
    if len(data) == 0:
        continue
    low_data = data[np.where(data[:,0] == 5)]
    global_reg_l.append(low_data)

# Medium
global_reg_m = []
for data in global_reg_init:
    if len(data) == 0:
        continue
    medium_data = data[np.where(data[:,0] == 8)]
    global_reg_m.append(medium_data)

# High
global_reg_h = []
for data in global_reg_init:
    if len(data) == 0:
        continue
    high_data = data[np.where(data[:,0] == 11)]
    global_reg_h.append(high_data)
    
global_reg = [global_reg_l, global_reg_m, global_reg_h]

# Get me average and peak flights in the air
for i, density in enumerate(densities): #Representing [5,8,11]
    for j in to_plot: # Representing N, NA, M, MA, O, OA, MC, MCA, OC, OCA
        method = methods[j]
        avg_flights = np.average(global_reg[i][j][:,2])

