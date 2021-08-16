#processes FLSTLOG
import numpy as np
import os

files = os.listdir('Attempt_3')
filtered_files = []

#Get needed files
for i, file in enumerate(files):
    if file[0:4] == 'FLST':
        filtered_files.append(file)

MVP_data = np.array([])
ORCA_data = np.array([])
VO_data = np.array([])
NONE_data = np.array([])
AIR_data = np.array([])

for file in filtered_files:
    # Position 22 is always the aimed-for flight density
    flight_den = int(file[22])
    # Extract data
    data = np.genfromtxt('Attempt_3/' + file, delimiter = ',', dtype = int)