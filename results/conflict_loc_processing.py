# %% Processes CONFLOG
import numpy as np
import os
import geopandas as gpd
import pandas as pd
import osmnx as ox
from shapely.geometry import MultiPoint

nm2m = 1852

results_folder = 'Attempt_8'

files = os.listdir(results_folder)
filtered_files = []

# %%% Get needed files
for i, file in enumerate(files):
    if file[0:4] == 'CONF':
        filtered_files.append(file)


MVP_conf_gdf = MVPC_conf_gdf = ORCA_conf_gdf = ORCAC_conf_gdf = NONE_conf_gdf = \
MVP_air_conf_gdf = MVPC_air_conf_gdf = ORCA_air_conf_gdf = ORCAC_air_conf_gdf = NONE_air_conf_gdf = file =\
     gpd.GeoDataFrame(columns=['flight_density','lat1','lon1','h1','lat2','lon2','h2','geometry','time'], crs='epsg:4326')

for file in filtered_files:
    # Position 22 (and maybe 23) is always the aimed-for flight density
    if file[23] !='_':
        flight_den = int(file[22:24])
    else:
        flight_den = int(file[22])

    # Extract data
    data = np.genfromtxt(results_folder + '/' + file, delimiter = ',', dtype = float)
    if data.size == 0:
        #Skip empty files
        continue
    # get last 6 columns and remove anything without zero
    data = data[data[:,4] != 0]
    time = data[:,0]
    data = data[:,4:]
    n_row, n_col = data.shape

    # add flight density at first column
    data = np.c_[np.ones((n_row,1))*flight_den, data[:,0:2], data[:,2]*3.2804, data[:,3:5], data[:,5]*3.2804, time]

    # translate to dataframe and add geometry as multipoint
    df = pd.DataFrame(data=data, columns=['flight_density','lat1','lon1','h1','lat2','lon2','h2','time'])        
    df['geometry'] = df.apply(lambda x: MultiPoint([ (x['lon1'],x['lat1']), (x['lon2'],x['lat2'])]), axis=1)
    
    # add file name
    df['file'] = file

    # Append data to appropriate place
    if '_MVP_' in file and 'air' not in file:
        MVP_conf_gdf = MVP_conf_gdf.append(df)
    elif '_MVPC_' in file and 'air' not in file:
        MVPC_conf_gdf = MVP_conf_gdf.append(df)
    elif '_ORCA_' in file and 'air' not in file:
        ORCA_conf_gdf = ORCA_conf_gdf.append(df)
    elif '_ORCAC_' in file and 'air' not in file:
        ORCAC_conf_gdf = ORCAC_conf_gdf.append(df)
    elif '_NONE_' in file and 'air' not in file:
        NONE_conf_gdf = NONE_conf_gdf.append(df)
    elif '_MVP_air' in file:
        MVP_air_conf_gdf = MVP_air_conf_gdf.append(df)
    elif '_MVPC_air' in file:
        MVPC_air_conf_gdf = MVP_air_conf_gdf.append(df)
    elif '_ORCA_air' in file:
        ORCA_air_conf_gdf = ORCA_air_conf_gdf.append(df)
    elif '_ORCAC_air' in file:
        ORCAC_air_conf_gdf = ORCAC_air_conf_gdf.append(df)
    elif '_NONE_air' in file:
        NONE_air_conf_gdf = NONE_air_conf_gdf.append(df)        
    else:
        print('Something is wrong.')

# save to geopackages
MVP_conf_gdf.to_file(os.path.join('conflict_locs', 'MVP.gpkg'), driver='GPKG')
MVPC_conf_gdf.to_file(os.path.join('conflict_locs', 'MVPC.gpkg'), driver='GPKG')
ORCA_conf_gdf.to_file(os.path.join('conflict_locs', 'ORCA.gpkg'), driver='GPKG')
ORCAC_conf_gdf.to_file(os.path.join('conflict_locs', 'ORCAC.gpkg'), driver='GPKG')
NONE_conf_gdf.to_file(os.path.join('conflict_locs', 'NONE.gpkg'), driver='GPKG')

MVP_air_conf_gdf.to_file(os.path.join('conflict_locs', 'MVP_air.gpkg'), driver='GPKG')
MVPC_air_conf_gdf.to_file(os.path.join('conflict_locs', 'MVPC_air.gpkg'), driver='GPKG')
ORCA_air_conf_gdf.to_file(os.path.join('conflict_locs', 'ORCA_air.gpkg'), driver='GPKG')
ORCAC_air_conf_gdf.to_file(os.path.join('conflict_locs', 'ORCAC_air.gpkg'), driver='GPKG')
NONE_air_conf_gdf.to_file(os.path.join('conflict_locs', 'NONE_air.gpkg'), driver='GPKG')
