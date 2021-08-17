#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Michael Windisch
"""


from osgeo import gdal
import numpy as np

"""
PATH
"""

outpath = "SPECIFY_OUTPUT_PATH"
inpath = "SPECIFY_INPUT_PATH"
biome_names = ["ebf", "enf", "dbf", "dnf","otf"]



"""
LOAD OPPORTUNITY MASK
"""

rasterload = gdal.Open(f"{inpath}/nonforest_mask_GLC2000_gte10lte19.tif") # external dataset, see readme
opportunity_mask = np.array(rasterload.GetRasterBand(1).ReadAsArray())

rasterload = None


"""
LOAD BASE BIOMASS
"""

rasterload = gdal.Open(f"{inpath}/20181023_IPCC_Biomass.tif") # external dataset, see readme
base_biomass = np.array(rasterload.GetRasterBand(1).ReadAsArray())

rasterload = None



"""
LOAD FOREST DATASETS
"""

biomass_list = []

for i in range(0,len(biome_names)):
    rasterload = gdal.Open(inpath + biome_names[i] + "_biomass.tif") # external dataset, see readme
    biomass_list.append(np.array(rasterload.GetRasterBand(1).ReadAsArray(),dtype=float))
    
    rasterload = None #CLEAN MEMORY
    
    
biomass_nan = np.array(biomass_list,dtype=float)
biomass_nan[biomass_nan == 0.] = np.nan
biomass_nan[biomass_nan < 1.] = np.nan

biomass_list = None #Delete biomass_list to free memory

#Marker
print("Datasets Loaded")


"""
WINDOW CALC
"""

opportunity_location = np.argwhere(opportunity_mask==1) #Get indices of where calculations are necessary because gra/cro/shr exists here
opportunity_basevalue = np.zeros(len(opportunity_location))
opportunity_basevalue[:] = base_biomass[opportunity_location[:,0],opportunity_location[:,1]]

carbon_potential_grid = np.zeros(np.shape(biomass_nan[0,:,:]))
transition_grid = np.zeros(np.shape(biomass_nan[0,:,:]))+np.nan

windows_lat_min = np.int32(np.zeros(len(opportunity_location)))
windows_lat_max = np.int32(np.zeros(len(opportunity_location)))
windows_lon_min = np.int32(np.zeros(len(opportunity_location)))
windows_lon_max = np.int32(np.zeros(len(opportunity_location)))

ebf_windows = np.zeros(len(opportunity_location))
enf_windows = np.zeros(len(opportunity_location))
dbf_windows = np.zeros(len(opportunity_location))
dnf_windows = np.zeros(len(opportunity_location))
otf_windows = np.zeros(len(opportunity_location))




window = np.int16(30/2)


 
    

windows_lat_min[:] = opportunity_location[:,0]-window
windows_lat_max[:] = opportunity_location[:,0]+window
windows_lon_min[:] = opportunity_location[:,1]-window 
windows_lon_max[:] = opportunity_location[:,1]+window

#longitudes close to zero are possible
windows_lon_min[windows_lon_min<0]=0
windows_lon_max[windows_lon_max<0]=0


#Marker
print("Windows Set")


#Slow calc. Be aware.
n = len(opportunity_location)

for i in range(0,len(opportunity_location)):
    ebf_windows[i] = np.nanmean(biomass_nan[0,windows_lat_min[i]:windows_lat_max[i],windows_lon_min[i]:windows_lon_max[i]])
    enf_windows[i] = np.nanmean(biomass_nan[1,windows_lat_min[i]:windows_lat_max[i],windows_lon_min[i]:windows_lon_max[i]])
    dbf_windows[i] = np.nanmean(biomass_nan[2,windows_lat_min[i]:windows_lat_max[i],windows_lon_min[i]:windows_lon_max[i]])
    dnf_windows[i] = np.nanmean(biomass_nan[3,windows_lat_min[i]:windows_lat_max[i],windows_lon_min[i]:windows_lon_max[i]])
    otf_windows[i] = np.nanmean(biomass_nan[4,windows_lat_min[i]:windows_lat_max[i],windows_lon_min[i]:windows_lon_max[i]])
    
    if(i%100==0):
        print("\r{0}".format((float(i)/n)*100), end='\r'),


biomass_nan = None #Delete biomass_nan to free up memory

#Marker
print("Windows Done")

#Gather all tree biomes to a single array
tree_biomes = []
tree_biomes.append(ebf_windows)
tree_biomes.append(enf_windows)
tree_biomes.append(dbf_windows)
tree_biomes.append(dnf_windows)
tree_biomes.append(otf_windows)

ebf_windows= None
enf_windows=None
dbf_windows=None
dnf_windows=None
otf_windows=None


#FIND POTENTIAL
carbon_potential=[]
for i in range(0,len(tree_biomes)):
    intermediate_step = tree_biomes[i] - opportunity_basevalue
    intermediate_step[intermediate_step<=0]=0
    
    carbon_potential.append(intermediate_step)

    intermediate_step = None

#Marker
print("Potential Done")

#FIND TRANSITION


#transition = np.zeros(len(opportunity_location))+np.nan
#tree_biomes = np.array(tree_biomes)
#
#d0 = np.nanmax(tree_biomes, axis=0)
#nonnan_indices = np.argwhere(~np.isnan(d0))[:,0]
#
#transition[nonnan_indices] = np.nanargmax(tree_biomes[:,nonnan_indices], axis=0)
#
#transition = transition + 1 #so first transition is 1 and not 0

##Marker
#print("Transition Done")




#ASSIGN THE CARBON POTENTIAL BACK TO THE ORIGINAL GRID


for i in range(0,len(biome_names)):
    carbon_potential_grid = np.zeros((21121,43201))
    carbon_potential_grid[opportunity_location[:,0],opportunity_location[:,1]] = carbon_potential[i][:]
    
    ds = gdal.Open(outpath + biome_names[i] + "_biomass_potential.tif", gdal.GA_Update) 
    ds.GetRasterBand(1).WriteArray(carbon_potential_grid)
    ds = None
    carbon_potential_grid=None

    






