#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Michael Windisch
"""
import numpy as np
from osgeo import gdal
import operator
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/


"""
Weighting / Decision
"""
def scenario(weight_c,weight_t,biomass_nan,LTS_nan,CRO_mask,GRA_mask,SHR_mask):
    weight_carbon = weight_c
    weight_temp = weight_t

    biomass_pre_weight = np.copy(biomass_nan)
    LTS_pre_weight = np.copy(LTS_nan)
    
    #Go through the relevant pixels only to speed up processing

    nonan_3darray = ~np.isnan(biomass_nan)
    nonan_indices = np.any(nonan_3darray,axis=0)
    nonan_3darray_LTS = ~np.isnan(LTS_nan)
    nonan_indices_LTS = np.any(nonan_3darray_LTS,axis=0)


    biomass_nan = biomass_nan*weight_carbon
    LTS_nan = LTS_nan*weight_temp





    transitions = []

    transitions.append(LTS_nan[0]+biomass_nan[0]) #norm_cro_ebf T0
    transitions.append(LTS_nan[1]+biomass_nan[0]) #norm_gra_ebf T1
    transitions.append(LTS_nan[2]+biomass_nan[0]) #norm_shr_ebf T2
    transitions.append(LTS_nan[3]+biomass_nan[1]) #norm_cro_enf T3
    transitions.append(LTS_nan[4]+biomass_nan[1]) #norm_gra_enf T4
    transitions.append(LTS_nan[5]+biomass_nan[1]) #norm_shr_enf T5
    transitions.append(LTS_nan[6]+biomass_nan[2]) #norm_cro_dbf T6
    transitions.append(LTS_nan[7]+biomass_nan[2]) #norm_gra_dbf T7
    transitions.append(LTS_nan[8]+biomass_nan[2]) #norm_shr_dbf T8
    transitions.append(LTS_nan[9]+biomass_nan[3]) #norm_cro_dnf T9
    transitions.append(LTS_nan[10]+biomass_nan[3]) #norm_gra_dnf T10
    transitions.append(LTS_nan[11]+biomass_nan[3]) #norm_shr_dnf T11
    
    

    transitions.append(np.nanmean(LTS_nan[0:12:3],axis=0)+biomass_nan[4]) #norm_cro_otf T12
    transitions.append(np.nanmean(LTS_nan[1:12:3],axis=0)+biomass_nan[4]) #norm_gra_otf T13
    transitions.append(np.nanmean(LTS_nan[2:12:3],axis=0)+biomass_nan[4]) #norm_shr_otf T14


    #CONSIDER THE ACTUAL LANDCOVER
    transitions[0:15:3] = transitions[0:15:3]*CRO_mask
    transitions[1:15:3] = transitions[1:15:3]*GRA_mask
    transitions[2:15:3] = transitions[2:15:3]*SHR_mask




    value_best_trans = np.nanmax(transitions,axis=0)
    value_best_trans[np.logical_not(nonan_indices)]=np.nan

    #mask transitions to eliviate problem with nan slices
    for i in range(0,len(transitions)):

        transitions[i] = np.ma.masked_invalid(transitions[i])



    index_best_trans = np.ma.argmax(transitions,axis=0,fill_value=-999999) #nanargmax does not work on full nan slices!


    logical_zero = (index_best_trans==0)
    logical_ebf = np.isnan(biomass_nan[0])
    logical_false_zero = np.logical_and(logical_ebf,logical_zero)

    index_best_trans[logical_false_zero==True]=33
    index_best_trans[np.logical_not(nonan_indices)]=33
    index_best_trans[np.logical_not(nonan_indices_LTS)]=33
    
    index_best_trans = np.ma.masked_where(index_best_trans==33,index_best_trans)
    


    """
    Reallocate the individual dT and dC values for display
    """



    dC_grid = np.zeros(np.shape(biomass_pre_weight[0,:,:]))+np.nan
    dT_grid = np.zeros(np.shape(LTS_pre_weight[0,:,:]))+np.nan

    for i in range(0,12): #must include all transitions but exclude 33

        location = np.argwhere(index_best_trans==i)

        dC_intervalue = np.zeros(len(location))+np.nan
        dC_intervalue[:] = biomass_pre_weight[np.int16(i/3)][location[:,0],location[:,1]]

        dT_intervalue = np.zeros(len(location))+np.nan
        dT_intervalue[:] = LTS_pre_weight[i][location[:,0],location[:,1]]

        #Assign back onto grid

        dC_grid[location[:,0],location[:,1]] = dC_intervalue[:]
        dT_grid[location[:,0],location[:,1]] = dT_intervalue[:]


    #For OTFs / Transition 8 and 9
    location = np.argwhere(index_best_trans==12)
    dC_intervalue = np.zeros(len(location))+np.nan
    dC_intervalue[:] = biomass_pre_weight[4][location[:,0],location[:,1]]
    dC_grid[location[:,0],location[:,1]] = dC_intervalue[:]

    dT_intervalue = np.zeros(len(location))+np.nan
    dT_intervalue[:] = np.nanmean(LTS_nan[0:12:3],axis=0)[location[:,0],location[:,1]]
    dT_grid[location[:,0],location[:,1]] = dT_intervalue[:]


    location = np.argwhere(index_best_trans==13)
    dC_intervalue = np.zeros(len(location))+np.nan
    dC_intervalue[:] = biomass_pre_weight[4][location[:,0],location[:,1]]
    dC_grid[location[:,0],location[:,1]] = dC_intervalue[:]

    dT_intervalue = np.zeros(len(location))+np.nan
    dT_intervalue[:] = np.nanmean(LTS_nan[1:12:3],axis=0)[location[:,0],location[:,1]]
    dT_grid[location[:,0],location[:,1]] = dT_intervalue[:]
    
    
    location = np.argwhere(index_best_trans==14)
    dC_intervalue = np.zeros(len(location))+np.nan
    dC_intervalue[:] = biomass_pre_weight[4][location[:,0],location[:,1]]
    dC_grid[location[:,0],location[:,1]] = dC_intervalue[:]

    dT_intervalue = np.zeros(len(location))+np.nan
    dT_intervalue[:] = np.nanmean(LTS_nan[2:12:3],axis=0)[location[:,0],location[:,1]]
    dT_grid[location[:,0],location[:,1]] = dT_intervalue[:]

    local_value = dC_grid+dT_grid
    
    
    
    


    
    
    return(local_value,index_best_trans,dT_grid,dC_grid)













"""
LOAD REFOR CARBON DATASETS
"""
def load_refor_carbon():

    in_path = "SPECIFY_INPUT_PATH" #
    biome_names = ["ebf", "enf", "dbf", "dnf","otf"]
    
    biomass_list = []
    add_zeros = np.zeros((24,4320))
    
    for i in range(0,len(biome_names)):
        rasterload = gdal.Open(in_path + biome_names[i] + "_biomass_potential.tif") # output of 1_biomass_potential.py
        biomass_list.append(np.concatenate((add_zeros,np.array(rasterload.GetRasterBand(1).ReadAsArray(),dtype=float),add_zeros),axis=0))
            
        rasterload = None #CLEAN MEMORY
        
    
        
    biomass_nan = np.array(biomass_list,dtype=float)
    biomass_nan[biomass_nan == 0.] = np.nan
    biomass_nan[biomass_nan < 1.] = np.nan
    
    biomass_list = None #Delete biomass_list to free memory
    
    return(biomass_nan)
    
    
"""
LOAD REFOR MIN CARBON DATASETS
"""
def load_refor_carbon_min():

    in_path = "SPECIFY_INPUT_PATH"
    biome_names = ["ebf", "enf", "dbf", "dnf","otf"]
    
    biomass_list = []
    add_zeros = np.zeros((24,4320))
    
    for i in range(0,len(biome_names)):
        rasterload = gdal.Open(in_path + biome_names[i] + "_biomass_potential_min_1SD.tif") # output of 1_biomass_potential.py
        biomass_list.append(np.concatenate((add_zeros,np.array(rasterload.GetRasterBand(1).ReadAsArray(),dtype=float),add_zeros),axis=0))
            
        rasterload = None #CLEAN MEMORY
        
    
        
    biomass_nan = np.array(biomass_list,dtype=float)
    biomass_nan[biomass_nan == 0.] = np.nan
    biomass_nan[biomass_nan < 1.] = np.nan
    
    biomass_list = None #Delete biomass_list to free memory
    
    return(biomass_nan)



"""
LOAD REFOR MAX CARBON DATASETS
"""
def load_refor_carbon_max():

    in_path = "SPECIFY_INPUT_PATH"
    biome_names = ["ebf", "enf", "dbf", "dnf","otf"]
    
    biomass_list = []
    add_zeros = np.zeros((24,4320))
    
    for i in range(0,len(biome_names)):
        rasterload = gdal.Open(in_path + biome_names[i] + "_biomass_potential_max_1SD.tif") # output of 1_biomass_potential.py
        biomass_list.append(np.concatenate((add_zeros,np.array(rasterload.GetRasterBand(1).ReadAsArray(),dtype=float),add_zeros),axis=0))
            
        rasterload = None #CLEAN MEMORY
        
    
        
    biomass_nan = np.array(biomass_list,dtype=float)
    biomass_nan[biomass_nan == 0.] = np.nan
    biomass_nan[biomass_nan < 1.] = np.nan
    
    biomass_list = None #Delete biomass_list to free memory
    
    return(biomass_nan)





"""
LOAD DEFOR CARBON DATASETS
"""
def load_defor_carbon():
    in_path = "SPECIFY_INPUT_PATH"
    biome_names = ["ebf", "enf", "dbf", "dnf","otf"]
    
    biomass_list = []
    add_zeros = np.zeros((24,4320))
    
    for i in range(0,len(biome_names)):
        rasterload = gdal.Open(in_path + biome_names[i] + "_biomass.tif") # external dataset, see readme
        biomass_list.append(np.concatenate((add_zeros,np.array(rasterload.GetRasterBand(1).ReadAsArray(),dtype=float),add_zeros),axis=0))
            
        rasterload = None #CLEAN MEMORY
        
    
        
    biomass_nan = np.array(biomass_list,dtype=float)
    biomass_nan = biomass_nan-500 #MINUS GLOBAL CULTIVATED LAND CARBON VALUE OF THE T1 ASSESSMENT
    biomass_nan[biomass_nan == 0.] = np.nan
    biomass_nan[biomass_nan < 1.] = np.nan
    
    biomass_list = None #Delete biomass_list to free memory
    
    return(biomass_nan)
    
    

"""
LOAD DEFOR MIN CARBON DATASETS
"""
def load_defor_carbon_min() :
    in_path = "SPECIFY_INPUT_PATH"
    biome_names = ["ebf", "enf", "dbf", "dnf","otf"]
    
    biomass_list = []
    add_zeros = np.zeros((24,4320))
    
    for i in range(0,len(biome_names)):
        rasterload = gdal.Open(in_path + biome_names[i] + "biomass_min_1SD.tif") # external dataset, see readme
        biomass_list.append(np.concatenate((add_zeros,np.array(rasterload.GetRasterBand(1).ReadAsArray(),dtype=float),add_zeros),axis=0))
            
        rasterload = None #CLEAN MEMORY
        
    
        
    biomass_nan = np.array(biomass_list,dtype=float)
    biomass_nan = biomass_nan-500 #MINUS GLOBAL CULTIVATED LAND CARBON VALUE OF THE T1 ASSESSMENT
    biomass_nan[biomass_nan == 0.] = np.nan
    biomass_nan[biomass_nan < 1.] = np.nan
    
    biomass_list = None #Delete biomass_list to free memory
    
    return(biomass_nan)    
    
    
    
"""
LOAD DEFOR MAX CARBON DATASETS
"""
def load_defor_carbon_max() :
    in_path = "SPECIFY_INPUT_PATH"
    biome_names = ["ebf", "enf", "dbf", "dnf","otf"]
    
    biomass_list = []
    add_zeros = np.zeros((24,4320))
    
    for i in range(0,len(biome_names)):
        rasterload = gdal.Open(in_path + biome_names[i] + "biomass_max_1SD.tif") # external dataset, see readme
        biomass_list.append(np.concatenate((add_zeros,np.array(rasterload.GetRasterBand(1).ReadAsArray(),dtype=float),add_zeros),axis=0))
            
        rasterload = None #CLEAN MEMORY
        
    
        
    biomass_nan = np.array(biomass_list,dtype=float)
    biomass_nan = biomass_nan-500 #MINUS GLOBAL CULTIVATED LAND CARBON VALUE OF THE T1 ASSESSMENT
    biomass_nan[biomass_nan == 0.] = np.nan
    biomass_nan[biomass_nan < 1.] = np.nan
    
    biomass_list = None #Delete biomass_list to free memory
    
    return(biomass_nan)  



"""
LOAD WINTER TEMPERATURE DATASETS
"""
def load_DJF_LTS(biomass_nan):
    in_path = "SPECIFY_INPUT_PATH"
    nc_D18 = f"{in_path}/D18_LST_IGBPdet.nc"  # external dataset, see readme
    nc_D18id = Dataset(nc_D18, "r")  # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class
    
    lats = nc_D18id.variables["lat"][:]
    lons = nc_D18id.variables["lon"][:]
    
    LTS_day_all = []
    LTS_night_all=[]
    
    
    for trans in range(0,45):
        for mon in range(0,12):
            LTS_day_all.append(nc_D18id.variables["Delta_LSTday"][mon,trans,:,:])
            LTS_night_all.append(nc_D18id.variables["Delta_LSTnight"][mon,trans,:,:])
            
            
    LTS_day_DJF = []        
    LTS_night_DJF = [] 
            
    for x in range(0,45):     
            LTS_day_array = np.ma.asarray(LTS_day_all)
            intermediate_day_array = LTS_day_array[x*12-1:(x+1)*12-10].mean(axis=0)
            LTS_day_DJF.append(intermediate_day_array)
            
            LTS_night_array = np.ma.asarray(LTS_night_all)
            intermediate_night_array = LTS_night_array[x*12-1:(x+1)*12-10].mean(axis=0)
            LTS_night_DJF.append(intermediate_night_array)
    
    LTS_DJF = (np.ma.asarray(LTS_day_DJF) + np.ma.asarray(LTS_night_DJF))/2   
    
    LTS_DJF_reverse = LTS_DJF*-1
    LTS_day_DJF_reverse = np.ma.asarray(LTS_day_DJF)*-1
    LTS_night_DJF_reverse = np.ma.asarray(LTS_night_DJF)*-1     
            
    
    LTS_D18_nan = []
    LTS_D18_nan.append(LTS_DJF_reverse[7,:,:]) #cro_ebf
    LTS_D18_nan.append(LTS_DJF_reverse[6,:,:]) #gra_ebf
    LTS_D18_nan.append(LTS_DJF_reverse[5,:,:]) #shr_ebf
    LTS_D18_nan.append(LTS_DJF_reverse[22,:,:])#cro_enf
    LTS_D18_nan.append(LTS_DJF_reverse[21,:,:])#gra_enf
    LTS_D18_nan.append(LTS_DJF_reverse[20,:,:])#shr_enf
    LTS_D18_nan.append(LTS_DJF_reverse[15,:,:])#cro_dbf
    LTS_D18_nan.append(LTS_DJF_reverse[14,:,:])#gra_dbf
    LTS_D18_nan.append(LTS_DJF_reverse[13,:,:])#shr_dbf
    LTS_D18_nan.append(LTS_DJF_reverse[28,:,:])#cro_dnf
    LTS_D18_nan.append(LTS_DJF_reverse[27,:,:])#gra_dnf
    LTS_D18_nan.append(LTS_DJF_reverse[26,:,:])#shr_dnf
    
    
    for z in range(0,len(LTS_D18_nan)):
        LTS_D18_nan[z] = np.ma.filled(LTS_D18_nan[z], np.nan)
    
    ###
    
    nc_f = f"{in_path}/B17_DTS.nc" # external dataset, see readme
    nc_fid = Dataset(nc_f, "r")  # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class
    
    lats = nc_fid.variables["lat"][:]
    lons = nc_fid.variables["lon"][:]
    
    dTs_ann = []
    
    for trans in range(0, 6):
    
        dTs_ann.append(nc_fid.variables["dTs_DJF"][trans,:,:]) #Transition, lat, lon
    
    
    cro_enf = dTs_ann[0]
    cro_dbf = dTs_ann[1]
    cro_ebf = dTs_ann[2]
    gra_enf = dTs_ann[3]
    gra_dbf = dTs_ann[4]
    gra_ebf = dTs_ann[5]
    
    LTS_B17 = []
    LTS_B17.append(cro_ebf)
    LTS_B17.append(gra_ebf)
    LTS_B17.append(LTS_D18_nan[2])
    LTS_B17.append(cro_enf)
    LTS_B17.append(gra_enf)
    LTS_B17.append(LTS_D18_nan[5])
    LTS_B17.append(cro_dbf)
    LTS_B17.append(gra_dbf)
    LTS_B17.append(LTS_D18_nan[8])
    LTS_B17.append(LTS_D18_nan[9])
    LTS_B17.append(LTS_D18_nan[10])
    LTS_B17.append(LTS_D18_nan[11])
    
    for z in range(0,len(LTS_B17)):
        LTS_B17[z] = np.ma.filled(LTS_B17[z], np.nan)
    
    
      
    
    
###MEAN
    LTS_out=[]
    LTS_nan=[]
    for i in range(0,len(LTS_B17)):
        LTS_intermediate_step = []
        LTS_intermediate_step.append(LTS_D18_nan[i])
        LTS_intermediate_step.append(LTS_B17[i])
        LTS_nan.append(np.nanmean(LTS_intermediate_step,axis=0))  
    
    
    #UPSCALE TO FIT THE BIOMASS ARRAY
    multi_lat = np.int16(np.shape(biomass_nan[0,:,:])[0]/np.shape(LTS_nan[0])[0])
    multi_lon = np.int16(np.shape(biomass_nan[0,:,:])[1]/np.shape(LTS_nan[0])[1])
    
    for i in range(0,len(LTS_nan)):
        LTS_nan[i] = np.kron(LTS_nan[i],np.ones((multi_lat,multi_lon)))
    
    LTS_nan = np.array(LTS_nan,dtype=float)
    
    LTS_out.append(LTS_nan)
    
###LOW
    LTS_nan=[]
    for i in range(0,len(LTS_B17)):
        LTS_intermediate_step = []
        LTS_intermediate_step.append(LTS_D18_nan[i])
        LTS_intermediate_step.append(LTS_B17[i])
        LTS_nan.append(np.nanmin(LTS_intermediate_step,axis=0))  
    
    
    #UPSCALE TO FIT THE BIOMASS ARRAY
    multi_lat = np.int16(np.shape(biomass_nan[0,:,:])[0]/np.shape(LTS_nan[0])[0])
    multi_lon = np.int16(np.shape(biomass_nan[0,:,:])[1]/np.shape(LTS_nan[0])[1])
    
    for i in range(0,len(LTS_nan)):
        LTS_nan[i] = np.kron(LTS_nan[i],np.ones((multi_lat,multi_lon)))
    
    LTS_nan = np.array(LTS_nan,dtype=float)
    
    LTS_out.append(LTS_nan)
    
    
###HIGH
    LTS_nan=[]
    for i in range(0,len(LTS_B17)):
        LTS_intermediate_step = []
        LTS_intermediate_step.append(LTS_D18_nan[i])
        LTS_intermediate_step.append(LTS_B17[i])
        LTS_nan.append(np.nanmax(LTS_intermediate_step,axis=0))  
    
    
    #UPSCALE TO FIT THE BIOMASS ARRAY
    multi_lat = np.int16(np.shape(biomass_nan[0,:,:])[0]/np.shape(LTS_nan[0])[0])
    multi_lon = np.int16(np.shape(biomass_nan[0,:,:])[1]/np.shape(LTS_nan[0])[1])
    
    for i in range(0,len(LTS_nan)):
        LTS_nan[i] = np.kron(LTS_nan[i],np.ones((multi_lat,multi_lon)))
    
    LTS_nan = np.array(LTS_nan,dtype=float)
    
    LTS_out.append(LTS_nan)
    
    
    
    
    return(LTS_out)





"""
LOAD SUMMER TEMPERATURE DATASETS
"""
def load_JJA_LTS(biomass_nan):
    in_path = "SPECIFY_INPUT_PATH"
    nc_D18 = f"{in_path}/D18_LST_IGBPdet.nc"  # external dataset, see readme
    nc_D18id = Dataset(nc_D18, "r")  # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class
    
    lats = nc_D18id.variables["lat"][:]
    lons = nc_D18id.variables["lon"][:]
    
    LTS_day_all = []
    LTS_night_all=[]
    
    
    for trans in range(0,45):
        for mon in range(0,12):
            LTS_day_all.append(nc_D18id.variables["Delta_LSTday"][mon,trans,:,:])
            LTS_night_all.append(nc_D18id.variables["Delta_LSTnight"][mon,trans,:,:])
            
            
    LTS_day_JJA = []        
    LTS_night_JJA = [] 
            
    for x in range(0,45):     
            LTS_day_array = np.ma.asarray(LTS_day_all)
            intermediate_day_array = LTS_day_array[x*12+5:(x+1)*12-4].mean(axis=0)
            LTS_day_JJA.append(intermediate_day_array)
            
            LTS_night_array = np.ma.asarray(LTS_night_all)
            intermediate_night_array = LTS_night_array[x*12+5:(x+1)*12-4].mean(axis=0)
            LTS_night_JJA.append(intermediate_night_array)
    
    LTS_JJA = (np.ma.asarray(LTS_day_JJA) + np.ma.asarray(LTS_night_JJA))/2   
    
    LTS_JJA_reverse = LTS_JJA*-1
    LTS_day_JJA_reverse = np.ma.asarray(LTS_day_JJA)*-1
    LTS_night_JJA_reverse = np.ma.asarray(LTS_night_JJA)*-1     
            
    
    LTS_D18_nan = []
    LTS_D18_nan.append(LTS_JJA_reverse[7,:,:]) #cro_ebf
    LTS_D18_nan.append(LTS_JJA_reverse[6,:,:]) #gra_ebf
    LTS_D18_nan.append(LTS_JJA_reverse[5,:,:]) #shr_ebf
    LTS_D18_nan.append(LTS_JJA_reverse[22,:,:])#cro_enf
    LTS_D18_nan.append(LTS_JJA_reverse[21,:,:])#gra_enf
    LTS_D18_nan.append(LTS_JJA_reverse[20,:,:]) #shr_enf
    LTS_D18_nan.append(LTS_JJA_reverse[15,:,:])#cro_dbf
    LTS_D18_nan.append(LTS_JJA_reverse[14,:,:])#gra_dbf
    LTS_D18_nan.append(LTS_JJA_reverse[13,:,:]) #shr_dbf
    LTS_D18_nan.append(LTS_JJA_reverse[28,:,:])#cro_dnf
    LTS_D18_nan.append(LTS_JJA_reverse[27,:,:])#gra_dnf
    LTS_D18_nan.append(LTS_JJA_reverse[26,:,:])#shr_dnf
    
    
    
    for z in range(0,len(LTS_D18_nan)):
        LTS_D18_nan[z] = np.ma.filled(LTS_D18_nan[z], np.nan)
    
    ###
    
    nc_f = f"{in_path}/B17_DTS.nc"  # external dataset, see readme
    nc_fid = Dataset(nc_f, "r")  # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class
    
    lats = nc_fid.variables["lat"][:]
    lons = nc_fid.variables["lon"][:]
    
    dTs_ann = []
    
    for trans in range(0, 6):
    
        dTs_ann.append(nc_fid.variables["dTs_JJA"][trans,:,:]) #Transition, lat, lon
    
    
    cro_enf = dTs_ann[0]
    cro_dbf = dTs_ann[1]
    cro_ebf = dTs_ann[2]
    gra_enf = dTs_ann[3]
    gra_dbf = dTs_ann[4]
    gra_ebf = dTs_ann[5]
    
    LTS_B17 = []
    LTS_B17.append(cro_ebf)
    LTS_B17.append(gra_ebf)
    LTS_B17.append(LTS_D18_nan[2])
    LTS_B17.append(cro_enf)
    LTS_B17.append(gra_enf)
    LTS_B17.append(LTS_D18_nan[5])
    LTS_B17.append(cro_dbf)
    LTS_B17.append(gra_dbf)
    LTS_B17.append(LTS_D18_nan[8])
    LTS_B17.append(LTS_D18_nan[9])
    LTS_B17.append(LTS_D18_nan[10])
    LTS_B17.append(LTS_D18_nan[11])
    
    for z in range(0,len(LTS_B17)):
        LTS_B17[z] = np.ma.filled(LTS_B17[z], np.nan)
    
    
    
    
    
###MEAN
    LTS_out=[]
    LTS_nan=[]
    for i in range(0,len(LTS_B17)):
        LTS_intermediate_step = []
        LTS_intermediate_step.append(LTS_D18_nan[i])
        LTS_intermediate_step.append(LTS_B17[i])
        LTS_nan.append(np.nanmean(LTS_intermediate_step,axis=0))  
    
    
    #UPSCALE TO FIT THE BIOMASS ARRAY
    multi_lat = np.int16(np.shape(biomass_nan[0,:,:])[0]/np.shape(LTS_nan[0])[0])
    multi_lon = np.int16(np.shape(biomass_nan[0,:,:])[1]/np.shape(LTS_nan[0])[1])
    
    for i in range(0,len(LTS_nan)):
        LTS_nan[i] = np.kron(LTS_nan[i],np.ones((multi_lat,multi_lon)))
    
    LTS_nan = np.array(LTS_nan,dtype=float)
    
    LTS_out.append(LTS_nan)
    
###LOW
    LTS_nan=[]
    for i in range(0,len(LTS_B17)):
        LTS_intermediate_step = []
        LTS_intermediate_step.append(LTS_D18_nan[i])
        LTS_intermediate_step.append(LTS_B17[i])
        LTS_nan.append(np.nanmin(LTS_intermediate_step,axis=0))  
    
    
    #UPSCALE TO FIT THE BIOMASS ARRAY
    multi_lat = np.int16(np.shape(biomass_nan[0,:,:])[0]/np.shape(LTS_nan[0])[0])
    multi_lon = np.int16(np.shape(biomass_nan[0,:,:])[1]/np.shape(LTS_nan[0])[1])
    
    for i in range(0,len(LTS_nan)):
        LTS_nan[i] = np.kron(LTS_nan[i],np.ones((multi_lat,multi_lon)))
    
    LTS_nan = np.array(LTS_nan,dtype=float)
    
    LTS_out.append(LTS_nan)
    
    
###HIGH
    LTS_nan=[]
    for i in range(0,len(LTS_B17)):
        LTS_intermediate_step = []
        LTS_intermediate_step.append(LTS_D18_nan[i])
        LTS_intermediate_step.append(LTS_B17[i])
        LTS_nan.append(np.nanmax(LTS_intermediate_step,axis=0))  
    
    
    #UPSCALE TO FIT THE BIOMASS ARRAY
    multi_lat = np.int16(np.shape(biomass_nan[0,:,:])[0]/np.shape(LTS_nan[0])[0])
    multi_lon = np.int16(np.shape(biomass_nan[0,:,:])[1]/np.shape(LTS_nan[0])[1])
    
    for i in range(0,len(LTS_nan)):
        LTS_nan[i] = np.kron(LTS_nan[i],np.ones((multi_lat,multi_lon)))
    
    LTS_nan = np.array(LTS_nan,dtype=float)
    
    LTS_out.append(LTS_nan)
    
    
    
    
    return(LTS_out)
    
    
"""
LOAD ANNUAL TEMPERATURE DATASETS
"""
def load_ANN_LTS(biomass_nan):
    in_path = "SPECIFY_INPUT_PATH"
    nc_D18 = f"{in_path}/D18_LST_IGBPdet.nc"  # external dataset, see readme
    nc_D18id = Dataset(nc_D18, "r")  # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class
    
    lats = nc_D18id.variables["lat"][:]
    lons = nc_D18id.variables["lon"][:]
    
    LTS_day_all = []
    LTS_night_all=[]
    
    
    for trans in range(0,45):
        for mon in range(0,12):
            LTS_day_all.append(nc_D18id.variables["Delta_LSTday"][mon,trans,:,:])
            LTS_night_all.append(nc_D18id.variables["Delta_LSTnight"][mon,trans,:,:])
            
            
    
    
    
    
    
    LTS_day_annual = []        
    LTS_night_annual = [] 
            
    for x in range(0,45):     
            LTS_day_array = np.ma.asarray(LTS_day_all)
            intermediate_day_array = LTS_day_array[x*12:(x+1)*12].mean(axis=0)
            LTS_day_annual.append(intermediate_day_array)
            
            LTS_night_array = np.ma.asarray(LTS_night_all)
            intermediate_night_array = LTS_night_array[x*12:(x+1)*12].mean(axis=0)
            LTS_night_annual.append(intermediate_night_array)
    
    LTS_annual = (np.ma.asarray(LTS_day_annual) + np.ma.asarray(LTS_night_annual))/2   
    
    LTS_annual_reverse = LTS_annual*-1    
    
    
    
    
            
    
    LTS_D18_nan = []
    LTS_D18_nan.append(LTS_annual_reverse[7,:,:]) #cro_ebf
    LTS_D18_nan.append(LTS_annual_reverse[6,:,:]) #gra_ebf
    LTS_D18_nan.append(LTS_annual_reverse[5,:,:])#shr_ebf
    LTS_D18_nan.append(LTS_annual_reverse[22,:,:])#cro_enf
    LTS_D18_nan.append(LTS_annual_reverse[21,:,:])#gra_enf
    LTS_D18_nan.append(LTS_annual_reverse[20,:,:])#shr_enf
    LTS_D18_nan.append(LTS_annual_reverse[15,:,:])#cro_dbf
    LTS_D18_nan.append(LTS_annual_reverse[14,:,:])#gra_dbf
    LTS_D18_nan.append(LTS_annual_reverse[13,:,:])#shr_dbf
    LTS_D18_nan.append(LTS_annual_reverse[28,:,:])#cro_dnf
    LTS_D18_nan.append(LTS_annual_reverse[27,:,:])#gra_dnf
    LTS_D18_nan.append(LTS_annual_reverse[26,:,:])#shr_dnf
    
    
    for z in range(0,len(LTS_D18_nan)):
        LTS_D18_nan[z] = np.ma.filled(LTS_D18_nan[z], np.nan)
    
    ###
    
    nc_f = f"{in_path}\B17_DTS.nc"  # external dataset, see readme
    nc_fid = Dataset(nc_f, "r")  # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class
    
    lats = nc_fid.variables["lat"][:]
    lons = nc_fid.variables["lon"][:]
    
    dTs_ann = []
    
    for trans in range(0, 6):
    
        dTs_ann.append(nc_fid.variables["dTs_Annual"][trans,:,:]) #Transition, lat, lon
    
    
    cro_enf = dTs_ann[0]
    cro_dbf = dTs_ann[1]
    cro_ebf = dTs_ann[2]
    gra_enf = dTs_ann[3]
    gra_dbf = dTs_ann[4]
    gra_ebf = dTs_ann[5]
    
    LTS_B17 = []
    LTS_B17.append(cro_ebf)
    LTS_B17.append(gra_ebf)
    LTS_B17.append(LTS_D18_nan[2])
    LTS_B17.append(cro_enf)
    LTS_B17.append(gra_enf)
    LTS_B17.append(LTS_D18_nan[5])
    LTS_B17.append(cro_dbf)
    LTS_B17.append(gra_dbf)
    LTS_B17.append(LTS_D18_nan[8])
    LTS_B17.append(LTS_D18_nan[9])
    LTS_B17.append(LTS_D18_nan[10])
    LTS_B17.append(LTS_D18_nan[11])
    
    for z in range(0,len(LTS_B17)):
        LTS_B17[z] = np.ma.filled(LTS_B17[z], np.nan)
    
    
    
    
    
###MEAN
    LTS_out=[]
    LTS_nan=[]
    for i in range(0,len(LTS_B17)):
        LTS_intermediate_step = []
        LTS_intermediate_step.append(LTS_D18_nan[i])
        LTS_intermediate_step.append(LTS_B17[i])
        LTS_nan.append(np.nanmean(LTS_intermediate_step,axis=0))  
    
    
    #UPSCALE TO FIT THE BIOMASS ARRAY
    multi_lat = np.int16(np.shape(biomass_nan[0,:,:])[0]/np.shape(LTS_nan[0])[0])
    multi_lon = np.int16(np.shape(biomass_nan[0,:,:])[1]/np.shape(LTS_nan[0])[1])
    
    for i in range(0,len(LTS_nan)):
        LTS_nan[i] = np.kron(LTS_nan[i],np.ones((multi_lat,multi_lon)))
    
    LTS_nan = np.array(LTS_nan,dtype=float)
    
    LTS_out.append(LTS_nan)
    
###LOW
    LTS_nan=[]
    for i in range(0,len(LTS_B17)):
        LTS_intermediate_step = []
        LTS_intermediate_step.append(LTS_D18_nan[i])
        LTS_intermediate_step.append(LTS_B17[i])
        LTS_nan.append(np.nanmin(LTS_intermediate_step,axis=0))  
    
    
    #UPSCALE TO FIT THE BIOMASS ARRAY
    multi_lat = np.int16(np.shape(biomass_nan[0,:,:])[0]/np.shape(LTS_nan[0])[0])
    multi_lon = np.int16(np.shape(biomass_nan[0,:,:])[1]/np.shape(LTS_nan[0])[1])
    
    for i in range(0,len(LTS_nan)):
        LTS_nan[i] = np.kron(LTS_nan[i],np.ones((multi_lat,multi_lon)))
    
    LTS_nan = np.array(LTS_nan,dtype=float)
    
    LTS_out.append(LTS_nan)
    
    
###HIGH
    LTS_nan=[]
    for i in range(0,len(LTS_B17)):
        LTS_intermediate_step = []
        LTS_intermediate_step.append(LTS_D18_nan[i])
        LTS_intermediate_step.append(LTS_B17[i])
        LTS_nan.append(np.nanmax(LTS_intermediate_step,axis=0))  
    
    
    #UPSCALE TO FIT THE BIOMASS ARRAY
    multi_lat = np.int16(np.shape(biomass_nan[0,:,:])[0]/np.shape(LTS_nan[0])[0])
    multi_lon = np.int16(np.shape(biomass_nan[0,:,:])[1]/np.shape(LTS_nan[0])[1])
    
    for i in range(0,len(LTS_nan)):
        LTS_nan[i] = np.kron(LTS_nan[i],np.ones((multi_lat,multi_lon)))
    
    LTS_nan = np.array(LTS_nan,dtype=float)
    
    LTS_out.append(LTS_nan)
    
    
    
    
    return(LTS_out)
    
    
    
"""
LOAD ANNUAL CMIP5 DATASET
"""
def load_ANN_cmip5():
    carbon_emitted = 2*(572-286)*2.13*10**9 / (510*10**6) #2* Additional CO2 in atm. to get emissions before sinks (GCP 2018) | *2.31 to get GtC out of ppm | per km2 |
    in_path = "SPECIFY_INPUT_PATH"
    nc_cmip = f"{in_path}/mean_ensemble.nc" # external dataset, see readme
    nc_cmipid = Dataset(nc_cmip, 'r')
    
    #January[0] to December[11]
    cmip_temp = nc_cmipid['ts'][:]
    
    cmip_temp_ann_mean = np.mean(cmip_temp,axis=0)
    
    #FLIP TO MATCH LTS AND BIOMASS DATASETS
    cmip_per_tC = np.flip(cmip_temp_ann_mean/carbon_emitted,axis=0)
    
    return(cmip_per_tC)
    
    
"""
LOAD ANNUAL CMIP5 STD DATASET
"""
def load_ANN_cmip5_std():
    carbon_emitted = 2*(572-286)*2.13*10**9 / (510*10**6) #2* Additional CO2 in atm. to get emissions before sinks (GCP 2018) | *2.31 to get GtC out of ppm | per km2 |
    in_path = "SPECIFY_INPUT_PATH"
    nc_cmip = f"{in_path}/std_ensemble.nc" # external dataset, see readme
    nc_cmipid = Dataset(nc_cmip, 'r')
    
    #January[0] to December[11]
    cmip_temp = nc_cmipid['ts'][:]
    
    cmip_temp_ann_mean = np.mean(cmip_temp,axis=0)
    
    cmip_per_tC = np.flip(cmip_temp_ann_mean/carbon_emitted,axis=0)
    
    return(cmip_per_tC)
    
    
"""
LOAD SUMMER CMIP5 DATASET
"""
def load_JJA_cmip5():
    carbon_emitted = 2*(572-286)*2.13*10**9 / (510*10**6) #2* Additional CO2 in atm. to get emissions before sinks (GCP 2018) | *2.31 to get GtC out of ppm | per km2 |
    in_path = "SPECIFY_INPUT_PATH"
    nc_cmip = f"{in_path}/mean_ensemble.nc" # external dataset, see readme
    nc_cmipid = Dataset(nc_cmip, 'r')
    
    #January[0] to December[11]
    cmip_temp = nc_cmipid['ts'][5:8] #0 = january
    
    cmip_temp_ann_mean = np.mean(cmip_temp,axis=0) 
    
    cmip_per_tC = np.flip(cmip_temp_ann_mean/carbon_emitted,axis=0)
    
    return(cmip_per_tC)
    
    
"""
LOAD SUMMER CMIP5 STD DATASET
"""
def load_JJA_cmip5_std():
    carbon_emitted = 2*(572-286)*2.13*10**9 / (510*10**6) #2* Additional CO2 in atm. to get emissions before sinks (GCP 2018) | *2.31 to get GtC out of ppm | per km2 |
    in_path = "SPECIFY_INPUT_PATH"
    nc_cmip = f"{in_path}/std_ensemble.nc" # external dataset, see readme
    nc_cmipid = Dataset(nc_cmip, 'r')
    
    #January[0] to December[11]
    cmip_temp = nc_cmipid['ts'][5:8] #0 = january
    
    cmip_temp_ann_mean = np.mean(cmip_temp,axis=0) 
    
    cmip_per_tC = np.flip(cmip_temp_ann_mean/carbon_emitted,axis=0)
    
    return(cmip_per_tC)
    
    
"""
LOAD WINTER CMIP5 DATASET
"""
def load_DJF_cmip5():
    carbon_emitted = 2*(572-286)*2.13*10**9 / (510*10**6) #2* Additional CO2 in atm. to get emissions before sinks (GCP 2018) | *2.31 to get GtC out of ppm | per km2 |
    in_path = "SPECIFY_INPUT_PATH"
    nc_cmip = f"{in_path}/mean_ensemble.nc" # external dataset, see readme
    nc_cmipid = Dataset(nc_cmip, 'r')
    
    #January[0] to December[11]
    cmip_temp = [nc_cmipid['ts'][-1],nc_cmipid['ts'][0],nc_cmipid['ts'][1]] #0 = january
    
    cmip_temp_ann_mean = np.mean(cmip_temp,axis=0) 
    
    cmip_per_tC = np.flip(cmip_temp_ann_mean/carbon_emitted,axis=0)
    
    return(cmip_per_tC)
    
    
"""
LOAD WINTER CMIP5 STD DATASET
"""
def load_DJF_cmip5_std():
    carbon_emitted = 2*(572-286)*2.13*10**9 / (510*10**6) #2* Additional CO2 in atm. to get emissions before sinks (GCP 2018) | *2.31 to get GtC out of ppm | per km2 |
    in_path = "SPECIFY_INPUT_PATH"
    nc_cmip = f"{in_path}/std_ensemble.nc" # external dataset, see readme
    nc_cmipid = Dataset(nc_cmip, 'r')
    
    #January[0] to December[11]
    cmip_temp = [nc_cmipid['ts'][-1],nc_cmipid['ts'][0],nc_cmipid['ts'][1]] #0 = january
    
    cmip_temp_ann_mean = np.mean(cmip_temp,axis=0) 
    
    cmip_per_tC = np.flip(cmip_temp_ann_mean/carbon_emitted,axis=0)
    
    return(cmip_per_tC)