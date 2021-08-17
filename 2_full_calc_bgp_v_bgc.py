#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Michael Windisch
"""

from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import sys
sys.path.insert(0, 'SPECIFY_functions_PATH')
from functions import *
import time
import numpy as np
from osgeo import gdal


start_time = time.time()

"""
PATH
"""

outpath = "SPECIFY_OUTPUT_PATH"
inpath = "SPECIFY_INPUT_PATH"


"""
LOAD LANDCOVER DATASET
"""
add_zeros = np.zeros((24,4320))

rasterload = gdal.Open(f"{inpath}/gra_mask_GLC2000_eq_13.tif") # external dataset, see readme
GRA_mask = np.concatenate((add_zeros,np.array(rasterload.GetRasterBand(1).ReadAsArray()),add_zeros),axis=0)
GRA_mask[GRA_mask>0]=1
GRA_mask[GRA_mask<=0]=np.nan


rasterload = None

rasterload = gdal.Open(f"{inpath}/cro_mask_GLC2000_eq_16_17_18.tif") # external dataset, see readme
CRO_mask = np.concatenate((add_zeros,np.array(rasterload.GetRasterBand(1).ReadAsArray()),add_zeros),axis=0)
CRO_mask[CRO_mask>0]=1
CRO_mask[CRO_mask<=0]=np.nan

rasterload = None

add_zeros = np.zeros((24,4320))

rasterload = gdal.Open(f"{inpath}/shr_mask_GLC2000_eq_11_12_15.tif") # external dataset, see readme
SHR_mask = np.concatenate((add_zeros,np.array(rasterload.GetRasterBand(1).ReadAsArray()),add_zeros),axis=0)
SHR_mask[SHR_mask>0]=1
SHR_mask[SHR_mask<=0]=np.nan


rasterload = None



apf_emission = 1.0
apf_uptake = 1.0

"""
SOC Sanderman 2017
"""
rasterload = gdal.Open(f"{inpath}/SOCS_0_30cm_year_NoLU_10km.tif") # external dataset, see readme
soc_tcha = np.array(rasterload.GetRasterBand(1).ReadAsArray(),dtype='d')
rasterload = None

soc_tcha[soc_tcha==-32767]=np.nan

soc_mean_loss_gra = 0.180
soc_mean_loss_cro = 0.266

soc_SD_loss_gra = 0.247
soc_SD_loss_cro = 0.287

#*100 to get /km2 out of /ha
defor_mean_soc = soc_tcha*soc_mean_loss_cro*100
defor_high_soc = soc_tcha*(soc_mean_loss_cro+soc_SD_loss_cro)*100
defor_low_soc = soc_tcha*(soc_mean_loss_cro-soc_SD_loss_cro)*100   

refor_mean_gra_soc = ~np.isnan(GRA_mask)*soc_tcha*soc_mean_loss_gra*100
refor_high_gra_soc = ~np.isnan(GRA_mask)*soc_tcha*(soc_mean_loss_gra+soc_SD_loss_gra)*100
refor_low_gra_soc = ~np.isnan(GRA_mask)*soc_tcha*(soc_mean_loss_gra-soc_SD_loss_gra)*100

refor_mean_cro_soc = ~np.isnan(CRO_mask)*soc_tcha*soc_mean_loss_cro*100 
refor_high_cro_soc = ~np.isnan(CRO_mask)*soc_tcha*(soc_mean_loss_cro+soc_SD_loss_cro)*100
refor_low_cro_soc = ~np.isnan(CRO_mask)*soc_tcha*(soc_mean_loss_cro-soc_SD_loss_cro)*100


"""
BIOMASS
"""
##NEED TO CREATE A LIST OF 3 ESTIMATES FOR EACH DATASET FOR THE ANALYSIS LOOP
defor_biomass_list=[]
defor_biomass_list.append(load_defor_carbon()[:] + defor_mean_soc)
defor_biomass_list.append(load_defor_carbon_min()[:] + defor_low_soc)
defor_biomass_list.append(load_defor_carbon_max()[:] + defor_high_soc)

refor_biomass_list=[]
refor_biomass_list.append(load_refor_carbon()[:] + refor_mean_cro_soc + refor_mean_gra_soc)
refor_biomass_list.append(load_refor_carbon_min()[:] + refor_low_cro_soc + refor_low_gra_soc)
refor_biomass_list.append(load_refor_carbon_max()[:] + refor_high_cro_soc + refor_high_gra_soc)



"""
ANNUAL
"""
#NEED TO CREATE A LIST OF 3 ESTIMATES FOR EACH DATASET FOR THE ANALYSIS LOOP

CMIP_list=[]
CMIP_list.append(load_ANN_cmip5())
CMIP_list.append(load_ANN_cmip5()-load_ANN_cmip5_std())
CMIP_list.append(load_ANN_cmip5()+load_ANN_cmip5_std())

# #CMIP_list FOR GLOBAL TCR
# carbon_emitted = 2*(572-286)*2.13*10**9 / (510*10**6)
# CMIP_list=[]
# CMIP_list.append(np.ones((2160,4320))*1.7928/carbon_emitted)
# CMIP_list.append(np.ones((2160,4320))*(1.7928-0.9382)/carbon_emitted)
# CMIP_list.append(np.ones((2160,4320))*(1.7928+0.9382)/carbon_emitted)


temperature_list = load_ANN_LTS(defor_biomass_list[0])


# """
# SUMMER
# """
##NEED TO CREATE A LIST OF 3 ESTIMATES FOR EACH DATASET FOR THE ANALYSIS LOOP


# CMIP_list=[]
# CMIP_list.append(load_JJA_cmip5())
# CMIP_list.append(load_JJA_cmip5()-load_JJA_cmip5_std())
# CMIP_list.append(load_JJA_cmip5()+load_JJA_cmip5_std())

# #CMIP_list FOR GLOBAL JJA TCR
# carbon_emitted = 2*(572-286)*2.13*10**9 / (510*10**6)
# CMIP_list=[]
# CMIP_list.append(np.ones((2160,4320))*1.7/carbon_emitted)
# CMIP_list.append(np.ones((2160,4320))*(1.7-0.8)/carbon_emitted)
# CMIP_list.append(np.ones((2160,4320))*(1.7+0.8)/carbon_emitted)


# temperature_list = load_JJA_LTS(defor_biomass_list[0])



# """
# WINTER
# """
##NEED TO CREATE A LIST OF 3 ESTIMATES FOR EACH DATASET FOR THE ANALYSIS LOOP

# CMIP_list=[]
# CMIP_list.append(load_DJF_cmip5())
# CMIP_list.append(load_DJF_cmip5()-load_DJF_cmip5_std())
# CMIP_list.append(load_DJF_cmip5()+load_DJF_cmip5_std())

# #CMIP_list FOR GLOBAL DJF TCR
# carbon_emitted = 2*(572-286)*2.13*10**9 / (510*10**6)
# CMIP_list=[]
# CMIP_list.append(np.ones((2160,4320))*1.87/carbon_emitted)
# CMIP_list.append(np.ones((2160,4320))*(1.87-0.9)/carbon_emitted)
# CMIP_list.append(np.ones((2160,4320))*(1.87+0.9)/carbon_emitted)


# temperature_list = load_DJF_LTS(defor_biomass_list[0])






"""
UNCERTAINTY ANALYSIS
""" 

defor_range_analysis = []
refor_range_analysis = []

defor_dT_analysis = []
refor_dT_analysis = []

defor_dC_analysis = []
refor_dC_analysis = []

for a in range(0,3):
    for b in range(0,3):
        for c in range(0,3):
            intermediate = None; dT = None; dC = None 
            intermediate,index,dT,dC = scenario(1,1,defor_biomass_list[a],(temperature_list[b]/CMIP_list[c]*-1),np.ones(np.shape(CRO_mask)),np.ones(np.shape(CRO_mask)),np.ones(np.shape(CRO_mask)))
            defor_range_analysis.append(intermediate)
            defor_dT_analysis.append(dT)
            defor_dC_analysis.append(dC)
            elapsed_time = time.time() - start_time
            print(a,b,c,elapsed_time/60)
            
            
for a in range(0,3):
    for b in range(0,3):
        for c in range(0,3):
            intermediate = None; dT = None; dC = None
            intermediate,index,dT,dC = scenario(1,1,refor_biomass_list[a],(temperature_list[b]/CMIP_list[c]*-1),CRO_mask,GRA_mask,SHR_mask)
            refor_range_analysis.append(intermediate)
            refor_dT_analysis.append(dT)
            refor_dC_analysis.append(dC)
            elapsed_time = time.time() - start_time
            print(a,b,c,elapsed_time/60)
            
elapsed_time = time.time() - start_time
print('before median',elapsed_time/60)

###MEDIAN
defor_range_mean = np.nanmedian(defor_range_analysis,axis=0)
refor_range_mean = np.nanmedian(refor_range_analysis,axis=0)

defor_dT_mean = np.nanmedian(defor_dT_analysis,axis=0)
refor_dT_mean = np.nanmedian(refor_dT_analysis,axis=0)

defor_dC_mean = np.nanmedian(defor_dC_analysis,axis=0)
refor_dC_mean = np.nanmedian(refor_dC_analysis,axis=0)

np.save(f"{outpath}/full_defor_median.npy",defor_range_mean)
np.save(f"{outpath}/full_refor_median.npy",refor_range_mean)

np.save(f"{outpath}/dT_defor_median.npy",defor_dT_mean)
np.save(f"{outpath}/dT_refor_median.npy",refor_dT_mean)

np.save(f"{outpath}/dC_defor_median.npy",defor_dC_mean)
np.save(f"{outpath}/dC_refor_median.npy",refor_dC_mean)

###STD
defor_range_std = np.nanstd(defor_range_analysis,axis=0)
refor_range_std = np.nanstd(refor_range_analysis,axis=0)

defor_dT_std = np.nanstd(defor_dT_analysis,axis=0)
refor_dT_std = np.nanstd(refor_dT_analysis,axis=0)

defor_dC_std = np.nanstd(defor_dC_analysis,axis=0)
refor_dC_std = np.nanstd(refor_dC_analysis,axis=0)

elapsed_time = time.time() - start_time
print('after median',elapsed_time/60)

np.save(f"{outpath}/full_defor_STD.npy",defor_range_std)
np.save(f"{outpath}/full_refor_STD.npy",refor_range_std)

np.save(f"{outpath}/dT_defor_STD.npy",defor_dT_std)
np.save(f"{outpath}/dT_refor_STD.npy",refor_dT_std)

np.save(f"{outpath}/dC_defor_STD.npy",defor_dC_std)
np.save(f"{outpath}/dC_refor_STD.npy",refor_dC_std)

