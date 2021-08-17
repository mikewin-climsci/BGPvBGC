# Prioritizing forestation based on biogeochemical and local biogeophysical impacts, Windisch et al.

## Contact
michael.gregory.windisch@alumni.ethz.ch

## Citation
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5211680.svg)](https://doi.org/10.5281/zenodo.5211680)

## How to run
### Create output data
Always specify in and output paths first as indicated by "SPECIFY IN/OUTPUT PATH" in the code. The main output calculation is performed by `2_full_calc_bgp_v_bgc.py` which imports functions from `functions.py` and depends on output from `1_biomass_potential.py`.

### Analyse output data
Jupyter Notebook scripts for figure creation and result analysis are prepared in the folder /output_analysis. For proper function the scripts need to be executed from top to bottom as they are separated into loading datasets and figure creation / analysis.

## Necessary external datasets

### GLC2000 land cover map:
<https://cdiac.ess-dive.lbl.gov/epubs/ndp/global_carbon/carbon_documentation.html>

As [0,1] masks of:
```
nonforest: 10<=x<=19 # in code call "nonforest_mask_GLC2000_gte10lte19"
grassland: x=13 # in code call "gra_mask_GLC2000_eq_13"
cropland: 16<=x<=18 # in code call "cro_mask_GLC2000_eq_16_17_18"
shrubland: x=11;x=12;x=15 # in code call "shr_mask_GLC2000_eq_11_12_15"
```

### IPCC Tier 1 biomass:
https://cdiac.ess-dive.lbl.gov/epubs/ndp/global_carbon/carbon_documentation.html

As full dataset and forest specific biomass prepared by multiplying with GLC2000 land cover classes (numbers provided below)

```
full dataset # in code call "20181023_IPCC_Biomass"

forest specific biomass only
evergreen broadleaf: x=1 # in code call "ebf_biomass"
decideous broadleaf: x=2;x=3 # in code call "dbf_biomass"
evergreen needleleaf: x=4 # in code call "enf_biomass"
decideous needleleaf: x=5 # in code call "dnf_biomass"
other forest types: 6<=x<=9 # in code call "otf_biomass"
```

### Standard deviation of biomass, Erb et al. 2017:
https://www.nature.com/articles/nature25138.

As forest specific biomass prepared by multiplying with GLC2000 land cover classes (numbers provided below)

```
forest specific biomass only
evergreen broadleaf: x=1 # in code call "ebf_biomass_max/min_1SD"
decideous broadleaf: x=2;x=3 # in code call "dbf_biomass_max/min_1SD"
evergreen needleleaf: x=4 # in code call "enf_biomass_max/min_1SD"
decideous needleleaf: x=5 # in code call "dnf_biomass_max/min_1SD"
other forest types: 6<=x<=9 # in code call "otf_biomass_max/min_1SD"
```

### Soil organic carbon, Sanderman et al. 2017:
https://www.pnas.org/content/114/36/9575.short

```
full dataset # in code call "SOCS_0_30cm_year_NoLU_10km"
```

### BGP temperature change, Duveiller et al. 2018; Bright et al. 2017:
<https://www.nature.com/articles/s41467-017-02810-8> \
<https://www.nature.com/articles/nclimate3250>

```
Duveiller # in code call "D18_LST_IGBPdet"
Bright # in code call "B17_DTS"
```

### Transient (after doubling of CO2) surface temperature ("ts") response to cumulative emissions of +1% CO2 ESM experiments:

Please see method section of the manuscript for detailed information about the experiment setup.

```
monthly ensemble mean # in code call "mean_ensemble"
monthly ensemble standard deviation # in code call "std_ensemble"
```
