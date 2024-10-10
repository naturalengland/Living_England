## ############################## 
##
## Script Name: LE Zonal Statistics Workflow script
##
## Description: Script for generating segment-by-segment zonal stats for a list  
## of rasters. Calls functions in zonalStats_exactExtract_function.R
##
## Author: Miles Clement, Becky Tripper, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2023-06-07
## 
## Date Last Modified: 2024-07-22
## 
## Versioning: v1.3
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## lubridate_1.9.2     forcats_1.0.0      
## dplyr_1.1.3         purrr_1.0.2        
## readr_2.1.4         tidyr_1.3.0        
## tibble_3.2.1        ggplot2_3.4.4      
## tidyverse_2.0.0     terra_1.7-29       
## stringr_1.5.0       sf_1.0-12          
## furrr_0.3.1         future_1.32.0      
## exactextractr_0.9.1 data.table_1.14.8 
## tidyterra_0.5.2           
## ############################## ##

#### Set up ####
source("./2-Functions/common_functions_library.r")

# Load in packages
required_packages <- c("sf", "terra", "tibble", "tidyverse","future", "data.table")
install_and_load(required_packages)

# Functions for calculating zonal stats
source("2-Functions/4_zonalStats_exactExtract_function.R")
source("2-Functions/4_calc_indices.R")
source("2-Functions/4_cloud_removal_NAs.R")

## Set user variables ####
# Set working directory to folder containing variable datasets
setwd("./LE2223_Publication/Input_Layers/")

# Create list of bgzs
bgzList <- str_pad(rep(1:13), 2, pad = "0")

# Set LE version name
version <- "LE2223"

# Set output folder to write cloud removed images to
s2_folder <- "./LE2223_Publication/temp/S2_noCloud/"

# Folder path to segmentation
seg_folder <- "./LE2223_Publication/Segmentation/"

# Set output folder to write zonal stats to
out_folder <- "./LE2223_Publication/zonalStats_cookie/"

# note - you will also need to edit the data input tables below to set the 'path' locations to where these data files are stored on your local machine 

# -----------------------------

## Cloud removal stage ##
for (bgz in bgzList) {
  print(paste("BGZ", bgz, " running zonal stats..."))

  ## Create data input table ##
  # var - Name of variable dataset to be used in output Zonal Stats
  # path - location of dataset on local filesystem
  # stat - which statisics to be calculated (common separated within string)
  # layers - when bands within dataset to be analysed
  rasterDF <- tribble(
    ~var, ~path, ~stat, ~layers,
    # Sentinel-2
    "S2_spr", paste0("S2/Spring/s2_Spring23_BGZ", bgz, ".img"), c("mean", "min", "max", "stdev"), c(2, 3, 4, 5, 6, 7, 8, 9, 11, 12),
    "S2_sum", paste0("S2/Summer/s2_Summer23_BGZ", bgz, ".img"), c("mean", "min", "max", "stdev"), c(2, 3, 4, 5, 6, 7, 8, 9, 11, 12),
    "S2_aut", paste0("S2/Autumn/s2_Autumn22_BGZ", bgz, ".img"), c("mean", "min", "max", "stdev"), c(2, 3, 4, 5, 6, 7, 8, 9, 11, 12)
  )
  # run function
  removeCloud(rasterDF, outFolder = s2_folder)
}

## Calculate NDVI and NDWI
for (bgz in bgzList[-1]) {
  print(paste("BGZ", bgz, " running indices calculator..."))
  imgList <- c(
    paste0("S2/Spring/s2_Spring23_BGZ", bgz, ".img"),
    paste0("S2/Summer/s2_Summer23_BGZ", bgz, ".img"),
    paste0("S2/Autumn/s2_Autumn22_BGZ", bgz, ".img")
  )
  if (!dir.exists("S2_Indices")) {
    dir.create("S2_indices")
  }
  # run through imageList calculating indices layers
  indicesCalc(imgList,
    indices = c("NDVI", "NDWI"),
    outfolder = "S2_indices/",
    numCores = availableCores() - 2
  )
}

## run zonal stats ##
# iterate through all BGZs
for (bgz in bgzList){
  print(paste('BGZ',bgz, ' running zonal stats...'))
  
## Create data input table ##
  # var - Name of variable dataset to be used in output Zonal Stats
  # path - location of dataset on local filesystem
  # stat - which statisics to be calculated (common separated within string)
  # layers - when bands within dataset to be analysed
  rasterDF <- tribble(~var,~path,~stat,~layers,
                    # Sentinel-2
                    "S2_spr", paste0("S2_noCloud/s2_Spring23_BGZ",bgz,".img"), c("mean","min","max","stdev"),c(1:10),
                    "S2_sum", paste0("S2_noCloud/s2_Summer23_BGZ",bgz,".img"), c("mean","min","max","stdev"),c(1:10),
                    "S2_aut", paste0("S2_noCloud/s2_Autumn22_BGZ",bgz,".img"), c("mean","min","max","stdev"),c(1:10),
                    # Sentinel-1
                   "S1_spr_a", paste0("S1/GRD/Spring/s1_Asc_Spring23_BGZ",bgz,".tif"), c("mean","min","max","stdev"),c(1, 2),
                    "S1_spr_d", paste0("S1/GRD/Spring/s1_Desc_Spring23_BGZ",bgz,".tif"), c("mean","min","max","stdev"),c(1, 2),
                    "S1_sum_a", paste0("S1/GRD/Summer/s1_Asc_Summer23_BGZ",bgz,".tif"), c("mean","min","max","stdev"),c(1, 2),
                    "S1_sum_d", paste0("S1/GRD/Summer/s1_Desc_Summer23_BGZ",bgz,".tif"), c("mean","min","max","stdev"),c(1, 2),
                    "S1_aut_a", paste0("S1/GRD/Autumn/s1_Asc_Autumn22_BGZ",bgz,".tif"), c("mean","min","max","stdev"),c(1, 2),
                    "S1_aut_d", paste0("S1/GRD/Autumn/s1_Desc_Autumn22_BGZ",bgz,".tif"), c("mean","min","max","stdev"),c(1, 2),
                    # Sentinel-2 Indices
                     "NDVI_spr", paste0("S2_indices/NDVI/s2_Spring23_BGZ",bgz,"_NDVI.img"), c("mean","min","max","stdev"),1,
                     "NDVI_sum", paste0("S2_indices/NDVI/s2_Summer23_BGZ",bgz,"_NDVI.img"), c("mean","min","max","stdev"),1,
                     "NDVI_aut", paste0("S2_indices/NDVI/s2_Autumn22_BGZ",bgz,"_NDVI.img"), c("mean","min","max","stdev"),1,
                     "NDWI_spr", paste0("S2_indices/NDWI/s2_Spring23_BGZ",bgz,"_NDWI.img"), c("mean","min","max","stdev"),1,
                     "NDWI_sum", paste0("S2_indices/NDWI/s2_Summer23_BGZ",bgz,"_NDWI.img"), c("mean","min","max","stdev"),1,
                     "NDWI_aut", paste0("S2_indices/NDWI/s2_Autumn22_BGZ",bgz,"_NDWI.img"), c("mean","min","max","stdev"),1,
                    # Sentinel-1 Coherence
                    "coh_spr_VH_a", paste0("S1/COH/Spring/bgz",bgz,"_spr2023_A_VH_cc.tif"), c("mean","min","max","stdev"),1,
                    "coh_spr_VV_a", paste0("S1/COH/Spring/bgz",bgz,"_spr2023_A_VV_cc.tif"), c("mean","min","max","stdev"),1,
                    "coh_sum_VH_a", paste0("S1/COH/Summer/bgz",bgz,"_sum2023_A_VH_cc.tif"), c("mean","min","max","stdev"),1,
                    "coh_sum_VV_a", paste0("S1/COH/Summer/bgz",bgz,"_sum2023_A_VV_cc.tif"), c("mean","min","max","stdev"),1,
                    "coh_aut_VH_a", paste0("S1/COH/Autumn/bgz",bgz,"_aut2022_A_VH_cc.tif"), c("mean","min","max","stdev"),1,
                    "coh_aut_VV_a", paste0("S1/COH/Autumn/bgz",bgz,"_aut2022_A_VV_cc.tif"), c("mean","min","max","stdev"),1,
                    "coh_spr_VH_d", paste0("S1/COH/Spring/bgz",bgz,"_spr2023_D_VH_cc.tif"), c("mean","min","max","stdev"),1,
                    "coh_spr_VV_d", paste0("S1/COH/Spring/bgz",bgz,"_spr2023_D_VV_cc.tif"), c("mean","min","max","stdev"),1,
                    "coh_sum_VH_d", paste0("S1/COH/Summer/bgz",bgz,"_sum2023_D_VH_cc.tif"), c("mean","min","max","stdev"),1,
                    "coh_sum_VV_d", paste0("S1/COH/Summer/bgz",bgz,"_sum2023_D_VV_cc.tif"), c("mean","min","max","stdev"),1,
                    "coh_aut_VH_d", paste0("S1/COH/Autumn/bgz",bgz,"_aut2022_D_VH_cc.tif"), c("mean","min","max","stdev"),1,
                    "coh_aut_VV_d", paste0("S1/COH/Autumn/bgz",bgz,"_aut2022_D_VV_cc.tif"), c("mean","min","max","stdev"),1,
                    # Terrain
                     "dtm", paste0("Ancillary/DTM/LiDAR_10m_DTM_BGZ",bgz,".tif"), c("mean","min","max","stdev"),1,
                     "slope", paste0("Ancillary/Slope/LiDAR_10m_Slope_BGZ",bgz,".tif"), c("mean","min","max","stdev"),1,
                     "aspect", paste0("Ancillary/Aspect/LiDAR_10m_Aspect_BGZ",bgz,".tif"),c("mean","min","max","stdev"),1,
                     "twi", paste0("Ancillary/TWI/LiDAR_10m_TWI_BGZ",bgz,".tif"), c("mean","min","max","stdev"),1,
                     "chm", paste0("Ancillary/CHM/LiDAR_10m_CHM_BGZ",bgz,".img"), c("mean","min","max","stdev"),1,
                     "focal10", paste0("Ancillary/Focal/LiDAR_10m_Focal_BGZ",bgz,"_SegIn.img"), c("mean","min","max","stdev"),1,
                     "focal_1m", paste0("Ancillary/Focal/LiDAR_10m_Focal_BGZ",bgz,"_SegIn.img"), c("mean","min","max","stdev"),1,
                    # Climate - Annual
                     "rain2Y", "Ancillary/Climate/Annual/RAINFALL/RAINFALL_Average_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tas2Y", "Ancillary/Climate/Annual/TAS/TAS_Average_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasmax2Y", "Ancillary/Climate/Annual/TASMAX/TASMAX_Average_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasmin2Y", "Ancillary/Climate/Annual/TASMIN/TASMIN_Average_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "rain20Y", "Ancillary/Climate/Annual/RAINFALL/RAINFALL_Average_2003_2022.tif", c("mean","min","max","stdev"),1,
                     "tas20Y", "Ancillary/Climate/Annual/TAS/TAS_Average_2003_2022.tif", c("mean","min","max","stdev"),1,
                     "tasmax20Y", "Ancillary/Climate/Annual/TASMAX/TASMAX_Average_2003_2022.tif", c("mean","min","max","stdev"),1,
                     "tasmin20Y", "Ancillary/Climate/Annual/TASMIN/TASMIN_Average_2003_2022.tif", c("mean","min","max","stdev"),1,
                    # Climate - Seasonal
                     "rainDJF", "Ancillary/Climate/Seasonal/RAINFALL/RAINFALL_Average_DJF_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "rainMAM", "Ancillary/Climate/Seasonal/RAINFALL/RAINFALL_Average_MAM_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "rainJJA", "Ancillary/Climate/Seasonal/RAINFALL/RAINFALL_Average_JJA_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "rainSON", "Ancillary/Climate/Seasonal/RAINFALL/RAINFALL_Average_SON_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasDJF", "Ancillary/Climate/Seasonal/TAS/TAS_Average_DJF_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasMAM", "Ancillary/Climate/Seasonal/TAS/TAS_Average_MAM_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasJJA", "Ancillary/Climate/Seasonal/TAS/TAS_Average_JJA_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasSON", "Ancillary/Climate/Seasonal/TAS/TAS_Average_SON_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasmaxDJF", "Ancillary/Climate/Seasonal/TASMAX/TASMAX_Average_DJF_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasmaxMAM", "Ancillary/Climate/Seasonal/TASMAX/TASMAX_Average_MAM_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasmaxJJA", "Ancillary/Climate/Seasonal/TASMAX/TASMAX_Average_JJA_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasmaxSON", "Ancillary/Climate/Seasonal/TASMAX/TASMAX_Average_SON_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasminDJF", "Ancillary/Climate/Seasonal/TASMIN/TASMIN_Average_DJF_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasminMAM", "Ancillary/Climate/Seasonal/TASMIN/TASMIN_Average_MAM_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasminJJA", "Ancillary/Climate/Seasonal/TASMIN/TASMIN_Average_JJA_2021_2022.tif", c("mean","min","max","stdev"),1,
                     "tasminSON", "Ancillary/Climate/Seasonal/TASMIN/TASMIN_Average_SON_2021_2022.tif", c("mean","min","max","stdev"),1,
                     # Proximity
                     "coastalProx", paste0("Ancillary/Proximity/Coastal/MHW_Coastal_Proximity_BGZ",bgz,".tif"), c("mean","min","max","stdev"),1,
                     "urbanProx", paste0("Ancillary/Proximity/Urban//OS_Urban_Proximity_BGZ",bgz,".tif"), c("mean","min","max","stdev"),1,
                     "waterProx", paste0("Ancillary/Proximity/Water//OS_Water_Proximity_BGZ",bgz,".tif"), c("mean","min","max","stdev"),1,
                     # Soil/Geology
                     "spmCarb", "Ancillary/Soils/BGS_SPM/BGS_SPM_Carbonate_National.tif", "mode",1,
                     "spmDepth", "Ancillary/Soils/BGS_SPM/BGS_SPM_Depth_National.tif", "mode",1,
                     "spmESB", "Ancillary/Soils/BGS_SPM/BGS_SPM_ESB_National.tif", "mode",1,
                     "spmSize", "Ancillary/Soils/BGS_SPM/BGS_SPM_GrainSize_National.tif", "mode",1,
                     "spmGroup", "Ancillary/Soils/BGS_SPM/BGS_SPM_SoilGroup_National.tif", "mode",1,
                     "spmTex", "Ancillary/Soils/BGS_SPM/BGS_SPM_Texture_National.tif", "mode",1,
                     "natmapSS", paste0("Ancillary/Soils/NATMAP/Soilscape/NATMAP_Soilscape_BGZ",bgz,".tif"), "mode",1,
                     "natmapWet", paste0("Ancillary/Soils/NATMAP/Wetness/NATMAP_Wetness_BGZ",bgz,".tif"), "mode",1,
                     "geology", paste0("Ancillary/Geology/gb_50k_bedrock_rcsx2_bgz/gb_50k_bedrock_rcsx2_bgz",bgz,".tif"), "mode",1
                    )

  # Load in segmentation
  segmentation <- st_read(paste0(seg_folder, version, "_BGZ", bgz, "_Segmentation_withIDs.shp"))
  # Returns a table  zonal statistics for each polygon.
  seg_zonalStats <- zonalStats_exactExtract(segmentation, rasterDF, parallel = T)
  # Write the training data zonal stats to a txt file
  fwrite(seg_zonalStats, paste0(out_folder, "BGZ", bgz, "_zonalStats.csv"), row.names = FALSE)
}