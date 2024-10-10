## ############################## 
##
##
## Script Name: Workflow for extracting data for manual QA checklist and QC statistics and data comparisons
##
## Description: Extracting various checks of pulling out coastal habitats away from a given distance by the coast, habitats on unexpected soils, and assessing comparisons with other datasets. These checks include:
## Check 1 - assessment of coastal habitats being flagged inland
## Check 2 - assessment of bracken segments on non-acidic soils
## Check 3 - assessment of bog segments on non-peaty soils
## Check 4 - committing adjustments from QA steps and additional visual checks
## Final data compiling into geopackage
## 
## Author: Becky Trippier, Rebecca Mein, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2024-02-05
## 
## Date Last Modified: 2024-08-12
## 
## Versioning:
## R version 4.2.3 (2023-03-15 ucrt)
## Dependencies:
## sf_1.0-12       lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0  
## dplyr_1.1.3     purrr_1.0.2     readr_2.1.4     tidyr_1.3.0    
## tibble_3.2.1    ggplot2_3.4.4   tidyverse_2.0.0
## ############################## ##

## User defined functions 
source('2-Functions/common_functions_library.R')

## Load dependencies ####
libaries <- c("tidyverse", "sf", "readxl", "janitor", "openxlsx")
install_and_load(libaries)

## user variables ####

## user file and folder paths
# folder path with LE habitat map outputs per bgz
hab_path <- './LE2223_Publication/Final_Outputs/FinalMerge/'
# folder path to store outputs
outputdir <- './LE2223_Publication/Final_Outputs/FinalMerge/'
# folder path to zonal stats
zs_path <- './LE2223_Publication/zonalStats/'
# folder path to write manual QA outputs to
QA_path <- './LE2223_Publication/Beta_Outputs/QA_checks/'

## filepaths for coastal check datasets ##
# filepath to mean high water springs shapefile
MHW_path <-  './Data/Ancillary_Raw/MHW/merged/MHW_valid.shp'
MHW_out <- './LE2223_Publication/Beta_Outputs/VectorOutputs/Specific_Habs/MHWclip/'
#path to qa adjustments
QAadjust <- './LE2223_Publication/Beta_Outputs/QA_checks/LE2223_Hab_QA_Adjusted.xlsx'

## filepaths for soils checks ##
#filepath to peaty soils dataset - https://www.data.gov.uk/dataset/9d494f48-f0d7-4333-96f0-8b736ac8fb18/peaty-soils-location
peatySoils <- './Data/Phase_Archive/LE2223_archive/Soils_checks/Peaty_Soils_Location_(England)___BGS_&_NSRI.shp'
#filepath to epm beta dataset
epm_beta <- './Data/Phase_Archive/LE2223_archive/Soils_checks/Beta_Vegetation_240501/'

#filepaths for final assessment and corrections
QAvec <- QAvec_path
QAmod <- QAmod_path

## user objects
#create list of bgzs
bgzList <- str_pad(rep(1:13),2,pad='0')
# set coastal threshold - 3.5km. See technical user guide for reference
coastThresh <- 3500
#dataset version name
vers = 'LE2223'

#habitat lookup table
habLookup <- tibble::tribble(
  ~HabCode,  ~Class,
  'AH','Arable and Horticultural',
  'BAG','Built-up Areas and Gardens',
  'BMYW','Broadleaved, Mixed and Yew Woodland',
  'BOG','Bog',
  'BRA','Bracken',
  'BS','Bare Sand',
  'BG','Bare Ground',
 # 'BSSP','Bare Ground',
  'CS','Coastal Saltmarsh',
  'CSD','Coastal Sand Dunes',
  'CW','Coniferous Woodland',
  'DSH','Dwarf Shrub Heath',
  'FMS','Fen, Marsh and Swamp',
 # 'IR','Bare Ground',
  'ISIG','Improved and Semi-Improved Grassland',
  'SCR','Scrub',
  'UG','Unimproved Grassland',
  'WAT','Water', 
  'SF', 'Solar Farms')


#--------------------------------#
#### LE Data Manual Checklist ####
#--------------------------------#
#### Check 1 - any coastal sand dunes, coastal saltmarsh, coastal vegetated shingle <3.5km from MHW ####

##create 3.5km line from MHW 
# create MHW polygon at buffer distance
#read in MHW line and make valid
  MHW <- st_read(MHW_path) %>% st_make_valid()
   # get inland boundary - 3.5km from MHW - crop to bgz
  MHW_3.5km <- MHW %>% st_buffer(dist = -coastThresh) 
  st_write(MHW_3.5km, paste0(QA_path, 'MHW_3_5km_coastalCheck.shp'))
  
#optional - read back in 3.5km line
#MHW_3.5km <- st_read(paste(QA_path, 'MHW_3_5km_coastalCheck.shp'))
  
#iterate through to grab segments outside on 3.5km threshold
coastSegs <- NULL
for (bgz in bgzList){
  #read in habitat map
  LEdat <- st_read(paste0(hab_path,'LE2223_BGZ',bgz,'_HabitatProbabilityMap.shp'))
  inlandCoast <- LEdat %>% 
    #filter to coastal habs
    filter(Prmry_H %in% c('Coastal Sand Dunes','Bare Sand',"Coastal Saltmarsh")) %>%
    #find those within buffer range
    st_intersection(MHW_3.5km)
  #get original geometries for inland segs
  inlandSegs <- LEdat %>% filter(SegID %in% inlandCoast$SegID)
  coastSegs <- rbind(coastSegs,inlandSegs)
}

#get primary habitat habcode
coastSegs <- coastSegs %>% left_join(habLookup, by=c('Prmry_H'='Class'))
#get secondary habitat
coastSegs_sec <- coastSegs %>% mutate(secHab = str_split_i(Mdl_Hbs,',',2))

#write out for QA
st_write(coastSegs_sec, paste0(QA_path,'LE2223_BGZall_CoastalInlandCheck.shp'))
coastSegs_sec %>% st_drop_geometry() %>% 
  write.csv(paste0(QA_path,'LE2223_BGZall_CoastalInlandCheck_secHab.csv'))

## manual QA visual checking ##

## Post-QA adjustments ##
#read in those incorrect
QAadjusted <- read_excel(QAadjust) %>% dplyr::select(-HabCode) %>% 
  left_join(habLookup, by=c('Adjusted'='HabCode')) %>% 
  dplyr::select(SegID,QAAdjust=Class,QAreas=Notes,QArel=Reliability)

#iterate through to update
for (bgz in bgzList){
  #read in habitat map
  LEdat <- st_read(paste0(hab_path,'LE2223_BGZ',bgz,'_HabitatProbabilityMap.shp'))
  LEUpdate<- LEdat %>% 
    #join to qa data
    left_join(QAadjusted,by='SegID') %>%
  #update in reasoning attribute
  mutate(SorcRsn = ifelse(!is.na(QAAdjust), QAreas, SorcRsn)) %>%
  #update reliability where qa class is incorrect
  mutate(Relblty = ifelse(!is.na(QAAdjust), QArel, Relblty)) %>%
#update where qa class is incorrect
  mutate(Prmry_H = ifelse(!is.na(QAAdjust), QAAdjust, Prmry_H)) %>%
  dplyr::select(-QAreas,-QAAdjust,-QArel)
  #overwrite data
  st_write(LEUpdate,paste0(hab_path,'LE2223_BGZ',bgz,'_HabitatProbabilityMap.shp'),delete_layer = T)
}

#-----------------------------------
#### Check 2 - any bracken segments on non-acidic soils ####

## get total number of bracken segments ##
brackenAll <- 0
for (bgz in bgzList){
  #read in habitat map
  LEdat <- st_read(paste0(hab_path,'LE2223_BGZ',bgz,'_HabitatProbabilityMap.shp'))
  bracken<- LEdat %>% 
    #filter to bracken
    dplyr::filter(Prmry_H=='Bracken') 
  brackenAll <- brackenAll+ nrow(bracken)
}
cat(brackenAll) # return total number of bracken segments

## Bracken on non-acidic soil ##
# check against cranfield soilscapes dataset, see technical report for full reference#
#create object to write to
brackenCheckSegs <- NULL
#iterate through bgzs
for (bgz in bgzList){
  #read in habitat map
  LEdat <- st_read(paste0(hab_path,'LE2223_BGZ',bgz,'_HabitatProbabilityMap.shp'))
  #read in zonal stats and select soil type
  zs <- read.csv(paste0(zs_path,'BGZ',bgz,'_zonalStats.csv'))
  #extract just the soils data
  zs <- zs %>% dplyr::select(SegID,natmapSS_mode)
  # select the soilscapes groups which are non-acidic
  zs_nonAcidic <- zs %>% dplyr::filter(natmapSS_mode %in% c(5,11,12,13,19,21,22,29))
  #subset segments 
  nonAcidSegs<- LEdat %>% 
    #filter to bracken
    dplyr::filter(Prmry_H=='Bracken') %>% 
    #filter to select those on non-acid soil
    dplyr::filter(SegID %in% zs_nonAcidic$SegID)
  #join back to soilscapes groups
  nonAcidSegs<- nonAcidSegs %>% left_join(zs_nonAcidic,by="SegID")
  #join to top object
  brackenCheckSegs <- rbind(brackenCheckSegs,nonAcidSegs)
}
#report count
cat(nrow(brackenCheckSegs))
#write out shapefile
st_write(brackenCheckSegs, paste0(QA_path,'LE2223_BGZall_BrackenNonAcidCheck.shp'))
#write out csv
brackenCheckSegs %>% st_drop_geometry() %>% write.csv(paste0(QA_path,'LE2223_BGZall_BrackenNonAcidCheck.csv'))

## Manual QA step of checkpoint ## - manually view in QGIS and remove plausible bracken segments through assessment against the LE2223 seasonal mosaics. Any which look to be obviously incorrect note on the 'LE2223_Hab_QA_Adjusted.xlsx' file

#----------------------------------------------
#### Check 3 - any bog segments on non-peat soils ####

## get total number of bog segments ##
bogAll <- 0
for (bgz in bgzList){
  #read in habitat map
  LEdat <- st_read(paste0(hab_path,'LE2223_BGZ',bgz,'_HabitatProbabilityMap.shp'))
  bogs<- LEdat %>% 
    #filter to bracken
    dplyr::filter(Prmry_H=='Bog') 
  bogAll <- bogAll+ nrow(bogs)
}
cat(bogAll)

## Bog on non-peat soil ##
# check against cranfield soilscapes dataset, see technical report for full reference#
#create object to write to
bogCheckSegs <- NULL
#iterate through bgzs
for (bgz in bgzList){
  #read in habitat map
  LEdat <- st_read(paste0(hab_path,'LE2223_BGZ',bgz,'_HabitatProbabilityMap.shp'))
  #read in zonal stats and select soil type
  zs <- read.csv(paste0(zs_path,'BGZ',bgz,'_zonalStats.csv'))
  #extract just the soils data
  zs <- zs %>%  dplyr::select(SegID,natmapSS_mode)
  # select the soilscapes groups which are not peat
  zs_nonPeat <- zs %>%  dplyr::filter(natmapSS_mode %in% c(3,4,5,6,7,8,9,10,11,12,13,15,16,18,19,20,21,22,24,25,26,29))
  #subset segments 
  nonPeatSegs<- LEdat %>% 
    #filter to bog
    dplyr::filter(Prmry_H=='Bog') %>% 
    #filter to select those on non-peat soil
    dplyr::filter(SegID %in% zs_nonPeat$SegID)
  #join back to soilscapes groups
  nonPeatSegs <- nonPeatSegs %>% left_join(zs_nonPeat,by="SegID")
  #join to top object
  bogCheckSegs <- rbind(bogCheckSegs,nonPeatSegs)
}
#report count
cat(nrow(bogCheckSegs))

## check against peaty soils dataset ##
#read in peaty soils data and check validity and crs
peatySoilsValid <- st_read(peatySoils) %>% st_make_valid() %>% st_transform(crs=27700)
#intersect remaining segments with the peaty soils layer
bogCheckTwo <- st_intersection(peatySoilsValid,bogCheckSegs)
#filter to segments not in the peaty soils layer
bogCheckSegs <- dplyr::filter(bogCheckSegs,!SegID %in% bogCheckTwos$SegID)
#report count
print(nrow(bogCheckSegs))


#write out shapefile
st_write(bogCheckSegs, paste0(QA_path,'LE2223_BGZall_BogNonPeatCheck.shp'))

# check against EPM beta outputs ##
EPM_bogOnPeat <- NULL
for (bgz in 1:13){
  #read in habitat map
  EPMdat <- st_read(paste0(epm_beta,'EPM2223_Veg_BGZ',bgz,'.shp'))
  #intersect with the peaty soils layer
  bog_EPM <- st_intersection(EPMdat,bogCheckSegsTwo)
  #add list seg ids
  bog_df <- data.frame(SegID = bog_EPM$SegID)
  #join out
 EPM_bogOnPeat <- rbind(EPM_bogOnPeat,bog_df)
}  
#filter remaining segments to remove those on epm peat map
bogCheck_nopeat <- bogCheckSegs %>% dplyr::filter(!SegID %in% EPM_bogOnPeat$SegID)
#report numbers
print(nrow(bogCheck_nopeat))

st_write(bogCheck_nopeat,paste0(QA_path,paste0(QA_path,'LE2223_BGZall_bogCheckFinal.shp')))

## Manual QA step of checkpoint ## - manually view in QGIS and remove plausible bracken segments through assessment against the LE2223 seasonal mosaics. Any which look to be obviously incorrect note on the 'LE2223_Hab_QA_Adjusted.xlsx' file

#----------------------
#### Check 4 - beta - final draft - Visual checks and adjustments all habs ####
#read in model qa adjustments
QAvec <- read_excel(QAvec_path, sheet='vector_add') %>% select(SegID, Notes, Reliability, Adjusted)
QAmod <- read_excel(QAmod_path, sheet='model_add') 

# convert adjusted habcodes to habitat class names
QAmod <- QAmod %>%
  dplyr::select(SegID, Notes, Reliability, Adjusted) %>% 
  left_join(habLookup, by=c('Adjusted'='HabCode')) %>% dplyr::select(-Adjusted)

# convert adjusted habcodes to habitat class names
QAvec <- QAvec %>%
  dplyr::select(SegID, Notes, Reliability, Adjusted) %>% 
  left_join(habLookup, by=c('Adjusted'='HabCode')) %>% dplyr::select(-Adjusted)

#iterate through to update attributes to QA assessed result
for (bgz in bgzList){
  #read in habitat map
  LEdat <- st_read(paste0(hab_path,'LE2223_BGZ',bgz,'_HabitatProbabilityMap.shp'))
  # update model layer changes
  LEUpdate<- LEdat %>% 
    #join to qa data
    left_join(QAmod,by='SegID') %>%
    #update source to LE adjusted
    mutate(Source = ifelse(!is.na(Class), 'LE QA Adjusted', Source)) %>%
    #update in reasoning attribute
    mutate(SorcRsn = ifelse(!is.na(Class), Notes, SorcRsn)) %>%
    #update reliability where qa class is incorrect
    mutate(Relblty = ifelse(!is.na(Class), Reliability, Relblty)) %>%
    #update where qa class is incorrect
    mutate(Prmry_H = ifelse(!is.na(Class), Class, Prmry_H)) %>%
    dplyr::select(-Notes,-Reliability,-Class)
  
  # update vector layer changes
  LEUpdate<- LEUpdate %>% 
    #join to qa data
    left_join(QAvec,by='SegID') %>%
    #update source to LE adjusted
    mutate(Source = ifelse(!is.na(Class), 'LE QA Adjusted', Source)) %>%
    #update in reasoning attribute
    mutate(SorcRsn = ifelse(!is.na(Class), Notes, SorcRsn)) %>%
    #update reliability where qa class is incorrect
    mutate(Relblty = ifelse(!is.na(Class), Reliability, Relblty)) %>%
    #remove modelled habs and probs
    mutate(Mdl_Hbs = ifelse(!is.na(Class), NA, Mdl_Hbs)) %>%
    #remove modelled habs and probs
    mutate(Mdl_Prb = ifelse(!is.na(Class), NA, Mdl_Prb)) %>%
    #update where qa class is incorrect
    mutate(Prmry_H = ifelse(!is.na(Class), Class, Prmry_H)) %>%
    dplyr::select(-Notes,-Reliability,-Class)
  
  # last check no duplicates
  LEUpdate %>% get_dupes(SegID)
  #overwrite data file
  st_write(LEUpdate,paste0(outputdir,'LE2223_BGZ',bgz,'_HabitatProbabilityMap.shp'),delete_layer = T)
}

#----------------------
## Final version of data compiled into geopackage ####
for (bgz in bgzList){
  bgzSeg<- st_read(paste0(outputdir,vers,'_BGZ',bgz,'_HabitatProbabilityMap.shp'))
  st_write(bgzSeg,
           paste0(outputdir,vers,"_HabitatProbabilityMap.gpkg"),
           layer = paste0(vers,"_BGZ",bgz),
           driver="GPKG",
           append=TRUE)
}

#view geopackage
sf::st_layers(paste0(outputdir,vers,"_HabitatProbabilityMap.gpkg"))

