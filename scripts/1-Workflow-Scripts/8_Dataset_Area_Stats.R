
## ##############################
##
##
## Script Name: Workflow for extracting Area Statistics for habitat predictions and reliabilities
##
## Description: Script to automate extract of summary statistics of habitat predictions and reliability scores. This produces an excel workbox populated with: raw summary data, segment counts, area totals and percentages by habitat, bgz, and reliability class. This is calculated for both the full dataset extent (OS England boundary) and up to the terrestrial boundary (Mean High Water Springs OS boundary) with the exclusion of coastal saltmarsh which is mapped within its full extent in the intertidal zone.
## further statistics have been added to draw out analysis from the mixd_hab, model_hab and model_prob attributes, and to create individual habitat datasets, to make reliability maps per habitat.  
##
## Author: Becky Trippier, Natural England
##
## Licence: MIT Licence
## Date Created: 2024-08-05
##
## Date Last Modified: 2024-08-12
##
## Versioning:
## R version 4.2.3 (2023-03-15 ucrt)
## Dependencies:
## sf_1.0-12       lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0
## dplyr_1.1.3     purrr_1.0.2     readr_2.1.4     tidyr_1.3.0
## tibble_3.2.1    ggplot2_3.4.4   tidyverse_2.0.0 openxlsx_4.2.5.2
## readxl_1.4.2
## ############################## ##
#----------------------

#### User variables ####
# Load libraries
source("2-Functions/common_functions_library.R")
libaries <- c("sf", "dplyr", "janitor", "tidyr", "purrr", "readxl", "openxlsx")
install_and_load(libaries)

# Path to final dataset
hab_path <- "./LE2223_Publication/Final_Outputs/FinalMerge/"

# path to MHWS list of segments
MHW_out <- "./LE2223_Publication/Beta_Outputs/VectorOutputs/Specific_Habs/MHWclip/"

# path to raw OS boundary folder
OS_folder <- "./Data/Ancillary_Raw/MHW/"

# path to raw model predictions
mod_path <- "D:/Projects/Living_England/Data/LE2223_Publication/Final_Outputs/ModelOutputs/model_score/"

# path to save outputs
out_path <- "./LE2223_Publication/Final_Outputs/FinalMerge/stats/"

## User objects
# Create list of bgzs
bgzList <- str_pad(rep(1:13), 2, pad = "0")

# habitat lookup table
habLookup <- tribble(
  ~HabCode, ~Class,
  "AH", "Arable and Horticultural",
  "BAG", "Built-up Areas and Gardens",
  "BMYW", "Broadleaved, Mixed and Yew Woodland",
  "BOG", "Bog",
  "BRA", "Bracken",
  "BS", "Bare Sand",
  "BG", "Bare Ground",
  "CS", "Coastal Saltmarsh",
  "CSD", "Coastal Sand Dunes",
  "CW", "Coniferous Woodland",
  "DSH", "Dwarf Shrub Heath",
  "FMS", "Fen, Marsh and Swamp",
  "ISIG", "Improved and Semi-Improved Grassland",
  "SCR", "Scrub",
  "UG", "Unimproved Grassland",
  "WAT", "Water",
  "SF", "Solar Farms"
)

#------------------------#
#### Area stats  for whole dataset extent (OS England boundary) ####
## Create summary dataset ##
# Set up summary object
all_summ <- NULL

# Iterate through to update summary object
for (bgz in bgzList) {
  # Read in habitat map
  le_dat <- st_read(paste0(hab_path, "LE2223_BGZ", bgz, "_HabitatProbabilityMap.shp"))

  # Drop geometries and unneccesary fields
  le_hab_summ <- le_dat %>%
    st_drop_geometry() %>%
    dplyr::select(SegID, Prmry_H, Relblty, Shap_Ar) %>%
    # Get habitat stats
    dplyr::group_by(Prmry_H, Relblty) %>%
    # Total segment count and area
    dplyr::summarise(totalSeg = n(), Area_m2 = sum(Shap_Ar)) %>%
    # Area in km
    dplyr::mutate(Area_km2 = Area_m2 / 1000000, zone = bgz)

  # add to object
  all_summ <- rbind(all_summ, le_hab_summ)
}
# Write out to excel
out <- createWorkbook()

# Ceate sheet
addWorksheet(out, "raw")
writeData(out, sheet = "raw", x = all_summ)

## Summary of habitats by bgz ##
hab_summ_bgz <- all_summ %>%
  # Get habitat stats by zone
  dplyr::group_by(Prmry_H, zone) %>%
  # Total segment count and area
  dplyr::summarise(
    HabtotalSeg = sum(totalSeg),
    HabArea_km2 = sum(Area_km2)
  )

# Segment area as a percentage
hab_summ_bgz <- hab_summ_bgz %>%
  dplyr::mutate(habArea_perc = (HabArea_km2 / sum(hab_summ_bgz$HabArea_km2)) * 100)

# Write out
addWorksheet(out, "habitat_bgz")
writeData(out, sheet = "habitat_bgz", x = hab_summ_bgz)

## Summary of habitats for all zones ##
hab_summ_all <- all_summ %>%
  # Get habitat stats by zone
  dplyr::group_by(Prmry_H) %>%
  # Total segment count and area
  dplyr::summarise(
    HabtotalSeg = sum(totalSeg),
    HabArea_km2 = sum(Area_km2)
  )

# Segment area as a percentage
hab_summ_all <- hab_summ_all %>%
  dplyr::mutate(HabArea_perc = (HabArea_km2 / sum(hab_summ_all$HabArea_km2)) * 100)

# Write out
addWorksheet(out, "habitat_all")
writeData(out, sheet = "habitat_all", x = hab_summ_all)

# Summary of reliability by bgz
rel_summ_bgz <- all_summ %>%
  # get habitat stats by zone
  dplyr::group_by(Relblty, zone) %>%
  # total segment count and area
  dplyr::summarise(
    ReltotalSeg = sum(totalSeg),
    RelArea_km2 = sum(Area_km2)
  )

# Segment area as a percentage
rel_summ_bgz <- rel_summ_bgz %>%
  dplyr::mutate(habArea_perc = (RelArea_km2 / sum(rel_summ_bgz$RelArea_km2)) * 100)
# write out
addWorksheet(out, "reliability_bgz")
writeData(out, sheet = "reliability_bgz", x = rel_summ_bgz)

## Summary of reliability for all zones ##
rel_summ_all <- all_summ %>%
  # Get habitat stats by zone
  dplyr::group_by(Relblty) %>%
  # Total segment count and area
  dplyr::summarise(
    ReltotalSeg = sum(totalSeg),
    RelArea_km2 = sum(Area_km2)
  )

# Segment area as a percentage
rel_summ_all <- rel_summ_all %>%
  dplyr::mutate(RelArea_perc = (RelArea_km2 / sum(rel_summ_all$RelArea_km2)) * 100)

# Write out
addWorksheet(out, "reliability_all")
writeData(out, sheet = "reliability_all", x = rel_summ_all)

# Summary of reliability by habitat
rel_hab_summ_all <- all_summ %>% 
  #get habitat stats by zone
  group_by(Relblty, Prmry_H) %>% 
  #total segment count and area
  summarise(ReltotalSeg = sum(totalSeg), 
            RelArea_km2 = sum(Area_km2))
# segment area as a percentage
rel_hab_summ_all<- rel_hab_summ_all %>% 
  mutate(RelArea_perc = (RelArea_km2/sum(rel_hab_summ_all$RelArea_km2))*100) 
#write out
addWorksheet(out, "reliability_hab")
writeData(out, sheet = "reliability_hab", x = rel_hab_summ_all)

#write out workbook
saveWorkbook(out, paste0(out_path,'LE2223_Habitat_Reliability_summary.xlsx'))


##--------------------------##
#### Area stats for terrestrial extent (clipped to MHWS except Coastal Saltmarsh) ####

## Create summary dataset ##
# Set up summary object
all_summ <- NULL

# iterate through zones
for (bgz in bgzList) {
  # Read in habitat map
  le_dat <- st_read(paste0(hab_path, "LE2223_BGZ", bgz, "_HabitatProbabilityMap.shp"))

  # Read in terrestrial segments up to MHWS incl. full coastal saltmarsh extent
  mhw_segs <- read.csv(paste0(MHW_out, "LE2223_BGZ_", bgz, "_withinMHW.csv"))

  # Clip to MHW with exception of saltmarsh
  le_dat <- le_dat %>%
    # Excluding saltmarsh, get seg ids within MHW
    dplyr::filter(SegID %in% mhw_segs$SegID | Prmry_H %in% c("Coastal Saltmarsh"))

  # wrangle data
  le_hab_summ <- le_dat %>%
    st_drop_geometry() %>%
    dplyr::select(SegID, Prmry_H, Relblty, Shap_Ar) %>%
    # get habitat stats
    dplyr::group_by(Prmry_H, Relblty) %>%
    # total segment count and area
    dplyr::summarise(totalSeg = n(), Area_m2 = sum(Shap_Ar)) %>%
    # area in km
    dplyr::mutate(Area_km2 = Area_m2 / 1000000, zone = bgz)

  # Bind together and loop
  all_summ <- rbind(all_summ, le_hab_summ)
}
# Create new excel workbook
out <- createWorkbook()
# Create sheet for raw data
addWorksheet(out, "raw")
# Write out to raw summery data
writeData(out, sheet = "raw", x = all_summ)

## Summary of habitats by bgz ##
hab_summ_bgz <- all_summ %>%
  # Get habitat stats by zone
  dplyr::group_by(Prmry_H, zone) %>%
  # total segment count and area
  dplyr::summarise(
    HabtotalSeg = sum(totalSeg),
    HabArea_km2 = sum(Area_km2)
  )

# Segment area as a percentage
hab_summ_bgz <- hab_summ_bgz %>%
  dplyr::mutate(habArea_perc = (HabArea_km2 / sum(hab_summ_bgz$HabArea_km2)) * 100)

# Write out
addWorksheet(out, "habitat_bgz")
writeData(out, sheet = "habitat_bgz", x = hab_summ_bgz)

## Summary of habitats for all zones ##
hab_summ_all <- all_summ %>%
  # get habitat stats by zone
  dplyr::group_by(Prmry_H) %>%
  # total segment count and area
  dplyr::summarise(
    HabtotalSeg = sum(totalSeg),
    HabArea_km2 = sum(Area_km2)
  )

# Segment area as a percentage
hab_summ_all <- hab_summ_all %>%
  dplyr::mutate(HabArea_perc = (HabArea_km2 / sum(hab_summ_all$HabArea_km2)) * 100)

# write out
addWorksheet(out, "habitat_all")
writeData(out, sheet = "habitat_all", x = hab_summ_all)

## Summary of reliability by bgz ##
rel_summ_bgz <- all_summ %>%
  # get habitat stats by zone
  dplyr::group_by(Relblty, zone) %>%
  # total segment count and area
  dplyr::summarise(
    ReltotalSeg = sum(totalSeg),
    RelArea_km2 = sum(Area_km2)
  )

# Segment area as a percentage
rel_summ_bgz <- rel_summ_bgz %>%
  dplyr::mutate(habArea_perc = (RelArea_km2 / sum(rel_summ_bgz$RelArea_km2)) * 100)

# Write out
addWorksheet(out, "reliability_bgz")
writeData(out, sheet = "reliability_bgz", x = rel_summ_bgz)

## Summary of reliability for all zones ##
rel_summ_all <- all_summ %>%
  # get habitat stats by zone
  dplyr::group_by(Relblty) %>%
  # total segment count and area
  dplyr::summarise(
    ReltotalSeg = sum(totalSeg),
    RelArea_km2 = sum(Area_km2)
  )

# Segment area as a percentage
rel_summ_all <- rel_summ_all %>%
  dplyr::mutate(RelArea_perc = (RelArea_km2 / sum(Rerel_summ_alllSummall$RelArea_km2)) * 100)

# Write out
addWorksheet(out, "reliability_all")
writeData(out, sheet = "reliability_all", x = rel_summ_all)

## Summary of reliability by habitat ##
rel_hab_summ_all <- all_summ %>%
  # get habitat stats by zone
  dplyr::group_by(Relblty, Prmry_H) %>%
  # total segment count and area
  dplyr::summarise(
    ReltotalSeg = sum(totalSeg),
    RelArea_km2 = sum(Area_km2)
  )

# Segment area as a percentage
rel_hab_summ_all <- rel_hab_summ_all %>%
  dplyr::mutate(RelArea_perc = (RelArea_km2 / sum(rel_hab_summ_all$RelArea_km2)) * 100)

# Write out
addWorksheet(out, "reliability_hab")
writeData(out, sheet = "reliability_hab", x = rel_hab_summ_all)

saveWorkbook(out, paste0(out_path,'LE2223_Habitat_Reliability_MHW.xlsx'))

##--------------------------##
#### raw vector layer area stats for boundaries - will differ from LE due to pixelisation with 10m S2 pixels of segmentation ####

## OS England regional boundary
# Read in OS Eng boundary polygon
eng_os <- st_read(paste0(OS_folder, "english_region_region.shp"))

# Merge into one object
eng_os_all <- eng_os %>% st_union()

# Get area size (m2)
st_area(eng_os_all)

## Mean High Water Springs (MHWS) OS boundary
# Read in OS MHWS boundary polygon
mhw <- st_read(paste0(OS_folder, "MHW_valid.shp"))
# get area size (m2)
st_area(mhw)

##--------------------------##
#### Area stats for mixed segments ####

# Iterate through BGZs
all_summ <- NA
for (bgz in bgzList) {
  # Read in habitat map
  le_dat <- st_read(paste0(hab_path, "LE2223_BGZ", bgz, "_HabitatProbabilityMap.shp"))
  
  # Remove geometries and unneeded fields
  le_hab_summ <- le_dat %>%
    st_drop_geometry() %>%
    dplyr::select(SegID, Prmry_H, Mixd_Sg, Shap_Ar) %>%
    # get stats per habitat and mixed seg status
    dplyr::group_by(Prmry_H, Mixd_Sg) %>%
    # total segment count and area
    dplyr::summarise(totalSeg = n(), Area_m2 = sum(Shap_Ar)) %>%
    # area in km
    dplyr::mutate(Area_km2 = Area_m2 / 1000000, zone = bgz)
  all_summ <- rbind(all_summ, le_hab_summ)
}

# Summary across all zones by habitat and mixed status
hab_summ_all <- all_summ %>%
  # get habitat stats by mixed status
  dplyr::group_by(Prmry_H, Mixd_Sg) %>%
  # total segment count and area
  dplyr::summarise(
    HabtotalSeg = sum(totalSeg),
    HabArea_km2 = sum(Area_km2)
  )

# Calculate percentage of area totals
hab_summ_all <- hab_summ_all %>%
  dplyr::mutate(habArea_perc = (HabArea_km2 / sum(hab_summ_bgz$HabArea_km2)) * 100)

# Summary by zones by habitat and mixed status
hab_summ_bgz <- all_summ %>%
  # Get habitat stats by mixed status and zone
  dplyr::group_by(Prmry_H, Mixd_Sg, zone) %>%
  # Total segment count and area
  dplyr::summarise(HabtotalSeg = sum(totalSeg), HabArea_km2 = sum(Area_km2))

# Calculate percentage of area totals
hab_summ_all <- hab_summ_all %>%
  dplyr::mutate(habArea_perc = (HabArea_km2 / sum(hab_summ_bgz$HabArea_km2)) * 100)

##--------------------------##
#### Stats for modelled habs and probabilities ####

# get list of vector classed segs

## iterate through BGZs and collate model predictions
# Set up object to write to
all_dat <- NULL

# Iterate through zones
for (bgz in bgzList) {
  # Read in habitat map
  model_score <- read.csv(paste0(mod_path, "bgz_", bgz, "_model_score.csv"))

  # Read in final data
  le_dat <- st_read(paste0(hab_path, "LE2223_BGZ", bgz, "_HabitatProbabilityMap.shp")) %>% st_drop_geometry()

  # Filter to just modelled habitats in final data
  mod_segs <- le_dat %>% dplyr::filter(!is.na(Mdl_Hbs))

  # Data wrangle
  m_scores <- model_score %>%
    dplyr::filter(SegID %in% mod_segs$SegID) %>%
    dplyr::select(SegID, a_pred, a_prob, b_pred, b_prob) %>%
    dplyr::mutate(zone = bgz)
  
  # Write to out object
  all_dat <- rbind(all_dat, m_scores)
  print(paste(bgz, "done"))
}

## Calculate overall statistics
# total number of segments per habitat class
all_dat %>%
  dplyr::group_by(a_pred) %>%
  dplyr::summarise(count = n())

# a_prob range per habitat class
all_dat %>%
  dplyr::group_by(a_pred) %>%
  dplyr::summarise(aprobmean = mean(a_prob), aprobmin = min(a_prob), aprobmax = max(a_prob))

# b_prob range per habitat class
all_dat %>%
  dplyr::group_by(a_pred) %>%
  dplyr::summarise(bprobmean = mean(b_prob), bprobmin = min(b_prob), bprobmax = max(b_prob))

# mean aprob_bprob difference per habitat class
all_dat %>%
  dplyr::group_by(a_pred) %>%
  dplyr::summarise(bprobmean = mean(a_prob - b_prob))

# most likely pairing per habitat class
all_dat %>%
  dplyr::group_by(a_pred, b_pred) %>%
  dplyr::summarise(count = n(), aprobmean = mean(a_prob), bprobmean = mean(b_prob))

## calculate statistics by zone
# total number of segments per habitat class
all_dat %>%
  dplyr::group_by(zone, a_pred) %>%
  dplyr::summarise(count = n())
# a_prob range per habitat class
all_dat %>%
  dplyr::group_by(zone, a_pred) %>%
  dplyr::summarise(aprobmean = mean(a_prob), aprobmin = min(a_prob), aprobmax = max(a_prob))

# b_prob range per habitat class
all_dat %>%
  dplyr::group_by(zone, a_pred) %>%
  dplyr::summarise(bprobmean = mean(b_prob), bprobmin = min(b_prob), bprobmax = max(b_prob))

# mean aprob_bprob difference per habitat class
all_dat %>%
  dplyr::group_by(zone, a_pred) %>%
  dplyr::summarise(bprobmean = mean(a_prob - b_prob))

# most likely pairing per habitat class
all_dat %>%
  dplyr::group_by(zone, a_pred, b_pred) %>%
  dplyr::summarise(count = n(), aprobmean = mean(a_prob), bprobmean = mean(b_prob))

##---------------------##
#### make individual reliability maps ####

## collate individual habitat class results
# Iterate through habitats
for (hab in habLookup$Class) {
  # Iterate through BGZs
  all_gab_dat <- NA
  for (bgz in bgzList) {
    # read in final data
    le_dat <- st_read(paste0(hab_path, "LE2223_BGZ", bgz, "_HabitatProbabilityMap.shp"))
    # filter to habitat class
    le_sub <- le_dat %>% filter(Prmry_H == hab)
    all_gab_dat <- rbind(all_gab_dat, le_sub)
  }
  # Get habcode
  classHabCode <- habLookup %>%
    dplyr::filter(Class == hab) %>%
    dplyr::select(HabCode) %>%
    as.character()
  # write out
  st_write(alldat, paste0(out_path, "LE2223_", classHabCode, ".shp"))
}
