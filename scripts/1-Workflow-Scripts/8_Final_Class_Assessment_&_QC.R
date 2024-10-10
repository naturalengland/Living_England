## ############################## ##
##
## Script Name: Workflow script for merging the modelled and vector classifications 
##
## Description:  This script runs through merging the vector and modelled based outputs and clipping to MHW
## 
## Author: Becky Trippier, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2023-04-14
## 
## Date Last Modified: 2024-09-25
## 
## Versioning:
## R version 4.2.3 (2023-03-15 ucrt)
## Dependencies:
## units_0.8-1     janitor_2.2.0   sf_1.0-12       lubridate_1.9.2
## forcats_1.0.0   stringr_1.5.0   dplyr_1.1.3     purrr_1.0.2    
## readr_2.1.4     tidyr_1.3.0     tibble_3.2.1    ggplot2_3.4.4  
## tidyverse_2.0.0
## ############################## ##

## User defined functions
source('2-Functions/common_functions_library.R')

## Load dependencies ####
libaries <- c("tidyverse", "sf", "janitor", "purrr", "units", "stringr", "readxl")
install_and_load(libaries)

#### User defined variables ####

## Env variables
vers <- "LE2223" # Data version
prob_thresh <- 0.1 # Threshold to report model_hab results, e.g. x => prob_thresh
mixed_thresh <- 0.3 # Threshold to find if a segment is mixed where b_prob is within the bounds of a_prob - e.g. if a_prob - mixed_thresh <= b_prob then segment is 'Probable' mixed segment

## File paths ##
# Path to model score outputs
model_path <- "./LE2223_Publication/Final_Outputs/ModelOutputs/model_score/"

# Path to vector outputs
vector_path <- "./LE2223_Publication/Final_Outputs/VectorOutputs/Vector_Classed/"

# Path to reliability outputs
rel_path <- "./LE2223_Publication/Final_Outputs/ModelOutputs/combined/"

# Ground data
rf_ground <- "./LE2223_Publication/Ground_data/20240521_LE2223_QA/7.RFinput/20240424_FieldDataQA_RFinput.csv"
vec_ground <- "./LE2223_Publication/ground_data/20240521_LE2223_QA/8.Spare/20240424_FieldDataQA_nonRF_Habs.csv"

# Saltmarsh reliability
salt_path <- "./LE2223_Publication/Final_Outputs/VectorOutputs/Specific_Habs/saltmarsh/"

# Path to mean high water shapefile
mhw_path <- "./Data/Ancillary_Raw/MHW/merged/MHW_valid.shp"

# Folder to save MHW outputs
mhw_out <- "./LE2223_Publication/Final_Outputs/VectorOutputs/Specific_Habs/MHWclip/"

# File path to save out to
output_dir <- "./LE2223_Publication/Final_Outputs/FinalMerge/"

# Create list of bgzs
bgz_list <- str_pad(rep(1:13), 2, pad = "0")
#habitat class lookup
hab_lookup <- tibble::tribble(
  ~HabCode, ~Class,
  "AH", "Arable and Horticultural",
  "BAG", "Built-up Areas and Gardens",
  "BMYW", "Broadleaved, Mixed and Yew Woodland",
  "BOG", "Bog",
  "BRA", "Bracken",
  "BS", "Bare Sand",
  "BSSP", "Bare Ground",
  "CS", "Coastal Saltmarsh",
  "CSD", "Coastal Sand Dunes",
  "CW", "Coniferous Woodland",
  "DSH", "Dwarf Shrub Heath",
  "FMS", "Fen, Marsh and Swamp",
  "IR", "Bare Ground",
  "ISIG", "Improved and Semi-Improved Grassland",
  "SCR", "Scrub",
  "UG", "Unimproved Grassland",
  "WAT", "Water",
  "BG", "Bare Ground",
  "SF", "Solar Farms",
  "HR", "Scrub"
)
#vector lookup (where set for vector layers)
vec_lookup <- tibble::tribble(
  ~Primary_Hab, ~Source, ~vecReliability,
  "Built-up Areas and Gardens", "Vector OSMM Urban","High",
  "Allotment", "Vector OpenStreetMap",NA,
  "Water", "Vector OSMM Water","Very High", 
  "Arable and Horticultural","Vector RPA CROME, ALC grades 1-4", "Medium",
  "Solar Farms", "LE QA Adjusted","High",
  "Bare Ground", "Vector LE Bare Ground Analysis","Medium",
  "Coastal Saltmarsh","Vector EA saltmarsh, LE saltmarsh & QA",NA
)

#------------------------------------------------------####
#### Generate segment list of those within MHW boundary ####
#  Load MHW polygon
MHW <- st_read(mhw_path)

# Iterate through BGZs to get list of overlapping segments
for (bgz in bgz_list) {
  # Read in segmentation
  vec_classed <- st_read(paste0(vector_path, vers, "_BGZ", bgz, "_VectorClassed.shp"))

  # Find overlapping segments
  mhw_segs <- vec_classed %>% st_intersection(MHW)

  # Get segment list
  mhw_segs <- mhw_segs %>%
    st_drop_geometry() %>%
    dplyr::select(SegID, Clss_nm)

  # Save out
  write.csv(mhw_segs, paste0(mhw_out, "LE2223_BGZ_", bgz, "_withinMHW.csv"))
}

## for speed - this was performed in qgis then brought into r to convert to desired format ##
# R code below for reference to put intersected clip from qgis into correct format
for (bgz in bgz_list) {
  # read in clipped segmentation intersected with MHW shape
  clip_shape <- st_read(paste0(mhw_out, "shape/BGZ", bgz, ".shp")) 

  # Clip to shape and keep just segID list
  segs <- clip_shape %>%
    st_drop_geometry() %>%
    dplyr::select(SegID)
  
  # Write out results
  write.csv(segs, paste0(mhw_out, "LE2223_BGZ_", bgz, "_withinMHW.csv")) # write out
}

#------------------------------------------------------####

#### merge outputs to create final dataset #### 

##iterate through bgzs
for (bgz in bgz_list) {
  # Load in rf results
  rf_classed <- read.csv(paste0(model_path, "bgz_", bgz, "_model_score.csv"))

  # Load in vector results
  vec_classed <- st_read(paste0(vector_path, vers, "_BGZ", bgz, "_VectorClassed.shp"))

  ## Reshape vector results and set source infos
  merge_classed <- vec_classed %>%
    dplyr::select(SegID, Primary_Hab = Prmry_H) %>%
    #use lookup to update data source for specific habitat mapping
    dplyr::left_join(vec_lookup[,c("Primary_Hab","Source")], by='Primary_Hab') %>% 
    #update source on remaining habitats to LE modelled
    mutate(Source = ifelse(is.na(Source), "LE Modelled",Source)) %>%
    # Change allotments back to BAG class
    dplyr::mutate(Primary_Hab = ifelse(Primary_Hab == "Allotment", "Built-up Areas and Gardens", Primary_Hab))

  ## Add in modelled results
  rf_all <- vec_classed %>%
    st_drop_geometry() %>%
    # Prep for cleaning join rf results to vector to remove any non-modelled segs
    dplyr::left_join(rf_classed, by = c("SegID" = "SegID")) %>%
    dplyr::filter(is.na(Prmry_H))

  # Calculate if segment is mixed
  rf_all <- rf_all %>%
    dplyr::mutate(Mixed_Seg = ifelse((a_prob - mixed_thresh) <= b_prob, "Probable", "Unlikely"))

  # List attributes of modelled predictions and probabilities
  rf_listed <- rf_all %>%
    # round all probabilities to 2 d.p.
    dplyr::mutate_at(vars(contains("prob")), ~ round(., digits = 2)) %>%
    # always add apred and aprob to list
    dplyr::mutate(A_acc = a_prob) %>%
    # if prob is less than prob_thresh convert to NA
    dplyr::mutate_at(vars(contains("prob")), ~ replace(., . < prob_thresh, NA))

  # If A-L prob is NA, replace relevant pred with NA
  letter_string <- c("b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l")
  for (letter in letter_string) {
    var_prob <- paste0(letter, "_prob")
    var_pred <- paste0(letter, "_pred")
    rf_listed <- rf_listed %>%
      dplyr::mutate(!!sym(var_pred) := if_else(is.na(!!sym(var_prob)), NA, !!sym(var_pred)))
  }

  # Compile list attributes
  rf_tidy <- rf_listed %>%
    # Collapse list of secondary probabilities
    tidyr::unite("Model_Probs", c("A_acc", paste0(letter_string, "_prob")), sep = ", ", na.rm = T, remove = F) %>%
    # Collapse list of secondary predictions
    tidyr::unite("Model_Habs", c("a_pred", paste0(letter_string, "_pred")), sep = ", ", na.rm = T, remove = F) %>%
    # Replace IR with BG (Bare ground)
    dplyr::mutate(Model_Habs = str_replace(Model_Habs, "IR", "BG"))

  # Select just modelled habs to join
  rf_merge <- rf_tidy %>% dplyr::select(SegID, A_pred = a_pred, Model_Habs, Model_Probs, Mixed_Seg)

  ## Merge to vector dataset
  all_classed <- merge_classed %>% dplyr::left_join(rf_merge, by = "SegID") %>%
    # join with lookup
    dplyr::left_join(hab_lookup, by = c("A_pred" = "HabCode")) %>%
    # assign primary hab from modelled data
    dplyr::mutate(Primary_Hab = ifelse(is.na(Primary_Hab), Class, Primary_Hab)) %>%
    # remove unnecessary attributes
    dplyr::select(-A_pred, -Class) %>%
    # default vector classified to 'Unlikely' mixed segment
    dplyr::mutate(Mixed_Seg = ifelse(is.na(Mixed_Seg), "Unlikely", Mixed_Seg))

  ## Add in reliability metrics for modelled habs
  reliability <- read.csv(paste0(rel_path, "bgz_", bgz, "_combined_results.csv")) %>%
    dplyr::select(SegID, Reliability = a_reliability_score) %>%
    # Filter to update just modelled habs
    dplyr::filter(SegID %in% rf_all$SegID) %>%
    # Rename to capitalise classes
    dplyr::mutate(Reliability = case_when(
      Reliability == "very high" ~ "Very High",
      Reliability == "high" ~ "High",
      Reliability == "medium" ~ "Medium",
      Reliability == "low" ~ "Low",
      Reliability == "very low" ~ "Very Low",
      .default = NA
    ))

  # Add in saltmarsh reliability metric
  salt_rel <- read.csv(paste0(salt_path, "BGZ", bgz, "_saltmarsh_segments.csv"))
  # Set saltmarsh reliability scores (see Technical User Guide for justification)
  salt_rel <- salt_rel %>%
    dplyr::mutate(Rel = case_when(
      Check == 0 & tot == 1 ~ "Low",
      Check == 0 & tot == 2 ~ "High",
      Check == 0 & tot == 3 ~ "Very High",
      Check == 1 & tot == 0 ~ "Very Low",
      Check == 1 & tot == 1 ~ "Low",
      Check == 1 & tot == 2 ~ "Medium",
      Check == 1 & tot == 3 ~ "High",
      Check == 2 ~ "Very High",
      .default = NA
    )) %>%
    dplyr::select(SegID, Rel) %>%
    dplyr::mutate(SegID = as.character(SegID))

  ## Check for duplicates ##
  reliability %>% get_dupes(SegID)
  salt_rel %>% get_dupes(SegID)

  # Join in reliability scores
  all_classed <- all_classed %>%
    dplyr::left_join(reliability, by = "SegID") %>%
    dplyr::left_join(salt_rel, by = "SegID") %>%
    dplyr::left_join(vec_lookup[,c("Primary_Hab","vecReliability")], by="Primary_Hab") 
    #only update bare ground reliability where non-modelled
  all_classed <- all_classed %>% 
    dplyr::mutate(vecReliability = ifelse(Primary_Hab == "Bare Ground" & Source != "Vector LE Bare Ground Analysis", NA,vecReliability)) %>% 
    #update saltmarsh reliability 
    dplyr::mutate(Reliability = ifelse(Primary_Hab == "Coastal Saltmarsh", Rel, Reliability)) %>% dplyr::select(-Rel) %>% 
    #update remaining vector habs
    dplyr::mutate(Reliability = ifelse(!is.na(vecReliability), vecReliability, Reliability)) %>% dplyr::select(-vecReliability)

  ## Check for duplicates ##
  all_classed %>% get_dupes(SegID)

  # Where a ground point is present update reliability
  rf_ground_data <- read.csv(rf_ground) %>% dplyr::select(SegID, HabCode)
  all_ground_data <- read.csv(vec_ground) %>%
    dplyr::select(SegID, HabCode) %>%
    rbind(rf_ground_data) %>%
    dplyr::left_join(hab_lookup, by = "HabCode") %>%
    dplyr::select(-HabCode)

  # Subset those segments with ground points
  all_classed_grnd <- all_classed %>%
    dplyr::filter(SegID %in% all_ground_data$SegID) %>%
    dplyr::left_join(all_ground_data, by = "SegID") %>%
    # where a group point is correct update socre to Very High
    dplyr::mutate(Reliability = ifelse(!is.na(Class) & Class == Primary_Hab, "Very High", Reliability)) %>%
    # where a group point is incorrect update score to Very Low
    dplyr::mutate(Reliability = ifelse(!is.na(Class) & Class != Primary_Hab, "Very Low", Reliability))

  # Separate out where multiple habs per seg in ground data
  grnd_dupes <- all_classed_grnd %>%
    janitor::get_dupes(SegID) %>%
    dplyr::arrange(SegID, Reliability)

  multiSegs <- length(unique(grnd_dupes$SegID))

  all_corr <- grnd_dupes %>%
    group_by(SegID) %>%
    dplyr::filter(Reliability == first(Reliability)) %>%
    dplyr::select(-dupe_count, -Class) %>%
    unique()
  multiSegs_post <- length(unique(all_corr$SegID))

  # Check to see if any data lost
  if (multiSegs == multiSegs_post) {
    print(paste0(multiSegs_post, " segments where multiple habitats input.All processed"))
  } else {
    stop(print("some multisegments loss, check code."))
  }

  # Combine ground data reliability results
  all_grnd_results <- all_classed_grnd %>%
    dplyr::select(-Class) %>%
    dplyr::filter(!SegID %in% grnd_dupes$SegID) %>%
    rbind(all_corr)

  # Flag those very low where ground data indicates a different hab
  seg_flag <- all_grnd_results %>% dplyr::filter(Reliability == "Very Low")
  if (!dir.exists(paste0(output_dir, "incorrectGroundPoint"))) {
    dir.create(paste0(output_dir, "incorrectGroundPoint"))
  }
  st_write(seg_flag, paste0(output_dir, "incorrectGroundPoint/", vers, "_BGZ", bgz, "_incorrectGroundSeg.shp"), delete_layer = T)

  # Replace reliability in all segment data
  all_classed_done <- all_classed %>%
    dplyr::filter(!SegID %in% all_grnd_results$SegID) %>%
    rbind(all_grnd_results)

  # Check for duplicates
  all_classed_done %>% get_dupes(SegID)

  # Add in check for mixed ground data segments
  multi_seg <- all_ground_data %>% get_dupes(SegID)

  # Make sure Mixed seg attribute is probable for these segments
  all_classed_done <- all_classed_done %>% mutate(Mixed_Seg = ifelse(SegID %in% multi_seg$SegID, "Probable", Mixed_Seg))

  all_classed <- all_classed_done %>%
    ## Add in field for qa updates
    dplyr::mutate(SourceReason = "N/A") %>%
    # Update for any adjusted with QA
    dplyr::mutate(SourceReason = ifelse(Primary_Hab == "Solar Farms", "LE digitised solar farms", SourceReason))

  # Clip to MHW with exception of saltmarsh
  mhw_segs <- read.csv(paste0(mhw_out, "LE2223_BGZ_", bgz, "_withinMHW.csv"))

  outside_mhw <- all_classed %>%
    # excluding saltmarsh, get seg ids outside of MHW
    dplyr::filter(!SegID %in% mhw_segs$SegID & !Primary_Hab %in% c("Coastal Saltmarsh", "Water"))

  # Maintain boundary but reclassify those outside mhw to water (except saltmarsh)
  all_classed_mhw <- all_classed %>% mutate(
    Primary_Hab = ifelse(SegID %in% outside_mhw$SegID, "Water", Primary_Hab),
    Source = ifelse(SegID %in% outside_mhw$SegID, "LE MHW clip", Source),
    Model_Habs = ifelse(SegID %in% outside_mhw$SegID, NA, Model_Habs),
    Model_Probs = ifelse(SegID %in% outside_mhw$SegID, NA, Model_Probs),
    Mixed_Seg = ifelse(SegID %in% outside_mhw$SegID, "Unlikely", Mixed_Seg),
    Reliability = ifelse(SegID %in% outside_mhw$SegID, "Very High", Reliability)
  )

  # Check geometries are valid
  all_classed_mhw <- all_classed_mhw %>% st_make_valid()

  # Add st_area
  all_classed_mhw$Shape_Area <- st_area(all_classed_mhw)

  # Reorder attributes
  all_classed_mhw <- all_classed_mhw %>%
    dplyr::mutate(Shape_Area = drop_units(Shape_Area)) %>%
    dplyr::select(SegID, Primary_Hab, Reliability, Model_Habs, Model_Probs, Mixed_Seg, Source, SourceReason, Shape_Area, geometry)

  # Write out dataset
  st_write(all_classed_mhw, paste0(output_dir, "/", vers, "_BGZ", bgz, "_HabitatProbabilityMap.shp"), delete_layer = T)
  
  # Clean house
  rm(all_classed_mhw, all_classed, outside_MHW, mhw_segs, seg_flag, rf_ground_data, all_ground_data, salt_rel, reliability, vec_classed, rf_classed, merge_classed, rf_tidy, rf_all, rf_listed, rf_merge)
}
