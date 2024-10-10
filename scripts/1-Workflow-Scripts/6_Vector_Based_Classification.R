## ############################## 
##
## Script Name: Vector_Based_Classification
##
## Description: Script for generating vector class segments based on overlap with thematic datasets 
## Calls functions in vector_classification_segment_overlap
## for classifying Urban areas, quarries, arable, saltmarsh
##
## Author: Miles Clement,Becky Trippier, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2023-09-23
## 
## Date Last Modified: 2024-07-22
## 
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## ggplot2_3.4.4 
## osmdata_0.2.5 
## tidyr_1.3.0   
## stringr_1.5.0 
## sf_1.0-12    
## dplyr_1.1.3 
## data.table_1.14.8
##
## ############################## ##

#### Set up ####
source("./2-Functions/common_functions_library.r")
source("./2-Functions/6_vector_classification_segment_overlap.R")

# Load functions
'%!in%' = function(x,y)!('%in%'(x,y))

# Load in packages
required_packages <- c("dplyr", "sf", "stringr", "tidyr","osmdata", "ggplot2", "data.table")
install_and_load(required_packages)

## Set user variables ####
# Create list of bgzs
bgz_list <- str_pad(rep(1:13), 2, pad = "0")

## File and folder paths for input data
# File path to segmentation
seg_path <- "./LE2223_Publication/Segmentation/"

# File path to zonal statistics
zs_path <- "./LE2223_Publication/zonalStats/"

## File paths to national datasets informing specific habitats
# Quarries - OS Mastermap Topography Layer - https://www.ordnancesurvey.co.uk/products/os-mastermap-topography-layer
quarry_path <- "./LE2223_Publication/Input_Layers/OSMM/Mastermap/RawCartoBGZ.gdb"

# Arable - RPAs Crop Map of England - data.gov.uk
crome_path <- "./LE2223_Publication/Input_Layers/Vector_Habs/CROME/CROME_BGZ.gdb"

# Agricultural land classification - https://www.data.gov.uk/dataset/952421ec-da63-4569-817d-4d6399df40a1/provisional-agricultural-land-classification-alc
alc_path <- "./Ancillary_Raw/Agricultural_Land_Classification_Dec2023/data/Agricultural_Land_Classification_Provisional_England.shp"

# Saltmarsh - LE JBA output
saltmarsh_jba_path <- "./LE2223_Publication/Input_Layers/Vector_Habs/LE_Saltmarsh_JBA/"

# Saltmarsh - EA saltmarsh extent and zonation - https://www.data.gov.uk/dataset/0e9982d3-1fef-47de-9af0-4b1398330d88/saltmarsh-extent-zonation
saltmarsh_ea_path <- "./LE2223_Publication/Input_Layers/Vector_Habs/LE_Saltmarsh_EA/June2021_saltmarshExtentZonation/data/Saltmarsh_Zonation.shp"

# Saltmarsh QA outputs folder (LE created)
salt_qa_folder <- "./LE2223_Publication/Beta_Outputs/VectorOutputs/Specific_Habs/saltmarsh_manual_checking/saltmarsh_QA/"

# Solar farms - LE team digitised - see LE 2022-23 technical user guide for more details
solar <- "./LE2223_Publication/Input_Layers/Vector_Habs/SolarFarms/Solar_Farms_LE2223.shp"

# Path to QA segments sheet
qa_adjustments_final <- "./LE2223_Publication/Beta_Outputs/QA_checks/LE2223_Hab_QA_Adjusted_FINAL.xlsx"

## Folder paths to save outputs
# To save individual habitat csvs
output_path <- "./LE2223_Publication/Beta_Outputs/VectorOutputs/Specific_Habs/"

# To save classified segmentation with specific habitat mapping applied
vector_out_path <- "./LE2223_Publication/Beta_Outputs/VectorOutputs/Vector_Classed/"

#-----------------------------------------------------#
## Vector Classification Overlap ####

#### ARABLE ####
## Get arable overlap ##
for (bgz in bgz_list) {
  print(paste0("BGZ", bgz, ": Processing Arable overlap"))

  ## Load segmentation ###
  segmentation <- st_read(paste0(seg_path, "LE2223_BGZ", bgz, "_Segmentation_withIDs.shp"))
  segmentation$Area <- as.numeric(st_area(segmentation))

  # Read in crome with bgz query
  qry <- sprintf("SELECT * FROM CROME_2022_BGZ%s", bgz)
  crome <- st_read(crome_path, query = qry)

  # Remove codes relating to non-arable habitats - see LE guidance for rationale
  bad.codes <- list("AC100", "CA02", "LG14", "SR01", "FA01", "HE02", "PG01", "NA01", "WA00", "TC01", "NU01", "WO12", "AC00")
  crome <- crome[crome$LUCODE %!in% bad.codes, ]

  # Remove any urban and vector segs which get burnt in before arable
  segmentation <- segmentation %>% dplyr::filter(!Class_name %in% c("Urban", "Water"))

  # Run vector overlap
  segments_overlap <- vector_overlap(segmentation, crome, target_pct = 1)

  # Write out segs
  if (!dir.exists(paste0(output_path, "arable"))) {
    dir.create(paste0(output_path, "arable"))
  }

  fwrite(segments_overlap, paste0(output_path, "arable/BGZ", bgz, "_arable_segments.csv"), row.names = FALSE)
  print(paste(bgz, "Arable generated."))
}
## Sense check of arable against the agricultural land classification grade ####
alc <- st_read(alc_path)

# Iterate through zones to collect alc grades
for (bgz in bgz_list) {
  # Read in overlapped segs - 1% threshold
  arable <- read.csv(paste0(output_path, "arable/BGZ", bgz, "_arable_segments.csv"))

  # Read in segmentation
  seg <- st_read(paste0(seg_path, "LE2223_BGZ", bgz, "_Segmentation_withIDs.shp"))

  # Join to segments
  all <- seg %>% dplyr::left_join(arable, by = "SegID")

  # Intersect to find the alc grade
  seg_alc <- all %>% st_intersection(alc)

  # Clean up data to assign alc grade to segment geometries
  seg_alc <- seg_alc %>%
    st_drop_geometry() %>%
    dplyr::select(SegID, alc_grade)

  all_out <- all %>%
    dplyr::left_join(seg_alc, by = "SegID") %>%
    dplyr::select(SegID, Class_name = Class_name.x, Overlap_Area, Overlap_Pct, alc_grade)
 
  # Calculate km2 area and write out as csv
  arable_segs <- all_out %>%
    dplyr::mutate(area = st_area(all_out)) %>%
    dplyr::mutate(area_km2 = units::set_units(area, "km^2")) %>%
    st_drop_geometry()

  arable_segs <- arable_segs %>% dplyr::filter(!is.na(Overlap_Pct)) %>% # subset just those which overlap with arable
    dplyr::filter(alc_grade %in% c("Grade 1", "Grade 2", "Grade 3", "Grade 4")) %>% # subset to just those Grade 1-4
    dplyr::select(SegID:alc_grade)

  # Write out
  fwrite(arable_segs, paste0(output_path, "arable/BGZ", bgz, "_arable_segments.csv"), row.names = FALSE)
}

#### SALTMARSH (JBA) ####
for (bgz in bgz_list) {
  print(paste0("BGZ", bgz, ": Processing Saltmarsh overlap"))

  ## Load segmentation ###
  segmentation <- st_read(paste0(seg_path, bgz, "_Segmentation_withIDs.shp"))
  segmentation$Area <- as.numeric(st_area(segmentation))

  # Read in data
  saltmarsh <- st_read(paste0(saltmarsh_jba_path, "BGZ", bgz, "_Saltmarsh_v1.shp"))

  # Run vector overlap
  segments_overlap <- vector_overlap(segmentation,
    saltmarsh,
    target_pct = 50
  )

  # Write out segs
  if (!dir.exists(paste0(output_path, "saltmarsh_JBA"))) {
    dir.create(paste0(output_path, "saltmarsh_JBA"))
  }

  fwrite(segments_overlap, paste0(output_path, "saltmarsh_JBA/BGZ", bgz, "_saltmarsh_segments.csv"), row.names = FALSE)
}

#### SALTMARSH (EA) ####
for (bgz in bgz_list) {
  print(paste0("BGZ", bgz, ": Processing Saltmarsh overlap"))

  ## Load segmentation ###
  segmentation <- st_read(paste0(seg_path, bgz, "_Segmentation_withIDs.shp"))
  segmentation$Area <- as.numeric(st_area(segmentation))

  # Read in data
  saltmarsh <- st_read(paste0(saltmarsh_ea_path))

  # Run vector overlap
  segments_overlap <- vector_overlap(segmentation,
    saltmarsh,
    target_pct = 50
  )

  # Write out segs
  if (!dir.exists(paste0(output_path, "saltmarsh_EA"))) {
    dir.create(paste0(output_path, "saltmarsh_EA"))
  }
  fwrite(segments_overlap, paste0(output_path, "saltmarsh_EA/BGZ", bgz, "_saltmarsh_segments.csv"), row.names = FALSE)
}

#### Combining Saltmarsh post-QA step ####
## Pull out additional segments indicating new saltmarsh areas to add, these were pulled into a shapefile to examine with visual inspections ##
# Read in additional qa created segments
additions <- read.csv(paste0(saltQAfolder,'saltmarsh_additions.csv'))
out <- NULL
for(bgz in bgz_list){
  #read in segmentation
  segmentation <- st_read(paste0(seg_path,'LE2223_BGZ',bgz,"_Segmentation_withIDs.shp"))
  seg_sub <- segmentation %>% #extract the segment shapes for additions
    dplyr::filter(SegID %in% additions$SegID)
  out <- rbind(out,seg_sub)
}
saltmarsh_add <- out %>% left_join(additions,by='SegID')
#save out additions with geometries
st_write(saltmarsh_add,paste0(saltQAfolder,'saltmarsh_additions.shp'))   

## Additional segments - these were examined against S2 mosaics and APGB, where a second assessor agreed new polygons were created to encompass this larger area across indicator segments  ##
## Extract segments at manually created rough polygons ##
## Assessor 1+3
add_a1 <- st_read(paste0(saltQAfolder, "LE_saltmarsh_Beth__additions.shp"))
add_a2 <- st_read(paste0(saltQAfolder, "LE_saltmarsh_Alison__additions.shp"))

# Iterate through collating geometries
all_add_a1 <- NULL
all_add_a2 <- NULL

for (bgz in bgz_list) {
  ## load segmentation ###
  segmentation <- st_read(paste0(seg_path, "LE2223_BGZ", bgz, "_Segmentation_withIDs.shp"))
  a1_segs <- st_intersection(segmentation, add_a1) # intersect with segmentation to get segIDs
  if (nrow(a1_segs) > 0) {
    a1_seg_geom <- segmentation %>% filter(SegID %in% a1_segs$SegID) # select just intersecting SegIDs
    all_add_a1 <- rbind(all_add_a1, a1_seg_geom) # add to list
  }
  a2_segs <- st_intersection(segmentation, add_a2) # intersect with segmentation to get segIDs
  if (nrow(a2_segs) > 0) {
    a2_seg_geom <- segmentation %>% filter(SegID %in% a2_segs$SegID) # select just intersecting SegIDs
    all_add_a2 <- rbind(all_add_a2, a2_seg_geom) # add to list
  }
}

# Combine assessor data of additions to QA layer
all_add_a1 <- all_add_a1 %>% dplyr::mutate(a1 = 1, a2 = 0, a3 = 1)
all_add_a2 <- all_add_a2 %>% dplyr::mutate(a1 = 1, a2 = 1, a3 = 0)
all_additions <- all_add_a1 %>%
  rbind(all_add_a2) %>%
  dplyr::rename(Clss_nm = Class_name) %>%
  dplyr::mutate(EA_sltm = 0, JBA_slt = 0, Check = 2)

# Write out addition data
st_write(all_additions, paste0(saltQAfolder, "LE_saltmarsh_QA_all_additions.shp"), row.names = FALSE)

## Combine all data ##
## Segments from QA dataset layers
# Assessor 1
assessor1 <- st_read(paste0(saltQAfolder, "LE_saltmarsh_manualCheck_Anne.shp")) %>%
  dplyr::mutate(a1 = 1) %>%
  st_drop_geometry()

# Assessor 2
assessor2 <- st_read(paste0(saltQAfolder, "LE_saltmarsh__manualCheck_Alison.shp")) %>%
  dplyr::mutate(a2 = 1) %>%
  st_drop_geometry()

# Assessor 3
assessor3 <- st_read(paste0(saltQAfolder, "LE_saltmarsh_manualCheck_Beth.shp")) %>%
  dplyr::mutate(a3 = 1) %>%
  st_drop_geometry()

## Read in qa dataset orig
saltmarsh_qa <- st_read(saltQAfolder, "LE_saltmarsh_manualCheck.shp")

# Join data
saltmarsh_combo <- saltmarsh_qa %>%
  # join columns per assessor
  dplyr::left_join(assessor1[, c("SegID", "a1")], by = "SegID") %>%
  dplyr::left_join(assessor2[, c("SegID", "a2")], by = "SegID") %>%
  dplyr::left_join(assessor3[, c("SegID", "a3")], by = "SegID") %>%
  # replace NAs with 0s
  dplyr::mutate(across(a1:a3, ~ replace_na(.x, 0))) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(total = sum(a1, a2, a3, na.rm = T)) %>% # total across assessors
  dplyr::mutate(pcnt = (total / 3) * 100) %>% # calculate as percentage
  dplyr::select(-verdict, -notes)

# Write out assessor report
saltmarsh_combo %>% st_drop_geometry() %>% write.csv(paste0(saltQAfolder,'LE_saltmarsh_QA_results.csv'))

## Add in user added segments
all_additions <- all_additions %>%
  rowwise() %>%
  dplyr::mutate(total = sum(a1, a2, a3, na.rm = T)) %>% # total across assessors
  dplyr::mutate(pcnt = (total / 3) * 100)

# Merge
saltmarsh_combo_all <- saltmarsh_combo %>%
  rbind(all_additions)
st_write(saltmarsh_combo_all, paste0(saltQAfolder, "LE_saltmarsh_QA_all_merged.shp"), delete_layer = T)

# Write as csv
saltmarsh_conf <- saltmarsh_combo_all %>%
  st_drop_geometry() %>%
  dplyr::rename(tot = total) %>%
  dplyr::filter(tot >= 2) # only keep as saltmarsh those where 2 or more assessors agree presence
fwrite(saltmarsh_conf, paste0(saltQAfolder, "LE_saltmarsh_QA_all_merged.csv"), row.names = FALSE)

# Write out as bgz csvs for merging
saltmarsh_conf <- read.csv(paste0(saltQAfolder, "LE_saltmarsh_QA_all_merged.csv"))
saltmarsh_bgz <- saltmarsh_conf %>% mutate(zone = str_extract(str_remove(SegID, "LE2223_BGZ"), "\\d{2}"))
for (bgz in bgz_list) {
  saltmarsh_sub <- saltmarsh_bgz %>%
    dplyr::filter(zone == bgz) %>%
    dplyr::select(SegID, Clss_nm, Area, Check, tot)
  fwrite(saltmarsh_sub, paste0(output_path, "saltmarsh/BGZ", bgz, "_saltmarsh_segments.csv"), row.names = FALSE)
}

#### Solar Farms ####
for (bgz in bgz_list) {
  print(paste0("BGZ", bgz, ": Processing Solar Farm overlap"))
  
  ## Load segmentation ###
  segmentation <- st_read(paste0(seg_path, "LE2223_BGZ", bgz, "_Segmentation_withIDs.shp"))
  segmentation$Area <- as.numeric(st_area(segmentation))
  
  # Read in solar farms
  solar_farms <- st_read(solar) %>%
    st_transform(27700) %>%
    st_make_valid()
  
  # Remove any urban and vector segs which get burnt in before arable
  segmentation <- segmentation %>% dplyr::filter(!Class_name %in% c("Urban", "Water"))
  
  # Run vector overlap
  segments_overlap <- vector_overlap(segmentation, solar_farms, target_pct = 50)
  segments_overlap <- segments_overlap %>% mutate(HabCode = "SF")
  
  # Write out segs
  if (!dir.exists(paste0(output_path, "solar"))) {
    dir.create(paste0(output_path, "solar"))
  }
  fwrite(segments_overlap, paste0(output_path, "solar/BGZ", bgz, "_solar_segments.csv"), row.names = FALSE)
  print(paste(bgz, "Solar Farms generated."))
}

#### Bare ground ####
## Classifying bare ground segments based on bare soil index
# Iterate through zones
for (bgz in bgz_list) {
  # read in zonal stats
  zonalstat <- read.csv(paste0(zs_path, "BGZ", bgz, "_zonalStats.csv"))

  # get required band names
  bandNames <- tibble(fields = names(zonalstat)) %>%
    dplyr::filter(str_detect(fields, c("mean_band11|mean_band4|mean_band8|mean_band2")))

  # calculate normalised bare soil index for all seasons
  zonalIndices <- zonalstat %>%
    dplyr::select(SegID, all_of(bandNames$fields)) %>%
    # s2 NBSI: (swir(B11)+red(B4))-(VNIR(B8)+blue(B2))/(swir(B11)+red(B4)) + (VNIR(B8)+blue(B2))
    dplyr::mutate(S2_spr_mean_ndbsi = ((S2_spr_mean_band11 + S2_spr_mean_band4) - (S2_spr_mean_band8 + S2_spr_mean_band2)) / ((S2_spr_mean_band11 + S2_spr_mean_band4) + (S2_spr_mean_band8 + S2_spr_mean_band2))) %>%
    dplyr::mutate(S2_sum_mean_ndbsi = ((S2_sum_mean_band11 + S2_sum_mean_band4) - (S2_sum_mean_band8 + S2_sum_mean_band2)) / ((S2_sum_mean_band11 + S2_sum_mean_band4) + (S2_sum_mean_band8 + S2_sum_mean_band2))) %>%
    dplyr::mutate(S2_aut_mean_ndbsi = ((S2_aut_mean_band11 + S2_aut_mean_band4) - (S2_aut_mean_band8 + S2_aut_mean_band2)) / ((S2_aut_mean_band11 + S2_aut_mean_band4) + (S2_aut_mean_band8 + S2_aut_mean_band2)))

  # select indices
  ZonalSeg <- zonalIndices %>% dplyr::select(SegID, S2_spr_mean_ndvi:S2_aut_mean_ndvi)

  # write out
  fwrite(ZonalSeg, paste0(out_path, "BGZ", bgz, "_bare_indices.csv"), row.names = FALSE)
}

## Assign to segments and threshold values
for (bgz in bgz_list) {
  zonal_seg <- read.csv(paste0(out_path, "BGZ", bgz, "_bare_indices.csv"))

  # Read in segmentation
  seg <- st_read(paste0(seg_path, "LE2223_BGZ", bgz, "_Segmentation_withIDs.shp")) %>% st_drop_geometry()
  bare_seg <- seg %>%
    dplyr::left_join(zonal_seg, by = "SegID") %>% # join to segments
    dplyr::filter(Class_name == "Segment") %>% # remove any water or urban
    # filter to keep just those where NBSI is over 0.2 in summer and autumn and 0.1 in spring.
    dplyr::mutate(bare = ifelse(S2_spr_mean_ndvi >= 0.1 & S2_sum_mean_ndvi >= 0.2 & S2_aut_mean_ndvi >= 0.2, 1, 0)) %>%
    dplyr::filter(bare == 1)

  # Write out
  fwrite(bare_seg, paste0(output_path, "bareGround/BGZ", bgz, "_bareGround_segments.csv"), row.names = FALSE)
  print(paste(bgz, " done."))
}


#### Quarries  ####
for (bgz in bgz_list) {
  print(paste0("BGZ", bgz, ": Processing Quarries overlap"))

  ## load segmentation ###
  segmentation <- st_read(paste0(seg_path, "LE2223_BGZ", bgz, "_Segmentation_withIDs.shp"))
  segmentation$Area <- as.numeric(st_area(segmentation))

  # read in data using database query for quarry features
  qry <- sprintf("SELECT FID, DESCRIPTIVEGROUP, DESCRIPTIVETERM, MAKE, THEME, SHAPE
               FROM MM_Carto_BGZ%s
               WHERE DESCRIPTIVETERM LIKE 'Mineral Workings' OR DESCRIPTIVETERM LIKE 'Spoil Heap'", bgz)
  quarry <- st_read(quarry_path, query = qry)

  # run vector overlap
  segments_overlap <- vector_overlap(segmentation, quarry, target_pct = 50)
  if (!dir.exists(paste0(output_path, "quarry"))) {
    dir.create(paste0(output_path, "quarry"))
  }

  # subset urban segments labelled as quarries
  if (!dir.exists(paste0(output_path, "quarry/urban_quarry"))) {
    dir.create(paste0(output_path, "quarry/urban_quarry"))
  }
  urbQuarries <- segments_overlap %>% filter(Class_name == "Urban")
  write.csv(urbQuarries, paste0(output_path, "quarry/urban_quarry/BGZ", bgz, "_urbanquarry_segments.csv"))

  # write out for qa
  if (!dir.exists(paste0(output_path, "quarry/quarry_qa"))) {
    dir.create(paste0(output_path, "quarry/quarry_qa"))
  }

  # join back to shapes
  seg_shp <- segments_overlap %>% left_join(segmentation, by = "SegID")

  # write out
  st_write(seg_shp, paste0(output_path, "quarry/quarry_qa/BGZ", bgz, "_quarry_segments.shp"))
}

# Quarry qa step
NBSIthresh <- 0.2
for (bgz in bgz_list) {
  quarry <- st_read(paste0(output_path, "quarry/quarry_qa/BGZ", bgz, "_quarry_segments.shp"))
  zonal <- read.csv(paste0(zs_path, "BGZ", bgz, "_zonalStats.csv"))
  # Get required band names
  bandNames <- tibble(fields = names(zonal)) %>%
    dplyr::filter(str_detect(fields, c("max_band4|max_band8")))
  
  # Calculate normalised bare soil index for all seasons
  zonalIndices <- zonal %>%
    dplyr::filter(SegID %in% quarry$SegID) %>%
    dplyr::select(SegID, all_of(bandNames$fields)) %>%
    # s2 NBSI: (nir(B8)-red(B4))/nir(B8)+red(B4)))
    dplyr::mutate(S2_spr_max_nbsi = (S2_spr_max_band8 - S2_spr_max_band4) / (S2_spr_max_band8 + S2_spr_max_band4)) %>%
    dplyr::mutate(S2_sum_max_nbsi = (S2_sum_max_band8 - S2_sum_max_band4) / (S2_sum_max_band8 + S2_sum_max_band4)) %>%
    dplyr::mutate(S2_aut_max_nbsi = (S2_aut_max_band8 - S2_aut_max_band4) / (S2_aut_max_band8 + S2_aut_max_band4)) %>%
    dplyr::select(SegID, S2_spr_max_nbsi:S2_aut_max_nbsi)
  
  # Select indices
  bare_quarry <- quarry %>%
    dplyr::left_join(zonalIndices, by = "SegID") %>% # join to segments
    # filter to keep just those where NBSI is over 0.2 in summer and autumn and 0.1 in spring.
    dplyr::mutate(bare = ifelse(S2_spr_max_nbsi <= NBSIthresh & S2_sum_max_nbsi <= NBSIthresh & S2_aut_max_nbsi <= NBSIthresh, 1, 0)) %>%
    dplyr::filter(bare == 1)

  # take a look
  # st_write(bare_quarry,paste0(output_path,'quarry/quarry_qa/BGZ',bgz,"_quarryBare_segments.shp"),delete_layer = T)

  # post QA join back in to write in standard format
  quarry <- bare_quarry %>% st_drop_geometry()
  fwrite(quarry, paste0(output_path, "quarry/BGZ", bgz, "_quarry_segments.csv"), row.names = FALSE)
  print(paste(bgz, "quarries done."))
}

#### allotments ####

# get data with osm api for all of england
# set england boundary
england_bb <- getbb("England")

# Call data download with api
q <- england_bb %>%
  opq() %>%
  add_osm_feature(key = "landuse", value = "allotments") %>%
  osmdata_sf()

# Plot to check
ggplot() +
  geom_sf(data = q$osm_polygons)

# Clean cols - get selected fields, transform to BNG, make valid
allotments <- q$osm_polygons %>%
  dplyr::select(id = osm_id, name, landuse, geometry) %>%
  st_transform(27700) %>%
  st_make_valid()

# Iterate through bgzs to allocate to segments
for (bgz in bgz_list) {
  segmentation <- st_read(paste0(seg_path, "LE2223_BGZ", bgz, "_Segmentation_withIDs.shp"))
  
  # remove any urban and vector segs which get burnt in before arable
  segmentation <- segmentation %>% dplyr::filter(!Class_name %in% c("Urban", "Water"))
  
  # Run vector overlap
  segments_overlap <- vector_overlap(segmentation, allotments, target_pct = 1)
  segments_overlap <- segments_overlap %>% dplyr::mutate(HabCode = "ALL")
  
  # Write out segs
  if (!dir.exists(paste0(output_path, "allotment"))) {
    dir.create(paste0(output_path, "allotment"))
  }
  
  fwrite(segments_overlap, paste0(output_path, "allotment/BGZ", bgz, "_allotment_segments.csv"), row.names = FALSE)
  print(paste0("BGZ", bgz, " allotments generated."))
}



#### Merging all vector layers ####

## Corrections prior to merge ##
# Correct urban quarries where QA has shown these are now decommissioned
for (bgz in bgz_list) {
  # read in urban
  urb_seg <- read.csv(paste0(output_path, "urban/BGZ", bgz, "_urban_segments.csv"))

  # reclass any urban quarries back to segment
  urb_quarries <- read.csv(paste0(output_path, "quarry/urban_quarry/BGZ", bgz, "_urbanquarry_segments.csv"))
  urb_seg <- urb_seg %>% dplyr::filter(!SegID %in% urb_quarries$SegID)
  fwrite(urb_seg, paste0(output_path, "urban/BGZ", bgz, "_urban_segments.csv"), row.names = FALSE)
}

# Correct arable segs where QA has shown these are incorrect
ah_qa <- read_excel(qa_adjustments_final, sheet = "AH_remove")

# Iterate through arable segs
for (bgz in bgz_list) {
  # read in arable
  ah_seg <- read.csv(paste0(output_path, "arable/BGZ", bgz, "_arable_segments.csv"))
  # remove any AH in qa doc back to segment
  ah_seg <- ah_seg %>% filter(!SegID %in% ah_qa$SegID)
  write.csv(ah_seg, paste0(output_path, "arable/BGZ", bgz, "_arable_segments.csv"), row.names = FALSE)
}

# iterate through bgzs merging layers by hierarchy
for (bgz in bgz_list) {
  vec_classes <- tribble(
    ~layer, ~Primary_Hab,
    "water", "Water",
    "urban", "Built-up Areas and Gardens",
    "saltmarsh", "Coastal Saltmarsh",
    "allotment", "Allotment",
    "solar", "Solar Farms",
    "arable", "Arable and Horticultural",
    "bareGround", "Bare Ground",
    "quarry", "Bare Ground"
  )
  # assign vector classes in heirarchical order
  class_out <- assign_habitat(vec_classes, bgz, output_path, seg_path)
  # write out
  st_write(class_out, paste0(vector_out_path, "LE2223_BGZ", bgz, "_VectorClassed.shp"), delete_layer = T)
}
