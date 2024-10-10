## Living England data ingestion workflow
##
## Description: Function for cleaning LE ground collected data pre-2024 app updates.
##
## Authors: Becky Trippier, Natural England
##
## Licence: MIT Licence
## Date created: 2024-10-02
##
## Date Last Modified: 2024-10-10
##
## Versioning:
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## readODS_2.1.0   DBI_1.1.3       purrr_1.0.2     readxl_1.4.2
## arrow_13.0.0    stringr_1.5.0   tidyr_1.3.0     lubridate_1.9.2
## janitor_2.2.0   dplyr_1.1.3     sf_1.0-12       openxlsx_4.2.5.2
## rnrfa_2.1.0     
#---------------------------------------------------

# run prep function for extracted LE ground data (pre-update to LE app 26-07-2024)
#' Title
#'
#' @param dat data table
#' @param x field with x coordinates
#' @param y field with y coordinates
#' @param source_crs epsg of data source coordinate reference system
#' @param target_crs epsg of data target coordinate reference system
#' @param date_field field with capture date 
#'
#' @return
#' @export
#'
#' @examples
le_prep_pre2024 <- function(dat,x_cord,y_cord,source_crs=4326,target_crs=27700, date_field, seg_path, out_folder){
  
  #### Set up ####
  ## load required functions
  source("./2-Functions/common_functions_library.r")
  source("./2-Functions/1_connect_to_access_database.R")
  source("./2-Functions/1_find_translation_path.R")
  source("./2-Functions/1_translate_survey_data.R")
  source("./2-Functions/1_assign_segments_to_survey_data.R")
  source("./2-Functions/1_clean_coord_check.R")
  
  # load in packages
  required_packages <- c("sf", "dplyr", "janitor", "lubridate","tidyr", "stringr", "arrow", "readxl", "purrr", "DBI", "dplyr", "readODS", "rnrfa")
  install_and_load(required_packages)
  
  #create folder structure for storing outputs
  if(!dir.exists(paste0(out_folder,'2_surveyData_translated'))){
    dir.create(paste0(out_folder,'2_surveyData_translated'))
  }
  if(!dir.exists(paste0(out_folder,'3_surveyData_segment_assigned'))){
    dir.create(paste0(out_folder,'3_surveyData_segment_assigned'))
  }
  if(!dir.exists(paste0(out_folder,'4_priorityHabitats'))){
    dir.create(paste0(out_folder,'4_priorityHabitats'))
  }
  ## 1. Data cleaning ##
  ## Check 1 - Coordinates
  # rename coord cols
  dat <- dat %>% dplyr::rename(x = x_cord,y = y_cord)
  
  # check if any NAs in coords
  if(any(is.na(dat$x)| is.na(dat$y))){
    start_length <- nrow(dat)
    dat <- dat %>% 
      dplyr::filter(!is.na(x)) %>%
      dplyr::filter(!is.na(y))
    end_length <- nrow(dat)
    #report to user
    print(paste(start_length-end_length, "record removed with invalid coordinates."))
  }
  
  #convert to spatial object
  lec_sp <- dat %>% st_as_sf(coords = c("x", "y"), crs = paste0("EPSG:",source_crs))
  #run through clean coord checks for crs, valid geometries and points/polygons provided
  lec_bng <- lec_sp %>% clean_coords(crs = target_crs)

  ## Check 2 - survey date
  lec_date_clean <- lec_bng %>% 
    dplyr::mutate(surveyDate = as.Date(get(date_field)))
  if(any(is.na(lec_date_clean$surveyDate))){
      start_length <- nrow(lec_date_clean)
      lec_date_clean<- lec_date_clean %>% 
        dplyr::filter(!is.na(surveyDate)) #remove any missing capture date
      start_length <- nrow(lec_date_clean)
      #report to user
      print(paste(start_length-end_length, "records removed with invalid capture dates."))
  }

  ## Check 3 - habitats
  # Compile the highest broad and primary habitats and covers
lec_clean <- lec_date_clean %>%
  dplyr::rename(
    primaryHab = "Primary LE Habitat", 
    primaryCov = "Primary LE Habitat % Cover",
    primarySub = "Primary Sub-class (inc. UKBAP Priority Habitat - PH)",
    secHab = "Secondary LE Habitat", 
    secCov = "Secondary LE Habitat % Cover",
    secSub = "Secondary Sub-class (inc. UKBAP Priority Habitat - PH)",
  ) %>%
  # remove any classed as multiple or NA
  dplyr::filter(primaryHab != 999) %>%
  ## compile highest class and %column between primary and secondary classes
  # replace  % cover with 0 where NAs
  dplyr::mutate(primaryCov = ifelse(is.na(primaryCov), 0, primaryCov)) %>%
  dplyr::mutate(secCov = ifelse(is.na(secCov), 0, secCov)) %>%
  # find highest cover between primary and secondary and replace primary habitat where this is the case
  dplyr::mutate(highestCov = ifelse((secCov > primaryCov & !is.na(secHab)), secCov, primaryCov)) %>%
  dplyr::mutate(highestHab = ifelse((secCov > primaryCov & !is.na(secHab)), secHab, primaryHab)) %>%
  dplyr::mutate(highestSub = ifelse(highestHab == primaryHab, primarySub, secSub)) %>%
  dplyr::mutate(habNotes = ifelse(is.na(Notes), paste("%Cov:", highestCov), paste0("%Cov: ", highestCov, ". ", Notes)))

## Clean up the data based on % cover
# Remove any below 40% cover
lec_clean <- lec_clean %>% dplyr::filter(highestCov >= 40) %>%
  # remove where 40% highest but primary and secondary don't match
  dplyr::filter(!(highestCov == 40 & primaryHab != secHab)) %>%
  # remove where 60% in both main and second cov but primary and secondary don't match
  dplyr::filter(!(primaryCov >= 60 & secCov >= 60 & primaryHab != secHab)) %>%
  # remove any remaining multiple classes
  dplyr::filter(highestHab != 999)

## Clean priority names ##
# Get priority hab names
# Connect to the habitat translations database
con <- connect_to_access_dbi(db_path)

# Query database
lebapp_lookup <- dbGetQuery(con, "SELECT * from Habitat_class WHERE Scheme_ID=3 ")

# Disconnect from database
dbDisconnect(con)

# Split string
lec_pstring <- lec_clean %>%
  dplyr::select(OBJECTID, highestSub) %>%
  st_drop_geometry()
lec_pstring[, 3:8] <- str_split_fixed(lec_pstring$highestSub, ",", 6)

# Bring into long format
lec_pstring <- lec_pstring %>%
  pivot_longer(V1:V6, values_to = "priority") %>%
  dplyr::select(OBJECTID, priority) %>%
  dplyr::filter(priority %in% lebapp_lookup$hab_code) %>%
  unique()

# Translation function wont handle multiple rows for OBJECTID
lec_pstring_dupes <- lec_pstring %>%
  get_dupes(OBJECTID) %>%
  dplyr::select(-dupe_count)

# Join back to broad data
lec_all <- lec_clean %>%
  dplyr::left_join(lec_pstring, by = "OBJECTID") %>%
  dplyr::select(OBJECTID, surveyDate, highestHab, priority, habNotes)

## 2. Habitat translations ##
# Find translation pathway - returns a list of available pathways to translate
sourceName <- "LE_alias"
pathways <- translation_path(sourceName, targetname, db_path = db_path)

# Find pathway for detailed hab
detTarget <- "LE_UKBAPpriority"
detpathways <- translation_path(sourceName = "LE_UKBAPpriority", detTarget, db_path)

# translate survey data
lec_cleaned <- translate_survey(surveyDat = lec_all,
                                Survey_Code = "LEC",
                                Survey_Source = "Living England Field Maps",
                                Source_ID = "OBJECTID",
                                Source_Date = "2023-10-17",
                                Survey_Date = "surveyDate",
                                Source_Broad_Habitat = "highestHab",
                                Source_Detailed_Habitat = "priority",
                                Source_Hab_Notes = "habNotes",
                                hab_codes = T,
                                sourceName = sourceName,
                                targetname = targetname,
                                pathway = pathways[[1]],
                                db_path = db_path,
                                detTranslate = T, 
                                detSourceName = "LE_UKBAPpriority",
                                detTarget = "LE_UKBAPpriority", 
                                detPath = detpathways[[1]])

# Select required fields
lec_broad <- lec_cleaned$broad %>% dplyr::select(Survey_ID, Source_ID, Survey_Code, Survey_Source, Source_Date, Survey_Date, Source_Scheme_ID, Source_Habitat_ID, Source_Broad_Name, Source_Detailed_Habitat, Source_Hab_Notes, LE_Broad_Habitat_ID = Target_Habitat_ID, LE_Broad_Habitat = Target_Broad_Name, LE_Habitat_Code = Target_Habitat_Code)

# Remove duplicates for source_id
lec_broad_clean <- lec_broad %>% dplyr::distinct(Source_ID, .keep_all = TRUE)
# Save cleaned survey
st_write(lec_broad_clean, paste0(out_folder, "2_surveyData_translated/LEC_20231110.shp"), delete_layer = T)

## 3. Segment association ##
lec_segs <- survey_seg_process(surveyDat = lec_broad_clean, seg_path = seg_path, bgz_path = bgz_path)
# Save cleaned survey
write.csv(lec_segs, paste0(out_folder, "3_surveyData_segment_assigned/LEC_20231110.csv"))

## 4. Collate priority habs ##
# Select required fields
lec_priority <- lec_priority %>%
  dplyr::filter(!is.na(Target_Broad_Name)) %>%
  dplyr::select(Survey_ID, LE_Priority_Habitat_ID = Target_Habitat_ID, LE_Priority_Habitat = Target_Broad_Name)
# Join to broad data
write.csv(lec_priority, paste0(out_folder, "4_priorityHabitats/LEC_20231110.csv"))
}
