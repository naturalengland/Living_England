## Living England data ingestion workflow
##
## Description: Ingestion script for cleaning and translating a habitat survey collected using the LE App. This is using the MS Access Habitat translations database. Decisions for steps taken should be noted in the data ingestion log.
##
## This will run through 3 steps for each survey:
##  1. Dataset cleaning
##  2. Translation of the habitat class
##  3. Association to a segmentation layer and extracting of an associated point (if polygons supplied)
##
## Authors: Becky Trippier, Natural England
##
## Licence: MIT Licence
## Date created: 2024-10-02
##
## Date Last Modified: 2024-10-02
##
## Versioning:
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## readODS_2.1.0   DBI_1.1.3       purrr_1.0.2     readxl_1.4.2
## arrow_13.0.0    stringr_1.5.0   tidyr_1.3.0     lubridate_1.9.2
## janitor_2.2.0   dplyr_1.1.3     sf_1.0-12       openxlsx_4.2.5.2
## rnrfa_2.1.0     
#---------------------------------------------------

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

## Set User Variables ####
# Target scheme LE_broad
targetname <- "LE_alias"
# Summary data year - year to subset in summary table, does not remove older data from the dataset
sum_year <- 2021

## File and Folder paths ####
# Target folder to save outputs
out_folder <- paste0(gd_path, "GroundSurveys_Processed/")
# Point to ground data folder location (may need to edit for your machine)
gd_path <- "./Data/Ground_Data/"
# Point to database location
db_path <- paste0(gd_path, "20240813_Habitat_Translation_DB_v3/Habitat_Translation_DB_v3.accdb")
# Point to segmentation data folder location
seg_path <- "./Data/LE2223_Publication/Segmentation/"
# Point to bgz boundaries shapefile data location
bgz_path <- "./Data/Ancillary_Raw/BGZ/LivingEngland_BioGeographicZones_Phase4.shp"
# LE data file
le_dat_extract <-paste0(gd_path, "LE_Field_Surveys/LEC/20231110_Field_Maps_2020_2023/LE_GRD_LIVE4_10112023.xlsx") 
#-------------------------
### LE Data cleaning and translation steps ####

### Old app extract format pre July-24 app update - 
### 09/10/2020 - 10/11/2023 dataset  ####
## read in dataset
lec <- read_excel(le_dat_extract) 

## look at LE scheme lookup
# connect to translation database
con <- connect_to_access_dbi(db_path)
#retrieve LE habitat broad scheme lookup
dbGetQuery(con,"SELECT * from Habitat_class WHERE Scheme_ID=1 ") %>% select(habitat_name,hab_code) %>% mutate(hab_code=as.numeric(hab_code))

# run LE data prep (function for pre July-24 app update)
le_prep_pre2024(dat=lec,
                x_cord="x2",y_cord="y2",
                date_field="CreationDate", seg_path=seg_path, out_folder=out_folder)