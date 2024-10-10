## ############################## 
##
## Script Name: Habitat Data QA
##
## Description: Script is intended to remove any habitat records that are no 
## longer representative of the habitat assigned to it.
## It has two main steps:
## a) Filters habitat records to remove duplicate segment/habitat combinations
## b) For each habitat, filters remaining records using the spectral reflectance
## of recently acquired, high-confidence national ground data.
## The function will create a folder structure to store outputs in in the specified output_folder directory.
## 
## Author: Miles Clement,Becky Trippier, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2023-10-04
## 
## Date Last Modified: 2024-08-12
## 
## Versioning:
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## data.table 1.14.8
## dplyr 1.1.2
## lwgeom 0.2-11
## sf 1.0-13
## tidyr_1.3.0 
## janitor_2.2.0 
## purrr_1.0.2 
## stringr_1.5.0
## ############################## ##

# Load libraries
source('2-Functions/common_functions_library.R')
source('2-Functions/1_habitat_data_QA_checks.R')
libaries <- c("data.table", "dplyr", "lwgeom", "sf", "janitor", "tidyr", "purrr")
install_and_load(libaries)

## Set User Variables ####
# Set date hab data ingested
data_date <- "20240424"
# Ground data anomolies to remove 
outliers <- 'Ground_Data_QA_outliers_20240424.csv'

# date to filter from - last 3 years
filter_date <- "2021-01-01"
vers <- "LE2223"
# filepath to habitat ground data
hab_survey_dat <- paste0('./', vers, '_Publication/Ground_data/', data_date, '_', vers, '_QA/', data_date, "_GroundSurveys_data.csv")

# folder path to zonal stats
zonal_stats <- './LE2223_Publication/zonalStats/'

# folder path to store outputs
output_folder <- './LE2223_Publication/Ground_data/20240521_LE2223_QA/'

# vector habitats to split out from final csv
vec_hab_list <- c("AH", "BAG", "WAT", "CS", "BSSP")

#high quality survey list
quality_surveys <- c("NES", "EPM", "LTM", "LEC")

# list of variables to create filter on
bands_df <- tribble(~var, ~stat, ~bands,
  "S2_spr", c("min", "mean", "max"), c(2:9, 11, 12),
  "S2_sum", c("min", "mean", "max"), c(2:9, 11, 12),
  "S2_aut", c("min", "mean", "max"), c(2:9, 11, 12),
  "S1_spr_a", c("min", "mean", "max"), 1:2,
  "S1_spr_d", c("min", "mean", "max"), 1:2,
  "S1_sum_a", c("min", "mean", "max"), 1:2,
  "S1_sum_d", c("min", "mean", "max"), 1:2,
  "S1_aut_a", c("min", "mean", "max"), 1:2,
  "S1_aut_d", c("min", "mean", "max"), 1:2
)

# create list of band names
band_list <- NULL
for (i in 1:nrow(bands_df)) {
  row <- bands_df[i, ]
  for (j in unlist(row$stat)) {
    for (k in unlist(row$bands)) {
      varname <- paste0(row$var, "_", j, "_band", k)
      band_list <- c(band_list, varname)
    }
  }
}

#-------------------------
## Run QA analysis ##

# Read in csv outputted from habitat data database
hab_input <- read.csv(hab_survey_dat)

# Recode Hab names for easier to read export filenames & graphics
hab_input <- hab_input %>%
  dplyr::mutate(HabCode = dplyr::recode(LE_Broad_Habitat,
    "Arable and horticultural" = "AH",
    "Bare sand" = "BS",
    "Bare soil/silt/peat" = "BSSP",
    "Bog" = "BOG",
    "Bracken" = "BRA",
    "Broadleaved, mixed and yew woodland" = "BMYW",
    "Built-up areas and gardens" = "BAG",
    "Coastal saltmarsh" = "CS",
    "Coastal sand dunes" = "CSD",
    "Coniferous woodland" = "CW",
    "Dwarf shrub heath" = "DSH",
    "Fen, marsh and swamp" = "FMS",
    "Improved and semi-improved grassland" = "ISIG",
    "Inland rock" = "IR",
    "Scrub" = "SCR",
    "Unimproved grassland" = "UG",
    "Water" = "WAT"
  )) %>%
  # select only required attributes
  dplyr::select(SegID, bgz, HabCode, Survey_Date, Survey_Code, easting, northing)

#remove ground data anomalies
point_rm <- read.csv(outliers) %>%
  dplyr::mutate(flag = 1) %>%
  dplyr::select(-easting, -northing)

# join in lookup and remove anomolies
hab_input_edit <- hab_input %>%
  dplyr::left_join(point_rm, by = c("SegID", "HabCode", "Survey_Date", "Survey_Code")) %>%
  dplyr::filter(is.na(flag)) %>%
  dplyr::select(-flag)

# Run QA filter function
habQAfilter(hab_Input=hab_input_edit, # input data object
            output_folder, # folder to save the outputs
            data_date, # ingestion date
            zonal_stats,  # folder path to zonal statistics
            quality_surveys, # list of quality surveys to filter on
            filter_date, # date to filter representative filter data to
            band_list, #list of bands to filter with
            vec_hab_list # list of vector habitats to split output into RF habs and vector habs
            )
