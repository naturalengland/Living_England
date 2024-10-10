## ############################## 
##
## Script Name: HadUK_Climate_PreProcess
##
## Description: Script for averaging gridded climate data from NETCDF inputs
## Includes 2 year and 20 year annual averages, and 2 year seasonal averages
## 
## Author: Miles Clement, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2023-08-10
## 
## Date Last Modified: 2024-04-18
## 
## Versioning:
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## stringr_1.5.0
## terra_1.7-29
## 
## ############################## ##

#load functions
source("HadUK_Climate_PreProcess_Functions.R")

## set user variables ##
## Climate data from: Met Office, Hollis, D., McCarthy, M., Kendon, M., Legg, T., Simpson, I. (2018) HadUK-Grid gridded and regional average climate observations for the UK. Centre for Environmental Data Analysis, date of citation. http://catalogue.ceda.ac.uk/uuid/4dc8450d889a491ebb20e724debe2dfb. ##
#folder path to climate variable netcdf data files, should be stored in separate 'Annual' and Seasonal folders
folderPath = paste0("RawData/")
  
#--------------

#list HadUK variables
# RAINFALL = Precipitation, TAS = Temperature mean, TASMAX = Temperature max, TASMIN = Temperature min
varList <- list("RAINFALL", "TAS", "TASMAX", "TASMIN")

#iterate through variable list
for (climVar in varList){
  #report variable
  print(climVar)
  #run annual climatic variable pre-processing
  climateAnnual(climVar,folderPath)
  #run seasonal climatic variable pre-processing
  # DJF = Dec, Jan, Feb
  # MAM = Mar, Apr, May
  # JJA = Jun, Jul, Aug
  # SON = Sep, Oct, Nov
  climateSeasonal(climVar,folderPath)
  #report functions completed
  print(paste(climVar, ' done.'))
}
