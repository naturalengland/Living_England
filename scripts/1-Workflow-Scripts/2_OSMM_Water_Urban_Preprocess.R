## ############################## 
##
##
## Script Name: OSMM_Water_Urban_Preprocess
##
## Description: Script to extract urban and water features from Ordnance Survey MasterMap data
## and preprocess these for input into Living England segmentation process
##
## NOTE: Due to the size of national Mastermap, data is split into BGZ subsets using ArcPro first
## 
## Author: Miles Clement, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2023-06-27
## 
## Date Last Modified: 2024-04-18
## 
## Versioning:
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## data.table_1.14.8
## dplyr_2.3.2
## sf_1.0-13
## smoothr_1.0.1
## stringr_1.5.0
## terra_1.7-29
## tidyterra_0.4.0
## 
## ############################## ##

#### Set up ####
# load in packages
library(data.table)
library(dplyr)
library(sf)
library(smoothr)
library(stringr)
library(terra)
library(tidyterra)

# load in functions
source("OSMM_Water_Urban_Preprocess_Functions.R")

## set user variables ##
bgzList <- str_pad(rep(1:13),2,pad='0')

# path to s2 template folder - created in 2_S2_Mosaic_creation.ipynb
templatePath <- './LE2223_Publication/Input_Layers/S2/Spring/'

##Ordnance Survey MasterMap Topography product - https://www.ordnancesurvey.co.uk/products/os-mastermap-topography-layer ##
# path to OSMM topography geodatabase file 
OSpath <- "./LE2223_Publication/Input_Layers/OSMM/MasterMap/OSMM_Raw.gdb"
#------------------

# Iterate through BGZs
for (bgz in bgzList){
  #set time report
  startTime <- Sys.time()
  print(paste0("Processing BGZ",bgz,"....."))
  # call OSMM processing function
  mm_preProcess(OSpath,bgz,refImage = paste0(templatePath,"/S2_Spring23_BGZ",bgz,".img"))
  #report timing
  runtime <- difftime(Sys.time(), startTime, units='mins')
  print(paste0("BGZ ",bgz," Runtime: ", format(round(runtime,2),nsmall=2)))
}
