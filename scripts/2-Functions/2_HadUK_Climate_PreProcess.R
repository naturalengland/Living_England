## ############################## 
##
## Script Name: HadUK_Climate_PreProcess_Functions
##
## Description: Functions for averaging gridded climate data from NETCDF inputs
## Includes 2 year and 20 year annual averages, and 2 year seasonal averages
## 
## Author: Miles Clement, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2023-08-10
## 
## Date Last Modified: 2024-04-22
## 
## Versioning:
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## stringr_1.5.0
## terra_1.7-29
## 
## ############################## ##

#' Averaging annual gridded climate data
#'
#' @param climVar climate variable name
#' @param folderPath folder path to climate variable data
#'
#' @return
#' @export
#'
#' @examples
climateAnnual <- function(climVar, folderPath){
  
  #load libraries
  library(stringr)
  library(terra)
  
  # List pre-downloaded climate .nc files
  fileList <- list.files(paste0(folderPath,climVar,"/Annual/"),full.names = TRUE)
  
  # Subset final two files for 2 year dataset
  files2Yr <- terra::rast(tail(fileList,2))
  files20Yr <- terra::rast(fileList)
  
  # Generate year range string for file output names
  numFiles <- length(fileList)
  outname_2Yr <- paste0(str_sub(fileList[numFiles-1],-9,-6),"_",str_sub(fileList[numFiles],-9,-6))
  outname_20Yr <- paste0(str_sub(fileList[1],-9,-6),"_",str_sub(fileList[numFiles],-9,-6))
  
  # Average data
  ave2Yr <- mean(files2Yr)
  ave20Yr <- mean(files20Yr)
  
  #write out rasters
  writeRaster(ave2Yr,paste0("Annual/",climVar,"/",climVar,"_Average_",outname_2Yr,".tif"))
  writeRaster(ave20Yr,paste0("Annual/",climVar,"/",climVar,"_Average_",outname_20Yr,".tif"))
  
}


#' Averaging annual gridded climate data
#'
#' @param climVar climate variable name
#' @param folderPath folder path to climate variable data
#'
#' @return
#' @export
#'
#' @examples
climateSeasonal <- function(climVar,folderPath){
  
  #load libraries
  library(stringr)
  library(terra)
  
  # List pre-downloaded climate .nc files
  fileList <- list.files(paste0(folderPath,climVar,"/Seasonal/"),full.names = TRUE)
  
  # Read in data (2 year average, so two files loaded)
  fileDate1 <- terra::rast(fileList[1])
  fileDate2 <- terra::rast(fileList[2])
  
  # Average data for each season/band in netcdf
  seasonDJF <- mean(fileDate1[[1]],fileDate2[[1]])
  seasonMAM <- mean(fileDate1[[2]],fileDate2[[2]])
  seasonJJA <- mean(fileDate1[[3]],fileDate2[[3]])
  seasonSON <- mean(fileDate1[[4]],fileDate2[[4]])
  
  #set output file name
  outname_2Yr <- paste0(str_sub(fileList[1],-9,-6),"_",str_sub(fileList[2],-9,-6))
  
  #write outputs
  writeRaster(seasonDJF,paste0("Seasonal/",climVar,"/",climVar,"_Average_DJF_",outname_2Yr,".tif"))
  writeRaster(seasonMAM,paste0("Seasonal/",climVar,"/",climVar,"_Average_MAM_",outname_2Yr,".tif"))
  writeRaster(seasonJJA,paste0("Seasonal/",climVar,"/",climVar,"_Average_JJA_",outname_2Yr,".tif"))
  writeRaster(seasonSON,paste0("Seasonal/",climVar,"/",climVar,"_Average_SON_",outname_2Yr,".tif"))
  
}

