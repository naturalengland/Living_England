## ############################## ##
##
## Script Name: Cloud removal function 
##
## Description: Function for removing cloud before the zonal stats function - converts values in GEE exported files where cloud is present from 0s to NAs, where a pixel has 0 across all bands.  
##
## Author: Becky Trippier, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2024-02-16
## 
## Date Last Modified: 2024-04-23
## 
## Versioning: v1.3
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## sf_1.0-13
## stringr_1.5.0
## terra_1.7-29
## tidyverse_2.0.0
## tidyterra_0.5.2
##
## ############################## ##

#' function for removing cloud from GEE exported mosaics
#'
#' @param rasterDF a data frame containing var, path, stat, layers
#' @param outFolder where to store outputs
#'
#' @return
#' @export
#'
#' @examples
removeCloud <- function(rasterDF,outFolder){
  
  # edit options for memory issues
  terraOptions(memfrac=0.8,todisk =T)
  
  #load packages
  library(sf)
  library(stringr)
  library(terra)
  library(tidyverse)
  library(tidyterra)
  
  # Abort programme function
  exit <- function(){invokeRestart("abort") }   
  
  # Check that list of raster/band combinations exist prior to processing using terra
  # Break loop if error
  for (r in 1:nrow(rasterDF)){
    rastRow <- rasterDF[r, ]
    errorCheck <- tryCatch (rast(rastRow$path, lyrs=rastRow$layers), error=function(e) e)
    if (inherits(errorCheck, "error")){
      print(paste0("Invalid Raster/Band combination: ", rastRow$path,": Band ", rastRow$layers))
      print("Process aborting...")
      exit()
    } 
  }
  rm(errorCheck,rastRow)
  print("Data I/O checks complete: no issues")
  
  #iterate through rasters 
  for (r in 1:nrow(rasterDF)){
    #select first variable
    rastRow <- rasterDF[r, ]
    print(paste('Processing', rastRow$var,'...'))
    # Separate out filepath, function and layers/bands 
    file <- rastRow$path
    stats <- unlist(rastRow$stat)
    bands <- rastRow$layers
    # Load data (specific bands only)
    rasterData <- rast(file, lyrs = bands)
    # Identify the S2 bands with 0s that should be NA and convert these to NAs
    # covert to dataframe
    rast_tbl <- tidyterra::as_tibble(rasterData, xy = TRUE)
    gc()     
    #pull out cell xy
    rast_xy <- rast_tbl %>% dplyr::select(x,y) 
    #calculate sum across band values
    rast_tbl<- rast_tbl %>% dplyr::select(-x,-y) %>% dplyr::mutate(newSum = rowSums(.))
    rm(rasterData)
    gc() 
    #make temp folder to write to
    if(!dir.exists(paste0(outFolder,'temp'))){
      dir.create(paste0(outFolder,'temp'))
    }
    #iterate through bands
    for(i in 1:length(bands[[1]])){#iterate through bands
      layerName <- names(rast_tbl)[i] #get dynamic band name
      #replace 0 values where row sum is 0 and band value is 0
      layer_tbl<- rast_tbl %>% dplyr::mutate(!!sym(layerName) := ifelse(!!sym(layerName)==0 & newSum==0,NA,!!sym(layerName))) %>% select(!!sym(layerName))
      #convert band to raster
      layer_out <- rast_xy %>% cbind(layer_tbl) 
      rm(layer_tbl,layerName)
      gc()
      layer_out <- layer_out %>% 
        tidyterra::as_spatraster(crs = "EPSG:27700",xycols=1:2)
      #join to outside object
      writeRaster(layer_out, paste0(outFolder,'temp/',i,'_', basename(rastRow$path)))
      rm(layer_out)
      gc()
    }
    rm(rast_tbl,rast_xy)
    gc()
    #list all
    files<- data.frame(file = list.files(paste0(outFolder,'temp'), pattern='img')) %>% filter(!str_detect(file,'img.aux.xml')) %>% mutate(no = as.numeric(str_remove(str_extract(file,"^.{2}"),'_'))) %>% arrange(no)
    #read back in
    rOut<- NULL
    for(tif in files$file){
      r <- rast(paste0(outFolder,'temp/',tif))
      rOut<- c(rOut,r)
    }
    #convert to raster brick
    rasterData_out<- rast(rOut)
    rm(rOut,r,rasterData)
    gc()
    #write out
    writeRaster(rasterData_out, paste0(outFolder, basename(rastRow$path)))
   # unlink(paste0(outFolder,'temp/'), recursive=TRUE)
    print(paste('S2 cloud processing done:', rastRow$var))
  }
}
