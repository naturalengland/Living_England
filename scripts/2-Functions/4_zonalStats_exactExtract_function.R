## ############################## ##
##
## Script Name: zonalStats_exactExtract_function.R
##
## Description: Function for generating segment-by-segment zonal stats for a 
## list of rasters 
##
## Author: Miles Clement, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2023-06-07
## 
## Date Last Modified: 2024-04-23
## 
## Versioning: v1.3
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## lubridate_1.9.2     forcats_1.0.0      
## dplyr_1.1.3         purrr_1.0.2        
## readr_2.1.4         tidyr_1.3.0        
## tibble_3.2.1        ggplot2_3.4.4      
## tidyverse_2.0.0     terra_1.7-29       
## stringr_1.5.0       sf_1.0-12          
## furrr_0.3.1         future_1.32.0      
## exactextractr_0.9.1
##
## ############################## ##


#' Function for extracting zonal statistics using the exactExtract function
#'
#' @param segmentation Polygon layer (opened with st_read) for which zonal statistics are to be calculated.
#                  The layer must have a column named ID containing unique values.
#' @param rasterDF a data frame containing var, path, stat, layers
#' @param parallel optional argument to determine if to be used in parallel. Default FALSE
#' @param numCores if parallel is TRUE, then number of cores to use. Default is availableCores() - 2
#'
#' @return Returns a table  zonal statistics for each polygon. 
#' @export
#'
#' @examples
zonalStats_exactExtract <- function(segmentation, rasterDF, 
                                    parallel = FALSE,
                                    numCores = availableCores() - 2){
  #load packages
  library(exactextractr)
  library(furrr)
  library(sf)
  library(stringr)
  library(terra)
  library(tidyverse)
  
  # Abort programme function
  exit <- function(){invokeRestart("abort") }   
  startTime <- Sys.time()
  
  # Check segmentation is sf object format
  if (class(segmentation)[1] != "sf"){
    stop("Segmented polygons must be a polygon layer opened using 'st_read')")
  }
  
  # Check that list of raster/band combinations exist prior to processing using terra
  # Break loop if error
  for (r in 1:nrow(rasterDF)){
    #select first row
    rastRow <- rasterDF[r,] 
    #check band combinations exist, otherwise return user error
    errorCheck <- tryCatch (rast(rastRow$path, lyrs=rastRow$layers),error=function(e) e)
    if (inherits(errorCheck, "error")){
      print(paste0("Invalid Raster/Band combination: ",rastRow$path,": Band ",rastRow$layers))
      print("Process aborting...")
      exit()
    } 
  }
  print("Data I/O checks complete: no issues")
  
  # Check that dataset CRS match
  # Filter listRasters so remove any duplicate checks
  for (r in 1:nrow(rasterDF)){
    # get raster
    rastRow <- rasterDF[r,] 
    #check it opens
    rastCheck <- rast(rastRow$path)
    #check coordinate reference systems match with segmentation
    rastCRS <- crs(rastCheck,describe=TRUE)
    vectCRS <- st_crs(segmentation) 
    # The format of the CRS in vector and raster datasets vary, so a simple == doesn't work
    # The function extracts the written version of the coordinate system:
    # In the testing example, SpatRaster will produce "British_National_Grid"
    # and sfc_polygon will produce "OSGB36 / British National Grid"
    # Checks are the made to see if the raster PROJ name (minus the '_') is within the vector PROJ
    if (grepl(str_replace_all(rastCRS$name, "_", " "),vectCRS$input) == FALSE){
      print("Raster CRS does not match Segmentation CRS")
      print(paste0(rastRow$path," CRS: ",rastCRS$name))
      print(paste0("Segmentation CRS: ",vectCRS$input))
      print("Process aborting...")
      exit()
    }
  }
  print("Data CRS checks complete: no issues")
  
  # Create output df with segment IDs
  zonalStats_df <- data.frame(SegID=segmentation$SegID)
  
  # if not running processing in parallel
  if (parallel == FALSE){
    print("Process not running in parallel...")
    #iterate through raster layers
    for (r in 1:nrow(rasterDF)){
      rastRow <- rasterDF[r,]   
      
      # Separate out filepath, function and layers/bands 
      file <- rastRow$path
      stats <- unlist(rastRow$stat)
      bands <- rastRow$layers
      
      # Load data (specific bands only)
      rasterData <- rast(file,lyrs=bands)
      # run for each statistic in reference table
      for (LEstat in stats){
        print(paste0("Processing ",file," ",LEstat))
        # Use exact_extract to undertake zonal stats
        zonalStats <- exact_extract(rasterData,
                                    segmentation,
                                    fun=LEstat,
                                    force_df=TRUE,
                                    full_colnames=TRUE,
                                    progress=FALSE)
        # Rename columns in standard format
        for (layername in names(rasterData)){
           name <- paste0(rastRow$var,"_",LEstat)
           if (length(unlist(bands)) > 1){ 
             name <- paste0(name, "_band", str_remove(as.character(layername),'Layer_'))
           }
           names(zonalStats)[names(zonalStats) == paste0(LEstat,'.',layername)] <- name
         }
        
        # Combine with segments IDs
        zonalStats_df <- cbind(zonalStats_df, zonalStats)
      }
    }
  }
  
  # for processing in parallel...
  if (parallel == TRUE){
    print("Attempting to run process in parallel...")
    #check on number of cores
    if (numCores > availableCores()){
      print("Specified cores greater than those available")
      print(paste0("Available cores: ",availableCores()))
      print("Process aborting...")
      exit()
    }
    #set up parallel
    plan(multisession, workers = numCores)
    # run zonal stats using future package via furrr
    zonalStats_furrr <- furrr::future_map_dfc(1:nrow(rasterDF),.f=function(r){
      #set seed
      furrr_options(seed=42)
      
      #open raster
      rastRow <- rasterDF[r,]   
      
      # Separate out filepath, function and layers/bands 
      file <- rastRow$path
      stats <- unlist(rastRow$stat)
      bands <- rastRow$layers
      
      # Load data (specific bands only)
      rasterData <- rast(file,lyrs=bands)
      
      # Create temp dataset for binding results within parallel processing
      zonalStatsTemp <- zonalStats_df
      # iterate through statistics from reference table
      for (LEstat in stats){
        print(paste0("Processing ",file," ",LEstat))
        # Use exact_extract to undertake zonal stats
        zonalStats <- exact_extract(rasterData,
                                    segmentation,
                                    fun=LEstat,
                                    force_df=TRUE,
                                    full_colnames=TRUE,
                                    progress=FALSE)
        
        # Rename columns
        for (layername in names(rasterData)){
          name <- paste0(rastRow$var,"_",LEstat)
          if (length(unlist(bands)) > 1){ 
            name <- paste0(name, "_band", str_remove(as.character(layername),'Layer_'))
          }
          names(zonalStats)[names(zonalStats) == paste0(LEstat,'.',layername)] <- name
        }
        # Bind stats to segment ID
        zonalStatsTemp <- cbind(zonalStatsTemp, zonalStats)
      }
      # Remove duplicate ID columns
      zonalStatsTemp %>% dplyr::select(-SegID)
    })
    
    # Combine with segments IDs
    zonalStats_df <- cbind(zonalStats_df, zonalStats_furrr)
  }
  # return reporting of runtimes
  runtime <- difftime(Sys.time(), startTime, units='hours')
  print("Zonal Stats processing complete")
  print(paste0("Runtime: ", format(round(runtime,2),nsmall=2)))
  #return zonal statistics dataframe
  return(zonalStats_df)
}
