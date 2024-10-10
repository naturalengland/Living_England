## Clean coordinates function
##
## Description: Function to make sure coodinates are in the correct geometries and crs required
##
## Authors: Becky Trippier, Natural England
##
## Licence: MIT Licence
## Date created: 2023-01-24
##
## Date Last Modified: 2024-10-02
## 
## Versioning:
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## dplyr_1.1.3 
## sf_1.0-12 
#---------------------------------------------------


#' clean coordinates
#'
#' @param spatialDat spatial dataframe
#' @param crs the EPSG code for the coordinate reference system, e.g. British National Grid EPSG:27700
#'
#' @return
#' @export
#'
#' @examples
  #' cleanCoords(BBT,27700)
#' 
clean_coords <- function(spatialDat, crs= 27700){
  
  #load packages
  library(sf)
  library(dplyr)
  
  # data check 1 - correct geometries
  if (all(unique(as.character(st_geometry_type(spatialDat, by_geometry = T))) %in% c("POINT", "MULTIPOINT", "POLYGON","MULTIPOLYGON"))==FALSE){stop("Data needs to be either point or polygon class.")}
  #report if mix of geometries
  if (all(unique(as.character(st_geometry_type(spatialDat, by_geometry = T))) %in% c("POINT", "MULTIPOINT"))){
    print("All geometries spatial points")
  } else if(all(unique(as.character(st_geometry_type(spatialDat, by_geometry = T))) %in% c("POLYGON", "MULTIPOLYGON"))){
    print("All geometries spatial polygons") 
  } else{
      print("Geometries mix of spatial points and polygons.")
    }
  
  #data check 2 - all in the correct crs
  if(st_crs(spatialDat)$epsg !=crs){
    spatialDat <- spatialDat %>% st_transform(st_crs(crs))
    print(paste0("reprojected to ",crs))
  }else{print("Correct crs")}
  
  #data check 3 - check all valid geometries
  if(all(st_is_valid(spatialDat)) == FALSE){
    habNotValid<- spatialDat %>% dplyr::filter(st_is_valid(spatialDat)==FALSE) %>% st_buffer(0.00001)
    spatialDat <- spatialDat %>% dplyr::filter(st_is_valid(spatialDat)==TRUE) %>% rbind(habNotValid)
    print(paste("not all valid geometries:", nrow(spatialDat), 'contained geometries not valid, buffered by 0.0001'))
  }else{print("All valid geometries")}
  #return spatial object
  return(spatialDat)
}
