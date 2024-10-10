## ############################## 
##
##
## Script Name: OSMM_Water_Urban_Preprocess_Functions
##
## Description: Functions for preprocessing MasterMap for input into Living England segmentation process
## 
## Author: Miles Clement, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2023-07-05
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
## units_0.8-1 
## 
## ############################## ##


#' OSMM preprocessing function - main function call
#'
#'@param OSpath filepath to OSMM file
#' @param bgz LE biogeographic zone
#' @param refImage reference raster image to align OSMM features to
#'
#' @return
#' @export
#'
#' @examples
mm_preProcess <- function(OSpath, bgz, refImage){
  
  #load packages
  library(data.table)
  library(dplyr)
  library(sf)
  library(smoothr)
  library(stringr)
  library(terra)
  library(tidyterra)
  
  print("Reading in data...")
  # Load in S2 image for BGZ to build grid to align data too
  raster_template <- rast(refImage)
  # Generate raster grid
  raster_template <- rast(ext(raster_template), 
                          resolution = 10,
                          crs = crs(raster_template))
  
  # Read in Mastermap for BGZ, using SQL to restrict to features of interest
  # For Urban, 'MAKE' column, 'Manmade' and 'Multiple' features
  # For Water, 'DESCRIPTIVEGROUP' column, 'Inland Water' and 'Tidal Water' features
  #create query
  qry <- sprintf("SELECT FID, DESCRIPTIVEGROUP, MAKE, THEME, SHAPE 
               FROM OSMM_Carto_BGZ_%s_subset
               WHERE MAKE = 'Manmade' OR MAKE = 'Multiple' OR THEME LIKE '%s'", bgz, '%Water%')
  #read in with query
  mm <- vect(OSpath, query = qry)
  
  # Split Data into separate Urban and Water subsets
  # Additional filters for urban to remove manmade water and landform features
  mm_Urban <- mm %>% 
    tidyterra::filter(MAKE %in% c("Manmade", "Multiple")) %>% 
    tidyterra::filter(!(DESCRIPTIVEGROUP %like% "Water" | DESCRIPTIVEGROUP %like% "Landform"))
  mm_Water <- mm %>% 
    tidyterra::filter(str_detect(THEME, "^Water"))
  
  # Run alignment functions and export
  mm_Water <- align_to_S2(mm_Water, "Water", raster_template)
  st_write(mm_Water,paste0("OUTPUTS/MM_Water_SegIN_BGZ",bgz,".shp"))
  
  mm_Urban <- align_to_S2(mm_Urban, "Urban", raster_template)
  st_write(mm_Urban,paste0("OUTPUTS/MM_Urban_SegIN_BGZ",bgz,".shp"))
}

#' Remove polygons below MMU (Minimum Mapping Unit)
#'
#' @param polygons OSMM polygon shapes
#' @param mmu minimum mapping unit, m^2. default is 300m^2
#'
#' @return
#' @export
#'
#' @examples
remove_below_mmu <- function(polygons,mmu=300){
  #load packages
  library(smoothr)
  library(units)
  # set 
  fill_params <- units::set_units(mmu, m^2)
  print("Dropping crumbs & Filling holes...")
  polygons <- drop_crumbs(polygons, threshold = fill_params) #remove polygons below threshold
  polygons <- fill_holes(polygons, threshold = fill_params) # fill polygons below threshold
  #return polygon object
  return(polygons)
}

#' Align polygons to S2 grid
#'
#' @param polygons OSMM polygon shapes
#' @param hab habitat feature name, 'Water' or 'Urban'
#' @param raster_template spatRaster object of S2 grid
#'
#' @return
#' @export
#'
#' @examples
align_to_S2 <- function(polygons, hab, raster_template){
  
  #load packages
  library(data.table)
  library(dplyr)
  library(sf)
  library(smoothr)
  library(stringr)
  library(terra)
  library(tidyterra)
  
  #reporting
  print(paste0("BGZ",bgz,": Processing Start, ",hab))
  print(Sys.time())
  
  # Aggregate polygons that are joined, and disaggregate to split multipolygons
  print("Aggregating/Disaggregating polygons...")
  polygons <- terra::aggregate(polygons, dissolve=TRUE) 
  polygons <- disagg(polygons) 
  
  # Convert from spatvector to sf
  polygons <- sf::st_as_sf(polygons)
  polygons <- remove_below_mmu(polygons)
  
  #process water features
  if (hab == "Water"){
    # Dilute and Expand water features to remove those no likely to be visible in S2 imagery 
    # Remove any empty features
    print("Buffering water polygons...")
    polygons <- st_buffer(polygons, -7)
    polygons <- polygons[!st_is_empty(polygons),]
    polygons <- st_buffer(polygons, 7)
  }
  
  #process urban features
  if (hab == "Urban"){
    # Dilute and Expand urban features to remove those no likely to be visible in S2 imagery 
    # Remove any empty features
    print("Buffering urban polygons...")
    polygons <- st_buffer(polygons, -3)
    polygons <- polygons[!st_is_empty(polygons),]
    polygons <- st_buffer(polygons, 3)
  }
  
  # Rasterise Data to S2 Grid, then reconvert to polygons and disaggregate
  # Convert sf to spatvector format
  polygons <- vect(polygons)
  print("Rasterising polygons...")
  polygons <- rasterize(polygons, raster_template,field=1)
  print("Vectorising raster...")
  polygons <- as.polygons(polygons, dissolve=TRUE) 
  polygons <- disagg(polygons) 
  
  # convert spatvector to sf
  polygons <- sf::st_as_sf(polygons)
  # Remove polygons and holes below min mapping unit
  polygons <- remove_below_mmu(polygons)
  
  #reporting
  print(paste0("BGZ",bgz,": Processing End, ",hab))
  print(Sys.time())
  
  #return polygon object
  return(polygons)
}
