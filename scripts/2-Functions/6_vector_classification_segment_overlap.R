## ############################## 
##
## Script Name: vector_classification_segment_overlap
##
## Description: Function for classifying segments based on overlap with thematic datasets 
## Returns a dataframe of segments with overlap above user-defined threshold
##
## Author: Miles Clement, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2023-09-23
## 
## Date Last Modified: 2024-07-22
## 
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## dplyr      v.2.3.2
## sf         v.1.0-13
##
## ############################## ##

# Abort programme function
exit <- function(){invokeRestart("abort") }   

#' Vector Overlap Function - For each segment, checks overlap with classifying dataset returns dataframe containing segments with overlap above user-defined threshold
#'
#' @param segIn segmentation
#' @param thematic classifying dataset
#' @param target_pct % overlap required for classification
#'
#' @return
#' @export
#'
#' @examples
vector_overlap <- function(segIn,thematic,target_pct){
  
  #read in packages
  library(dplyr)
  library(sf)
  
  #report timing
  startTime <- Sys.time()
  
  # Check data CRS match, abort if they don't
  segCRS <- st_crs(segIn) 
  thematicCRS <- st_crs(thematic) 
  if (segCRS != thematicCRS){
    print("Data CRS does not match")
    print(paste0("Segmentation CRS: ",segCRS))
    print(paste0("Thematic CRS: ",thematicCRS))
    print("Process aborting...")
    exit()
  }
  print("Data CRS checks complete: no issues")
  
  #set up area overlap reporting
  segIn$Overlap_Area = 0
  segIn$Overlap_Pct = 0
  
  # For each segment calculate overlap with each CROME hexagon
  overlap_pct = st_intersection(segIn, thematic) %>% 
    mutate(intersect_area = st_area(.)) %>%     
    dplyr::select(SegID, intersect_area) %>%   
    st_drop_geometry()
  
  # Combine overlap areas for each segment to get total
  overlap_unique_ids = unique(overlap_pct$SegID)
  for(id in overlap_unique_ids){
    id_overlap = overlap_pct[overlap_pct$SegID==id,]
    area = sum(id_overlap$intersect_area)
    segIn <- within(segIn, Overlap_Area[SegID == id] <- area)
  }
  
  # Calculate as %, threshold based on user input, drop geometry to convert to dataframe
  segIn$Overlap_Pct = as.numeric(format(round((segIn$Overlap_Area/segIn$Area)*100,2),nsmall=2))
  segInHab <- segIn[segIn$Overlap_Pct >= target_pct,]
  segInHab <- st_drop_geometry(segInHab)
  
  #reporting script run times
  runtime <- difftime(Sys.time(), startTime, units='hours')
  print("Vector Overlap processing complete")
  print(paste0("Runtime: ", format(round(runtime,2),nsmall=2)))
  
  #return output object
  return(segInHab)
}