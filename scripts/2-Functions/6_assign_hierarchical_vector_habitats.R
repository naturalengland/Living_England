## ############################## 
##
## Script Name: Hierarchical layer assignment to segments
##
## Description: Script for assigning vector layers to segmentation in a given hierarchical order
##
## Author: Becky Trippier, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2024-01-26
## 
## Date Last Modified: 2024-07-22
## 
## R version 4.2.3 (2024-01-26 ucrt)
## Dependencies:
## tidyr_1.3.0   
## stringr_1.5.0 sf_1.0-12    
## dplyr_1.1.3 
##
## ############################## ##


#' Hierarchical layer assignment to segments
#'
#' @param vec_classes data frame listing layer names and the LE habitat class name. The layer names should reference folder named with the same convention, containing identified segments for that habitat class
#' @param bgz biogeographic zone
#' @param output_path file path to where the vector habitat data is stored
#' @param seg_path file path to segmentation shapefiles
#' 
#'
#' @return
#' @export
#'
#' @examples
assign_habitat <- function(vec_classes,bgz,output_path,seg_path){
  
  #read in packages
  library(tidyr)
  library(stringr)
  library(dplyr)
  
  #read in segmentation
  segments <- st_read(paste0(seg_path,"LE2223_BGZ",bgz,"_Segmentation_withIDs.shp"))
  
  #create new field to assign classes
  classification <-segments %>% mutate(Primary_Hab=NA) #create new class attribute to populate
  
  #iterate through layers in the supplied data frame
  for(i in 1:length(vec_classes$layer)){
    hab_layer <- vec_classes[i,] #select layer
    #read in layer data for bgz
    habCSV <- read.csv(paste0(output_path,hab_layer$layer,"/BGZ",bgz,"_",hab_layer$layer,"_segments.csv"))
    
    # where class attribute is empty and segid is in habitat folder, reclassify class to given habitat
    classification <- classification %>% mutate(Primary_Hab = ifelse(is.na(Primary_Hab) & SegID %in% habCSV$SegID,hab_layer$Primary_Hab, Primary_Hab))
    print(paste('Assigned', hab_layer$layer))  #report progress to user
  }
  #return output object
  return(classification)
}
