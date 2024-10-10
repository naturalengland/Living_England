## ############################## 
##
##
## Script Name: Generate Segment IDs
##
## Description: Update segmentation attributes to include unique segment IDs
## ID Format: LEXXXX_BGZYY_ZZZZZZZ
## X = Imagery Years
## Y = BioGeographic Zone
## Z = 7 digit unique number
## Example: LE2223_BGZ06_1234567
## 
## Author: Clement NE, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2023-08-24
## 
## Date Last Modified: 2024-04-18
## 
## Versioning:
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## sf_1.0-13
## 
## ############################## ##


#' Function to generate segment ids
#'
#' @param bgz biogeographic zone
#' @param year LE iteration year e.g. 2223
#' @param folderPath path to ecog raw output
#' @param outputPath path to save outputs
#'
#' @return
#' @export
#'
#' @examples
generateIDs = function(bgz, year, folderPath,outputPath){
  
  #load packages
  library(sf)
  
  #read in segmentation bgz output 
  segmentation <- st_read(paste0(folderPath, "LE",year,"_BGZ",bgz,"_Segmentation_eCogOut.shp"),stringsAsFactors=FALSE)
  
  # Generate Unique IDs
  segmentation$SegID = paste0("LE",year,"_BGZ",bgz,"_",sprintf("%07d", as.numeric(rownames(segmentation))))
  # Reorder Columns
  segmentation = segmentation[,c("SegID","Class_name","geometry")]
  #set standard output file name
  outname = paste0(outputPath,"LE",year,"_BGZ",bgz,"_Segmentation_withIDs.shp")
  #write out features
  st_write(segmentation,outname,append=FALSE)
  
}