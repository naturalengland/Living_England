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
## Author: Miles Clement, Natural England
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
## stringr_1.5.0 
## ############################## ##

## load libraries
library(stringr)

## set user variables ##
bgzList <- str_pad(rep(1:13),2,pad='0')

# Update years of input imagery each iteration
year <- "2223"
#path to ecog output files
folderPath <- "./LE2223_Publication/temp/Segmentation/"
#folder to save output to
outputPath <- "./LE2223_Publication/Segmentation/"
#--------------------

#iterate through biogeographic zones
for (bgz in bgzList){
  #run segment id function 
  generateIDs(bgz, year,folderPath,outputPath)
  #report back
  print(paste('Segment ids generated',bgz))
}


