## Living England summarise all ground data
##
## Description: Script to summarise all the ground data ingested for Living England. This is prior to any qc steps and follows on after the data ingestion pipeline.
##
## Authors: Becky Trippier, Natural England
##
## Licence: MIT Licence  
##
## Date Created: 2023-01-18
## 
## Date Last Modified: 2024-10-02
## 
## Versioning:
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
##xlsx_0.6.5       openxlsx_4.2.5.2
## forcats_1.0.0    tidyr_1.3.0     
## tibble_3.2.1     ggplot2_3.4.4   
## tidyverse_2.0.0  purrr_1.0.2     
## readr_2.1.4      stringr_1.5.0   
## lubridate_1.9.2  dplyr_1.1.3     
## sf_1.0-12  
#---------------------------------------------------

#' Run summary function
#'
#' @param datfolder required, character string of the folder location the processed data are stored in
#' @param selectYear optional, numeric year to subset data and output on a separate tab a summary of the points per habitat for that cut-off year to present.
#'
#' @return
#' @export
#'
#' @examples
#' runSumm(datfolder = 'D:/Projects/Living England/Data Layers/Ground Datasets/GroundSurveys_Processed',selectYear=2020)
run_summ <- function(datfolder = './Ground Datasets/GroundSurveys_Processed',
                    selectYear = NA){
  
  # Load libraries
  source('2-Functions/common_functions_library.R')
  libaries <- c("sf", "dplyr", "lubridate", "stringr", "readr", "purrr", "tidyverse","openxlsx","lubridate","xlsx")
  install_and_load(libaries)
  
  #create file list
  all_list <- data.frame(files=list.files(paste0(datfolder,'/3_surveyData_segment_assigned'), 
                                         pattern='csv', full.names = T))
  
  #iterate through files and join all datasets
  all_dat <- map_df(all_list$files, .f=function(survey){
    shp_post <- read.csv(survey) %>% 
      dplyr::select(-X) %>%
      dplyr::mutate(Source_ID = as.character(Source_ID))%>%
      dplyr::mutate(Source_Broad_Name = as.character(Source_Broad_Name))%>%
      dplyr::mutate(Source_Detailed_Habitat = as.character(Source_Detailed_Habitat))
  })
  
  #add unique identifier
  all_dat$point_id <- 1:nrow(all_dat)
  
  #get current date
  date <- str_remove_all(Sys.Date(),'-')
  
  #write out all compiled data
  write.csv(all_dat,paste0(datfolder,'/', date,'_GroundSurveys_data.csv'))
  
  # summarise results by survey
  broad_post <- all_dat %>% 
    group_by(LE_Broad_Habitat,Survey_Code) %>% 
    summarise(total=n()) %>%
    pivot_wider(names_from=LE_Broad_Habitat,values_from = total,values_fill=0)
  #write out summary by survey tab
  xlsx::write.xlsx(broad_post,file=paste0(datfolder,'/', date,'_GroundSurveys_summary.xlsx'), sheetName = 'summary_hab_survey')
  
  # summarise results by bgz
  broadbgz_post <- all_dat %>% 
    group_by(LE_Broad_Habitat,bgz) %>% 
    summarise(total=n()) %>%
    pivot_wider(names_from=LE_Broad_Habitat,values_from = total,values_fill=0) 
  # write out summery by bgz tab
  xlsx::write.xlsx(broadbgz_post, file=paste0(datfolder,'/',date,'_GroundSurveys_summary.xlsx'), 
                   sheetName = 'summary_hab_bgz',append=TRUE)
  
  # if cutoff year provided
  if(is.na(selectYear)==FALSE){
    if(is.numeric(selectYear)==FALSE){selectYear <- as.numeric(selectYear)}
    # summarise results by selected year cutoff e.g. 2020
    broad_recent <- all_dat %>% 
    mutate(Survey_Date = as.Date(Survey_Date)) %>% 
    mutate(year = year(Survey_Date)) %>%
    dplyr::filter(year>=selectYear) %>%
    group_by(LE_Broad_Habitat,bgz) %>% 
    summarise(total=n()) %>%
    pivot_wider(names_from=LE_Broad_Habitat,values_from = total,values_fill=0) 
  #write out summary by select cutoff year
  xlsx::write.xlsx(broad_recent, file=paste0(datfolder,'/',date,'_GroundSurveys_summary.xlsx'), 
                   sheetName = paste0('summary_since_',selectYear),append=TRUE)
  }
  print(paste0('file written: ',date,'_GroundSurveys_summary.xlsx'))
  
}


