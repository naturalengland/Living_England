## ##############################
##
## Script Name: Habitat Data QA Functions
##
## Description: Function to run through habitat data qa process. This will remove any habitat records that are no longer representative of the habitat assigned to it. Determines expected spectral ranges for each habitat, and filters remaining records based on these
##
## Author: Miles Clement, Becky Trippier, Natural England
##
## Licence: MIT Licence
## Date Created: 2023-11-07
##
## Date Last Modified: 2024-04-23
##
## Versioning:
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## data.table 1.14.8
## dplyr 1.1.2
## lwgeom 0.2-11
## sf 1.0-13
## tidyr_1.3.0
## janitor_2.2.0
## purrr_1.0.2
## stringr_1.5.0
## ############################## ##

# Main qa function for preparing habitat data
#'
#' @param hab_Input input data object
#' @param output_folder folder to save the outputs
#' @param data_date ingestion date
#' @param zonal_stats folder path to zonal statistics
#' @param qualitySurveys list of quality surveys to filter on
#' @param filterDate date to filter representative filter data to
#' @param band_list list of bands to filter with
#' @param vecHabList list of vector habitats to split output into RF habs and vector habs
#'
#' @return
#' @export
#'
#' @examples
habQAfilter <- function(hab_Input, output_folder, data_date, zonal_stats ,qualitySurveys, filterDate, band_list, vecHabList){

  #load libraries
  library(data.table)
  library(dplyr)
  library(lwgeom)
  library(sf)
  library(janitor)
  library(tidyr)
  library(purrr)
  library(stringr)
  
  ## create folder structure
  dirs <- c('2.OPS','ZonalStats_QA','3.ZonalMerge','4.CloudRemoved','5.Filters','6.Filtered','7.RFinput','8.Spare')
  for (i in dirs){
    if(!dir.exists(paste0(output_folder, i))){
      dir.create(paste0(output_folder, i))
    }
  }
  
  #### Run One per segment check ####
  print('processing One-per-segment...')
  
  # sort by latest date
  hab_OPS <- hab_Input %>%
    dplyr::arrange(desc(Survey_Date)) %>%
    dplyr::distinct(SegID, HabCode, .keep_all = TRUE) # take  unique segment ID and habitats - taking the first row (latest date) where duplicates
  
  # store One per segment output
  write.csv(hab_OPS, paste0(output_folder, "2.OPS/", data_date, "_FieldData_OPS.csv"))
  
  # report summary
  OPSsumm <- hab_OPS %>%
    dplyr::group_by(HabCode, bgz) %>%
    dplyr::summarise(count = n()) %>%
    tidyr::spread(key = bgz, value = count)
  
  # store
  write.csv(OPSsumm, paste0(output_folder, "2.OPS/", data_date, "_FieldData_OPS_summ.csv"))
  
  ## collating zonal stats ####
  cat('Collating zonal stats...\n')
  
  # Iterate through Zonal Stats csvs, and reduce columns to those needed
  for (bgz in sort(unique(hab_OPS$bgz))) {
    cat("Processing bgz", bgz, "zonal stats\n")

    # read in zonal stats file for bgz
    zs_temp <- read.csv(paste0(zonal_stats, bgz, "_zonalStats.csv"))

    ## Check all bands required in ZS file
    if (any(!band_list %in% names(zs_temp)) == TRUE) {
      stop("Some required band_list variables are not present in zonal stats file, please check names.")
    }

    # collects required variables
    zs_temp <- zs_temp %>% dplyr::select(SegID, all_of(band_list))

    # writes out subset of zs
    write.csv(zs_temp, paste0(output_folder, "ZonalStats_QA/", bgz, "_zonalStats_QA.csv"))
  }
    
    ## Load in Zonal Stats csv's, and merge with habitat records using Segment ID ####
    hab_zs <- purrr::map_df(sort(unique(hab_OPS$bgz)), .f = function(bgz) {
      # read in zs subset above
      zs_temp <- read.csv(paste0(output_folder, "ZonalStats_QA/", bgz, "_zonalStats_QA.csv"))
      # filter to just those within habitat dataset
      hab_zs <- hab_OPS %>%
        dplyr::filter(SegID %in% zs_temp$SegID) %>%
        dplyr::left_join(zs_temp, by = "SegID")
      hab_zs
    })
    # Write out to save rerunning above if needed
    write.csv(hab_zs, paste0(output_folder, "3.ZonalMerge/", data_date, "_FieldData_zonalMerge.csv"))
    
    ## Optional reread
    # hab_zs <- read.csv(paste0(output_folder,"3.ZonalMerge/",data_date,"_FieldData_zonalMerge.csv"))
    # subset to just necessary fields
    hab_zs <- hab_zs %>% dplyr::select(SegID, bgz, Survey_Code, Survey_Date, HabCode, easting, northing, all_of(band_list))
    
    # Create list of bgzs from zonal stats
    bgz_list <- sort(unique(hab_zs$bgz))
    
    # Add Flag column for Later Filtering Steps
    hab_zs$Flag <- 0
    
    ### Filter out points that fall in S2 cloud regions (0s across all three seasonal bands) ####
    cat('Finding cloudy segments...\n')
    
    # extract segments where cloud present in a season
    spring <- hab_zs %>%
      dplyr::select(S2_spr_mean_band2:S2_spr_mean_band12) %>%
      mutate(rowSum = rowSums(.))
    summer <- hab_zs %>%
      dplyr::select(S2_sum_mean_band2:S2_sum_mean_band12) %>%
      mutate(rowSum = rowSums(.))
    autumn <- hab_zs %>%
      dplyr::select(S2_aut_mean_band2:S2_aut_mean_band12) %>%
      mutate(rowSum = rowSums(.))

    # flag cloudy segs
    hab_zs_cloud <- hab_zs %>% mutate(
      springSum = spring$rowSum,
      summerSum = summer$rowSum,
      autumnSum = autumn$rowSum
    )
    
    # report number with cloudy seg in a season
    cloudy_segs <- hab_zs_cloud %>% dplyr::filter(is.na(springSum) | is.na(summerSum) | is.na(autumnSum))

    # report number with cloudy segs in all seasons
    all_cloudy_segs <- hab_zs_cloud %>% dplyr::filter(is.na(springSum) & is.na(summerSum) & is.na(autumnSum))

    # report numbers of cloudy
    cat(nrow(cloudy_segs), "segments with some cloud present.", nrow(all_cloudy_segs), "segments where all seasons have cloud present\n")
    
    # remove any with no cloud across all three
    if (nrow(all_cloudy_segs) > 0) {
      hab_zs_cloud_na <- hab_zs_cloud %>% dplyr::filter(SegID %in% all_cloudy_segs$SegID)
      print(paste(nrow(all_cloudy_segs), "segments removed where all mosaics with cloud present."))
    } else {
      hab_zs_cloud_na <- hab_zs_cloud
    }
    # write out reporting
    cat(c(nrow(cloudy_segs), "segments with some cloud present.", nrow(all_cloudy_segs), "segments where all seasons have cloud present."), file = file.path(output_folder, "4.CloudRemoved/cloud_reporting.txt"), sep = "\n")

    # Write Out - dataset with cloud converted to NAs
    write.csv(hab_zs_cloud_na, paste0(output_folder, "4.CloudRemoved/", data_date, "_FieldData_cloudNA.csv"))
    
    # Create Subset for high reliability baseline points and filter against ####
    cat('create filter...\n')

    # Points collected since cut off, and high-confidence surveys, with no cloud
    hab_yr <- hab_zs_cloud_na %>%
      dplyr::filter(Survey_Date >= filterDate) %>%
      dplyr::filter(Survey_Code %in% qualitySurveys)

    # Write Out
    write.csv(hab_yr, paste0(output_folder, "5.Filters/", data_date, "_FieldData_Filter.csv"))

    # reporting filter
    hab_yr_summ <- hab_yr %>%
      dplyr::group_by(HabCode, bgz) %>%
      dplyr::summarise(count = n()) %>%
      tidyr::spread(key = bgz, value = count)
    
    # export to csv
    write.csv(hab_yr_summ, paste0(output_folder, "5.Filters/", data_date, "_FieldData_Filter_summ.csv"))

    # Generate List of Habitats
    hab_list <- unique(hab_zs_cloud_na$HabCode)

    ## Check High quality survey data and remove outliers ####
    outliers <- quality_check(band_list, hab_list, hab_yr, out_path = paste0(output_folder, "5.Filters/"))
    outliers <- outliers %>%
      dplyr::select(SegID, Survey_Code, Survey_Date, HabCode) %>%
      mutate(flag2 = 1)
    
    ### remove outliers here ##
    hab_yr <- hab_yr %>%
      left_join(outliers, by = c("SegID", "Survey_Code", "Survey_Date", "HabCode")) %>%
      dplyr::filter(is.na(flag2)) %>%
      dplyr::select(-flag2)
    cat(nrow(outliers), "outliers removed\n")

    # Run function to build spectral boundaries ####
    cat("Building spectral ranges\n")
    calculate_ranges(band_list, hab_list, hab_yr, out_path = paste0(output_folder, "5.Filters/"))

    ### Filter data according to filter range ####
    # Create output dataset with same colnames as hab_zs_cloud
    hab_filtered <- hab_zs_cloud_na[1, ]
    hab_filtered <- hab_filtered[-1, ]
    
    # Call filter function - filter those with no cloud by subset ranges
    hab_filtered <- filter_data(hab_zs_cloud_na, hab_filtered, band_list = band_list, out_path = paste0(output_folder, "5.Filters/"))

    ## QC check - for cloud in output - may not be present if lucky, but likely some
    hab_filtered %>%
      dplyr::filter(is.na(S2_spr_mean_band2)) %>%
      nrow()

    # Output all Data filtered data
    write.csv(hab_filtered, paste0(output_folder, "6.Filtered/", data_date, "_FieldData_Filtered.csv"))

    ### Summarise and Tidy output for modelling ####
    # Remove unneeded columns
    hab_filtered <- hab_filtered %>% dplyr::select(SegID:northing)

    # report summary
    filter_summ <- hab_filtered %>%
      dplyr::group_by(HabCode, bgz) %>%
      dplyr::summarise(count = n()) %>%
      tidyr::spread(key = bgz, value = count)
    
    # store
    write.csv(filter_summ, paste0(output_folder, "6.Filtered/", data_date, "_FieldData_Filtered_summ.csv"))
    
    ## create final exported datasets ####
    # Final RF habs
    hab_final <- hab_filtered %>% dplyr::filter(!HabCode %in% vecHabList)
    write.csv(hab_final, paste0(output_folder, "7.RFinput/", data_date, "_FieldDataQA_RFinput.csv"))

    # Separate out habitats not classified by RF
    hab_spare <- hab_filtered %>% dplyr::filter(HabCode %in% vecHabList)
    write.csv(hab_spare, paste0(output_folder, "8.Spare/", data_date, "_FieldDataQA_nonRF_Habs.csv"))
    
    # report to user
    cat("RF and non-RF filtered datasets generated. hab QA complete\n")
  }
  
#' QA check of high quality data
#'
#' @param band_list List of spectral bands to filter on
#' @param hab_list List of habitat classes to filter
#' @param hab_yr input data filtered by year and survey
#' @param out_path filepath to save to
#'
#' @return
#' @export
#'
#' @examples
quality_check <- function(band_list, hab_list, hab_yr, out_path){
  
  # create folder for storing check images
  if (!dir.exists(paste0(out_path, "HighQualCheck"))) {
    dir.create(paste0(out_path, "HighQualCheck"))
  }
  
  # iterate through habitats in list
  flag_outlier_all <- NULL
  for (hab in hab_list) {
    # Subset points to Hab
    points_hab <- hab_yr %>% dplyr::filter(HabCode == hab)
    band_check <- data.frame(band = band_list) %>%
      filter(str_detect(band, "S2")) %>%
      filter(str_detect(band, "min|max"))

    # iterate through S2 band
    for (bandno in paste0("band", 2:12)) {
      band_check_it <- band_check %>% filter(str_detect(band, bandno))

      # iterate through min and max
      for (minmax in c("min", "max")) {
        min_max_check_it <- band_check_it %>% filter(str_detect(band, minmax))
        flag_outlier <- NULL

        # iterate through season
        for (season in c("spr", "sum", "aut")) {
          band <- min_max_check_it %>% filter(str_detect(band, season))
          point_hab_band <- points_hab %>%
            dplyr::select(SegID, Survey_Code, Survey_Date, HabCode, easting, northing, all_of(band$band)) %>%
            rename(bandvar = all_of(band$band))

          # find Q1, Q3, and interquartile range
          Q1 <- quantile(point_hab_band$bandvar, .25, na.rm = T)
          Q3 <- quantile(point_hab_band$bandvar, .75, na.rm = T)
          IQR <- IQR(point_hab_band$bandvar, na.rm = T)

          # subset data where points value is outside 1.5*IQR of Q1 and Q3
          outliers <- subset(point_hab_band, point_hab_band$bandvar < (Q1 - 1.5 * IQR) | point_hab_band$bandvar > (Q3 + 1.5 * IQR))
          outliers <- outliers %>%
            dplyr::select(SegID, Survey_Code, Survey_Date, HabCode, easting, northing) %>%
            mutate(seas = season)

          # add to all outliers
          flag_outlier <- rbind(flag_outlier, outliers)
        }

        # make wider format to flag across seasons
        flag_outlier <- flag_outlier %>%
          mutate(val = 1) %>%
          pivot_wider(names_from = "seas", values_from = "val")

        # add any missing season headers
        if (!all(c("spr", "sum", "aut") %in% names(flag_outlier))) {
          if (!"spr" %in% names(flag_outlier)) {
            flag_outlier <- flag_outlier %>% mutate(spr = NA)
          }
          if (!"sum" %in% names(flag_outlier)) {
            flag_outlier <- flag_outlier %>% mutate(sum = NA)
          }
          if (!"aut" %in% names(flag_outlier)) {
            flag_outlier <- flag_outlier %>% mutate(aut = NA)
          }
        } # else not all seasons flagged, no outliers to add
        flag_outlier_all = rbind(flag_outlier_all, flag_outlier)
      } # close minmax iterate
    }
  }
  # collate all outliers and summarise number of seasons flagged as outlier
  flag_outlier_all<- flag_outlier_all %>% rowwise() %>% mutate(allseas = sum(spr,sum,aut, na.rm=T))
  
  #remove duplicate rows
  flag_outlier_all<- flag_outlier_all %>% dplyr::distinct(SegID,HabCode,Survey_Code, .keep_all=TRUE)
  write.csv(flag_outlier_all,paste0(out_path,"all_outliers_reasons.csv"))
  
  # pull out where outlier in spring, summer and autumn s2 bands
  out_all_seas <- flag_outlier_all %>% dplyr::filter(allseas == 3) 
  
  #report how many
  print(paste(nrow(out_all_seas), 'points flagged as outliers.', 'Total high quality data points =', nrow(hab_yr)))
  out_all_seas %>% group_by(Survey_Code, HabCode) %>% summarise(nPoints = n()) %>% print()
  write.csv(out_all_seas,paste0(out_path,"outliers_all_seasons.csv"))
  return(out_all_seas)
    
}

#' Calculating the band range function
#'
#' @param band_list List of spectral bands to filter on
#' @param hab_list List of habitat classes to filter
#' @param hab_yr input data filtered by year and survey
#' @param out_path filepath to save to
#'
#' @return
#' @export
#'
#' @examples
calculate_ranges <- function(band_list, hab_list, hab_yr, out_path) {
  # iterate through habitats in list
  for (hab in hab_list) {
    # Subset points to Hab
    points_hab <- hab_yr %>% dplyr::filter(HabCode == hab)

    # Create Output DF with Min/Max for Each Band
    outlier_df <- data.frame(matrix(nrow = length(band_list), ncol = 3))
    colnames(outlier_df) <- c("Band", "Min", "Max")
    outlier_df$Band <- band_list

    # Find min/max spectral data for each band
    for (band in band_list) {
      band_min <- min(points_hab[[band]], na.rm = T)
      band_max <- max(points_hab[[band]], na.rm = T)

      outlier_df$Min[outlier_df$Band == band] <- band_min
      outlier_df$Max[outlier_df$Band == band] <- band_max
    }

    # Export Collector app data for specific habitat
    write.csv(outlier_df, paste0(out_path, hab, "_ranges.csv"), row.names = FALSE)
  }
}

#' function to filter the remaining data by range values
#'
#' @param hab Input field data DF
#' @param hab_filtered output DF
#' @param band_list list of bands to assess
#' @param out_path  filepath to save to
#'
#' @return
#' @export
#'
#' @examples
filter_data <- function(hab, hab_filtered, band_list, out_path){
  
  # set up flag list
  flagged_segs <- NULL
  
  #Iterate through habitats
  for (h in sort(unique(hab$HabCode))) {
    # Create data subset
    hab_subset <- hab %>% dplyr::filter(HabCode == h)
    # iterate through band
    for (band in band_list) {
      # If a limits file for that dataset exists, read in
      fname <- paste0(out_path, h, "_ranges.csv")
      if (file.exists(fname)) {
        limits <- read.csv(fname)
        # get band min and max
        band_min <- limits$Min[limits$Band == band]
        band_max <- limits$Max[limits$Band == band]
        # filter subset
        habFlagged <- hab_subset %>%
          dplyr::select(SegID, HabCode, all_of(band)) %>%
          dplyr::filter(!is.na(get(band))) %>% # remove any NAs - treated as unflagged
          dplyr::filter(get(band) < band_min | get(band) > band_max) %>%
          dplyr::select(SegID, HabCode) # condition less than band min and greater than band max, keep to flag
        flagged_segs <- rbind(flagged_segs, habFlagged) # add to flag list
      }
    }
    cat(h, "filtered\n")
  }
  
  # get unique segIDs flagged
  flag_out <- distinct(flagged_segs)

  # report number of duplicate ground data points
  dupe_habs <- flag_out %>% janitor::get_dupes(SegID)

  # filter habitat dataset to remove those flagged
  hab_out <- hab %>% dplyr::filter(!SegID %in% flag_out$SegID)
  
  # report
  print(paste(
    nrow(hab_out), "ground data rows passed filters representing", length(unique(hab_out$SegID)), "segments.", nrow(flag_out), "rows removed during filtering stage.",
    nrow(dupe_habs), "duplicate segment rows with", length(unique(dupe_habs$SegID)), "unique segments duplicated."
  ))
  
  cat(
    c(nrow(hab_out), "ground data rows passed filters representing", length(unique(hab_out$SegID)), "segments.", nrow(flag_out), "rows removed during filtering stage.",
      nrow(dupe_habs), "duplicate segment rows with", length(unique(dupe_habs$SegID)), "unique segments duplicated."),
    file = file.path(out_path, "QA_output.txt"), sep = "\n"
  )
  return(hab_out)
}