
## Living England aligning habitat records to segmentation
##
## Description: Function to add LE segment ID to survey data. Where survey data is a polygon this will find overlapping segments and create centroid points where the overlap meets a defined threshold.
##
## Authors: Miles Clement, Becky Trippier, Natural England
##
## Licence: MIT Licence
## Date created: 2022-10-24
##
## Date Last Modified: 2024-10-02
## 
## Versioning:
## R version 4.3.0 (2023-04-21 ucrt)
## Dependencies:
## furrr_0.3.1       future_1.32.0     purrr_1.0.2      
## units_0.8-1       janitor_2.2.0     sf_1.0-12        
## lwgeom_0.2-11     dplyr_1.1.3       data.table_1.14.8
## stringr_1.5.0
#---------------------------------------------------

# Living England has two types of training data, point and polygon
## Point data - segments points fall within put forward as training segments
## Polygon data - overlap between polygons and segments calculated, segments 
### over user defined overlap put forward as training segments

#' Function to assign segments to survey data
#'
#' @param surveyDat survey datatable
#' @param segPath path to segmentation file
#' @param bgzPath path to biogeographic zone file
#' @param overlapThreshold % threshold of overlap with segment in order to assign
#' @param cores optional - number of cores
#'
#' @return
#' @export
#'
#' @examples
survey_seg_process <- function(surveyDat, segPath, bgzPath, overlapThreshold= 60,cores=detectCores()/2){
  
  #load libraries
  dplyr::source('./2-Functions/1_connect_to_access_database.R')
  libaries <- c("data.table", "dplyr", "lwgeom", "data.table","sf", "janitor", "units", "purrr","furrr","parallel", "stringr")
  install_and_load(libaries)
  
  ## data checks ##
  # check is correct spatial format, point or polygon
  if (all(unique(as.character(st_geometry_type(surveyDat, by_geometry = T))) %in% c("POINT", "MULTIPOINT", "POLYGON","MULTIPOLYGON"))==FALSE){stop("Data needs to be either point or polygon class.")}
  #read shapefile of bgz outlines
  bgz_map <- st_read(bgzPath, quiet=T)
  # if BGZ 14 noted then join bgz 13 and 14 where required (Isle of Scilly zone 14)
  if(any(str_detect(bgz_map$Zone,'14'))){
    bgz_map <- bgz_map %>% 
      dplyr::mutate(ZoneID = c(1:12,13,13)) %>% #edit zone ids
      group_by(ZoneID) %>%
      dplyr::mutate(geometry=st_union(geometry)) %>% ungroup() %>% #join 13 and 14
      dplyr::filter(Zone != 14) %>% dplyr::select(-ZoneID)
  }
  #check crs matches
  if(st_crs(surveyDat)$input!="EPSG:27700"){ #check crs is British National Grid
    habitatDataInput <- surveyDat %>% st_transform(st_crs(bgz_map)) #if not transform to 227700
  } else {habitatDataInput = surveyDat}
  
  ## check which bgz's are included in the data
  #get list of bgzs where survey data is covering
  zones <- bgz_map %>% st_intersection(habitatDataInput) %>% suppressWarnings() %>% 
    st_drop_geometry %>% 
    dplyr::select(Zone) %>% unique() %>% 
    mutate(Zone = paste('BGZ',str_pad(Zone,2,pad='0'),sep=""))
  
  # list relevant segmentation files
  segFiles <- data.frame(files=list.files(segPath, pattern='.shp',full.names = T)) %>% 
    dplyr::filter(str_detect(files,pattern=paste(zones$Zone,collapse='|'))) %>%
    dplyr::mutate(zone = sort(zones$Zone))
  
  ## data processing ##
  # separate points and polygons
  habDataPnt <- habitatDataInput %>% dplyr::filter(st_geometry_type(habitatDataInput, by_geometry = T) %in% c("POINT", "MULTIPOINT"))
  habDataPoly <- habitatDataInput %>% dplyr::filter(st_geometry_type(habitatDataInput, by_geometry = T) %in% c("POLYGON", "MULTIPOLYGON"))
  
  #set up parallel processing
  plan(multisession, workers = cores)
  
  #### POINT DATA ####
  if (nrow(habDataPnt)>0){
    # iterate through bgzs - Merge Training Data with Segmentation
    survey_segs <-  furrr::future_map_dfr(segFiles$files, .options = furrr_options(seed = T),
                                          .f=function(bgz){
      #read in bgz data
      seg_bgz <- st_read(bgz, quiet=T) %>% filter(Class_name=='Segment')
      if(any(c('Shp_Area','Shp_Len') %in% names(seg_bgz))){
        seg_bgz <- seg_bgz %>% dplyr::select(-Shp_Area, -Shp_Len)
      }
      #get zone
      zone <- segFiles %>% dplyr::filter(files == bgz) %>% dplyr::select(zone) %>% as.character()
      # Combine training data and segments when point found within polygon
      habitatDataIDs <- habDataPnt %>% st_intersection(seg_bgz) %>% suppressWarnings() %>%
        mutate(bgz = zone) 
      #test if any points fall on boundaries (i.e. assigned to more than 1 segment)
      hab_dupes <- habitatDataIDs %>% get_dupes(Survey_ID)
      #clean up any duplicates where needed
      if(nrow(hab_dupes)>0){
        # where a point falls between two segments - take segment where point is closest to the centroid
        hab_dupes_clean <- map_df(unique(hab_dupes$Survey_ID), .f=function(dupID){
          #filter to selected survey id
          ID_select <- hab_dupes %>% dplyr::filter(Survey_ID ==dupID)
          # filter to adjacent segments and get centroid points
          adjacentPolyPoint <- seg_bgz %>% dplyr::filter(SegID %in% ID_select$SegID) %>% 
            st_centroid() %>% suppressWarnings()
          #calculate distance from point to adjacent segment centroids
          adjacentPolyPoint$distance <- st_distance(ID_select,adjacentPolyPoint, by_element = TRUE)
          #select segment with the smallest centroid distance
          closestSeg <- adjacentPolyPoint %>% dplyr::filter(distance == min(adjacentPolyPoint$distance)) 
          #select segment with smallest centroid distance
          adjacentPoly <- ID_select %>% dplyr::filter(SegID %in% closestSeg$SegID) %>% dplyr::select(-dupe_count)
        })
        #compile with clean duplicate rows
        habitatDataIDs <- habitatDataIDs %>% 
          dplyr::filter(!Survey_ID %in% hab_dupes$Survey_ID) %>%
          rbind(hab_dupes_clean)
      }
      habitatDataIDs
    }) #close iteration
    habDataPntOut <-survey_segs
    #split easting northing
    habDataPntOut <- habDataPntOut %>%
      dplyr::mutate(easting = unlist(map(habDataPntOut$geometry,1)),
             northing = unlist(map(habDataPntOut$geometry,2))) %>%
      st_drop_geometry()
  } else {habDataPntOut <- NULL} #close point 
  
  #### POLYGON DATA #### 
  if (nrow(habDataPoly)>0){
    # iterate through bgzs - Merge Training Data with Segmentation
    survey_segs <- furrr::future_map_dfr(segFiles$files, .options = furrr_options(seed = T),
                                         .f=function(bgz){
    #survey_segs <- map_df(segFiles$files, .f=function(bgz){
      #read in segmentation data for bgz 
      seg_bgz <- st_read(bgz, quiet=T) %>% dplyr::filter(Class_name=='Segment')
      #get zone name
      zone <- segFiles %>% dplyr::filter(files == bgz) %>% dplyr::select(zone) %>% as.character()
      #remove additional spatial attributes
      seg_bgz <- seg_bgz %>% dplyr::mutate(area = drop_units(st_area(seg_bgz)))
      if(any(c('Shp_Area','Shp_Len') %in% names(seg_bgz))){
        seg_bgz <- seg_bgz %>% dplyr::select(-Shp_Area, -Shp_Len)
      }
      #extract segment area of intersecting segments
      overlap_pct <- habDataPoly %>% 
        st_intersection(seg_bgz) %>% suppressWarnings() %>% #itersect 
        dplyr::mutate(intersect_area = drop_units(st_area(.))) %>%
        # Calculate overlap percentage for each segment
        dplyr::mutate(Overlap_Pct = round((intersect_area/area)*100,2)) %>% 
        st_drop_geometry()  
      
      # Filter to keep segments over used defined overlap
      segOverlap <- overlap_pct %>% dplyr::filter(Overlap_Pct >= overlapThreshold)
      #get segment centroids
      seg_centroids <- seg_bgz %>% dplyr::filter(SegID %in% overlap_pct$SegID) %>% st_centroid() %>% suppressWarnings()
      #join Seg Ids to centroid geometry 
      overlap_segs <- segOverlap %>% left_join(seg_centroids[,c('SegID','geometry')],by='SegID') %>% 
        dplyr::select(-area,-intersect_area,-Overlap_Pct) %>%
        dplyr::mutate(bgz = zone)
      overlap_segs
    })
    habDataPolyOut <- survey_segs
    #split easting northing
    habDataPolyOut <- habDataPolyOut %>%
      dplyr::mutate(easting = unlist(map(habDataPolyOut$geometry,1)),
             northing = unlist(map(habDataPolyOut$geometry,2))) %>%
      dplyr::select(-geometry)
  } else{habDataPolyOut <-NULL} #close polygon 

  #bring point and polygon data together
  survey_segs <- rbind(habDataPntOut,habDataPolyOut)
  
  #report and write out
  print("Segmentation assignment complete.")
  return(survey_segs)
} #close function
      
