## Function to find neighbouring translations to translate between classifications. 
##
## Description: Function for finding the translation path between two habitat classifications using the MS Access Habitat translations database. This will return a list of possible translation relationships.
##
## Authors: Becky Trippier, Natural England
##
## Licence: MIT Licence
## Date created: 2022-09-26
##
## Date Last Modified: 2024-10-02
## 
## Versioning:
## R version 4.2.3 (2023-03-15 ucrt)
## Dependencies:
## dplyr_1.1.3 
## DBI_1.1.3  
#---------------------------------------------------

#' Finding the relationship path to translate between habitat classifications
#'
#' @param sourceName Scheme_ID you are translating from
#' @param targetname Scheme_ID you are translating to
#' @param db_path database file path for habitat translations database
#'
#' @return
#' @export
#'
#' @examples
translation_path <- function(sourceName, targetname,db_path){
  
  #load libraries
  require(DBI)
  library(dplyr)
  
  #load functions
  source('./2-Functions/connect_to_access_database.R')

  #make a connection
  #connect to the habitat translations database
  if(file.exists(db_path )==F){stop("Unable to find db_path. Please check.")}
  con<- connect_to_access_dbi(db_path )
  
  #find source and target IDs
  Schemes <- dbGetQuery(con,paste0("SELECT * from Habitat_Scheme"))
  sourceID <- Schemes %>% dplyr::filter(scheme_name_short==sourceName) %>% dplyr::select(scheme_id) %>% as.numeric()
  targetID <- Schemes %>% dplyr::filter(scheme_name_short==targetname) %>% dplyr::select(scheme_id) %>% as.numeric()
  
  #read in relevant table
  relationshipsAll <- dbGetQuery(con,paste0("SELECT * from Relationships_log"))

  # does direct relationship exist? 
  #filter relationships to where both target and source present
  direct <- relationshipsAll %>% dplyr::filter(scheme_id_from==sourceID & scheme_id_to==targetID | 
                                                 scheme_id_from==targetID & scheme_id_to==sourceID)
  # if direct relationship - stop and return relationship
  if(nrow(direct)>0){
    print("Target found! direct relationship available.")
    return(direct)
  }
  
  #if no direct relationship...
  # find source scheme relationships
  sourceRels <- relationshipsAll %>% dplyr::filter(scheme_id_from==sourceID | scheme_id_to==sourceID)
  #find neighbour schemes
  neighbours <- c(sourceRels$scheme_id_from[sourceRels$scheme_id_from!=sourceID],
                  sourceRels$scheme_id_to[sourceRels$scheme_id_to!=sourceID])
  # get source relationships
  nhbr1_rels <- relationshipsAll %>% dplyr::filter(scheme_id_from %in% neighbours | scheme_id_to %in% neighbours)
  # get neighbours of source neighbours (1st neighbours)
  nhbr1_neighbours <- unique(c(nhbr1_rels$scheme_id_from, nhbr1_rels$scheme_id_to))
  
  #is target in 1st neighbours?
  if(targetID %in% nhbr1_neighbours){ 
    # get target to neighbour relationship
    targetRel <- nhbr1_rels %>% dplyr::filter(scheme_id_from==targetID | scheme_id_to==targetID)
    #iterate through listing the relationship combinations
    relshpOut <- list()
    for (j in targetRel$log_id){
      #find selected relationship
      nhbr_target <- targetRel %>% dplyr::filter(log_id==j)
      # find good neighbour
      if(nhbr_target$scheme_id_from==targetID){
        good_nhbr <- nhbr_target %>% dplyr::select(int_ID=scheme_id_to,int_name=scheme_name_to)
        } else{
          good_nhbr <- nhbr_target %>% dplyr::select(int_ID=scheme_id_from,int_name=scheme_name_from)
          }
      # get good neighbour to source relationship
      nhbr1_source <- sourceRels %>% dplyr::filter(scheme_id_from==good_nhbr$int_ID | 
                                              scheme_id_to==good_nhbr$int_ID )
    #join all relationships and list
    all_df <- rbind(nhbr1_source, nhbr_target)
    relshpOut <- append(relshpOut, list(all_df))
    }
    print(paste("Target found! 1 step neighbour"))
    return(relshpOut)
  }
  
  # no 1st neighbour relationship...
  # get neighbours of 1st neighbours relationships
  nhbr2_rels <- relationshipsAll %>% dplyr::filter(scheme_id_from %in% nhbr1_neighbours | scheme_id_to %in% nhbr1_neighbours)
  # get al 2nd neighbours
  nhbr2_neighbours <- unique(c(nhbr2_rels$scheme_id_from, nhbr2_rels$scheme_id_to))
    
  #is target in 2nd neighbours?
  if(targetID %in% nhbr2_neighbours){ 

    # get target to 2nd neighbour relationship
    target_nhbr2 <- nhbr2_rels %>% dplyr::filter(scheme_id_from==targetID | scheme_id_to==targetID)
    # find good 2nd neighbour
    from_target <- target_nhbr2 %>% dplyr::select(int_ID=scheme_id_to,int_name=scheme_name_to)
    good_nhbr_2 <- target_nhbr2 %>% dplyr::select(int_ID=scheme_id_to,int_name=scheme_name_to) %>% 
      rbind(from_target)%>% dplyr::filter(int_ID!=targetID) %>% unique()
    
    #get source to 1st neighbour relationship
    source_nhbr1 <- nhbr1_rels %>% dplyr::filter(scheme_id_from==sourceID | scheme_id_to==sourceID)
    # find good 2nd neighbour
    from_source <- source_nhbr1 %>% dplyr::select(int_ID=scheme_id_from,int_name=scheme_name_from)
    good_nhbr_1 <- source_nhbr1 %>% dplyr::select(int_ID=scheme_id_to,int_name=scheme_name_to) %>% 
      rbind(from_source)%>% dplyr::filter(int_ID!=sourceID) %>% unique()
    
    #find neighbour relationship
    nhbr_nhbr <- nhbr2_rels %>% filter(scheme_id_from %in% good_nhbr_1$int_ID & scheme_id_to %in% good_nhbr_2$int_ID |
                                       scheme_id_from %in% good_nhbr_2$int_ID & scheme_id_to %in% good_nhbr_1$int_ID)
    
    relshpOut <- list()
    #iterate through listing the relationship combinations
    for (i in nhbr_nhbr$log_id){
      #filter to relationship
      log_rel <- nhbr_nhbr %>% dplyr::filter(log_id==i)
      #find source to nhbr 1 relationship
      nhbr1_df <- source_nhbr1 %>% dplyr::filter(scheme_id_from==log_rel$scheme_id_from | scheme_id_to==log_rel$scheme_id_from| scheme_id_from==log_rel$scheme_id_to| scheme_id_to==log_rel$scheme_id_to)
      #find source to nhbr 1 relationship
      nhbr2_df <- target_nhbr2 %>% dplyr::filter(scheme_id_from==log_rel$scheme_id_from | scheme_id_to==log_rel$scheme_id_from| scheme_id_from==log_rel$scheme_id_to| scheme_id_to==log_rel$scheme_id_to)
      #join all relationships and list
      all_df <- rbind(log_rel, nhbr1_df, nhbr2_df)
      relshpOut <- append(relshpOut, list(all_df))
      }
    print(paste("Target found! 2 step neighbours."))
    return(relshpOut)
  } else{
    stop("No relationship found within two intermediate classifications.")
  }
  
  }
  
    
  
