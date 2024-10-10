## Living England translate survey data function 
##
## Description: Functions for cleaning and translating a habitat survey data using the MS Access Habitat translations database
##
## Authors: Becky Trippier, Natural England
##
## Licence: MIT Licence
## Date created: 2022-09-20
##
## Date Last Modified: 2024-10-02
## 
## Versioning:
## R version 4.2.3 (2023-03-15 ucrt)
## Dependencies:
## odbc_1.3.5    stringr_1.5.0 readxl_1.4.2  units_0.8-1  
## sf_1.0-12     dplyr_1.1.3   DBI_1.1.3 
#---------------------------------------------------

#' Extract and join translation function
#'
#' @param surveyData data frame to be translated
#' @param from_habs_table Habitat table with source classes
#' @param to_habs_table Habitat table with target classes
#' @param from_is_source logical, if source scheme is the scheme_id_fromin the relationships table set to "T", or if it is listed as the to_ID then put as "F".
#' @param db_path filepath to database
#'
#' @return
#' @export
#'
#' @examples
  extract_translation <- function(surveyData=habDat,
                                 from_habs_table = Schemehabs,
                                 to_habs_table=targetSchemehabs,
                                 from_is_source = T,
                                 db_path){
    #load  packages
    library(DBI)
    library(dplyr)
    
    #connect to the habitat translations database
    con<- connect_to_access_dbi(db_path)
    
    #if data is listed as the source scheme (check extracted pathway or habitat translation database relationships table)
    if(from_is_source == T){ 
      #pull data from database
      relateTable <- dbGetQuery(con,paste0("SELECT * from Habitat_Relationships WHERE habitat_id_from IN (SELECT habitat_id FROM Habitat_Class WHERE scheme_id = ", unique(from_habs_table$scheme_id), ")")) %>% 
        dplyr::filter(habitat_id_to %in% to_habs_table[["habitat_id"]])
      #join the data together from the scheme tables
      habDatTranslated <- surveyData %>% 
        dplyr::select(Survey_Code, Survey_Source, Source_ID, Source_Date, Survey_Date,
               habitat_name = Source_Broad_Habitat,
               Source_Detailed_Habitat,Source_Hab_Notes,Survey_ID, habitat_id, scheme_id) %>% 
        left_join(relateTable, by = c('habitat_id' = 'habitat_id_from')) %>% 
        left_join(to_habs_table, by = c('habitat_id_to' ='habitat_id')) %>%
        dplyr::select(Survey_ID, Source_ID, Survey_Code, Survey_Source, Source_Date, Survey_Date, 
               Source_Scheme_ID = scheme_id.x, Source_Habitat_ID = habitat_id,
               Source_Broad_Name = habitat_name.x, Source_Detailed_Habitat, Source_Hab_Notes,
               Target_Scheme_ID = scheme_id.y, Target_Habitat_ID = habitat_id_to,
               Target_Broad_Name = habitat_name.y,Target_Habitat_Code = hab_code, geometry)
    } else { #if data is listed as the to scheme
      #pull data from database
      relateTable <- dbGetQuery(con,paste0("SELECT * from Habitat_Relationships WHERE habitat_id_from IN (SELECT habitat_id from Habitat_Class WHERE scheme_id = ", unique(to_habs_table$scheme_id), ")")) %>% 
        filter(habitat_id_to %in% from_habs_table[["habitat_id"]])
      #join the data together from the scheme tables
      habDatTranslated <- surveyData %>% 
        dplyr::select(Survey_Code,Survey_Source, Source_ID, Source_Date, Survey_Date,
                      habitat_name = Source_Broad_Habitat,
               Source_Detailed_Habitat,Source_Hab_Notes,Survey_ID, habitat_id, scheme_id)  %>%
        left_join(relateTable, by = c('habitat_id' = 'habitat_id_to'),relationship="many-to-many") %>% 
        left_join(to_habs_table, by = c('habitat_id_from' ='habitat_id'),relationship="many-to-many") %>%
        dplyr::select(Survey_ID, Source_ID, Survey_Code,Survey_Source, Source_Date, Survey_Date, 
               Source_Scheme_ID = scheme_id.x,Source_Habitat_ID = habitat_id,
               Source_Broad_Name = habitat_name.x,Source_Detailed_Habitat, Source_Hab_Notes,
               Target_Scheme_ID = scheme_id.y, Target_Habitat_ID = habitat_id_from,
               Target_Broad_Name = habitat_name.y, Target_Habitat_Code = hab_code,geometry)
    }  
    # return out the translated data
    return(habDatTranslated)
  }
  
  ## FUNCTION - Translate across supplied translation pathway
  translate_classes <- function(cleanDat, hab_codes, db_path, sourceName, targetname, pathway){
    ## Extract source habitat classes ##
    #connect to the habitat translations database
    con<- connect_to_access_dbi(db_path)
    
    # extract scheme data from Habitat_Scheme table associated with the classification
    ## SQL: select all from the Habitat_Scheme table where the scheme_name_short equals sourceName
    Scheme <- dbGetQuery(con,paste0("SELECT * from Habitat_Scheme WHERE scheme_name_short = '",sourceName ,"'"))
    
    # extract habitat classes associated with the habitat scheme
    ## SQL: select all from the Habitat_Class table where scheme_id equals the scheme id from above query
    Schemehabs <- dbGetQuery(con,paste0("SELECT * from Habitat_Class WHERE scheme_id = ",Scheme$scheme_id)) %>%
      dplyr::mutate(Hab_Code_lower = tolower(as.character(hab_code)), 
             Name_lower = tolower(habitat_name))
    
    # join to cleaned data to pull in db class names and ids
    ## if hab_codes provided join with hab codes
    if(hab_codes==T){
      habDat <- cleanDat %>% 
        dplyr::mutate(Source_Broad_Habitat_lower = tolower(as.character(Source_Broad_Habitat))) %>% 
        left_join(Schemehabs, by=c('Source_Broad_Habitat_lower' = 'Hab_Code_lower'), keep=T) %>% 
        dplyr::select(-Source_Broad_Habitat_lower, -Hab_Code_lower)
    } else{
      habDat <- cleanDat %>% 
        dplyr::mutate(Source_Broad_Habitat_lower = tolower(Source_Broad_Habitat)) %>% 
        left_join(Schemehabs, by=c('Source_Broad_Habitat_lower' = 'Name_lower')) %>% 
        dplyr::select(-Source_Broad_Habitat_lower, -Hab_Code_lower)
    }
    
    ## Extract target habitat classes ##
    # extract target scheme data
    ## SQL: select all from the Habitat_Scheme table where the scheme_name_short equals targetname
    targetScheme <- dbGetQuery(con,paste0("SELECT * from Habitat_Scheme WHERE scheme_name_short = '",targetname ,"'"))
    
    # extract target habitat classes associated with the target habitat scheme
    ## SQL: select all from the Habitat_Class table where scheme_id equals the scheme id from above query
    targetSchemehabs <- dbGetQuery(con,paste0("SELECT * from Habitat_Class WHERE scheme_id = ",targetScheme$scheme_id)) %>% 
      dplyr::select(habitat_id:hab_code)
    ## Translate source classes to target classes ##
    if (nrow(pathway)==1){ 
      #### Direct Translation ####
      # Source scheme = Scheme$scheme_id , Target Scheme = targetScheme$scheme_id
      if(pathway$scheme_id_from==Scheme$scheme_id){ # if the habitat_id_from  = Source scheme
        habDatTranslated <- extract_translation(surveyData = habDat,
                                               from_habs_table = Schemehabs,
                                               to_habs_table=targetSchemehabs,
                                               from_is_source = T,
                                               db_path=db_path)
      } else{ # if the habitat_id_from = target scheme
        habDatTranslated <- extract_translation(surveyData = habDat,
                                               from_habs_table = Schemehabs,
                                               to_habs_table=targetSchemehabs,
                                               from_is_source = F,
                                               db_path=db_path)
      }
      print("Direct translation.")
      #check if translation relationship is 1:1 and print warning if not
      duple <- habDatTranslated %>% group_by(Survey_ID) %>% filter(n()>1)
      if(nrow(duple)>0){print('contains duplicate data points, relationship not 1:1.')}
      #disconnect from database
      dbDisconnect(con)
      return(habDatTranslated)
      
    } else { 
      
      #### indirect translation ####
      # first translation - source to first neighbour
      sourceRel <-pathway %>% dplyr::filter(scheme_name_from==sourceName|scheme_name_to ==sourceName)
      if(sourceRel$scheme_id_from ==Scheme$scheme_id){ # if the scheme_id_from = Source scheme
        intHabs <- dbGetQuery(con,paste0("SELECT * from Habitat_Class WHERE scheme_id = ",sourceRel$scheme_id_to)) %>% 
          dplyr::select(habitat_id:hab_code)
        habDatTranslated <- extract_translation(surveyData = habDat,
                                               from_habs_table = Schemehabs,
                                               to_habs_table=intHabs,
                                               from_is_source = T,
                                               db_path=db_path)
      } else{
        intHabs <- dbGetQuery(con,paste0("SELECT * from Habitat_Class WHERE scheme_id = ",sourceRel$scheme_id_from)) %>% dplyr::select(habitat_id:hab_code)
        habDatTranslated <- extract_translation(surveyData = habDat,
                                               from_habs_table = Schemehabs,
                                               to_habs_table=intHabs,
                                               from_is_source = F,
                                               db_path=db_path)
      }
      if(nrow(pathway)==2){ 
        #### translation via 1 neighbour ####
        #second translation
        intPrep <- habDatTranslated %>% dplyr::select(Survey_ID:Survey_Date, scheme_id = Target_Scheme_ID, habitat_id =Target_Habitat_ID, Source_Broad_Habitat = Target_Broad_Name, Source_Detailed_Habitat, Source_Hab_Notes,geometry)
        secondRel <-pathway %>% filter(log_id!= sourceRel$log_id)
        if(secondRel$scheme_id_from==targetScheme$scheme_id){ # if the scheme_id_from = target scheme
          intTranslated <- extract_translation(surveyData = intPrep,
                                              from_habs_table = intHabs,
                                              to_habs_table=targetSchemehabs,
                                              from_is_source = F,
                                              db_path=db_path) %>% 
            st_drop_geometry()
        } else{
          intTranslated <- extract_translation(surveyData = intPrep,
                                              from_habs_table = intHabs,
                                              to_habs_table=targetSchemehabs,
                                              from_is_source = T,
                                              db_path=db_path) %>% 
            st_drop_geometry()
        }
        intOut <-  intTranslated %>% 
          dplyr::select(Survey_ID,Source_Scheme_ID:Source_Broad_Name, 
                        Target_Scheme_ID:Target_Habitat_Code) %>%
          rename(Int_Scheme_ID =Source_Scheme_ID,
                 Int_Habitat_ID=Source_Habitat_ID,
                 Int_Broad_Name = Source_Broad_Name)
        habTranslatedOut <- habDatTranslated %>% 
          rename(Int_Scheme_ID =Target_Scheme_ID,
                 Int_Habitat_ID=Target_Habitat_ID,
                 Int_Broad_Name = Target_Broad_Name) %>% 
          left_join(intOut, by=c('Survey_ID', "Int_Scheme_ID", "Int_Habitat_ID", "Int_Broad_Name"))
        print(paste("Translate via Scheme:", unique(intHabs$scheme_id)))
        #check if translation relationship is 1:1 and print warning if not
        duple <- habTranslatedOut %>% group_by(Survey_ID) %>% filter(n()>1)
        if(nrow(duple)>0){print('contains duplicate data points, relationship not 1:1.')}
        #disconnect from database
        dbDisconnect(con)
        return(habTranslatedOut)
      } else { 
        #### translation via 2 neighbours ####
        #prepare neighbour 1 data for translation
        intPrep <- habDatTranslated %>% dplyr::select(Survey_ID:Survey_Date, scheme_id = Target_Scheme_ID, habitat_id =Target_Habitat_ID, Source_Broad_Habitat = Target_Broad_Name, Source_Detailed_Habitat, Source_Hab_Notes,geometry)
        ##find first to second neighbour relationship
        # filter relationships where neighbour 1 is present but not the source scheme ID
        nhbr_nhbr <- pathway %>% filter(scheme_id_from == unique(intHabs$scheme_id) & 
                                          scheme_id_to != unique(Schemehabs$scheme_id)|
                                          scheme_id_from != unique(Schemehabs$scheme_id) & 
                                          scheme_id_to == unique(intHabs$scheme_id))
        # translate neighbour 1 to neighbour 2
        if(nhbr_nhbr$scheme_id_from== unique(intHabs$scheme_id)){ #where neighbour scheme is the from_ID
          intSecHabs <- dbGetQuery(con,paste0("SELECT * from Habitat_Class WHERE scheme_id = ",nhbr_nhbr$scheme_id_to)) %>% 
            dplyr::select(habitat_id:hab_code)
          intTranslated <- extract_translation(surveyData = intPrep,
                                              from_habs_table = intHabs,
                                              to_habs_table=intSecHabs,
                                              from_is_source = T,
                                              db_path=db_path)
        } else { #where neighbour scheme is the to_ID
          intSecHabs <- dbGetQuery(con,paste0("SELECT * from Habitat_Class WHERE scheme_id = ",nhbr_nhbr$scheme_id_from )) %>% 
            dplyr::select(habitat_id:hab_code)
          intTranslated <- extract_translation(surveyData = intPrep,
                                              from_habs_table = intHabs,
                                              to_habs_table=intSecHabs,
                                              from_is_source = F,
                                              db_path=db_path) 
        }
        #prepare neighbour 2 data for translation
        intSecPrep <- intTranslated %>% dplyr::select(Survey_ID:Survey_Date, scheme_id = Target_Scheme_ID, habitat_id =Target_Habitat_ID, Source_Broad_Habitat = Target_Broad_Name, Source_Detailed_Habitat,Source_Hab_Notes)
        ##find second neighbour to target relationship
        # filter relationships where neighbour 1 is present but not the source scheme ID
        nhbr_target <- pathway %>% filter(scheme_id_from== unique(intSecHabs$scheme_id) & 
                                            scheme_id_to == unique(targetSchemehabs$scheme_id)|
                                            scheme_id_from == unique(targetSchemehabs$scheme_id) & 
                                            scheme_id_to == unique(intSecHabs$scheme_id))
        # translate neighbour 1 to neighbour 2
        if(nhbr_target$scheme_id_from== unique(intSecHabs$scheme_id)){ #where neighbour scheme is the from_ID
          targetTranslated <- extract_translation(surveyData = intSecPrep,
                                                 from_habs_table = intSecHabs,
                                                 to_habs_table=targetSchemehabs,
                                                 from_is_source = T,
                                                 db_path=db_path) 
        } else { #where neighbour scheme is the to_ID
          targetTranslated <- extract_translation(surveyData = intSecPrep,
                                                 from_habs_table = intSecHabs,
                                                 to_habs_table=targetSchemehabs,
                                                 from_is_source = F,
                                                 db_path=db_path) 
        }
        
        #habDatTranslated #source - nhbr1
        #intTranslated # nhbr1 - nhbr 2
        #targetTranslated # nhbr 2 - target
        
        #rename translated columns
        int1Out <-  intTranslated %>% st_drop_geometry() %>%
          dplyr::select(Survey_ID, 
                        Int_1_Scheme_ID =Source_Scheme_ID,
                        Int_1_Habitat_ID=Source_Habitat_ID,
                        Int_1_Broad_Name = Source_Broad_Name,
                        Int_2_Scheme_ID =Target_Scheme_ID,
                        Int_2_Habitat_ID=Target_Habitat_ID,
                        Int_2_Broad_Name = Target_Broad_Name) 
        int2Out <-  targetTranslated %>% st_drop_geometry() %>%
          dplyr::select(Survey_ID, 
                        Int_2_Scheme_ID =Source_Scheme_ID,
                        Int_2_Habitat_ID=Source_Habitat_ID,
                        Int_2_Broad_Name = Source_Broad_Name,
                        Target_Scheme_ID:Target_Habitat_Code) 
        # join translated data
        habTranslatedOut <- habDatTranslated %>% 
          rename(Int_1_Scheme_ID =Target_Scheme_ID,
                 Int_1_Habitat_ID=Target_Habitat_ID,
                 Int_1_Broad_Name = Target_Broad_Name) %>% 
          left_join(int1Out, by=c('Survey_ID', "Int_1_Scheme_ID", "Int_1_Habitat_ID", "Int_1_Broad_Name")) %>%
          left_join(int2Out, by=c('Survey_ID', "Int_2_Scheme_ID", "Int_2_Habitat_ID", "Int_2_Broad_Name"))
        
        #out
        print(paste("Translate via 2 Schemes:", unique(intHabs$Scheme_ID),",", 
                    unique(intSecHabs$Scheme_ID)))
        #check if translation relationship is 1:1 and print warning if not
        duple <- habTranslatedOut %>% group_by(Survey_ID) %>% filter(n()>1)
        if(nrow(duple)>0){print('contains duplicate data points, relationship not 1:1.')}
        #disconnect from database
        dbDisconnect(con)
        return(habTranslatedOut)
      }
      
    }
  }
  
  direct_assignment <- function(cleanDat,sourceName,hab_codes,db_path){
    #pull in scheme data 
    #connect to the habitat translations database
    con<- connect_to_access_dbi(db_path)
    # extract scheme data from Habitat_Scheme table associated with the classification
    ## SQL: select all from the Habitat_Scheme table where the scheme_name_short equals sourceName
    Scheme <- dbGetQuery(con,paste0("SELECT * from Habitat_Scheme WHERE scheme_name_short = '",sourceName ,"'"))
    
    # extract habitat classes associated with the habitat scheme
    ## SQL: select all from the Habitat_Class table where scheme_id equals the scheme id from above query
    Schemehabs <- dbGetQuery(con,paste0("SELECT * from Habitat_Class WHERE scheme_id = ",Scheme$scheme_id)) %>%
      mutate(Hab_Code_lower = tolower(as.character(hab_code)), 
             Name_lower = tolower(habitat_name))
    
    # join to cleaned data to pull in db class names and ids
    ## if hab_codes provided join with hab codes
    if(hab_codes==T){
      habDat <- cleanDat %>% 
        mutate(Source_Broad_Habitat_lower = tolower(as.character(Source_Broad_Habitat))) %>% 
        left_join(Schemehabs, by=c('Source_Broad_Habitat_lower' = 'Hab_Code_lower'), keep=T) %>% 
        dplyr::select(-Source_Broad_Habitat_lower, -Hab_Code_lower)
    } else{
      habDat <- cleanDat %>% 
        mutate(Source_Broad_Habitat_lower = tolower(Source_Broad_Habitat)) %>% 
        left_join(Schemehabs, by=c('Source_Broad_Habitat_lower' = 'Name_lower')) %>% 
        dplyr::select(-Source_Broad_Habitat_lower, -Hab_Code_lower)
    }
    assignIDsOut <- habDat %>% mutate(Target_Scheme_ID=scheme_id,Target_Habitat_ID=habitat_id,Target_Broad_Name=habitat_name,Target_Habitat_Code=hab_code) %>%
      dplyr::select(Survey_ID,Source_ID,Survey_Code,Survey_Source,Source_Date,Survey_Date,Source_Scheme_ID = scheme_id,Source_Habitat_ID=habitat_id,Source_Broad_Name=habitat_name,Source_Detailed_Habitat,Source_Hab_Notes,Target_Scheme_ID,Target_Habitat_ID,Target_Broad_Name,Target_Habitat_Code)
    dbDisconnect(con)
    return(assignIDsOut)
  }
  
#' Clean and Translate survey into LE broad habitats
#'
#' @param surveyDat survey dataframe
#' @param Survey_Source 
#' @param Source_ID ID field within the source data per polygon/point 
#' @param Source_Date date data was provided
#' @param Survey_Date field name with the survey date
#' @param Source_Habitat_Classification habitat classification name, corresponding to the database 'Habitat_Scheme' table field 'Scheme_name_short' the data is provided in
#' @param Source_Broad_Habitat field name containing the broad habitat classes, see hab_codes below.
#' @param Source_Detailed_Habitat field name corresponding to a more detailed habitat class or extra information to store with the record. This will not be translated but retained for detail.
#' @param db_path file path to habitat translations database file
#' @param hab_codes logical, if Source_Broad_Habitat corresponds to the 'Hab_Codes' field of the 'Habitat_Class' table in the database. Default is False and will assume the data corresponds to the 'Name' field.
#' @param Target_Habitat_Classification habitat classification name, corresponding to the database 'Habitat_Scheme' table field 'Scheme_name_short' you wish to translate to
#'
#' @return
#' @export
#'
#' @examples
  translate_survey <- function(surveyDat,
                              Survey_Code,
                              Survey_Source,
                              Source_ID,
                              Source_Date = "2019-11-12",
                              Survey_Date,
                              Source_Broad_Habitat,
                              Source_Detailed_Habitat = NA,
                              Source_Hab_Notes = NA,
                              db_path,
                              hab_codes = F,
                              sourceName = 'LE_alias',
                              targetname  = 'LE_alias',
                              pathway,
                              detTranslate = F,
                              detSourceName=sourceName,
                              detTarget = 'LE_UKBAPpriority',
                              detPath=NA){
    
    ## Set up and Checks ####
    # load packages
    source('./2-Functions/1_find_translation_path.R')
    libaries <- c("sf", "dplyr", "units", "readxl","DBI", "stringr", "odbc")
    install_and_load(libaries)
    
    ## Add check for essential fields ##
    if(any(is.na(Survey_Code), is.na(Survey_Source), is.na(Source_ID), is.na(Source_Date), is.na(Survey_Date), is.na(Source_Broad_Habitat))){
      stop("Essential data missing from survey, check inputs")
    }
    ## Add check for detailed hab fields ##
    if(detTranslate==T & any(is.na(detTarget), is.na(detPath))){
      stop("Detailed hab data missing from survey, check inputs")
    }
    
    ## Clean Dataset and add data fields####
    # collate required data fields
    cleanDat <- surveyDat
    essenFields<-data.frame('Survey_Code'=Survey_Code,'Survey_Source'=Survey_Source,'Source_ID'=Source_ID,'Source_Date'=Source_Date,'Survey_Date'=Survey_Date,'Source_Broad_Habitat'=Source_Broad_Habitat)
    for(i in names(essenFields)){
      fieldName<- i
      datName <- essenFields[,i]
      if(datName %in% names(surveyDat)){
        cleanDat <- cleanDat %>% dplyr::mutate(!!fieldName := get(!!datName))
      } else{cleanDat <- cleanDat %>% dplyr::mutate(!!fieldName := datName)}
    }
    #create detailed hab field - handle NAs
    if(is.na(Source_Detailed_Habitat)){
      cleanDat <- cleanDat %>% dplyr::mutate(Source_Detailed_Habitat = NA)
    }else{
      cleanDat <- cleanDat %>% dplyr::mutate(Source_Detailed_Habitat = get(!!Source_Detailed_Habitat))
    }
    #create detailed hab field - handle NAs
    if(is.na(Source_Hab_Notes)){
      cleanDat <- cleanDat %>% dplyr::mutate(Source_Hab_Notes = NA)
    }else{
      cleanDat <- cleanDat %>% dplyr::mutate(Source_Hab_Notes = get(!!Source_Hab_Notes))
    }
    cleanDat <- cleanDat %>% dplyr::select(Survey_Code,Survey_Source,Source_ID,Source_Date,Survey_Date,Source_Broad_Habitat,Source_Detailed_Habitat,Source_Hab_Notes) 
    cleanDat$Survey_ID <- 1:nrow(cleanDat) # add unique ID per row, row_number() not working as expected for some
    ## Add check if target classification == source classification
    if(targetname==sourceName){
      broadOut<- direct_assignment(cleanDat,sourceName=sourceName,hab_codes,db_path)
    } else {
      ## translate broad habitat classes via supplied pathway
      broadOut <- translate_classes(cleanDat, hab_codes, db_path, sourceName, targetname, pathway)
    }
    
    # Optional - detailed habitat translation
    if(detTranslate==T){
      detailedHabDat <- cleanDat %>% dplyr::mutate(Source_Broad_Habitat=Source_Detailed_Habitat)
      if(detTarget==detSourceName){
        detailedOut<- direct_assignment(detailedHabDat,sourceName=detSourceName,hab_codes,db_path)
      } else{
        detailedOut <- translate_classes(detailedHabDat, hab_codes, db_path, 
                                        sourceName=detSourceName, targetname=detTarget, pathway=detPath)
        detailedOut<- detailedOut %>% st_drop_geometry() %>% 
          dplyr::select(Survey_ID,Source_Scheme_ID,Source_Habitat_ID, Source_Detailed_Habitat,
                        Source_Hab_Notes:Target_Broad_Name)
      }
      return(list(broad=broadOut,priority=detailedOut))
    } else{
      return(broadOut)
    }
    
  }


