## Living England translate survey data function
##
## Description: Function for establishing a connection with the MS Access Habitat translations database
##
## Authors: Becky Trippier, Natural England
##
## Licence: MIT Licence
## Date created: 2022-09-20
##
## Date Last Modified: 2024-08-13
## 
## Versioning:
## R version 4.2.3 (2023-04-21 ucrt)
## Dependencies:
## odbc_1.3.5 
## DBI_1.1.3 
#---------------------------------------------------

## Functions

#' Create connection string to connect to the habitat translation database
#'
#' @param db_file_path file path to the MS Access database
#'
#' @return
#' @export
#'
#' @examples
connect_to_access_dbi <- function(db_file_path)  {
  #load packages
  require(DBI)
  require(odbc)
  #check file exists
  if (!file.exists(db_file_path)) {
    stop("DB file does not exist at ", db_file_path)
  }
  # Assemble connection strings
  dbq_string <- paste0("DBQ=", db_file_path)
  driver_string <- "Driver={Microsoft Access Driver (*.mdb, *.accdb)};"
  db_connect_string <- paste0(driver_string, dbq_string)
  #connect to database
  myconn <- dbConnect(odbc::odbc(), .connection_string = db_connect_string)
  #return connection
  return(myconn)
}
