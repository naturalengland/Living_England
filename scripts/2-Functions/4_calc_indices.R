## ############################## ##
##
## Script Name: indices calculator
##
## Description: Function for generating band indices from imagery
##
## Author: Becky Trippier, Natural England
## 
## Licence: MIT Licence  
## Date Created: 2024-06-03
## 
## Date Last Modified: 2024-06-03
## 
## R version 4.2.3 (2023-04-21 ucrt)
## Dependencies:
## future_1.32.0   purrr_1.0.2     
## terra_1.7-29   stringr_1.5.0
## dplyr_1.1.3 
##
## ############################## ##


#' Function for calculating spectral indicies
#'
#' @param path imagery list to calculate indices layers with
#' @param indices list of indices to calculate
#' @param outfolder folder path to save outputs to
#' @param numCores if parallel is TRUE, then number of cores to use. Default is availableCores() - 2
#'
#' @return Returns a table zonal statistics for each polygon.
#' @export
#'
#' @examples
indicesCalc <- function(imgList, indices = c('NDVI','NDWI'), outfolder = 'S2_indices/',numCores = 2, nirband = 8, redband = 4, greenband = 3) {
  # Load common functions library
  source("./2-Functions/common_functions_library.r")

  # Load in packages
  required_packages <- c("stringr", "terra", "dplyr", "future", "purrr", "tidyverse")
  install_and_load(required_packages, verbose = FALSE)

  # Check if output folder to write to exists
  if (!dir.exists("S2_indices/")) {
    stop("outfolder not found.")
  }
  
  ## Set up ##
  # Define bands numbering
  nir <- nirband
  r <- redband
  g <- greenband

  # Define indices equations
  ndvi_fun <- function(x) {
    (x[[8]] - x[[4]]) / (x[[8]] + x[[4]])
  }

  ndwi_fun <- function(x) {
    (x[[3]] - x[[8]]) / (x[[3]] + x[[8]])
  }
  
  # Create indices equation lookup
  eqLookup <- tribble(
    ~name, ~eq,
    "NDVI", ndvi_fun,
    "NDWI", ndwi_fun
  )

  # Filter to specified indices
  indices_tbl <- eqLookup %>% dplyr::filter(name %in% indices)
  
  # Create grid to run through
  inputgrid <- expand.grid(x = imgList, y = indices_tbl$name)
  
  ## iterate through files and indices calculations
  for (i in 1:nrow(inputgrid)){
    x <- as.character(inputgrid$x[i])
    y <- as.character(inputgrid$y[i])

    # Get file and break up into raster bands
    filename <- x
    filebasename <- gsub(basename(filename), pattern = ".img", replacement = "")

    # Read in raster layer
    imgrast <- terra::rast(x)

    # create folder for indices outputs
    if(!dir.exists(paste0(outfolder, y))){
      dir.create(paste0(outfolder, y))
    }

    # Get raster function
    indi_row <- indices_tbl %>% dplyr::filter(name == y)

    # Run indices function over raster
    out <- terra::app(imgrast, fun = indi_row$eq[[1]], filename = paste0(outfolder, y, "/", basename(filebasename), "_", y, ".img"), cores = numCores, overwrite = TRUE)
    
    # Report back
    print(paste(y, 'calculated for', basename(filebasename)))
    }
  
}