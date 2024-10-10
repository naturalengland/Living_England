############### Automated object-based classification using RF ##############################
#
# Title: Random Forest Classification, Evaluation and Prediction - Function calls
#
# Version: 0.1 (dev)
#
# Author: Max Fancourt. Natural England, Evidence Earth Observation Service
#
# Description: Calls the functions needed by the '4.2_Random_Forest_Classification.R' workflow. This trains a Random Forest classifcation model for each biogeographic region and evaluates performance. Then uses the random forest to predict for the rest of the biogeographic zones. 
# 
## Licence: MIT Licence
## Dat Created: 2024-01-05
#
# Last modified 29/04/2023
#
#############################################################################################

# Load in common functions library
source('2-Functions/common_functions_library.r')

check_directories_and_files <- function(primary_folder_data_path){
    install_and_load(c("Dict"), verbose = FALSE)

    cat("Checking folder structure\n")
    # Build paths to directories relative to creating subdirectories if needed
    lookup_dictonary <- Dict$new("primary_folder_data_path" = primary_folder_data_path,
    "input_data_directory" = validate_dir(paste0(primary_folder_data_path, "/input_data")), # Directory to the input folder
    "segmentation_directory" = validate_dir(paste0(primary_folder_data_path, "/input_data", "/segmentation")), # Directory to the segmentation data
    "zonal_stats_directory" = validate_dir(paste0(primary_folder_data_path, "/input_data", "/zonal_stats")), # Directory to the zonal stats data
    "training_data_directory" = validate_dir(paste0(primary_folder_data_path, "/input_data", "/training_data")), # Directory to the training data
    "output_path" = validate_dir(paste0(primary_folder_data_path, "/outputs")), # Directory to where outputs will be stored
    "model_output_path" = validate_dir(paste0(primary_folder_data_path, "/outputs", "/model_output")), # Directory to where models will be saved
    "model_evaluation_output_path" = validate_dir(paste0(primary_folder_data_path, "/outputs", "/model_evaluation")), # Directory to where model evaluation statistics will be saved
    "segment_output_path" = validate_dir(paste0(primary_folder_data_path, "/outputs", "/all_segment_predictions")), # Directory where prediction for all segments will be predicted
    "segment_output_path_model_score" = validate_dir(paste0(primary_folder_data_path, "/outputs", "/all_segment_predictions", "/model_score")), # Directory where model results will be outputted
    "segment_output_path_calibrated_results" = validate_dir(paste0(primary_folder_data_path, "/outputs", "/all_segment_predictions", "/calibrated_probability")), # Directory where calibrated probabilities will be outputted
    "segment_output_path_combined_results" = validate_dir(paste0(primary_folder_data_path, "/outputs", "/all_segment_predictions", "/combined")) # Directory where combined results will be outputted 
    )
    cat("Folder structure checking complete\n")

    # Build paths to specific input data required 
    cat("Checking that required files have been provided\n")
    if (is.null(validate_file_path(paste0(lookup_dictonary$get("training_data_directory"), "/training_data.csv")))){
        cat(crayon::red("Unable to find training data, please make sure that it is located at:\n"))
        cat(crayon::red(lookup_dictonary$get("training_data_directory"), "/training_data.csv\n"))
    } else {
      cat(crayon::green(lookup_dictonary$get("training_data_directory"), "/training_data.csv found\n"))
      lookup_dictonary$add("training_data_path" = validate_file_path(paste0(lookup_dictonary$get("training_data_directory"), "/training_data.csv"))) # Path to training data
    }
    
    if (is.null(validate_file_path(paste0(lookup_dictonary$get("primary_folder_data_path"), "/input_data/preprocessed_data.csv")))){
      cat(crayon::red("Unable to find preprocessed data, if not provided then it will be rebuilt which will increase run time\n"))
      cat(crayon::red("If already created place in", lookup_dictonary$get("training_data_directory"), "/training_data.csv\n"))
    } else {
      cat(crayon::green(lookup_dictonary$get("primary_folder_data_path"), "/input_data/preprocessed_data.csv found\n"))
      lookup_dictonary$add("preprocessed_data_path" = validate_file_path(paste0(lookup_dictonary$get("primary_folder_data_path"), "/input_data/preprocessed_data.csv"))) # Path to training data
    }
        
    if (is.null(validate_file_path(paste0(lookup_dictonary$get("output_path"), "/mtry_tuning.csv")))){
        cat(crayon::red("Unable to find mtry tuning info, if not provided then make sure to either set a manual mtry value, or else default values will be used\n"))
        cat(crayon::red("If created place in", lookup_dictonary$get("training_data_directory"), "/training_data.csv\n"))
    } else {
       lookup_dictonary$add("mtry_tuning_data_path" = validate_file_path(paste0(lookup_dictonary$get("output_path"), "mtry_tuning.csv"))) # Path to mtry test files 
    }
    cat("Required file checking completed\n")
    # Return lookup dictionary with directories
    return(lookup_dictonary)
}

SetMtry <- function(tuned_mtry_selected, mtry_manual_override) {
  if (!is.null(tuned_mtry_selected)) {
    mtries <- tuned_mtry_selected
    cat(crayon::green("Mtry set to", mtries, "\n"))
  } else if (mtry_manual_override) {
    mtries <- mtry_manual_override
    cat(crayon::green("Mtry set to", mtries, "\n"))
  } else {
    mtries <- -1
    cat(crayon::green("Mtry set to H20ai default as the square root of the number of predictors\n"))
  }
  return(mtries)
}

# Checks if java is installed and if not returns instructions

# Checks to see if fold models already exists return list of paths to fold models if found
GetFoldModels <- function(model_directory){
  list_of_cross_validation_models <- list.files(model_directory, full.names = TRUE) %>% stringr::str_subset(pattern = "_cv_")
  if (length(list_of_cross_validation_models) == 0) {
    cat(crayon::red("No cross validation models found\n"))
    return(NULL)
  }else{
    cat(crayon::green("Cross validation models found returning paths\n"))
    return(list_of_cross_validation_models)
  }
}

CreatePreprocessedData <- function(training_data_path, zonal_stats_directory, id_column_name){
  # load the training data
    training_data <- data.table::fread(training_data_path, header = TRUE)

    # load a single zonal stats file in order to get meta data
    zonal_stats_seg <- data.table::fread(normalizePath(file.path(zonal_stats_directory, paste("BGZ", "01", "_zonalStats.csv", sep = ""))), header = TRUE)
    column_names <- append(colnames(training_data), colnames(zonal_stats_seg)[!(colnames(zonal_stats_seg) %in% c(id_column_name, "V1"))])

    # Create a blank dataframe to store output
    all_data <- data.frame(setNames(data.table(matrix(nrow = 0, ncol = length(column_names))), column_names))

    for (i in str_pad(rep(1:13), 2, pad = "0")){
        cat("Starting BGZ", i, "\n")
        # create file name using loop iteration
        file_name <- paste("BGZ", i, "_zonalStats.csv", sep = "")

        # filter points to just the zone of interest 
        training_data_subset <- training_data %>% dplyr::filter(bgz == paste0("BGZ", i))

        # read in zonal stats for zone 
        zonal_stats_seg <- data.table::fread(normalizePath(file.path(zonal_stats_directory, file_name)), header = TRUE)

        # left join training points with zonal stats, many to many as some segments will have multiple classes
        to_bind <- training_data_subset %>% dplyr::left_join(zonal_stats_seg, by=id_column_name, relationship = "many-to-many")
        
        all_data <- rbind(all_data, to_bind)
        cat("Finishing BGZ", i,"\n")
    }

    # ensure that all factor variables are factors
    all_data <- all_data %>%
      rowwise() %>%
      mutate(across(all_of(list_of_factor_variables), as.factor))
      print("factorisation of variables complete")

    return(all_data)
}

# Function to identify and remove highly correlated variables
RemoveHighlyCorrelatedVariables = function(cor_data, training_data_processed, correlation_threshold, output_path){
  cor_matrix <- h2o::cor(cor_data, use = "complete.obs")
  highly_correlated <- caret::findCorrelation(cor_matrix, cutoff = correlation_threshold, verbose = FALSE, exact = TRUE)

  # Return the names of the highly correlated variables to be removed
  highly_correlated_variables <- colnames(cor_data %>% dplyr::select(all_of(highly_correlated)))

  # Save out the variables to be removed 
  write.csv(as.data.frame(highly_correlated_variables), file.path(output_path, "highly_correlated_variables.csv"))

  # Remove highly correlated variables
  return(training_data_processed %>% dplyr::select(-all_of(highly_correlated_variables)))
}

# Calculate mode of a vector
Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

# Retrieve optimal mtry from file
GetOptimalMtryFromTesting <- function(path_to_mtry_testing){
  if (is.null(path_to_mtry_testing)) {
    cat(crayon::red("Mtry has not yet been tuned, either run tuning code to find optimal, ensure that a mtry_manual_override is set, or else defaults will be used\n"))
    return(NULL)
  } else {
    # Open mtry output file 
    mtry_data <- read.csv(path_to_mtry_testing)
    head(mtry_data)
    # Find mtry optima for
    logloss_mtry_optimum <- mtry_data %>% dplyr::arrange("logloss") %>% dplyr::slice(1) %>% dplyr::pull("run_id") # logloss
    mse_mtry_optimum <- mtry_data %>% dplyr::arrange("MSE") %>% dplyr::slice(1) %>% dplyr::pull("run_id")# mse
    rmse_mtry_optimum <- mtry_data %>% dplyr::arrange("RMSE") %>% dplyr::slice(1) %>% dplyr::pull("run_id")# rmse

    # If more than one method agrees, then chose the most common response, if not then select logloss variables 
    if (length(unique(c(logloss_mtry_optimum, mse_mtry_optimum, rmse_mtry_optimum))) != 3){
        tuned_mtry_selected <- Mode(c(logloss_mtry_optimum, mse_mtry_optimum, rmse_mtry_optimum))
    } else {
        tuned_mtry_selected <- logloss_mtry_optimum
    }
    return(tuned_mtry_selected)
    cat(crayon::green("Mtry tuned and optimal found to be:", tuned_mtry_selected))
  }
}

# Remove duplicate segments
RemoveDuplicateSegments <- function(full_dataset, fold_validation_predictions) {

}

