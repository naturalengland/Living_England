###############################
#
# Script Name: Model Based Classification
#
# Description: Trains a Random Forest classification model using data from all zones and evaluates performance. 
# This uses a distributional random forest, and codes in missing data using a missing in attribute methodology 
# that allows for prediction where data is missing.
#
# Author: Max Fancourt, Natural England
#
# Licence: MIT Licence
# Date Created: 2024-01-05
#
# Date Last Modified: 26/06/2024
#
# Versioning: 1.0
# R version 4.3.0
# Java version: Java jdk 17
# Dependencies:
# dplyr        v.1.1.2
# caret        v.6.0-94
# data.table   v.1.14.8
# h20          v.3.44.0.2
# stringr      v.1.5.0
# probably     v.1.0.2
# matrixStats  v.1.2.0
# lattice      v.0.21-8
# ggplot2      v.3.5.0
#
###############################

# 1. Setup environment ####
# Load custom functions
source('2-Functions/common_functions_library.r')
source('2-Functions/5_RF_modelling_functions.R')

# Install all required packages and load them 
required_packages <- c("dplyr", "caret", "data.table", "h2o", "stringr", "probably", "matrixStats", "crayon", "Dict")
install_and_load(required_packages)

# Setup java installation
# if you don't have java installed, then you need to download java version 17 from 
# https://jdk.java.net/java-se-ri/17
# download and unzip to a folder and set java_installation below to point to this folder
java_installation <- "D:/Program Files/Java/jdk-17/"
Sys.setenv(JAVA_HOME='D:/Program Files/Java/jdk-17/')

# initialise the H20 session
#h2o.init(max_mem_size = "20480000000") # If client failing then increase amount of avaiable memory like so
h2o.init()

# 2. Global variables, edit these to match the data you have ####
primary_folder_data_path <- "E:/NE Data/LE_2022_2023"
directories <- check_directories_and_files(primary_folder_data_path) ## Check file structure build paths, check  all required files present

cat("Setting Global variables\n")
global_variables <- Dict$new("id_column_name" = "SegID", # Name of the id column that matches training, segmentation etc
                            # Which of your variables should be treated as factor variables
                            "factor_variables" = c("bgz","spmCarb_mode", "spmDepth_mode", "spmESB_mode", "spmSize_mode", "spmGroup_mode", "spmTex_mode", "natmapSS_mode", "natmapWet_mode", "geology_mode"),
                            # Define the levels of the response variables
                            "class_codes" = read.csv(directories["training_data_path"]) %>% dplyr::distinct(HabCode) %>% dplyr::pull(HabCode),
                            "training_proportion" = 0.8,
                            "correlation_threshold" = 0.9,
                            "seed" = round(runif(1, min = 0, max = 9999999)),
                            "year" = "2023",
                            "version" = "v1",
                            # Random forest hyperparameter setting
                            "k_folds" = 5,
                            "num_of_trees" = 100,
                            "tuned_mtry_selected" = GetOptimalMtryFromTesting(directories["mtry_tuning_data_path"]),
                            "mtry_manual_override" = 60, # Set to not null if you want to use this specified version of mtry rather than the one found during tuning
                            "tune_mtry" = FALSE, # Set this to TRUE if you would like to tune MTRY by searching the full variables space (warning this can take a long time)
                            "boot_iterations" = 100 # Number of bootstrap iterations used for calibration suggested 100
                            )

# print summary of global variables used, and save to file
global_variables_summary <- c(c("id column:", global_variables["id_column_name"]),"",
                                c("factor variables:", global_variables["factor_variables"], ""),
                                c("class_codes:", global_variables["class_codes"], ""),
                                c("training proportion:", global_variables["training_proportion"], ""),
                                c("correlation threshold:", global_variables["correlation_threshold"], ""),
                                c("seed:", global_variables["seed"], ""),
                                c("year:", global_variables["year"], ""),
                                c("version:", global_variables["version"]),
                                c("Random forest hyperparameters", ""),
                                c("k_folds", global_variables["k_folds"], ""),
                                c("number of trees:", global_variables["num_of_trees"]),
                                c("mtries used:", SetMtry(global_variables["tuned_mtry_selected"], global_variables["tune_mtry"]), ""),
                                c("bootstrap iterations used during calibration", global_variables["boot_iterations"]))

cat(global_variables_summary, sep="\n")

# 3. Process training data for modelling ####
# Report starting the preprocessing steps
print("Starting processing of data")

# only runs if no processed data path is provided
if (is.null(directories["preprocessed_data_path"])) {
    # Report starting preprocessing
    cat("No exisiting data found, starting preprocessing\n")

    # Create preprocessed data 
    all_data <- CreatePreprocessedData(directories["training_data_path"], directories["zonal_stats_directory"], id_column_name)

    # Save out preprocessed data
    data.table::fwrite(all_data, file.path(directories["input_data_directory"], "preprocessed_data.csv"))

    # Check and set preprocessed data file path
    directories["preprocessed_data_path"] <- CheckFilePath(paste0(directories["primary_folder_data_path"], "/input_data/preprocessed_data.csv"))

    # Report finishing
    cat(crayon::green("Preprocessing of data complete\n"))
} else {
    cat(crayon::green("Training data has already been preprocessed and path confirmed\n"))
}

# 4. Random forest modelling ####
# Report starting the random forest modelling section
cat("Starting Random Forest modelling\n")

# load in the preprocessed training data
cat("Loading preprocessed data\n")
training_data_processed <- data.table::fread(directories["preprocessed_data_path"], header = TRUE)

# remove classes if required
class_codes <- training_data_processed %>% dplyr::distinct(HabCode) %>% dplyr::pull(HabCode) # get class codes from training data

# Select required columns 
training_data_processed  <- training_data_processed %>% dplyr::select(contains(c(global_variables["id_column_name"], "min", "max", "mean", "mode", "stdev", "TEMP", "HabCode", "bgz", "Area", "width", "height")))

# Identify and remove highly correlated variables
cat("Identifying and removing highlight correlated variables\n")

# Subset to just the variables to be tested 
cor_data <- training_data_processed %>% dplyr::select(!append(c(global_variables["id_column_name"], "HabCode"), global_variables["factor_variables"]))

# Identify and return the dataset without highlight correlated variables
training_data_all_hcv_removed <- RemoveHighlyCorrelatedVariables(cor_data, training_data_processed, global_variables["correlation_threshold"], directories["output_path"])

# Create validation and training sets balanced by class
# Generate class counts
class_counts <- table(training_data_all_hcv_removed$HabCode)
write.csv(as.data.frame(class_counts), file.path(directories["output_path"], "class_counts.csv"))
    
# calculate class weights from the training set
# calculate frequency of each class and convert to weight 1/squareroot(frequency)
cat("Calculating class weights\n")
class_weights <- 1 / sqrt(as.numeric(table(training_data_all_hcv_removed$HabCode)))
class_weights_export_raw <- setNames(class_weights, names(table(training_data_all_hcv_removed$HabCode)))
class_weights_export <- data.frame("class" = names(class_weights_export_raw), "class_weight" = unname(class_weights_export_raw))
write.csv(class_weights_export, file.path(directories["output_path"], "class_weights.csv"))

# Save out class weights to file
write.csv(as.data.frame(1 / sqrt(table(training_data_all_hcv_removed$HabCode))), file.path(directories["output_path"], "class_weights.csv"))

# convert to dataframe for ease of selection
frequencies <- as.data.frame(table(training_data_all_hcv_removed$HabCode))
frequencies$weights_column <- sum(frequencies$Freq) / (nrow(frequencies) * frequencies$Freq)

# clean up for left joining
to_join <- frequencies %>% dplyr::select(all_of(c("Var1", "weights_column")))
colnames(to_join) <- c("HabCode", "weights_column")

# left join to create class weights
training_data_all_hcv_removed <- training_data_all_hcv_removed %>% dplyr::left_join(to_join, by = "HabCode")
training_data_all_hcv_removed$HabCode <- as.factor(training_data_all_hcv_removed$HabCode)

# Convert all character columns that should be factors are converted to factors
training_data_all_hcv_removed <- training_data_all_hcv_removed %>% dplyr::mutate(across(all_of(global_variables["factor_variables"]), as.factor))
# factorise the dependent variable
training_data_all_hcv_removed$HabCode <- as.factor(training_data_all_hcv_removed$HabCode)

# convert to h20 frames for training
training_data_all_hcv_removed <- as.h2o(training_data_all_hcv_removed)

# define predictors and response variables for random forest
predictors <- colnames(training_data_all_hcv_removed)
predictors <- predictors[!predictors %in% c("HabCode", global_variables["id_column_name"])]
response <- "HabCode"

# If mtry is not defined in the global variables then presume that it needs to be tuned, WARNING: this takes a long time to run depending on the number of variables
if (is.null(global_variables["tuned_mtry_selected"]) && global_variables["tune_mtry"]){
        # tune for mtries
        cat("Beginning tuning of mtry\n")
        
        # create results table 
        tuning_data_frame <- data.frame("run_id" = integer(), "logloss" = double(), "MSE" = double(), "RMSE" = double())

    for (mtry in 1:(length(predictors)-2)){
        cat("Trying mtry of ", mtry, "\n")
        random_forest_model <- h2o.randomForest(x = predictors,
                                        y = response,
                                        mtries = SetMtry(global_variables["tuned_mtry_selected"], global_variables["mtry_manual_override"]),
                                        ntrees = 100,
                                        training_frame = training_data_all_hcv_removed,
                                        weights_column = "weights_column",
                                        balance_classes = FALSE,
                                        seed = seed,
                                        fold_assignment = "Stratified",
                                        nfolds = 5,
                                        keep_cross_validation_models = TRUE, 
                                        keep_cross_validation_fold_assignment = TRUE,
                                        keep_cross_validation_predictions = TRUE,
                                        stopping_metric = "logloss",
                                        stopping_tolerance = 1e-3,
                                        stopping_rounds = 100
                                        )
        cat("Scoring RF\n")
        model_performance <- h2o.performance(random_forest_model)

        temp_data <- data.frame("run_id" = mtry,
                                "logloss" = h2o.logloss(model_performance),
                                "MSE" = h2o.mse(model_performance), 
                                "RMSE" = h2o.rmse(model_performance))
        tuning_data_frame <- rbind(tuning_data_frame, temp_data)
        
        # save output so that it crashes exisitng runs aren't lost
        write.csv(tuning_data_frame, file.path(directories["output_path"], "mtry_tuning.csv"))
    }
    # Open mtry output file 
    mtry_data <- read.csv(file.path(directories["output_path"], "mtry_tuning.csv"))

    # Find mtry optima for
    logloss_mtry_optimum <- mtry_data %>% dplyr::arrange(logloss) %>% dplyr::slice(1) %>% dplyr::pull(run_id) # logloss
    mse_mtry_optimum <- mtry_data %>% dplyr::arrange(MSE) %>% dplyr::slice(1) %>% dplyr::pull(run_id)# mse
    rmse_mtry_optimum <- mtry_data %>% dplyr::arrange(RMSE) %>% dplyr::slice(1) %>% dplyr::pull(run_id)# rmse

    # If more than one method agrees, then chose the most common response, if not then select logloss variables 
    if (length(unique(c(logloss_mtry_optimum, mse_mtry_optimum, rmse_mtry_optimum))) != 3){
        tuned_mtry_selected <- Mode(c(logloss_mtry_optimum, mse_mtry_optimum, rmse_mtry_optimum))
    } else {
        tuned_mtry_selected <- logloss_mtry_optimum
    }
    cat(crayon::green("Mtry tuned and optimal found to be: ", tuned_mtry_selected, "\n"))
} else{
   cat(crayon::red("Mtry is not being tuned\n"))
}
global_variables
cat(crayon::green("Beginning Random Forest Training\n"))
# Create the final random forest model, using tuned mtry and k-fold validation
random_forest_model <- h2o.randomForest(x = predictors,
                    y = response,
                    mtries = SetMtry(global_variables["tuned_mtry_selected"], global_variables["mtry_manual_override"]),
                    ntrees = 10,
                    training_frame = training_data_all_hcv_removed,
                    weights_column = "weights_column",
                    balance_classes = FALSE,
                    seed = seed,
                    fold_assignment = "Stratified",
                    nfolds = 5,
                    keep_cross_validation_models = TRUE, 
                    keep_cross_validation_fold_assignment = TRUE,
                    keep_cross_validation_predictions = FALSE,
                    stopping_metric = "logloss",
                    min_rows = 10
                    )

cat(crayon::green("Random Forest succesfully trainined\n"))
# Save out overall model
cat("Saving Random Forest Model\n")
h2o.saveModel(random_forest_model, file.path(directories["model_output_path"]))

# Save out each fold model
list_of_cross_validation_models <- h2o.cross_validation_models(random_forest_model)
for (fold in 1:5){
    h2o.saveModel(list_of_cross_validation_models[[fold]], file.path(directories["model_output_path"]))
}

# save out the fold assignments
SegID <- as.data.frame(training_data_all_hcv_removed[global_variables["id_column_name"]])
fold_assignment <- as.data.frame(h2o.cross_validation_fold_assignment(random_forest_model))
fold_output <- data.frame(SegID, fold_assignment)
colnames(fold_output) <- c(global_variables["id_column_name"], "fold_assignment")
write.csv(fold_output, file.path(directories["output_path"], "fold_assignment.csv"))

# 5. Evaluate the performance of all random forest fold models #####
print("Evaluating Random Forest Model Performance")

# Check model output folder for fold models, and retrieve paths
list_of_cross_validation_models <- GetFoldModels(directories["model_output_path"])

# Load the processed data
data_processed <- data.table::fread(directories["preprocessed_data_path"], header = TRUE)

# Remove any extra columns that may exist
data_processed  <- data_processed %>% dplyr::select(contains(c(global_variables["id_column_name"],"min", "max", "mean", "mode", "stdev", "TEMP", "HabCode", "bgz")))

# Ensure that factor variables are factors 
data_processed <- data_processed %>% dplyr::mutate(across(all_of(global_variables["factor_variables"]), as.factor))

# Remove all highly correlated variables 
highly_correlated_variables <- read.csv(file.path(directories["output_path"], "highly_correlated_variables.csv"), header = TRUE)
data_processed <- data_processed %>% dplyr::select(-all_of(highly_correlated_variables$highly_correlated_variables))

# load the fold assigments from file
fold_assignment_dataframe <- read.csv(file.path(directories["output_path"], "fold_assignment.csv"), header = TRUE)

# Predict for each cross validation model
for (fold in 1:global_variables["k_folds"]){
    fold_output_path <- validate_dir(paste0(directories["model_evaluation_output_path"], paste0("/model_score_kfold_", fold)))
    
    # load the relevent model from storage
    fold_model <- h2o.loadModel(list_of_cross_validation_models[fold])

    # get the relevent cross validation model and predict on all training data
    validation_predictions <- as.data.frame(h2o.predict(fold_model, as.h2o(data_processed)))

    # create the confusion matrix and performance evaluation of each fold model 
    # Get the validation fold assigment from the model object to identify the correct validation data 
    validation_fold_ids <- fold_assignment_dataframe %>% dplyr::filter(fold_assignment == (fold - 1))
    
    # filter to just these seg ids in the training data 
    fold_validation_data <- data_processed %>% dplyr::filter(get(global_variables["id_column_name"]) %in% validation_fold_ids[,global_variables["id_column_name"]])
   
    # predict for all fold validation rows 
    fold_validation_predictions <- as.data.frame(h2o.predict(fold_model, as.h2o(fold_validation_data)))
    
    # add segID back
    fold_validation_predictions[[global_variables["id_column_name"]]] <- fold_validation_data %>% dplyr::select(all_of(global_variables["id_column_name"])) %>% dplyr::pull(global_variables["id_column_name"])

    # add Habcode back in 
    fold_validation_predictions$HabCode <- fold_validation_data$HabCode
    # add BGZ back in 
    fold_validation_predictions$bgz <- fold_validation_data$bgz

    # identify all segments with duplicates 
    segments_with_duplicates_ids <- data_processed %>% 
                                        dplyr::group_by(get(global_variables["id_column_name"])) %>% 
                                        dplyr::filter(n() > 1) %>% 
                                        dplyr::distinct(get(global_variables["id_column_name"])) %>%
                                        dplyr::pull()
        
    # get the duplicated rows to work out which to keep if any 
    segments_with_duplicates <- fold_validation_predictions %>% dplyr::filter(get(global_variables["id_column_name"]) %in% segments_with_duplicates_ids)

    # Create a results dataframe to store the results
    results <- segments_with_duplicates[FALSE, ]

    dupl_segment_id <- unique(segments_with_duplicates[,global_variables["id_column_name"]])[1]
    for (dupl_segment_id in unique(segments_with_duplicates[,global_variables["id_column_name"]])){
        # Get the habitat that is is predicted as
        predicted_habitat <- segments_with_duplicates %>%
                                    dplyr::filter(get(global_variables["id_column_name"]) == dupl_segment_id) %>% 
                                    dplyr::distinct(predict) %>% 
                                    dplyr::pull()

        # Get all the habitats that it could be from the validation dataset
        actual_habitat <- segments_with_duplicates %>% 
                                    dplyr::filter(get(global_variables["id_column_name"]) == dupl_segment_id) %>% 
                                    dplyr::filter(HabCode %in% predicted_habitat)

        # if length of zero then means both guess were wrong, therefore add both rows, else add just the matching row
        if (nrow(actual_habitat) == 0) {
            # Get the original data and add to results 
            results <- rbind(results, segments_with_duplicates %>% dplyr::filter(get(global_variables["id_column_name"]) == dupl_segment_id))
        } else {
            results <- rbind(results, actual_habitat)
        }
    }

    # remove all duplicate rows from the data
    fold_validation_predictions <- fold_validation_predictions %>% dplyr::filter(!get(global_variables["id_column_name"]) %in% segments_with_duplicates_ids)

    # add the updated results back to the validation predictions
    fold_validation_predictions <- rbind(fold_validation_predictions, results)

    # Convert prediction and habcode to factors 
    actual_class <- as.data.frame(fold_validation_predictions) %>% mutate_at("HabCode", as.factor) %>% dplyr::select("HabCode")
    predicted_class <- as.data.frame(fold_validation_predictions) %>% mutate_at("predict", as.factor) %>% dplyr::select("predict")

    actual_class <- actual_class %>% mutate(HabCode = case_when(HabCode == "HR" ~ "SCR", TRUE ~ HabCode))
    predicted_class <- predicted_class %>% mutate(predict = case_when(predict == "HR" ~ "SCR", TRUE ~ predict))

    # Extract the highest model score for reliability testing
    aprob <- apply(as.data.frame(fold_validation_predictions) %>% dplyr::select(!all_of(c("predict", "HabCode", "bgz", global_variables["id_column_name"]))), 1, max, na.rm=TRUE)
    export_a_prob <- cbind(as.data.frame(aprob), as.data.frame(fold_validation_predictions %>% dplyr::select(global_variables["id_column_name"]))) # join to the segID to allow joining later
    write.csv(export_a_prob, file.path(fold_output_path, "all_predictions_aprob.csv"))

    # Output validation predictions to file joined to id for spatial visualisation
    exportdata <- fold_validation_predictions %>% dplyr::select(all_of(c("predict", global_variables["id_column_name"], "HabCode"))) %>% as.data.frame()
    write.csv(exportdata, file.path(fold_output_path, "all_predictions.csv"))

    # Ensure that both predicted class factor, and validation class factor have the same levels 
    #fillinmissinglevels <- factor(validation_predictions_data$predict, levels = c(levels(validation_class$HabCode), levels(validation_predictions_data$predict)[!(levels(validation_predictions_data$predict) %in% levels(validation_class$HabCode))]))

    # Calculate model evaluation statistics
    # Calculate overall accuracy
    accuracy <- mean(actual_class$HabCode == predicted_class$predict, na.rm = TRUE)
    accuracy

    # Calculate confusion matrix
    cm <- caret::confusionMatrix(data = as.factor(predicted_class$predict), reference = as.factor(actual_class$HabCode))
    cm

    # write out confusion matrix
    write.csv(data.frame(cm$table), file.path(fold_output_path, "confusion_matrix.csv"))

    # Calculate mean balanced accuracy across all classes
    mean(cm$byClass[,"Balanced Accuracy"])

    # Calculate mean F1 score across all classes
    mean(cm$byClass[,"F1"], na.rm = TRUE)

    # print classwise statistics
    data.frame(cm$byClass)

    # save out classwise statistics
    write.csv(data.frame(cm$byClass), file.path(fold_output_path, "classwise_statistics.csv"))

    # macro f1 score
    macrof1 <- cm$byClass[, "F1"]
    macrof1[is.na(macrof1)] <- 0
    macrof1

    # save overall statistics to a single file 
    cat(c("Mean accuracy", accuracy, ""), 
        c("Mean Balanced accuracy", mean(cm$byClass[,"Balanced Accuracy"]), ""),
        c("Mean F1 score", mean(cm$byClass[,"F1"], na.rm = TRUE), ""),
        file=file.path(fold_output_path, "random_forest_overall_statistics.txt"), sep="\n")

    # classwise statsitcs
    cat(c("Mean accuracy", accuracy, ""), 
        c("Mean Balanced accuracy", mean(cm$byClass[,"Balanced Accuracy"]), ""),
        c("Mean F1 score", mean(cm$byClass[,"F1"], na.rm = TRUE), ""),
        file = file.path(fold_output_path, "random_forest_overall_statistics.txt"), sep="\n")

    # calculate and save out variable importance
    varimp <- h2o.varimp(fold_model)
    varimp
    write.csv(varimp, file.path(fold_output_path, "variable_importance.csv"))

    # Combine outputs for the validation set to create the output for reliability testing 
    SegIDs <- fold_validation_predictions[, global_variables["id_column_name"]]
    BGZ <- fold_validation_predictions$bgz
    Validation_HabCode <- fold_validation_predictions$HabCode
    
    #validation_predictions_data
    Prediction_Correct <- Validation_HabCode == fold_validation_predictions$predict

    combined_output <- data.frame(SegIDs, BGZ, Validation_HabCode, fold_validation_predictions, Prediction_Correct)
    
    write.csv(combined_output, file.path(fold_output_path, "reliability_output.csv"))

    # Calibrate the probabilities using a one-vs-all methodology, isotonic regression and scaling so that probabilies sum to 1
    print("Beginning model calibration ")
    # Create folder to store all calibration related files  
    # check if output path folder exists, if not then make it
    calibration_output_path <- validate_dir(paste0(directories["model_evaluation_output_path"], paste0("/calibrated_probability_kfold_", fold)))

    # Predict across all test/train points to provide data to train calibration model
    calibration_predictions_data <- validation_predictions
    
    # Deal with segments which have more than one entry as before 
    # add segID back 
    calibration_predictions_data[[global_variables["id_column_name"]]] <- data_processed %>% dplyr::select(all_of(global_variables["id_column_name"])) %>% dplyr::pull(global_variables["id_column_name"])

    # add Habcode back in 
    calibration_predictions_data$HabCode <- data_processed$HabCode
    # add BGZ back in 
    calibration_predictions_data$bgz <- data_processed$bgz

    # identify all segments with duplicates 
    segments_with_duplicates_ids <- data_processed %>%
                                        dplyr::group_by(get(global_variables["id_column_name"])) %>% 
                                        dplyr::filter(n() > 1) %>% 
                                        dplyr::distinct(get(global_variables["id_column_name"])) %>%
                                        dplyr::pull()

    # get the duplicated rows to work out which to keep if any 
    segments_with_duplicates <- calibration_predictions_data %>% dplyr::filter(get(global_variables["id_column_name"]) %in% segments_with_duplicates_ids)
    
    # Create a results dataframe to store the results 
    results <- segments_with_duplicates[FALSE, ]

    for (dupl_segment_id in unique(segments_with_duplicates[, global_variables["id_column_name"]])){
        # Get the habitat that is is predicted as 
        predicted_habitat <- segments_with_duplicates %>%
                                    dplyr::filter(get(global_variables["id_column_name"]) == dupl_segment_id) %>% 
                                    dplyr::distinct(predict) %>%
                                    dplyr::pull(predict)

        # Get all the habitats that it could be from the validation dataset
        actual_habitat <- segments_with_duplicates %>% 
                                    dplyr::filter(get(global_variables["id_column_name"]) == dupl_segment_id) %>%
                                    dplyr::filter(HabCode %in% predicted_habitat)

        # if length of zero then means both guess were wrong, therefore add both rows, else add just the matching row
        if (nrow(actual_habitat) == 0) {
            # Get the original data and add to results 
            results <- rbind(results, segments_with_duplicates %>% dplyr::filter(SegID == dupl_segment_id))
        } else {
            results <- rbind(results, actual_habitat)
        }
    }

    # remove all duplicate rows from the data 
    calibration_predictions_data <- calibration_predictions_data %>% dplyr::filter(!get(global_variables["id_column_name"]) %in% segments_with_duplicates_ids)

    # add the updated results back to the validation predictions
    calibration_predictions_data <- rbind(calibration_predictions_data, results)

    # Create a new subfolder to store calibration data 
    # check if output path folder exists, if not then make it
    calibration_data_path <- validate_dir(paste0(calibration_output_path, "/calibration_files"))

    # Create a results table to store all the data
    results <- calibration_predictions_data
    
    #current_class <- class_codes[1]
    # Create the calibration models for all classes
    for (current_class in global_variables["class_codes"]){
        # Report start
        print(paste0("Calibrating ", current_class))

        # From the calibration dataset extract the probability
        probability <- calibration_predictions_data %>% dplyr::select(all_of(current_class)) %>% dplyr::pull()
        
        # Create column that checks if prediction was correct or not 
        class <- calibration_predictions_data %>% dplyr::select(predict) %>% dplyr::mutate(class_predicted_correctly = predict == current_class) %>% dplyr::pull(class_predicted_correctly)

        # Combine together to create a binary matrix one vs all
        binary_matrix <- data.frame(class, probability)

        binary_matrix %>% probably::cal_plot_breaks(class, probability, include_rug = FALSE, num_breaks = 20)
        ggsave(file.path(calibration_data_path, paste0(current_class,"_before",".png")))

        # fit a bootstrapped isotonic regression to calibrate
        print(paste0("Fitting bootstrapped isotonic regression for ", current_class))
        
        # Create a frame to store results 
        results_itr <- matrix(ncol = global_variables["boot_iterations"], nrow = nrow(binary_matrix))
        
        for (bootstrap in 1:global_variables["boot_iterations"]){
            print(paste("Isotonic bootstrap", bootstrap))

            # for each bootstrap resample the data to constant sample size across bins 
            results_table <- data.frame(lower=numeric(), mid=numeric(), upper=numeric(), event_rate=numeric(), pos_count = numeric(), neg_count = numeric(), total_count = numeric())
            number_of_breaks <- 10
            break_width <- 1/number_of_breaks
            for (lower_limit in seq(0.00, 1 - break_width, break_width)){
                # get all model scores above a certain threshold
                test <- binary_matrix %>% dplyr::filter(probability >= lower_limit & probability <= (lower_limit + break_width))

                # count positives 
                positives <- test %>% dplyr::filter(class == "TRUE") %>% dplyr::summarise(n = n()) %>% dplyr::pull(n)

                # count negatives 
                negatives <- test %>% dplyr::filter(class == "FALSE") %>% dplyr::summarise(n = n()) %>% dplyr::pull(n)

                # calculate event rate 
                event_rate <- positives / (positives + negatives)

                # save to results dataframe
                results_table <- rbind(results_table, data.frame(lower=lower_limit, mid=(lower_limit+lower_limit+break_width)/2,upper=lower_limit+break_width, event_rate = event_rate, pos_count=positives, neg_count = negatives, total_count = positives+negatives))
            }

            resampled_binary_matrix <- data.frame(class = character(), probability = numeric())
            # find smallest bin
            min_bin_size <- min(results_table$total_count[results_table$total_count > 0])

            # if bin size is too small then set to 100 as default
            if (min_bin_size < 100){
                min_bin_size <- 100
            }

            for (lower_limit in seq(0.00, 1-break_width, break_width)){
                # get all model scores within threshold
                test <- binary_matrix %>% dplyr::filter(probability >= lower_limit & probability <= (lower_limit + break_width))

                # If greater than bin size available, then subsmaple down to min_bin_size
                if (nrow(test) >= min_bin_size){
                    # get sample these without replacement up to the min bin size to create new 
                    subsample <- sample(1:nrow(test), size = min_bin_size, replace = FALSE)
                    # select rows and append to new data frame 
                    resampled_binary_matrix <- rbind(resampled_binary_matrix, test[subsample, ])
                } else {
                    subsample <- test
                    # select rows and append to new data frame 
                    resampled_binary_matrix <- rbind(resampled_binary_matrix, test)
                }
            }
            
            # fit an isotonic regression
            iso_fit <- as.stepfun(isoreg(x = resampled_binary_matrix$probability, y = resampled_binary_matrix$class))
            print("Isotonic regression fitted")
            
            # predict on rest of samples
            isoreg_calibrated_probabilies <- iso_fit(calibration_predictions_data %>% dplyr::select(all_of(current_class)) %>% dplyr::pull())
            
            binary_matrix_post_calibration <- data.frame("class" = class, "probability" = isoreg_calibrated_probabilies)
            
            # Save the calibration model out for use later
            saveRDS(iso_fit, file.path(calibration_data_path, paste0(current_class,"_calibration_model","_",bootstrap)))
            print("Model saved to file")

            # save in results
            results_itr[,bootstrap] <- isoreg_calibrated_probabilies
        }

        # calculate mean bootstrap values
        isoreg_calibrated_probabilies_final <- rowMeans(results_itr)

        # store calibrated mean and upper and lower confidence interval 
        results[current_class] <- isoreg_calibrated_probabilies_final

        # convert to binary matrix for plotting 
        binary_matrix_post_calibration <- data.frame("class" = class, "probability" = isoreg_calibrated_probabilies_final)
        
        # Plot the calibration error after and save to file
        binary_matrix_post_calibration %>% probably::cal_plot_breaks(class, probability, include_rug = FALSE, num_breaks = 20)
        ggsave(file.path(calibration_data_path, paste0(current_class,"_after",".png")))
    }
    # calculate mean 
    # Get just the numeric columns 
    results_numeric <- results %>% dplyr::select(all_of(global_variables["class_codes"]))
    # Scale each row to sum to 1
    #results_scaled <- results_numeric / rowSums(results_numeric)
    results_scaled <- results_numeric

    # updated predict
    updated_predict_raw <- (colnames(results_scaled)[max.col(results_scaled, ties.method="first")])

    # find rows with all 0's for later
    zero_columns <- (rowSums(results_numeric) == 0)

    # add back the segID and other data
    results_scaled$predict <- as.character(results$predict)
    results_scaled <- results_scaled %>% mutate(SegID = results %>% dplyr::select(all_of(global_variables["id_column_name"])) %>% dplyr::pull())
    results_scaled$HabCode <- results$HabCode
    results_scaled$bgz <- results$bgz

    # add the updated prediction class for ones with 0's keep the original prediction 
    results_scaled <- results_scaled %>% dplyr::mutate(updated_predict = ifelse(zero_columns, results_scaled$predict, updated_predict_raw))

    # convert probability to a class confidence
    for (class_code in global_variables["class_codes"]){
        name <- paste0(class_code, "_confidence_category")
        
        temp <- results_scaled %>% dplyr::select({{class_code}}) %>%
        mutate(
           "{name}" := dplyr::case_when(
                .[[1]] <= 0.2 ~ "very low confidence",
                .[[1]] > 0.2 & .[[1]] <= 0.4 ~ "low confidence",
                .[[1]] > 0.4 & .[[1]] <= 0.6 ~ "medium confidence",
                .[[1]] > 0.6 & .[[1]] <= 0.8 ~ "high confidence",
                .[[1]] > 0.8 ~ "very high confidence", 
            )
        ) %>% dplyr::select({{name}})
        results_scaled <- cbind(results_scaled, temp)
    }

    # Export validation predictions 
    write.csv(results_scaled, file.path(calibration_output_path, "calibrated_validation_predictions.csv"))

    # select just the rows for validation for this fold 
    # Get the validation fold assigment from the model object to identify the correct validation data 
    
    # filter to just these seg ids in the training data 
    results_scaled_validation_data <- results_scaled %>% dplyr::filter(get(global_variables["id_column_name"]) %in% validation_fold_ids[[global_variables["id_column_name"]]])

    # Convert prediction and habcode to factors 
    actual_class <- as.data.frame(results_scaled_validation_data) %>% mutate_at("HabCode", as.factor) %>% dplyr::select("HabCode")
    predicted_class <- as.data.frame(results_scaled_validation_data) %>% mutate_at("updated_predict", as.factor) %>% dplyr::select("updated_predict")

    accuracy <- mean(actual_class$HabCode == predicted_class$updated_predict, na.rm = TRUE)
    accuracy
    
    # Calculate confusion matrix
    cm <- caret::confusionMatrix(data = as.factor(predicted_class$updated_predict), reference = as.factor(actual_class$HabCode))
    cm

    # write out confusion matrix
    write.csv(data.frame(cm$table), file.path(calibration_output_path, "confusion_matrix.csv"))

    # Calculate mean balanced accuracy across all classes
    mean(cm$byClass[,"Balanced Accuracy"])

    # Calculate mean F1 score across all classes
    mean(cm$byClass[,"F1"], na.rm = TRUE)

    # print classwise statistics
    data.frame(cm$byClass)

    # save out classwise statistics
    write.csv(data.frame(cm$byClass), file.path(calibration_output_path, "classwise_statistics.csv"))

    # macro f1 score
    macrof1 <- cm$byClass[, "F1"]
    macrof1[is.na(macrof1)] <- 0
    macrof1

    # save overall statistics to a single file 
    cat(c("Mean accuracy", accuracy, ""), 
        c("Mean Balanced accuracy", mean(cm$byClass[,"Balanced Accuracy"]), ""),
        c("Mean F1 score", mean(cm$byClass[,"F1"], na.rm = TRUE), ""),
        file=file.path(calibration_output_path, "random_forest_overall_statistics.txt"), sep="\n")

    # classwise statsitcs
    cat(c("Mean accuracy", accuracy, ""), 
        c("Mean Balanced accuracy", mean(cm$byClass[,"Balanced Accuracy"]), ""),
        c("Mean F1 score", mean(cm$byClass[,"F1"], na.rm = TRUE), ""),
        file = file.path(calibration_output_path, "random_forest_overall_statistics.txt"), sep="\n")

    # Combine outputs for the validation set to create the output for reliability testing 
    SegIDs <- results_scaled_validation_data %>% dplyr::select(all_of(global_variables["id_column_name"])) %>% dplyr::pull()
    BGZ <- results_scaled_validation_data$bgz
    Validation_HabCode <- results_scaled_validation_data$HabCode
    
    #validation_predictions_data
    Prediction_Correct <- Validation_HabCode == results_scaled_validation_data$updated_predict
    
    combined_output <- data.frame(SegIDs, BGZ, Validation_HabCode, results_scaled_validation_data %>% dplyr::select(all_of(global_variables["class_codes"])), Prediction_Correct)
    write.csv(combined_output, file.path(calibration_output_path, "reliability_output.csv"))
}

# 6. Predict for all segments using a specified model ####
model_to_use <- 4 # provide the number of the model to use
list_of_cross_validation_models <- GetFoldModels(directories["model_output_path"])

#c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13")
for (zone in c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13")){
    # Load the model identified as being the best 
    best_model <- h2o.loadModel(list_of_cross_validation_models[model_to_use])

    # Load zonal stats
    cat("Starting BGZ", zone, "\n")
    file_name <- paste("BGZ", zone, "_zonalStats.csv", sep = "")
    zonal_stats_seg <- data.table::fread(file.path(directories["zonal_stats_directory"], file_name), header = TRUE)
    zonal_stats_seg$bgz <- paste("BGZ", zone, sep = "")

    # save out the seg ids to add back later
    segid_to_join <- zonal_stats_seg[[global_variables["id_column_name"]]]

    # convert required variables to factor variables
    zonal_stats_seg <- zonal_stats_seg %>%
        dplyr::mutate(across(all_of(global_variables["factor_variables"]), as.factor))

    # Check if a small number of points to predict on, and if so then do them in one go, else chunk    
    if (nrow(zonal_stats_seg) <= 100000){
        data <- as.h2o(zonal_stats_seg)
        results <- as.data.frame(matrix(ncol = length(global_variables["class_codes"]) + 1, nrow = 0))
        output_predict <- as.data.frame(predict(best_model, newdata = data))
        results <- rbind(results, output_predict)
        cat("Prediction complete", "\n")
    } else {
        # Split the zone into chunks to predict on
        print(paste("Creating partitions to predict on"))
        seq_1 <- floor(seq(from = 1, to = nrow(zonal_stats_seg), length.out = ceiling(nrow(zonal_stats_seg)/100000)))
        seq_2 <- seq_1[2:length(seq_1)]
        seq_1 <- c(1, seq_1[-c(1, length(seq_1))] + 1)

        # create results dataframe to store
        print(paste("Create dataframe to store predictions"))
        results <- as.data.frame(matrix(ncol = length(global_variables["class_codes"]) + 1, nrow = 0))
        colnames(results) <- append("predict", global_variables["class_codes"])
        for (i in 1:length(seq_1)){
            cat("Predicting for rows", seq_1[i], "to", seq_2[i], "\n")
            # subset and convert to h20
            data <- as.h2o(zonal_stats_seg[seq_1[i]:seq_2[i], ])
            # predict using model 
            output_predict <- as.data.frame(predict(best_model, newdata = data))
            results <- rbind(results, output_predict)
            cat("Chunk", i, "of", length(seq_1), "complete", "\n")
        }
    }

    # Format results for export
    # drop predict column
    results_for_model_score <- results %>% dplyr::select(!predict)

    ## Read a class probability prediction matrix and return the n class with highest probabilities for each id
    nProbClass <- function(row, n) {
        colnames(results_for_model_score)[order(row, decreasing = TRUE)[1:n]]
    }

    ## Read a class probability prediction matrix and return the n highest probabilities for each id
    nProbProb <- function(row,n) {
        unname(row[order(row, decreasing = TRUE)[1:n]])
    }

    print("Formatting results")
    class_predictions <- data.frame(t(apply(results_for_model_score, 1, nProbClass,n=length(global_variables["class_codes"]))), row.names = 1:nrow(results_for_model_score))
    colnames(class_predictions) <- paste(letters[1:length(class_predictions)], "_pred", sep="")
 
    class_predictions_probabilities <- data.frame(t(apply(results_for_model_score, 1, nProbProb,n=length(global_variables["class_codes"]))),row.names = 1:nrow(results_for_model_score))
    colnames(class_predictions_probabilities) <- paste(letters[1:length(class_predictions)], "_prob", sep="")

    # Combine together into a single dataframe
    # zip colnames together 
    names_tuple <- rbind(paste(letters[1:length(class_predictions)], "_pred", sep=""),paste(letters[1:length(class_predictions)], "_prob", sep=""))
    
    # create blank dataframe to store data 
    results_for_model_score <- as.data.frame(matrix(data=NA, nrow = nrow(class_predictions), ncol = length(names_tuple) + 1))
    
    # assign colnames from the names tuple and add SegID at the start 
    colnames(results_for_model_score) <- c(global_variables["id_column_name"], names_tuple)
    
    # Fill dataframe with data
    for (column_name in colnames(results_for_model_score)){
        # check if pred column
        if (grepl("pred", column_name, fixed = TRUE)){
            results_for_model_score[column_name] <- class_predictions[column_name]
        }
        # check if prob column
        if (grepl("prob", column_name, fixed = TRUE)){
            results_for_model_score[column_name] <- class_predictions_probabilities[column_name]
        } 
        # if not then it is the id column 
        if (grepl(global_variables["id_column_name"], column_name, fixed = TRUE)){
            results_for_model_score[column_name] <- segid_to_join
        }
    }

    # Get a prob and pred and join with 
    model_score_lookup <- read.csv(file.path(paste0(directories["input_data_directory"], "/parametric_reliability_model.csv")))

    # get a pred and prob columns 
    aprob_pred <- results_for_model_score %>% dplyr::select(all_of(c("a_pred", "a_prob", global_variables["id_column_name"])))

    # round prob column nto accuracy of model score lookup
    aprob_pred$a_prob <- round(aprob_pred$a_prob, 2)

    joined_data <- aprob_pred %>%
                dplyr::left_join(model_score_lookup, by = c("a_pred" = "pred_habitat", "a_prob" = "A_prob")) %>%
                dplyr::select(!all_of(c("X")))
    
    # convert the estimated accuracy to a categorical variable 
    joined_data_test <- joined_data %>%
        mutate(
           "estimated_accuracy_class" := dplyr::case_when(
                R_pred <= 0.2 ~ "very low",
                R_pred > 0.2 & R_pred <= 0.4 ~ "low",
                R_pred > 0.4 & R_pred <= 0.6 ~ "medium",
                R_pred > 0.6 & R_pred <= 0.8 ~ "high",
                R_pred > 0.8 ~ "very high",
            )
        ) %>%
        dplyr::select(all_of(c("R_lower", "R_pred", "R_upper", "estimated_accuracy_class")))

    # Join onto dataframe for export
    results_for_model_score <- cbind(results_for_model_score, joined_data_test)
 
    # save output
    cat("Saving output to disk\n")
    data.table::fwrite(results_for_model_score, paste0(directories["segment_output_path_model_score"],"/", "bgz_", zone, "_model_score.csv"))

    # begin model calibration
    cat("Calibrating probabilies")
    
    # create copy of the results for calibrating
    calibration_results <- results
    #class_code <- class_codes[1] 
    for (class_code in global_variables["class_codes"][order(global_variables["class_codes"])]){
        cat("Calibrating ", class_code, "\n")
        # build the file path to the correct calibration folder
        path_to_models <- file.path(paste0(directories["model_evaluation_output_path"], paste0("/calibrated_probability_kfold_", model_to_use, "/calibration_files")))

        # iterate over all bootstrap models
        calibration_bootstrap_results_store <- data.frame(matrix(nrow = nrow(results), ncol = global_variables["boot_iterations"]))
        for (bootstrap_iteration in 1:global_variables["boot_iterations"]){
            # load the model
            iso_fit <- readRDS(file.path(paste0(path_to_models, "/",class_code ,"_calibration_model_", bootstrap_iteration)))
            
            # predict for all segments 
            calibration_bootstrap_results_store[ ,bootstrap_iteration] <- iso_fit(results[ ,class_code])
        }
        # convert bootstrap to averages across rows
        isoreg_calibrated_probabilies_final <- rowMeans(calibration_bootstrap_results_store)

        # store calibrated probabilities
        calibration_results[class_code] <- isoreg_calibrated_probabilies_final
    }
    
    # Get just the numeric columns 
    results_numeric <- calibration_results %>% dplyr::select(all_of(global_variables["class_codes"][order(global_variables["class_codes"])]))

    # Scale each row to sum to 1
    #results_scaled <- results_numeric / rowSums(results_numeric)
    # previous step introduces nans where row sum is 0, convert these back to 0
    #results_scaled <- results_scaled %>% replace(is.na(.), 0)
    results_scaled <- results_numeric

    # Read a class probability prediction matrix and return the n class with highest probabilities for each id
    nProbClass <- function(row, n) {
        colnames(results_scaled)[order(row, decreasing = TRUE)[1:n]]
    }

    # Read a class probability prediction matrix and return the n highest probabilities for each id
    nProbProb <- function(row,n) {
        unname(row[order(row, decreasing = TRUE)[1:n]])
    }

    print("Formatting calibrated results")
    class_predictions_calibrated <- data.frame(t(apply(results_scaled, 1, nProbClass, n=length(global_variables["class_codes"]))), row.names = 1:nrow(results_scaled))
    colnames(class_predictions_calibrated) <- paste(letters[1:length(class_predictions_calibrated)], "_pred", sep="")
 
    class_predictions_probabilities_calibrated <- data.frame(t(apply(results_scaled, 1, nProbProb,n=length(global_variables["class_codes"]))),row.names = 1:nrow(results_scaled))
    colnames(class_predictions_probabilities_calibrated) <- paste(letters[1:length(results_scaled)], "_prob", sep="")

    # find rows with all 0's
    zero_columns <- (rowSums(results_numeric) == 0)

    # for all rows with 0 replace the row with the order from the model score output
    class_predictions_calibrated[zero_columns,] <- class_predictions[zero_columns,]

    # categorise the probabilities for each class
    # make dataframe to store results 
    cal_prob_prob_class <- class_predictions_probabilities_calibrated
    colnames(cal_prob_prob_class) <- paste(letters[1:length(results_scaled)], "_reliability_score", sep="")

    for (prob_letter in letters[1:length(results_scaled)]){
        name <- paste(prob_letter, "_reliability_score", sep="")
        
        temp <- class_predictions_probabilities_calibrated %>% dplyr::select(paste(prob_letter, "_prob", sep="")) %>%
        mutate(
            "{name}" := dplyr::case_when(
                .[[1]] <= 0.2 ~ "very low",
                .[[1]] > 0.2 & .[[1]] <= 0.4 ~ "low",
                .[[1]] > 0.4 & .[[1]] <= 0.6 ~ "medium",
                .[[1]] > 0.6 & .[[1]] <= 0.8 ~ "high",
                .[[1]] > 0.8 ~ "very high", 
            )
        ) %>% dplyr::select({{name}})
        cal_prob_prob_class[, name] <- temp
    }

    # Combine together into a single dataframe
    # zip colnames together 
    names_tuple <- rbind(paste(letters[1:length(class_predictions)], "_pred", sep=""),paste(letters[1:length(class_predictions)], "_prob", sep=""), paste(letters[1:length(results_scaled)], "_reliability_score", sep=""))
    
    # create blank dataframe to store data 
    results_for_calibrated_probs <- as.data.frame(matrix(data=NA, nrow = nrow(class_predictions_calibrated), ncol = length(names_tuple) + 1))
    
    # assign colnames from the names tuple and add SegID at the start 
    colnames(results_for_calibrated_probs) <- c(global_variables["id_column_name"], names_tuple)
    
    # Fill dataframe with data
    for (column_name in colnames(results_for_calibrated_probs)){
        # check if pred column
        if (grepl("pred", column_name, fixed = TRUE)){
            results_for_calibrated_probs[column_name] <- class_predictions_calibrated[column_name]
        }
        # check if prob column 
        if (grepl("prob", column_name, fixed = TRUE)){
            results_for_calibrated_probs[column_name] <- class_predictions_probabilities_calibrated[column_name]
        }
        # check if prob column 
        if (grepl("_reliability_score", column_name, fixed = TRUE)){
            results_for_calibrated_probs[column_name] <- cal_prob_prob_class[column_name]
        } 
        # if not then it is the id column 
        if (grepl(global_variables["id_column_name"], column_name, fixed = TRUE)){
            results_for_calibrated_probs[column_name] <- segid_to_join
        }
    }
    cat(crayon::green("Calibration complete saving results to file\n"))

    # Export validation predictions 
    data.table::fwrite(results_for_calibrated_probs, paste0(directories["segment_output_path_calibrated_results"], "/bgz_", zone, "_calibrated_results.csv"))
    
    # Create a combo version just for a prob of both categories attached 
    orderCalibratedResults <- function(position, dataframe) {
        as.numeric(results_numeric[position, order(as.numeric(dataframe[position, ]), decreasing = TRUE)])
    }

    positions <- 1:nrow(results_numeric)
    calibrated_results_ordered_by_model_score <- lapply(positions, orderCalibratedResults, results %>% dplyr::select(!predict))
    
    df <- data.frame(matrix(unlist(calibrated_results_ordered_by_model_score), nrow=length(calibrated_results_ordered_by_model_score), byrow=TRUE))
    colnames(df) <- paste(letters[1:length(df)], "_reliability_score", sep="")
    
    # convert the probability to a category based on the score 
    for (column_name in paste(letters[1:length(df)], "_reliability_score", sep="")){
        cat(paste("Converting", column_name, "\n"))
        
        temp <- df %>% dplyr::select(all_of(column_name)) %>%
        mutate(
            "{column_name}" := dplyr::case_when(
                .[[1]] <= 0.2 ~ "very low",
                .[[1]] > 0.2 & .[[1]] <= 0.4 ~ "low",
                .[[1]] > 0.4 & .[[1]] <= 0.6 ~ "medium",
                .[[1]] > 0.6 & .[[1]] <= 0.8 ~ "high",
                .[[1]] > 0.8 ~ "very high", 
            )
        ) %>% dplyr::select(all_of({{column_name}}))
        df[,column_name] <- temp
    }

    # assemble the final dataframe from the constituent results dataframes s
    names_tuple <- rbind(paste(letters[1:length(class_predictions)], "_pred", sep=""), paste(letters[1:length(class_predictions)], "_prob", sep=""), paste(letters[1:length(df)], "_reliability_score", sep=""))

    # create blank dataframe to store data 
    results_for_combined_output <- as.data.frame(matrix(data=NA, nrow = nrow(df), ncol = length(names_tuple) + 1))
    
    # assign colnames from the names tuple and add SegID at the start 
    colnames(results_for_combined_output) <- c(global_variables["id_column_name"], names_tuple)

    # Fill dataframe with data
    for (column_name in colnames(results_for_calibrated_probs)){
        # check if pred column
        if (grepl("pred", column_name, fixed = TRUE)){
            results_for_combined_output[column_name] <- class_predictions[column_name]
        }
        # check if prob column 
        if (grepl("prob", column_name, fixed = TRUE)){
            results_for_combined_output[column_name] <- class_predictions_probabilities[column_name]
        }
        # check if prob column 
        if (grepl("_reliability_score", column_name, fixed = TRUE)){
            results_for_combined_output[column_name] <- df[column_name]
        } 
        # if not then it is the id column 
        if (grepl(global_variables["id_column_name"], column_name, fixed = TRUE)){
            results_for_combined_output[column_name] <- segid_to_join
        }
    }

    # Export combo version  
    cat(crayon::green("Saving Combined results\n"))
    data.table::fwrite(results_for_combined_output, paste0(directories["segment_output_path_combined_results"], "/bgz_",zone,"_combined_results.csv"))

    cat("Cleaning up\n")
    # clear up 
    rm(results_scaled, results_numeric, calibration_results, calibration_bootstrap_results_store, results_for_model_score, zonal_stats_seg, output_predict, data, results_for_calibrated_probs, cal_prob_prob_class, results_for_combined_output, df)
    h2o.removeAll()
}
# Shut down 
h2o.shutdown(prompt = FALSE)

