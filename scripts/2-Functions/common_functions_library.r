#' Function for install and loading a list of libraries 
#' @param list_of_libaries A character vector of libraries to be checked installed and loaded
#'
#' @param verbose default TRUE, prints status updates to terminal 
#'
#' @examples
#' # Install dplyr and tidyverse
#' install_and_load(c("dplyr", "tidyverse"))
install_and_load <- function(list_of_libaries, verbose = TRUE){
    if (verbose) {
        cat("Comparing requirements against system installed pacakges\n")
        # Find all packages that need to be installed
        to_install <- setdiff(list_of_libaries, rownames(installed.packages()))
        # If all installed then report else attempt to install missing packages
        if (length(to_install) == 0){
            cat("All required libraries installed\n")
        } else {
            install.packages(setdiff(list_of_libaries, rownames(installed.packages())))
        }

        # Check that all packages installed successfully and report if not 
        failed_to_install <- setdiff(list_of_libaries, rownames(installed.packages()))
        if (length(failed_to_install) == 0){
            cat("All libraries successfully installed\n")
        } else {
        cat("The following libraries failed to install and require manual installation\n", failed_to_install, "\n")
        }

        # Attempt to load all libraries
        cat("Loading required libraries into environment\n")
        for (ind_library in list_of_libaries){
            require(ind_library, character.only = TRUE)
        }

        # Check that all required libraries are loaded and report any that are not as they will require user intevention 
        failed_packages <- setdiff(list_of_libaries, (.packages()))
        if (length(failed_packages) == 0){
            cat(crayon::green("All libraries loaded successfully\n"))
        } else {
            cat("The following packages failed to install\n")
        }
    } else {
        # Find all packages that need to be installed
        to_install <- setdiff(list_of_libaries, rownames(installed.packages()))
        # If all installed then report else attempt to install missing packages
        if (length(to_install) != 0){
            install.packages(setdiff(list_of_libaries, rownames(installed.packages())))
        }

        # Attempt to load all libraries
        for (ind_library in list_of_libaries){
            require(ind_library, character.only = TRUE)
        }
    }
}

#' Function to check if the provided directory is valid and exists, can create if not found 
#' @param path A string containing the directory to check
#' @param create default True, will create the directory if directory is found to not exist
#'
#' @return A string containing the path to the validated directory
#' 
#' @examples
#' # Check if directory exists and create it if not 
#' validate_dir("C:/User/Desktop/test_dir")
validate_dir <- function(path, create = TRUE){
    # Make sure crayon is installed for coloured terminal output 
    install_and_load(c("crayon"), FALSE)

    # Check if folder exists at given path, if not and flagged to create then create 
    if (!dir.exists(file.path(path))) {
        if (create) {
            dir.create(file.path(path))
            cat(crayon::green(path, "didn't exist and has been created\n"))
            return(file.path(path))
        }
    cat(crayon::red(path, "doesn't exist, and has not created, change create_if_doesnt_exist to TRUE if you wish to create it\n"))
    return(file.path(path))
    } else {
        cat(crayon::green(path, "exists\n"))
        return(file.path(path))
    }
}

#' Function to check if the provided file directory is valid and exists
#' @param path A string containing the file directory to check
#'
#' @return A string containing the validated file path
#' 
#' @examples
#' # Check if directory exists and create it if not 
#' validate_file_path("C:/User/Desktop/test_dir/test.R")

# Check if file exists and if so then return checked path, if not raise error
validate_file_path <- function(path) {
    # Make sure crayon is installed for coloured terminal output 
    install_and_load(c("crayon"), FALSE)

    # Check if file path exists and reprot success
    if (file.exists(path)) {
        cat(crayon::green(path, "exists\n"))
        return(file.path(path))
    } else {
        cat(crayon::red(path, "doesn't exist\n"))
        return(NULL)
    }
}
