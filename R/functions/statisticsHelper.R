#***************************************#
# Functions used in statistics analysis #
#***************************************#

# ---- IMPORT SECTION ----

# statisticsHelper.R #

#' This is set of functions to help the
#' statistics analysis
#'
#' @author
#' Igor Brandão

# Import the necessary libraries

#*******************************************************************************************#

# ---- FUNCTIONS SECTION ----

#' Get all data from a folder and bind it together
#'
#' @param test_name_ The name of the test that is calling this function
#' @param folder_name_ The folder name that contains the necessary data
#' @param verbose_ Print every status message.
#'
#' @return This function returns a data frame containing the data from all dataSets inside a folder
#'
#' @examples
#' \dontrun{
#' dataSet <- generateDataSet()
#' }
#'
#' @author
#' Igor Brandão

generateDataSet <- function(test_name_ = '', folder_name_ = 'totalFrequency', verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage("GENERATING THE DATASET BASE")
  }

  # Get the list of files
  folder = paste0("./output/", folder_name_, "/")
  file_list <- list.files(path = folder, pattern = '*.RData')

  # Check if the folder contains files
  if (is.null(file_list) | length(file_list) == 0) {
    return(FALSE)
  }

  # Load all files at once
  big.list.of.data.frames <- lapply(file_list, function(file) {
    get(load(file = paste0(folder, file)))
  })

  # Combine multiple data frames in one
  enzymeList <- do.call(rbind, big.list.of.data.frames)

  # Remove temporaly variables
  rm(big.list.of.data.frames)

  # Handle empty graph
  if (is.null(enzymeList) | length(enzymeList) == 0) {
    # Status message
    if (verbose_) {
      printMessage("AN ERROR OCCURRED DURING THE DATASET GENERATION")
    }

    return(NULL)
  } else {
    # Remove unnecessary columns
    enzymeList <- enzymeList[, c('pathway','freq','total_species','percentage','is_bottleneck','bottleneck_classification')]

    # Status message
    if (verbose_) {
      printMessage(paste0(toupper(test_name_), " DATASET GENERATED WITH SUCCESS!"))
      printMessage("SAVING THE GENERATED DATASET...")
    }

    # Export the pathway data
    if (!dir.exists(file.path('./output/statistics/'))) {
      dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
      dir.create(file.path(paste0('./output/statistics/', test_name_, '/')), showWarnings = FALSE, mode = "0775")
    }

    if (dir.exists(file.path('./output/statistics/', test_name_, '/'))) {
      save(enzymeList, file = paste0('./output/statistics/', test_name_, '/', test_name_, '.RData'))
    }

    if (dir.exists(file.path('~/data3/'))) {
      save(enzymeList, file = paste0('~/data3/kegg-pathway-bottleneck/output/statistics/', test_name_, '/',
                                     test_name_, '.RData'))
    }

    # Return the generated dataSet
    return(enzymeList)
  }
}

#*******************************************************************************************#
