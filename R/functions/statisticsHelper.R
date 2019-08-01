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

# ---- SETTINGS SECTION ----

# Folder containing the necessary info
folder_name = 'totalFrequency'

#*******************************************************************************************#

# ---- FUNCTIONS SECTION ----


#' Get the reference KEGG pathway
#'
#' Given a KEGG pathway ID and a list of KOs/ECs, this function returns a
#' data.frame ready to create an igraph object with ECs.
#'
#' @param pathway A KEGG pathway ID.
#' @param ko_ec_dictionnaire_ Data.frame containing the KOs/ECs equivalence
#'
#' @return This function returns a data.frame containing the edges from a
#' KEGG pathway by its ECs.
#'
#' @examples
#' \dontrun{
#' df <- pathwayToDataframe("00010", KO2EC)
#' }
#'
#' @importFrom KEGGREST keggGet
#' @importFrom KEGGgraph parseKGML
#' @importFrom KEGGgraph KEGGpathway2Graph
#'
#' @author
#' Igor Brandão

generateDataSet <- function(verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage("GENERATING THE DATASET BASE")
  }

  # Get the list of files
  folder = paste0("./output/", folder_name, "/")
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
    enzymeList <- enzymeList[, c('pathway','freq','percentage','is_bottleneck','bottleneck_classification')]

    # Status message
    if (verbose_) {
      printMessage("DATASET GENERATED WITH SUCCESS!")
      printMessage("SAVING THE GENERATED DATASET...")
    }

    # Export the pathway data
    if (!dir.exists(file.path('./output/statistics/'))) {
      dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
      dir.create(file.path(paste0('./output/statistics/hypergeometric/')), showWarnings = FALSE, mode = "0775")
    }

    if (dir.exists(file.path('./output/statistics/hypergeometric/'))) {
      save(enzymeList, file = paste0('./output/statistics/hypergeometric/hypergeometric.RData'))
    }

    if (dir.exists(file.path('~/data3/'))) {
      save(enzymeList, file = paste0('~/data3/kegg-pathway-bottleneck/output/statistics/hypergeometric/hypergeometric.RData'))
    }

    # Return the generated dataSet
    return(enzymeList)
  }
}

#*******************************************************************************************#
