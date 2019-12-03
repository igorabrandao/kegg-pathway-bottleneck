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
#' @param testName_ The name of the test that is calling this function
#' @param folderName_ The folder name that contains the necessary data
#' @param filterColumns_ Flag to determine whether or not the dataSet should be filtered into default columns
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

generateDataSet <- function(testName_ = '', folderName_ = 'totalFrequency', filterColumns_ = TRUE, verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage("GENERATING THE DATASET BASE")
  }

  # Get the list of files
  folder = paste0("./output/", folderName_, "/")
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
  dataSet <- do.call(rbind, big.list.of.data.frames)

  # Remove temporaly variables
  rm(big.list.of.data.frames)

  # Handle empty graph
  if (is.null(dataSet) | length(dataSet) == 0) {
    # Status message
    if (verbose_) {
      printMessage("AN ERROR OCCURRED DURING THE DATASET GENERATION")
    }

    return(NULL)
  } else {
    # Remove unnecessary columns
    if (filterColumns_) {
      dataSet <- dataSet[, c('pathway','freq','total_species','percentage','is_bottleneck','bottleneck_classification')]
    }

    # Status message
    if (verbose_) {
      printMessage(paste0(toupper(testName_), " DATASET GENERATED WITH SUCCESS!"))
      printMessage("SAVING THE GENERATED DATASET...")
    }

    # Export the pathway data
    if (!dir.exists(file.path('./output/statistics/'))) {
      dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
    }

    if (!dir.exists(file.path(paste0('./output/statistics/', testName_, '/')))) {
      dir.create(file.path(paste0('./output/statistics/', testName_, '/')), showWarnings = FALSE, mode = "0775")
    }

    if (dir.exists(file.path('./output/statistics/', testName_, '/'))) {
      save(dataSet, file = paste0('./output/statistics/', testName_, '/dataSet.RData'))
    }

    if (dir.exists(file.path('~/data3/'))) {
      save(dataSet, file = paste0('~/data3/kegg-pathway-bottleneck/output/statistics/', testName_, '/dataSet.RData'))
    }

    # Return the generated dataSet
    return(dataSet)
  }
}

generateDataSetCSV <- function(testName_ = '', folderName_ = 'totalFrequency', filterColumns_ = TRUE, verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage("GENERATING THE DATASET BASE")
  }

  # Get the list of files
  folder = paste0("./output/", folderName_, "/")
  file_list <- list.files(path = folder, pattern = '*.csv')

  # Check if the folder contains files
  if (is.null(file_list) | length(file_list) == 0) {
    return(FALSE)
  }

  # Load all files at once
  big.list.of.data.frames <- lapply(file_list, function(file) {
    read.csv(file=paste0(folder, file), header=TRUE, sep=",", stringsAsFactors=FALSE)
  })

  # Combine multiple data frames in one
  dataSet <- do.call(rbind, big.list.of.data.frames)

  # Remove temporaly variables
  rm(big.list.of.data.frames)

  # Handle empty graph
  if (is.null(dataSet) | length(dataSet) == 0) {
    # Status message
    if (verbose_) {
      printMessage("AN ERROR OCCURRED DURING THE DATASET GENERATION")
    }

    return(NULL)
  } else {
    # Status message
    if (verbose_) {
      printMessage(paste0(toupper(testName_), " DATASET GENERATED WITH SUCCESS!"))
      printMessage("SAVING THE GENERATED DATASET...")
    }

    # Export the pathway data
    if (!dir.exists(file.path('./output/statistics/'))) {
      dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
    }

    if (!dir.exists(file.path(paste0('./output/statistics/', testName_, '/')))) {
      dir.create(file.path(paste0('./output/statistics/', testName_, '/')), showWarnings = FALSE, mode = "0775")
    }

    if (dir.exists(file.path('./output/statistics/', testName_, '/'))) {
      write.csv(dataSet, file = paste0('./output/statistics/', testName_, '/dataSet.csv'))
    }

    # Return the generated dataSet
    return(dataSet)
  }
}

# In order to avoid bias into the analysis, the bottlenecks with ZERO frequency should be removed
#'
#' @param dataSet_ Entrez number withou specie
#' @examples
#' \dontrun{
#' dataSet_ <- removeZeroBottlenecks(dataSet_, verbose_ = FALSE)
#' }
#' @author
#' Igor Brandão

removeZeroBottlenecks <- function(dataSet_, verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage(paste0("REMOVING BOTTLENECKS WITHOUT FREQUENCY..."))
  }

  # First of all, get the list of pathways that contains protein bottlenecks with ZERO frequency
  pathwaysWithZeroBottleneck <- unique(dataSet_[dataSet_$freq == 0 & dataSet_$is_bottleneck == 1, ]$pathway)

  # Dataframe to receive the zero bottleneck data
  zeroBottleneckDf <- data.frame(pathway = numeric(), zeroBottleneckPerc = numeric(), zeroBottleneck = numeric(),
                                 bottleneckNonZero = numeric(), nonBottleneckZero = numeric(), nonBottleneckNonZero = numeric(),
                                 allProteins = numeric())

  for (pathway in pathwaysWithZeroBottleneck) {
    # Get the data related to bottlenecks with ZERO frequency
    zeroBottleneck <- nrow(dataSet_[dataSet_$freq == 0 & dataSet_$is_bottleneck == 1 & dataSet_$pathway == pathway,])

    # Get the data related to bottlenecks with non ZERO frequency
    bottleneckNonZero <- nrow(dataSet_[dataSet_$freq != 0 & dataSet_$is_bottleneck == 1 & dataSet_$pathway == pathway,])

    # Get the data related to non bottlenecks with ZERO frequency
    nonBottleneckZero <- nrow(dataSet_[dataSet_$freq == 0 & dataSet_$is_bottleneck == 0 & dataSet_$pathway == pathway,])

    # Get the data related to non bottlenecks with non ZERO frequency
    nonBottleneckNonZero <- nrow(dataSet_[dataSet_$freq != 0 & dataSet_$is_bottleneck == 0 & dataSet_$pathway == pathway,])

    # Get the data related to all proteins
    allProteins <- nrow(dataSet_[dataSet_$pathway == pathway,])

    # Apply the values into zeroBottleneckDf
    zeroBottleneckDf[nrow(zeroBottleneckDf) + 1, "pathway"] <- pathway
    zeroBottleneckDf[nrow(zeroBottleneckDf), "zeroBottleneckPerc"] <- (zeroBottleneck / allProteins)
    zeroBottleneckDf[nrow(zeroBottleneckDf), "zeroBottleneck"] <- zeroBottleneck
    zeroBottleneckDf[nrow(zeroBottleneckDf), "bottleneckNonZero"] <- bottleneckNonZero
    zeroBottleneckDf[nrow(zeroBottleneckDf), "nonBottleneckZero"] <- nonBottleneckZero
    zeroBottleneckDf[nrow(zeroBottleneckDf), "nonBottleneckNonZero"] <- nonBottleneckNonZero
    zeroBottleneckDf[nrow(zeroBottleneckDf), "allProteins"] <- allProteins
  }

  # Export the zeroBottleneck data
  if (!dir.exists(file.path('./output/statistics/'))) {
    dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
    dir.create(file.path(paste0('./output/statistics/hypergeometric/')), showWarnings = FALSE, mode = "0775")
  }

  if (dir.exists(file.path('./output/statistics/hypergeometric/'))) {
    save(zeroBottleneckDf, file=paste0("./output/statistics/hypergeometric/zeroBottleneckPathways.RData"))
  }

  # Remove the pathways containing ZERO bottlenecks from the dataSet
  dataSet_ <- dataSet_[!(dataSet_$pathway%in%pathwaysWithZeroBottleneck),]

  # Return the result
  return(dataSet_)
}

fillPathwayCodeWithZeros <- function(dataSet_, verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage(paste0("FILLING THE PATHWAYS CODE WITH ZEROS..."))
  }

  # Count the dataSet rows
  rows <- nrow(dataSet_)

  # Loop over the reference pathway rows
  for (idx in 1:rows) {
    current_pathway <- as.integer(dataSet_[idx,]$pathway)

    if (current_pathway >= 1000) {
      dataSet_[idx,]$pathway <- paste0('0', dataSet_[idx,]$pathway)
    } else if (current_pathway >= 100 & current_pathway < 1000) {
      dataSet_[idx,]$pathway <- paste0('00', dataSet_[idx,]$pathway)
    } else {
      dataSet_[idx,]$pathway <- paste0('000', dataSet_[idx,]$pathway)
    }
  }

  # Return the result
  return(dataSet_)
}

#*******************************************************************************************#
