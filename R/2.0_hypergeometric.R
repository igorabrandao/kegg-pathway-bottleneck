#*********************************************#
# Pipeline to perform hypergeometric analysis #
#*********************************************#

# ---- IMPORT SECTION ----

# 2.0_hypergeometric #

#' This is the pipeline script to perform
#' the hypergeometric analysis
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library(ggplot2)
library(dplyr)
options(scipen = 999, digits = 2) # sig digits

#*******************************************************************************************#

# ---- SETTINGS SECTION ----

#*************************#
# Pipeline basic settings #
#*************************#

# Import the graphLoader functions
files.sources = NULL
files.sources[1] = paste0("./R/functions", "/", "helperFunctions.R")
sapply(files.sources, source)

# Folder containing the necessary info
folder_name = 'totalFrequency'

#*******************************************************************************************#

# ---- DATASET SECTION ----

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
    enzymeList <-
      enzymeList[, c('pathway',
                     'freq',
                     'percentage',
                     'is_bottleneck',
                     'bottleneck_classification')]

    # Status message
    if (verbose_) {
      printMessage("DATASET GENERATED WITH SUCCESS!")
      printMessage("SAVING THE GENERATED DATASET...")
    }

    # Export the pathway data
    if (!dir.exists(file.path('./output/statistics/'))) {
      dir.create(file.path(paste0('./output/statistics/')),
                 showWarnings = FALSE,
                 mode = "0775")
      dir.create(file.path(paste0(
        './output/statistics/hypergeometric/'
      )),
      showWarnings = FALSE,
      mode = "0775")
    }

    if (dir.exists(file.path('./output/statistics/hypergeometric/'))) {
      save(
        enzymeList,
        file = paste0(
          './output/statistics/hypergeometric/hypergeometric.RData'
        )
      )
    }

    if (dir.exists(file.path('~/data3/'))) {
      save(
        enzymeList,
        file = paste0(
          '~/data3/kegg-pathway-bottleneck/output/statistics/hypergeometric/hypergeometric.RData'
        )
      )
    }

    # Return the generated dataSet
    return(enzymeList)
  }
}

#*******************************************************************************************#

# ---- HYPERGEOMETRIC SECTION ----

hypergeometricAnalysis <- function(dataSet_, verbose_ = TRUE) {
  # Handle empty dataSet
  if (is.null(dataSet_) | length(dataSet_) == 0) {
    # Status message
    if (verbose_) {
      printMessage("THE DATASET DOES NOT CONTAIN DATA. SKIPPING IT...")
    }

    return(NULL)
  } else {
    #**********************************************************************************#
    # Test parameters
    # What is the probability of selecting x bottlenecks from a sample of k taken from an
    # dataSet containing m bottlenecks proteins and n non bottlenecks proteins?
    #**********************************************************************************#

    m = as.integer(nrow(dataSet_[dataSet_$is_bottleneck == 1, ]))  # Number of success states in the population
    n = as.integer(nrow(dataSet_[dataSet_$is_bottleneck == 0, ]))  # Number of unsuccessful states in the population
    k = as.integer(m)  # Number of draws (i.e. quantity drawn in each trial)
    x = as.integer(m / 2)  # Number of observed successes
    #**********************************************************************************#

    # Status message
    if (verbose_) {
      printMessage("CALCULATING THE H-PROBABILITY")
    }

    # Returns the cumulative probability (percentile) p at the specified value (quantile) q
    hProbability = dhyper(x = x, m = m, n = n, k = k)

    # Expected number of bottlenecks proteins
    expectedBottlenecks = as.integer(k * m / (m + n))

    # Variance
    var = k * m / (m + n) * (m + n - k) / (m + n) * n / (m + n - 1)

    #**********************************************************************************#

    # Status message
    if (verbose_) {
      printMessage("PLOTTING THE DISTRIBUTION")
    }

    # Determine the distribution range
    rangeDistance = 30
    range = (expectedBottlenecks-rangeDistance):(expectedBottlenecks+rangeDistance)

    # Calculates the distribution for a range
    density = dhyper(x = range, m = m, n = n, k = k)

    # Plot the bottlenecks distribution
    data.frame(numBottlenecks = range, density) %>%
      mutate(legend = ifelse(numBottlenecks == expectedBottlenecks, paste0("x = ", expectedBottlenecks), "others")) %>%
      ggplot(aes(x = factor(numBottlenecks), y = density, fill = legend)) +
      geom_col() +
      geom_text(
        aes(label = round(density,2), y = density + 0.01),
        position = position_dodge(0.9),
        size = 3,
        vjust = 0
      ) +
      labs(title = "Probability mass function [PMF] of X = x Bottlenecks Proteins",
           subtitle = paste0("Hypergeometric(k = ", k, ", M = ", m, ", N = ", n, ")"),
           x = "Number of bottlenecks proteins (x)",
           y = "Density")
  }
}

#*******************************************************************************************#

# ---- PIPELINE SECTION ----

#***************#
# Pipeline flow #
#***************#

#******************************#
# Step 1: Generate the dataSet #
#******************************#

# Data frame to receive the generated data
dataSet <- generateDataSet()

#*********************************************#
# Step 2: Perform the hypergeometric analysis #
#*********************************************#

hypergeometricAnalysis(dataSet)
