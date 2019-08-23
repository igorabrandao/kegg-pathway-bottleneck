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

# Import the graphLoader functions
files.sources = NULL
files.sources[1] = paste0("./R/functions", "/", "statisticsHelper.R")
files.sources[2] = paste0("./R/functions", "/", "helperFunctions.R")
sapply(files.sources, source)

#*******************************************************************************************#

# ---- HYPERGEOMETRIC SECTION ----

#' Function to generate hypergeometric analysis
#'
#' @param dataSet_ Dataframe containing the data to be analysed.
#' @param verbose_ Print every status message.
#'
#' @return This function does not return nothing, just export files
#' containing the analysis.
#'
#' @examples
#' \dontrun{
#' hypergeometricAnalysis(dataSet)
#' }
#'
#' @author
#' Igor Brandão

hypergeometricAnalysis <- function(dataSet_, verbose_ = TRUE) {
  # Handle empty dataSet
  if (is.null(dataSet_) | length(dataSet_) == 0) {
    # Status message
    if (verbose_) {
      printMessage("THE DATASET DOES NOT CONTAIN DATA. SKIPPING IT...")
    }

    return(NULL)
  } else {
    # Define the set of tests based on proteins frequencies
    # e.g: Most 20% frequents proteins
    testSet <- c(75, 50, 25, 20, 15, 10, 5, 2, 1, -2, -5, -10, -15, -20, -25, -50, -75)

    # Run each test instance
    for (idx in 1:length(testSet)) {
      # Define if it's TOP or LOW frequencies
      if (testSet[idx] > 0) {
        orderType = "TOP"
      } else {
        orderType = "LOW"
      }

      # Status message
      if (verbose_) {
        printMessage(paste0("RUNNING ", orderType, " ", abs(testSet[idx]), "%"))
      }

      # Select the most frequent proteins based on testSet parameter
      freq_percentual_rate_ <- abs(testSet[idx]/100)
      topFreq = NULL

      if (strcmp(orderType, "TOP")) {
        topFreq <- dataSet_[order(-dataSet_$freq),]
      } else {
        topFreq <- dataSet_[order(dataSet_$freq),]
      }

      topFreq <- topFreq[1:as.integer(nrow(topFreq) * freq_percentual_rate_),]

      #**********************************************************************************#
      # Test parameters
      # What is the probability of selecting x bottlenecks from a sample of k taken from an
      # dataSet containing m bottlenecks proteins and n non bottlenecks proteins?
      #**********************************************************************************#
      m = as.integer(nrow(topFreq[topFreq$is_bottleneck == 1, ]))  # Number of success states in the population
      n = as.integer(nrow(topFreq[topFreq$is_bottleneck == 0, ]))  # Number of unsuccessful states in the population
      k = as.integer(m)  # Number of draws (i.e. quantity drawn in each trial)
      x = as.integer(k / 2)  # Number of observed successes
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
      rangeDistance = 15
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
        labs(title = paste0("Probability mass function [PMF] of X = x Bottlenecks Proteins [", orderType, " ", abs(testSet[idx]), "% frequencies]"),
             subtitle = paste0("Hypergeometric(B = ", m, ", NB = ", n, ", Draws = ", k,
                               ", Obs. = ", x, ", Expected = ", expectedBottlenecks, ")\n\n",
                               "Prabability: ", hProbability, "\n\n",
                               "Var: ", format(round(var, 3), nsmall = 3)),
             x = "Number of bottlenecks proteins (x)",
             y = "Density")

      # Export the hypergeometric analysis
      if (!dir.exists(file.path('./output/statistics/'))) {
        dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
        dir.create(file.path(paste0('./output/statistics/hypergeometric/')), showWarnings = FALSE, mode = "0775")
      }

      if (dir.exists(file.path('./output/statistics/hypergeometric/'))) {
        ggsave(paste0("./output/statistics/hypergeometric/", tolower(orderType), abs(testSet[idx]), ".png"), width = 20, height = 15, units = "cm")
      }
    }
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
