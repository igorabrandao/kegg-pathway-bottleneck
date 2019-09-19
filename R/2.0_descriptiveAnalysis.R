#******************************************#
# Pipeline to perform descriptive analysis #
#******************************************#

# ---- IMPORT SECTION ----

# 2.0_descriptiveAnalysis.R #

#' This is the pipeline script to perform
#' the descriptive analysis
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library(ggplot2)
library(GGally)
library(dplyr)

# Import the basic functions
files.sources = NULL
files.sources[1] = paste0("./R/functions", "/", "statisticsHelper.R")
files.sources[2] = paste0("./R/functions", "/", "helperFunctions.R")
sapply(files.sources, source)

#*******************************************************************************************#

# ---- DESCRIPTIVE SECTION ----

#' Function to generate hypergeometric distribution with Descriptive ranges
#'
#' @param dataSet_ Dataframe containing the data to be analysed.
#' @param removeZeroBottlenecks_ Flag to determine whether or not the bottlenecks without frequency will be included into the analysis.
#' @param verbose_ Print every status message.
#'
#' @return This functions returns nothing.
#'
#' @examples
#' \dontrun{
#' descriptiveAnalysis(dataSet)
#' descriptiveAnalysis(dataSet, removeZeroBottlenecks_ = TRUE, verbose_ = TRUE)
#' }
#'
#' @author
#' Igor Brandão

descriptiveAnalysis <- function(dataSet_, removeZeroBottlenecks_ = FALSE, verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage(paste0("RUNNING THE DESCRIPTIVE ANALYSIS..."))
  }

  # First of all, check whether or not remove the ZERO bottlenecks
  if (removeZeroBottlenecks_) {
    dataSet_ <- removeZeroBottlenecks(dataSet_)

    # Filter dataSet from proteins with ZERO frequency
    dataSet_ <- dataSet_[!dataSet_$freq ==0,]

    exportFile <- "descriptiveAnalysisWithoutZeroBottleneck"
    plotTitle <- "Descriptive analysis Without Zero Bottleneck"
  } else {
    exportFile <- "descriptiveAnalysis"
    plotTitle <- "Descriptive analysis"
  }

  # Display or print the first 10 observations of our dataset
  print(head(dataSet_, n = 10))

  # To see the variable in the dataset; Use names(dataset) or ls(dataset)
  print(ls(dataSet_))

  # To see the number of columns and number of rows in the Prestige dataset; use ncol(dataset) and nrow(dataset)
  print(ncol(dataSet_))
  print(nrow(dataSet_))

  # A more advanced and complete way to see the structure of our dataset
  str(dataSet_)

  # Summarize the data
  print(summary(dataSet_))
  print(summary(as.factor(dataSet_$pathway)))
  print(summary(as.factor(dataSet_$is_bottleneck)))
  print(summary(as.factor(dataSet_$bottleneck_classification)))

  #**********************************************************************************#

  # Status message
  if (verbose_) {
    printMessage("PLOTTING THE DISTRIBUTION...")
  }

  # Drawing a scatterplot matrix of freq, total_species, percentage, and is_bottleneck using the pairs function
  plot1 <- ggpairs(dataSet_, columns = c("freq", "total_species", "percentage", "is_bottleneck", "bottleneck_classification"),
                   columnLabels = c("Frequency", "Processed Species", "Frequency (%)", "Is Bottleneck?", "Bottleneck Classification"),
                   title = plotTitle,
                   mapping = aes(color = bottleneck_classification),
                   lower = list(
                     continuous = "smooth",
                     combo = "facetdensity"
                   )) + theme_bw()

  # Recolor the matrix
  for(i in 1:plot1$nrow) {
    for (j in 1:plot1$ncol) {
      plot1[i, j] <- plot1[i, j] +
        scale_fill_manual(values = c("#ED553B", "#3CAEA3", "#20639B", "#173F5F")) +
        scale_color_manual(values = c("#ED553B", "#3CAEA3", "#20639B", "#173F5F"))
    }
  }

  # Export the hypergeometric Descriptive analysis
  if (!dir.exists(file.path('./output/statistics/'))) {
    dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
  }

  if (!dir.exists(file.path('./output/statistics/descriptive/'))) {
    dir.create(file.path(paste0('./output/statistics/descriptive/')), showWarnings = FALSE, mode = "0775")
  }

  if (dir.exists(file.path('./output/statistics/descriptive/'))) {
    print(plot1)
    ggsave(paste0("./output/statistics/descriptive/", exportFile, ".png"), width = 20, height = 15, units = "cm")
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

# Generate the dataSet
dataSet <- generateDataSet(test_name_ = 'descriptive', filter_columns_ = FALSE)

# OR

# Load the dataSet
dataSet <- get(load("./output/statistics/descriptive/descriptive.RData"))

#******************************************#
# Step 2: Perform the descriptive analysis #
#******************************************#

descriptiveAnalysis(dataSet, removeZeroBottlenecks_ = TRUE, verbose_ = TRUE)
descriptiveAnalysis(dataSet, removeZeroBottlenecks_ = FALSE, verbose_ = TRUE)
