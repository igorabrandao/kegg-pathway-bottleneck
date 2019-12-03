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
library(GGally)
library(corrplot)
library(dplyr)
library(pracma)
#options(scipen = 999, digits = 2) # sig digits

# Import the basic functions
files.sources = NULL
files.sources[1] = paste0("./R/functions", "/", "statisticsHelper.R")
files.sources[2] = paste0("./R/functions", "/", "helperFunctions.R")
sapply(files.sources, source)

#*******************************************************************************************#

# ---- HYPERGEOMETRIC SECTION ----

#' Function to generate hypergeometric distribution with discrete ranges
#'
#' @param dataSet_ Dataframe containing the data to be analysed.
#' @param p_value_ The probability to observe a statistic value higher than found.
#' @param rangeInterval_ The quantity of range values
#' @param cumulative_ Flag to determine if the ranges will accumulate the significance values
#' @param normalize_ Flag to determine if the proteins frequencies will be normalized by its pathway size
#' @param verbose_ Print every status message.
#'
#' @return This functions returns the hypergeometric distribution.
#'
#' @examples
#' \dontrun{
#' hypergeometricDistributionDiscrete(dataSet)
#' hypergeometricDistributionDiscrete(dataSet, p_value_ = 0.01)
#' hypergeometricDistributionDiscrete(dataSet, p_value_ = 0.05, verbose_ = FALSE)
#' hypergeometricDistributionDiscrete(dataSet, p_value_ = 0.05, rangeInterval_ = 10, verbose_ = FALSE)
#' hypergeometricDistributionDiscrete(dataSet, p_value_ = 0.05, rangeInterval_ = 10, cumulative_ = FALSE, verbose_ = FALSE)
#' }
#'
#' @author
#' Clóvis F. Reis / Igor Brandão

hypergeometricDistribution <- function(dataSet_, p_value_ = 0.05, rangeInterval_ = 20,
                                               cumulative_ = FALSE, normalize_ = FALSE, verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage(paste0("RUNNING THE DISCRETE HYPERGEOMETRIC DISTRIBUTION ANALYSIS..."))
  }

  # Apply the pathway size normalization
  if (normalize_) {
    # Get the unique pathways
    uniquePathways <- unique(dataSet_$pathway)

    # Calculates the normalized frequency
    for (item in uniquePathways) {
      # Max frequency in a pathway
      pathwayMaxFrequency <- max(dataSet_$percentage[dataSet_$pathway==item])

      # Min frequency in a pathway
      pathwayMinFrequency <- min(dataSet_$percentage[dataSet_$pathway==item])

      # Normalized frequency for each protein
      dataSet_$normalizedFrequency[dataSet_$pathway==item] <-
        (dataSet_$percentage[dataSet_$pathway==item]-pathwayMinFrequency)/
        (pathwayMaxFrequency-pathwayMinFrequency)
    }

    # Fix for NAN cases (when min and max frequency have the same values)
    dataSet_$normalizedFrequency[is.nan(dataSet_$normalizedFrequency)] <- dataSet_$percentage[is.nan(dataSet_$normalizedFrequency)]/100
  }

  # Filter dataSet from proteins with ZERO frequency
  dataSet_ <- dataSet_[!dataSet_$occurrences==0,]

  if (cumulative_) {
    exportFile <- "hypergeometricDistributionCumulative"
  } else {
    exportFile <- "hypergeometricDistribution"
  }

  # Order the dataSet
  dataSet_ <- dataSet_[order(dataSet_$normalizedFrequency, decreasing = T),]

  # Get the unique frequency percentages
  rangeVal <- c(seq(1, nrow(dataSet_), ceiling(nrow(dataSet_) / rangeInterval_)), nrow(dataSet_))

  # Count the bottlenecks and non-bottlenecks
  countsBase <- c(bottleneck=nrow(dataSet_[dataSet_$is_bottleneck ==1,]),
                  non_bottleneck=nrow(dataSet_[dataSet_$is_bottleneck !=1,]))
  countsBase[1]+countsBase[2]

  # Create a dataFrame for the result
  distribution <- data.frame(range=0,
                             bottleneck=0,
                             non_bottleneck=0,
                             drawn=0,
                             freq=0,
                             hyp=0.01)

  # Define the loop ranges
  range <- 1
  lastRange <- (length(rangeVal) - 1)

  for (range in 1:lastRange) {
    # Temporaly dataFrame for indexed results
    countsTop <- data.frame(range=numeric(),
                            bottleneck=numeric(),
                            non_bottleneck=numeric(),
                            drawn=numeric(),
                            freq=numeric(),
                            hyp=numeric())

    # Check the method to divide the range values
    if (cumulative_) {
      # Cumulative ranges 1:range
      initVal <- 1
    } else {
      # Non cumulative ranges previous range:range
      initVal <- rangeVal[range]
    }

    # Check if it is the last range of values
    if (range == (length(rangeVal) - 1)) {
      top <- dataSet_[initVal:(rangeVal[range + 1]), ]
    } else{
      top <- dataSet_[initVal:(rangeVal[range + 1] - 1), ]
    }

    # Number of draws
    drawn <- nrow(top)

    # The frequency itself
    freq <- nrow(top[top$is_bottleneck == 1,])

    countsTop[1,] <- t(c(range, c(countsBase), drawn, freq, NA))

    countsTop[1,"hyp"] <- phyper(countsTop[1,"freq"],
                                 countsTop[1,"bottleneck"],
                                 countsTop[1,"non_bottleneck"],
                                 countsTop[1,"drawn"],
                                 lower.tail = F)

    # Bind each result
    distribution <- rbind(distribution,countsTop)
  }

  # Adjust the p-value in order to compensate the accumulated error with
  # Benjamini-Hochberg method
  distribution$pCor <- p.adjust(distribution$hyp, method = "BH")
  nrow(distribution[distribution$pCor<=0.01,])

  # Status message
  if (verbose_) {
    printMessage(paste0("RESULT WITH ", nrow(distribution), " RANGES"))
  }

  #**********************************************************************************#

  # Export the hypergeometric discrete analysis
  if (!dir.exists(file.path('./output/statistics/'))) {
    dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
    dir.create(file.path(paste0('./output/statistics/hypergeometric/')), showWarnings = FALSE, mode = "0775")
  }

  # Status message
  if (verbose_) {
    printMessage("PLOTTING THE DISTRIBUTION...")
  }

  # Set colors properties
  distribution$cor <- "red"
  distribution$cor[distribution$pCor <= 0.01] <- "blue"
  distribution$cor[distribution$pCor == 0.0] <- "red"
  distribution$cor[distribution$range == 0] <- "gray"

  # Set the continuity param to the graph
  distribution$contin <- 'n'

  # Set the loop indexes
  idx <- 1
  lastIdx <- (nrow(distribution) - 1)

  for (idx in 2:lastIdx) {
    if (distribution$cor[idx] == distribution$cor[idx + 1] |
        distribution$cor[idx] == distribution$cor[idx - 1]) {
      distribution$contin[idx] <- 's'
    } else{
      distribution$contin[idx] <- 'n'
    }
  }

  # Set the continuity of the last range
  idx = nrow(distribution)
  if (distribution$cor[idx] == distribution$cor[idx - 1]) {
    distribution$contin[idx] <- 's'
  }

  g <- ggplot() + theme_bw() +
    xlab("Proteins") +
    ylab("Bottlenecks") +
    geom_point(data = distribution[1, ], color = "gray", pch = 3,
              aes(x = drawn, y = freq)) +
    geom_line(data = distribution[distribution$contin == 's', ],
              aes(x = drawn, y = freq, color = cor)) +
    geom_point(data = distribution[distribution$cor != "gray", ],
              aes(x = drawn, y = freq, color = cor)) +
    scale_color_manual(values = c("blue" = "blue", "red" = "red", "gray" = "gray"),
              labels = c("blue" = "Significant", "red" = "Non Significant", "gray" = NA ),
              name = "p-value")

  g

  if (dir.exists(file.path('./output/statistics/hypergeometric/'))) {
    ggsave(paste0("./output/statistics/hypergeometric/", exportFile, rangeInterval_, "bins.png"), width = 20, height = 15, units = "cm")
    write.csv(distribution, file=paste0("./output/statistics/hypergeometric/", exportFile, ".csv"))
  }

  #**********************************************************************************#

  inicio = 1
  range = 1

  result<-data.frame(range=numeric(),
                     count=numeric(),
                     isBtn=numeric(),
                     signif=numeric())

  while (inicio <= nrow(dataSet_)) {
    fim = inicio + (nrow(dataSet_)/rangeInterval_)

    if (fim > nrow(dataSet_)) {
      fim <- nrow(dataSet_)
    }

    tmp <- dataSet_[inicio:fim, ]
    btn = nrow(tmp[tmp$is_bottleneck == 1, ])
    nbtn = nrow(tmp[tmp$is_bottleneck == 0, ])
    signif <- ifelse(distribution$cor[distribution$range == range] == "blue", 1, 0)

    if (signif == 1) {
      b <- 1
      n <- 0
    } else{
      b <- 3
      n <- 2
    }
    result[nrow(result) + 1, ] <- c(range, btn, b, signif)

    inicio <- fim + 1
    range <- range + 1
  }

  p <- ggplot() +
    geom_col(data=result[result$signif==1&result$isBtn==1,],
             aes(x = range, y = count, fill = isBtn), width = 0.9) +
    geom_col(data=result[result$signif==0&result$isBtn==3,],
             aes(x = range, y = count, fill = isBtn), width = 0.9) +
    geom_line(data=distribution,
              aes(x = range, y = freq/drawn*50), col="red") +
    xlab("Ranges") +
    ylab("Bottlenecks count")
  p

  if (dir.exists(file.path('./output/statistics/hypergeometric/'))) {
    ggsave(paste0("./output/statistics/hypergeometric/hypergeometricDistribution", rangeInterval_, "bins.png"), width = 20, height = 15, units = "cm")
    write.csv(distribution, file=paste0("./output/statistics/hypergeometric/", exportFile, ".csv"))
  }

  # Return the result
  return(distribution)
}

#' Function to generate additional plots of hypergeomtric
#'
#' @param dataSet_ Dataframe containing the data to be analysed.
#' @param removeZeroBottlenecks_ Flag to determine whether or not the bottlenecks without frequency will be included into the analysis.
#' @param columns_ A vector of strings containing the dataSet columns name to be plotted.
#' @param columnLabels_ The label associated with the plotted columns.
#' @param labelAngle_ Define the angle of bottom labels.
#' @param title_ The plot title.
#' @param exportFile_ The exported plot filename.
#' @param verbose_ Print every status message.
#'
#' @return This functions returns nothing.
#'
#' @examples
#' \dontrun{
#' generateAdditionalPlots(dataSet)
#' generateAdditionalPlots(dataSet, removeZeroBottlenecks_ = TRUE, verbose_ = TRUE)
#' }
#'
#' @author
#' Igor Brandão

generateAdditionalPlots <- function(dataSet_, columns_ = NULL,
                                columnLabels_ = NULL, labelAngle_ = 45, title_ = NULL,
                                exportFile_ = NULL, verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage("GENERATING ADDITIONAL PLOTS...")
  }

  # Define which dataSet column will be displayed into the plot
  if (is.null(columns_) | length(columns_) == 0) {
    columns = c("range", "bottleneck", "non_bottleneck", "drawn", "freq", "hyp", "pCor")
  } else {
    columns = columns_
  }

  if (is.null(columnLabels_) | length(columnLabels_) == 0) {
    columnLabels = c("Range", "Qtd. of bottlenecks", "Qtd. of non bottlenecks",
                     "Drawns", "Bottleneck frequency", "Hypergeomtric", "Adjusted p-value")
  } else {
    columnLabels = columnLabels_
  }

  if (!is.null(title_)) {
    plotTitle <- title_
  }

  # Drawing a scatterplot matrix of freq, total_species, percentage, and is_bottleneck using the pairs function
  plot1 <- ggpairs(dataSet_, columns = columns, columnLabels = columnLabels, title = plotTitle,
                   mapping = aes(color = cor),
                   lower = list(
                     continuous = "smooth",
                     combo = "facetdensity"
                   ), cardinality_threshold = 1000) +
    theme(axis.text.x = element_text(angle = labelAngle_, hjust = 1)) + theme_bw()

  # Recolor the matrix
  for(i in 1:plot1$nrow) {
    for (j in 1:plot1$ncol) {
      plot1[i, j] <- plot1[i, j] +
        scale_fill_manual(values = c("#173F5F", "#4D4E4F", "#ED553B")) +
        scale_color_manual(values = c("#173F5F", "#4D4E4F", "#ED553B"))
    }
  }

  # Export the hypergeometric hypergeometric analysis
  if (!dir.exists(file.path('./output/statistics/'))) {
    dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
  }

  if (!dir.exists(file.path('./output/statistics/hypergeometric/'))) {
    dir.create(file.path(paste0('./output/statistics/hypergeometric/')), showWarnings = FALSE, mode = "0775")
  }

  if (dir.exists(file.path('./output/statistics/hypergeometric/'))) {
    print(plot1)
    if (!is.null(exportFile_)) {
      exportFile <- exportFile_
    }

    ggsave(paste0("./output/statistics/hypergeometric/", exportFile, ".png"), width = 25, height = 20, units = "cm")
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
dataSet <- generateDataSetCSV(testName_ = 'hypergeometric')
dataSet <- fillPathwayCodeWithZeros(dataSet)

#*********************************************************************#
# Step 2: Perform the hypergeometric distribution with discretization #
#*********************************************************************#

#' Call the hypergeometric distribution with the follwing parameters:
#' p-value = 0.01
#' intervals = 20 (each one representing 5% of the total data)
#' cumulative = TRUE (each range will accumulate the quantities from the previous ranges)
#' normalize = TRUE (the proteins frequency will be normalized by its pathway size)
#' verbose = TRUE (all status messages will be shown)
#'
distribution <- hypergeometricDistribution(dataSet, p_value_ = 0.01, rangeInterval_ = 100,
                                   cumulative_ = TRUE, normalize_ = TRUE, verbose_ = TRUE)

#***********************************#
# Step 3: Generate additional plots #
#***********************************#

# Default hypergeometric analysis
generateAdditionalPlots(distribution, verbose_ = TRUE,
                    columns_ = c("range", "bottleneck", "non_bottleneck", "drawn", "freq", "hyp", "pCor"),
                    columnLabels_ = c("Range", "Qtd. of bottlenecks", "Qtd. of non bottlenecks",
                                      "Drawns", "Bottleneck frequency", "Hypergeomtric", "Adjusted p-value"),
                    title_ = 'Hypergeometric Overview',
                    exportFile_ = 'hyperOverview')
