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
library(pracma)
#options(scipen = 999, digits = 2) # sig digits

# Import the basic functions
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
#' hypergeometricAnalysis(dataSet, FALSE)
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
                               ", Obs. = ", x, ", Expected = ", expectedBottlenecks, ")\n\n"),
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

#' Function to generate hypergeometric distribution
#'
#' @param dataSet_ Dataframe containing the data to be analysed.
#' @param p_value_ The probability to observe a statistic value higher than found.
#' @param removeZeroBottlenecks_ Flag to determine whether or not the bottlenecks without frequency will
#' be included into the analysis.
#' @param verbose_ Print every status message.
#'
#' @return This functions returns the hypergeometric distribution.
#'
#' @examples
#' \dontrun{
#' hypergeometricDistribution(dataSet)
#' hypergeometricDistribution(dataSet, p_value_ = 0.01)
#' hypergeometricDistribution(dataSet, p_value_ = 0.05, verbose_ = FALSE)
#' hypergeometricDistribution(dataSet, p_value_ = 0.05, removeZeroBottlenecks_ = TRUE, verbose_ = TRUE)
#' hypergeometricDistribution(dataSet, p_value_ = 0.01, removeZeroBottlenecks_ = FALSE, verbose_ = FALSE)
#' }
#'
#' @author
#' Clóvis F. Reis / Igor Brandão

hypergeometricDistribution <- function(dataSet_, p_value_ = 0.05, removeZeroBottlenecks_ = FALSE, verbose_ = TRUE) {

  # Status message
  if (verbose_) {
    printMessage(paste0("RUNNING THE HYPERGEOMETRIC DISTRIBUTION ANALYSIS..."))
  }

  # First of all, check whether or not remove the ZERO bottlenecks
  if (removeZeroBottlenecks_) {
    dataSet_ <- removeZeroBottlenecks(dataSet_)
    exportFile <- "distributionWithoutZeroBottleneck"
  } else {
    exportFile <- "distributionWithZeroBottleneck"
  }

  # Basic metrics
  mean <- mean(dataSet_$freq[dataSet_$freq>1])
  sd <- sd(dataSet_$freq[dataSet_$freq>1])

  # Order the dataSet
  dataSet_ <- dataSet_[order(dataSet_$percentage,decreasing = T),]

  # Get the unique frequency percentages
  percRange <- unique(round(dataSet_$percentage,1))

  # Count the bottlenecks and non-bottlenecks
  countsBase <- c(bottleneck=nrow(dataSet_[dataSet_$is_bottleneck ==1,]),
                non_bottleneck=nrow(dataSet_[dataSet_$is_bottleneck !=1,]))
  countsBase[1]+countsBase[2]

  # Create a dataFrame for the result
  distribution<-data.frame(perc=numeric(),
                     white=numeric(),
                     black=numeric(),
                     drawn=numeric(),
                     freq=numeric(),
                     hyp=numeric())

  # Loop over the unique percentages
  for (perc in percRange) {
    # Temporaly dataFrame for indexed results
    countsTop <- data.frame(perc=numeric(),
                            bottleneck=numeric(),
                            non_bottleneck=numeric(),
                            drawn=numeric(),
                            freq=numeric(),
                            hyp=numeric())

    # Get the enzymes with the highests frequencies
    top <- dataSet_[dataSet_$percentage>=perc,]

    # Number of draws
    drawn <- nrow(top)

    # The frequency itself
    freq <- nrow(top[top$is_bottleneck == 1,])

    countsTop[1,] <- t(c(perc,c(countsBase),drawn,freq,NA))

    countsTop[1,"hyp"] <- phyper(countsTop[1,"freq"],
                                 countsTop[1,"bottleneck"],
                                 countsTop[1,"non_bottleneck"],
                                 countsTop[1,"drawn"],
                                 lower.tail = F)

    # Bind each result
    distribution <- rbind(distribution,countsTop)
  }

  # Filter the result according to the p_value
  # result <- result[result$hyp<=p_value_&round(result$hyp,5)!=0,]

  # Status message
  if (verbose_) {
    printMessage(paste0("RESULT WITH ", nrow(distribution), " PERCENTAGES"))
  }

  #**********************************************************************************#

  # Status message
  if (verbose_) {
    printMessage("PLOTTING THE DISTRIBUTION...")
  }

  # Plot the hyper distribution
  g <- ggplot() + theme_bw() +
    xlab("Presence in species (%)") +
    ylab("Bottleneck (%)") +
    geom_line(data = distribution[distribution$perc > 0 & distribution$hyp <= 0.01, ],
              aes(x = perc, y = freq / drawn * 100), col = "#000099") +
    geom_line(data = distribution[distribution$perc > 0 & distribution$hyp > 0.01, ],
              aes(x = perc, y = freq / drawn * 100), col = "#b20000")

  # Export the hypergeometric analysis
  if (!dir.exists(file.path('./output/statistics/'))) {
    dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
    dir.create(file.path(paste0('./output/statistics/hypergeometric/')), showWarnings = FALSE, mode = "0775")
  }

  if (dir.exists(file.path('./output/statistics/hypergeometric/'))) {
    ggsave(paste0("./output/statistics/hypergeometric/", exportFile, ".png"), width = 20, height = 15, units = "cm")
    save(distribution, file=paste0("./output/statistics/hypergeometric/", exportFile, ".RData"))
  }

  # Return the result
  return(distribution)
}

#' Function to generate hypergeometric distribution with discrete ranges
#'
#' @param dataSet_ Dataframe containing the data to be analysed.
#' @param p_value_ The probability to observe a statistic value higher than found.
#' @param rangeInterval_ The quantity of range values
#' @param cumulative_ Flag to determine if the ranges will accumulate the significance values
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

hypergeometricDistributionDiscrete <- function(dataSet_, p_value_ = 0.05, rangeInterval_ = 10,
                                               cumulative_ = FALSE, verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage(paste0("RUNNING THE DISCRETE HYPERGEOMETRIC DISTRIBUTION ANALYSIS..."))
  }

  # First of all, remove the ZERO bottlenecks
  dataSet_ <- removeZeroBottlenecks(dataSet_)

  if (cumulative_) {
    exportFile <- "discreteDistributionCumulative"
  } else {
    exportFile <- "discreteDistribution"
  }

  # Filter dataSet from proteins with ZERO frequency
  dataSet_ <- dataSet_[!dataSet_$freq ==0,]

  # Order the dataSet
  dataSet_ <- dataSet_[order(dataSet_$percentage,decreasing = T),]

  # Get the unique frequency percentages
  rangeVal<- c(seq(1, nrow(dataSet_), ceiling(nrow(dataSet_) / rangeInterval_)), nrow(dataSet_))

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

  # Filter the result according to the p_value
  # distribution <- distribution[distribution$hyp<=p_value_&round(distribution$hyp,5)!=0,]

  # Adjust the p-value in order to compensate the accumulated error with
  # Benjamini-Hochberg method
  distribution$pCor <- p.adjust(distribution$hyp, method = "BH")
  nrow(distribution[distribution$pCor<=0.01,])

  # Status message
  if (verbose_) {
    printMessage(paste0("RESULT WITH ", nrow(distribution), " RANGES"))
  }

  #**********************************************************************************#

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

  # Export the hypergeometric discrete analysis
  if (!dir.exists(file.path('./output/statistics/'))) {
    dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
    dir.create(file.path(paste0('./output/statistics/hypergeometric/')), showWarnings = FALSE, mode = "0775")
  }

  if (dir.exists(file.path('./output/statistics/hypergeometric/'))) {
    ggsave(paste0("./output/statistics/hypergeometric/", exportFile, ".png"), width = 20, height = 15, units = "cm")
    save(distribution, file=paste0("./output/statistics/hypergeometric/", exportFile, ".RData"))
  }

  # Return the result
  return(distribution)
}

#*******************************************************************************************#

# ---- PIPELINE SECTION ----

#***************#
# Pipeline flow #
#***************#

#******************************#
# Step 1: Generate the dataSet #
#******************************#

# Load the dataSet
dataSet <- get(load("./output/statistics/hypergeometric/hypergeometric.RData"))

# OR

# Generate the dataSet
dataSet <- generateDataSet(testName_ = 'hypergeometric')

#*********************************************#
# Step 2: Perform the hypergeometric analysis #
#*********************************************#

hypergeometricAnalysis(dataSet)

#***************************************************#
# Step 2.1: Perform the hypergeometric distribution #
#***************************************************#

hypergeometricDistribution(dataSet, p_value_ = 0.01, removeZeroBottlenecks_ = TRUE)

#***********************************************************************#
# Step 2.2: Perform the hypergeometric distribution with discretization #
#***********************************************************************#

hypergeometricDistributionDiscrete(dataSet, p_value_ = 0.01, rangeInterval_ = 20,
                                   cumulative_ = TRUE, verbose_ = TRUE)
