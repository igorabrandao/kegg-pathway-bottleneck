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
library(ggpubr)
library(grid)
library(GGally)
library(corrplot)
library(dplyr)
library(pracma)
library(tidyr)
library(gghighlight)
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

hypergeometricDistribution <- function(dataSet_, p_value_ = 0.05, rangeInterval_ = 1, cumulative_ = TRUE, normalize_ = TRUE, verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage(paste0("RUNNING THE DISCRETE HYPERGEOMETRIC DISTRIBUTION ANALYSIS..."))
  }

  #**************************************************************************##
  # Apply the 2nd normalization:                                              #
  # 100% of occurrence in a pathway can be compared with 50% of other pathway #
  #**************************************************************************##

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

  #*****************************************##
  # Generate the hypergeometric distribution #
  #*****************************************##

  # Create a dataFrame for the result
  distribution <- data.frame(range=0, bottleneck=0, non_bottleneck=0, drawn=0, freq=0, percentage=0, hyp=0.01)

  # Get the unique frequency percentages
  rangeVal <- c(seq(1, nrow(dataSet_), ceiling(nrow(dataSet_) / rangeInterval_)), nrow(dataSet_))

  # Count the bottlenecks and non-bottlenecks
  countsBase <- c(bottleneck=nrow(dataSet_[dataSet_$is_bottleneck ==1,]), non_bottleneck=nrow(dataSet_[dataSet_$is_bottleneck !=1,]))

  # Verification
  if ((countsBase[1] +countsBase[2]) != nrow(dataSet_)) {
    # Save the log file
    printLog(message_='The sum of bottlenecks with non-bottlenecks is different from the dataSet total rows. Skipping it...',
             file_='hypergeometricDistribution')

    return(FALSE)
  }

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
                            percentage=numeric(),
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

    # Average frequencies
    percentage <- top[range,]$normalizedFrequency

    countsTop[1,] <- t(c(range, c(countsBase), drawn, freq, percentage, NA))

    countsTop[1,"hyp"] <- phyper(countsTop[1,"freq"],
                                 countsTop[1,"bottleneck"],
                                 countsTop[1,"non_bottleneck"],
                                 countsTop[1,"drawn"],
                                 lower.tail = F)

    # Bind each result
    distribution <- rbind(distribution,countsTop)
  }

  # Adjust the p-value in order to compensate the accumulated error with Benjamini-Hochberg method
  distribution$pCor <- p.adjust(distribution$hyp, method = "BH")

  # Set colors properties (significance)
  distribution$cor <- "red"
  distribution$cor[distribution$hyp <= p_value_] <- "blue"
  distribution$cor[distribution$hyp == 0.0] <- "red"
  distribution$cor[distribution$range == 0] <- "gray"

  # Set the continuity param to the graph
  distribution$contin <- 'n'

  # Status message
  if (verbose_) {
    printMessage(paste0("RESULT WITH ", nrow(distribution), " RANGES"))
  }

  #*************************************##
  # Plot the hypergeometric distribution #
  #*************************************##

  # Export the hypergeometric discrete analysis
  if (!dir.exists(file.path('./output/statistics/'))) {
    dir.create(file.path(paste0('./output/statistics/')), showWarnings = FALSE, mode = "0775")
    dir.create(file.path(paste0('./output/statistics/hypergeometric/')), showWarnings = FALSE, mode = "0775")
  }

  # Status message
  if (verbose_) {
    printMessage("PLOTTING THE DISTRIBUTION...")
  }

  # Set the loop indexes
  idx <- 1
  lastIdx <- (nrow(distribution) - 1)

  # Define whether or not the plot is continuous
  for (idx in 2:lastIdx) {
    if (distribution$cor[idx] == distribution$cor[idx + 1] | distribution$cor[idx] == distribution$cor[idx - 1]) {
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

  #*************##
  # Plot 1: Dots #
  #*************##

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

# Generate the dataSet
dataSet <- generateDataSetCSV(testName_ = 'hypergeometric')
dataSet <- fillPathwayCodeWithZeros(dataSet)
dataSet <- dataSet[!dataSet$occurrences==0,]

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
distribution <- hypergeometricDistribution(dataSet, p_value_ = 0.01, rangeInterval_ = nrow(dataSet),
                                            cumulative_ = TRUE, normalize_ = TRUE, verbose_ = TRUE)

#***********************************#
# Step 3: Generate additional plots #
#***********************************#

# Remove proteins with zero occurrence
enzymeList <- dataSet[!dataSet$occurrences==0,]

faixas <- data.frame(percent=numeric(), btn=numeric(), nbtn=numeric(), sig=numeric(), stringsAsFactors=F)
faixa = 0
step = 0.5

for (faixa in seq(0, 95, step)) {
  btn <- nrow(enzymeList[enzymeList$percentage > faixa & enzymeList$percentage <= faixa + 5 & enzymeList$is_bottleneck == 1, ])
  nbtn <- nrow(enzymeList[enzymeList$percentage > faixa & enzymeList$percentage <= faixa + 5 & enzymeList$is_bottleneck == 0, ])

  faixas[nrow(faixas) + 1, 1] <- faixa + 5
  faixas[nrow(faixas), 2] <- btn
  faixas[nrow(faixas), 3] <- nbtn

  faixas[nrow(faixas), 4] <- 0

  print(faixa)
}

faixas$btnN <- (faixas$btn-min(faixas$btn)) / (max(faixas$btn)-min(faixas$btn))
faixas$nbtnN <- (faixas$nbtn - min(faixas$nbtn)) / (max(faixas$nbtn)-min(faixas$nbtn))

# Classificacao dos APS que são significativos de acordo com a hipergeometrica
significativos <- distribution[distribution$cor=='blue',]
percentSig <- significativos[nrow(significativos),]$percentage * 100
faixas[faixas$percent >= percentSig,] $sig <- 1

tot <-sum(faixas$btn)+sum(faixas$nbtn)

p <- ggplot(data = faixas) +
  geom_line(aes(x=percent, y=btn/(btn+nbtn), col=sig)) +
  geom_line(aes(x=percent, y=nbtn/(btn+nbtn)), col="#ED553B") +
  geom_rect(aes(xmin = min(faixas[faixas$sig==1,]$percent), xmax = max(faixas[faixas$sig==1,]$percent), ymin = -Inf, ymax = Inf),
            fill = "pink", alpha = 0.002, linetype = "dashed") +
  geom_vline(xintercept=min(faixas[faixas$sig==1,]$percent), linetype='dashed', color="#E75480", size=0.5) +
  geom_vline(xintercept=max(faixas[faixas$sig==1,]$percent), linetype='dashed', color="#E75480", size=0.5) +

  # Chart visual properties
  xlab("") +
  ylab("") +
  ggtitle("") +
  scale_x_discrete(limits=seq(5, 100, by=5)) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", size=20, hjust = 0),
        axis.title.x = element_text(face="bold", size=20, margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=16),
        axis.title.y = element_text(face="bold", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(size=16),
        legend.title = element_text(face="bold", size=16),
        legend.text = element_text(size=16),
        legend.position = 'none')
p

ggsave(paste0("./output/statistics/hypergeometric/apXnap_dist.png"), width = 25, height = 20, units = "cm")

p <- ggplot(data = enzymeList) +
  # Non significant non APs
  geom_jitter(data = enzymeList[enzymeList$percentage < min(faixas[faixas$sig==1,]$percent) & enzymeList$is_bottleneck==0,],
              aes(y=percentage, x=1), color='#d0e4f4', fill='white', width=0.5, shape = 3, size = 3) +

  # Significant non APs
  geom_jitter(data = enzymeList[enzymeList$percentage >= min(faixas[faixas$sig==1,]$percent) & enzymeList$is_bottleneck==0,],
              aes(y=percentage, x=1), color='#81b6e0', fill='white', width=0.5, shape = 3, size = 3) +

  # Non Significant APs
  geom_jitter(data = enzymeList[enzymeList$percentage < min(faixas[faixas$sig==1,]$percent) & enzymeList$is_bottleneck==1,],
              aes(y=percentage, x=1), color='#2e7ebe', fill='white', width=0.5, shape = 23, size = 3) +

  # Significant APs
  geom_jitter(data = enzymeList[enzymeList$percentage >= min(faixas[faixas$sig==1,]$percent) & enzymeList$is_bottleneck==1,],
              aes(y=percentage, x=1), color='#173F5F', fill='#173F5F', width=0.5, shape = 23, size = 3) +

  geom_hline(yintercept=min(faixas[faixas$sig==1,]$percent), linetype='dashed', color="#E75480", size=1) +

  # Chart visual properties
  xlab("") +
  ylab("") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(face="bold", size=20, hjust = 0),
        axis.title.x = element_text(face="bold", size=20, margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=16),
        axis.title.y = element_text(face="bold", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(size=16),
        legend.title = element_text(face="bold", size=16),
        legend.text = element_text(size=16),
        legend.position = 'none')
p

ggsave(paste0("./output/statistics/hypergeometric/apXnap_scatter.png"), width = 25, height = 20, units = "cm")

p <- ggplot(data = faixas)+theme_bw()+
    geom_line(aes(x=percent,y=btn/(btn+nbtn)), col="green")+
    geom_line(aes(x=percent,y=nbtn/(btn+nbtn)), col="magenta")
p

p <- ggplot(data = faixas)+theme_bw()+
  geom_line(aes(x=percent,y=btnN/(btn+nbtn)), col="green")+
  geom_line(aes(x=percent,y=nbtnN/(btn+nbtn)), col="magenta")
p

p <- ggplot(data = faixas[faixas$percent>25,])+theme_bw()+
    geom_line(aes(x=percent,y=btnN), col="red")+
    geom_line(aes(x=percent,y=nbtnN), col="blue")
p

p <- ggplot(data = faixas[faixas$percent>25,])+theme_bw()+
  geom_line(aes(x=percent,y=btn), col="red")+
  geom_line(aes(x=percent,y=nbtn), col="blue")
p

ggsave(paste0("./output/statistics/hypergeometric/apXnap_tendency.png"), width = 25, height = 20, units = "cm")

p <- ggplot(data = faixas)+theme_bw() +
  geom_line(aes(x=percent,y=btn/(tot)), col="green") +
  geom_line(aes(x=percent,y=nbtn/(tot)), col="magenta")
p

faixas <- faixas[faixas$percent>15,]
faixas$btnN <- (faixas$btn-min(faixas$btn)) / (max(faixas$btn)-min(faixas$btn))
faixas$nbtnN <- (faixas$nbtn - min(faixas$nbtn)) / (max(faixas$nbtn)-min(faixas$nbtn))

p <- ggplot(data = faixas) +
  geom_tile(aes(y=percent, x=1, fill=(btnN))) +
  geom_tile(aes(y=percent, x=2), fill='white', color='white') +
  geom_tile(aes(y=percent, x=3, fill=(nbtnN))) +

  geom_hline(yintercept=min(faixas[faixas$sig==1,]$percent), linetype='dashed', color="#E75480", size=1) +

  # Chart visual properties
  xlab("") +
  ylab("") +
  ggtitle("") +
  theme_bw() +
  theme(plot.title = element_text(face="bold", size=20, hjust = 0),
        axis.title.x = element_text(face="bold", size=20, margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=16),
        axis.title.y = element_text(face="bold", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(size=16),
        legend.title = element_text(face="bold", size=16),
        legend.text = element_text(size=16),
        legend.position = 'none')
p

ggsave(paste0("./output/statistics/hypergeometric/apXnap_heatmap.png"), width = 25, height = 20, units = "cm")
