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
library(corrplot)
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
#' descriptiveAnalysis(dataSet)
#' descriptiveAnalysis(dataSet, removeZeroBottlenecks_ = TRUE, verbose_ = TRUE)
#' }
#'
#' @author
#' Igor Brandão

descriptiveAnalysis <- function(dataSet_, removeZeroBottlenecks_ = FALSE, columns_ = NULL,
                                columnLabels_ = NULL, labelAngle_ = 45, title_ = NULL,
                                exportFile_ = NULL, verbose_ = TRUE) {
  # Status message
  if (verbose_) {
    printMessage(paste0("RUNNING THE DESCRIPTIVE ANALYSIS..."))
  }

  # First of all, check whether or not remove the ZERO bottlenecks
  if (removeZeroBottlenecks_) {
    dataSet_ <- removeZeroBottlenecks(dataSet_)

    # Filter dataSet from proteins with ZERO frequency
    dataSet_ <- dataSet_[!dataSet_$occurrences ==0,]

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

  # Define which dataSet column will be displayed into the plot
  if (is.null(columns_) | length(columns_) == 0) {
    columns = c("occurrences", "totalSpecies", "percentage", "is_bottleneck", "bottleneck_classification")
  } else {
    columns = columns_
  }

  if (is.null(columnLabels_) | length(columnLabels_) == 0) {
    columnLabels = c("Frequency", "Processed Species", "Frequency (%)", "Is Bottleneck?", "Bottleneck Classification")
  } else {
    columnLabels = columnLabels_
  }

  if (!is.null(title_)) {
    plotTitle <- title_
  }

  # Drawing a scatterplot matrix of freq, totalSpecies, percentage, and is_bottleneck using the pairs function
  plot1 <- ggpairs(dataSet_, columns = columns, columnLabels = columnLabels, title = plotTitle,
                   mapping = aes(color = bottleneck_classification),
                   lower = list(
                     continuous = "smooth",
                     combo = "facetdensity"
                   ), cardinality_threshold = 1000) +
    theme(axis.text.x = element_text(angle = labelAngle_, hjust = 1)) + theme_bw()

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
    if (!is.null(exportFile_)) {
      exportFile <- exportFile_
    }

    ggsave(paste0("./output/statistics/descriptive/", exportFile, ".png"), width = 25, height = 20, units = "cm")
  }
}

#' Function to generate the correlation study
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
#' descriptiveAnalysis(dataSet, verbose_ = FALSE)
#' }
#'
#' @author
#' Igor Brandão

generateCorrelationStudy <- function(dataSet_, removeZeroBottlenecks_=TRUE, verbose_=TRUE) {
  # Status message
  if (verbose_) {
    printMessage(paste0("GENERATING CORRELATION STUDY..."))
  }

  # First of all, check whether or not remove the ZERO bottlenecks
  if (removeZeroBottlenecks_) {
    # Filter dataSet from proteins with ZERO frequency
    dataSet_ <- dataSet_[!dataSet_$occurrences ==0,]
  }

  # Backup the complete dataSet
  dataSetBkp_ <- dataSet_

  # Remove non numeric data
  dataSet_ <- dataSet_[ , -which(names(dataSet_) %in% c("bottleneck_classification","pathway"))]

  # Remove the pathwayData unnecessary columns
  dataSet_ <- dataSet_[,!(names(dataSet_) %in% c('X', 'entryID', 'name', 'type', 'link', 'reaction', 'x', 'y',
                                                 'graphicalType', 'width', 'height', 'fgcolor', 'bgcolor', 'reaction_type',
                                                 'org', 'bottleneck_classification', 'pathway'))]

  #*********************#
  # General correlation #
  #*********************#

  # Generate the correlation matrix
  correlationMatrix <- cor(dataSet_, method = c("spearman"), use = "complete.obs")

  # Define the color scale
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

  # Save the correlogram (plot1)
  png(file="./output/statistics/descriptive/correlationAll.png", res=300, width=3500, height=3500)
  corrplot(correlationMatrix, method="color",
           addCoef.col = "black", # Add coefficient of correlation
           tl.col = "black", tl.srt = 45, #Text label color and rotation
           # hide correlation coefficient on the principal diagonal
           diag=TRUE)
  dev.off()

  #**********************#
  # Correlation by group #
  #**********************#

  # Generate the correlation by protein classificaion
  correlationByGroup <- by(dataSetBkp_, dataSetBkp_$bottleneck_classification,
                           FUN = function(X) X[ , -which(names(X) %in% c('X', 'entryID', 'name', 'type', 'link', 'reaction', 'x', 'y',
                                                                         'graphicalType', 'width', 'height', 'fgcolor', 'bgcolor', 'reaction_type',
                                                                         'org', 'bottleneck_classification', 'pathway'))])

  generateCorrelationByGroup <- function(idx_) {
    # Define the filename
    filename <- names(correlationByGroup[idx_])

    # Save the group dataSet
    dataSetByGroup <- as.data.frame(correlationByGroup[idx_])
    save(dataSetByGroup, file = paste0('./output/statistics/descriptive/descriptive', filename, '.RData'))

    # Generate the correlation matrix
    correlationMatrixGroup <- cor(as.data.frame(correlationByGroup[idx_]), method = c("spearman"), use = "complete.obs")

    # Define the color scale
    col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

    # Save the correlogram (plot1)
    png(file=paste0("./output/statistics/descriptive/correlationGroup", filename, ".png"), res=300, width=3500, height=3500)
    corrplot(correlationMatrixGroup, method="color",
             addCoef.col = "black", # Add coefficient of correlation
             tl.col = "black", tl.srt = 45, #Text label color and rotation
             # hide correlation coefficient on the principal diagonal
             diag=TRUE)
    dev.off()
  }

  # Run the group correlation
  lapply(1:nrow(correlationByGroup), generateCorrelationByGroup)

  # Status message
  if (verbose_) {
    printMessage(paste0("CORRELATION STUDY GENERATED WITH SUCCESS!"))
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
dataSet <- generateDataSetCSV(testName_ = 'descriptive', filterColumns_ = FALSE)
dataSet <- fillPathwayCodeWithZeros(dataSet)

#******************************************#
# Step 2: Perform the descriptive analysis #
#******************************************#

# Default descriptive analysis
descriptiveAnalysis(dataSet, removeZeroBottlenecks_ = TRUE, verbose_ = TRUE,
                    columns_ = c("occurrences", "totalSpecies", "percentage", "is_bottleneck", "bottleneck_classification"),
                    columnLabels_ = c("Frequency", "Processed Species", "Frequency (%)", "Is Bottleneck", "Protein Classification"))

#******************************************#

# Protein classification
data <- dataSet[!dataSet$occurrences==0,]

ggplot(data, aes(fill=bottleneck_classification, x=bottleneck_classification)) +
  geom_bar(position="stack", stat="count") +
  xlab("Proteins Group") +
  ylab("Proteins Count") +
  ggtitle("Proteins Classification") + theme_bw() +
  guides(fill=guide_legend(title="Proteins Group")) +
  scale_fill_manual(values = c("#ED553B", "#3CAEA3", "#20639B", "#173F5F")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0("./output/statistics/descriptive/proteinClassification.png"), width = 30, height = 20, units = "cm")

#******************************************#

# Proteins by pathway without zero bottlenecks
data <- dataSet[!dataSet$occurrences==0,]

ggplot(data, aes(fill=bottleneck_classification, x=pathway)) +
    geom_bar(position="stack", stat="count") +
    xlab("Pathways") +
    ylab("Proteins Count") +
    ggtitle("Proteins by pathway") + theme_bw() +
    guides(fill=guide_legend(title="Proteins Classification")) +
    scale_fill_manual(values = c("#ED553B", "#3CAEA3", "#20639B", "#173F5F")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0("./output/statistics/descriptive/proteinByPathway.png"), width = 40, height = 20, units = "cm")

#******************************************#

# Bottleneck x Betweenness X Degree
descriptiveAnalysis(dataSet, removeZeroBottlenecks_ = TRUE, verbose_ = TRUE,
                    columns_ = c("is_bottleneck", "betweenness", "degree"),
                    columnLabels_ = c("Is Bottleneck", "Betweenness", "Degree"),
                    title_ = 'Bottleneck x Betweenness X Degree',
                    exportFile_ = 'bottleneckBetweennessDegree')

#******************************************#

# Graph metrics
descriptiveAnalysis(dataSet, removeZeroBottlenecks_ = TRUE, verbose_ = TRUE,
                    columns_ = c("authorityScore", "betweenness", "closenessCoef", "clusteringCoef", "community", "connectivity", "degree", "diameter", "eccentricity"
                                 ,"radius", "triangles"),
                    columnLabels_ = c("authorityScore", "betweenness", "closenessCoef", "clusteringCoef", "community", "connectivity", "degree", "diameter", "eccentricity"
                                      ,"radius", "triangles"),
                    labelAngle_ = 90, title_ = 'Graph Metrics',
                    exportFile_ = 'graphMetrics')

#******************************************#

# Proteins by pathway without zero bottlenecks
data <- dataSet[!dataSet$occurrences==0,]

# Aggregate pathways total species
orgByPath <- aggregate(data$totalSpecies, by=list(data$pathway), FUN=mean, stringsAsFactors=FALSE)

# Rename the group columns
names(orgByPath)[names(orgByPath) == "x"] <- "totalSpecies"
names(orgByPath)[names(orgByPath) == "Group.1"] <- "pathway"

## Use n equally spaced breaks to assign each value to n-1 equal sized bins
ii <- cut(orgByPath$totalSpecies, breaks = seq(min(orgByPath$totalSpecies), max(orgByPath$totalSpecies), len = 20),
          include.lowest = TRUE)
## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
colors <- colorRampPalette(c("#20639B", "#0c273e"))(19)[ii]

# Plot rganisms by pathways
ggplot(orgByPath) +
  geom_bar(aes(x=pathway, y=totalSpecies), color='#f6f6f6', fill=colors, stat="identity") +
  xlab("Pathways") +
  ylab("Organism Count") +
  ggtitle("Organism by pathway") + theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0("./output/statistics/descriptive/organismsByPathway.png"), width = 40, height = 20, units = "cm")

#***************************#
# Step 3: Correlation study #
#***************************#

# Generate the correlation study
generateCorrelationStudy(dataSet, removeZeroBottlenecks_ = TRUE, verbose_ = TRUE)
