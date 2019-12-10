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
library(ggpubr)
library(gghighlight)
library(grid)
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

# Rename the AP classification
dataSet[dataSet$bottleneck_classification=='HB',]$bottleneck_classification <- 'HAP'
dataSet[dataSet$bottleneck_classification=='HNB',]$bottleneck_classification <- 'HUB'
dataSet[dataSet$bottleneck_classification=='NHB',]$bottleneck_classification <- 'AP'
dataSet[dataSet$bottleneck_classification=='NHNB',]$bottleneck_classification <- 'Others'

#******************************************#
# Step 2: Perform the descriptive analysis #
#******************************************#

# Default descriptive analysis
descriptiveAnalysis(dataSet, removeZeroBottlenecks_ = TRUE, verbose_ = TRUE,
                    columns_ = c("occurrences", "totalSpecies", "percentage", "is_bottleneck", "bottleneckNormalizedImpact", "bottleneck_classification"),
                    columnLabels_ = c("Occurence", "Processed Species", "Frequency (%)", "Is AP", "AP Impact", "AP Classification"))

#******************************************#

#************************#
# Protein classification #
#************************#

data <- dataSet[!dataSet$occurrences==0,]

ggplot(data, aes(fill=bottleneck_classification, x=bottleneck_classification)) +
  # Add the bars
  geom_bar(position="stack", stat="count") +

  # Add the label text
  geom_text(stat='count', aes(label = paste0(..count.., ' (', format( (..count.. / nrow(data)) * 100, digits=3), '%)')),
            size=6, fontface="bold", vjust=-1) +

  # Chart visual properties
  xlab("Proteins Group") +
  ylab("Proteins Count") +
  ggtitle("") +
  guides(fill=guide_legend(title="Proteins Group")) +
  scale_fill_manual(values = c("#ED553B", "#3CAEA3", "#20639B", "#173F5F")) +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=20, margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=18),
        axis.title.y = element_text(face="bold", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(size=18),
        legend.title = element_text(face="bold", size=16),
        legend.text = element_text(size=16))

ggsave(paste0("./output/statistics/descriptive/proteinClassification.png"), width = 33, height = 25, units = "cm")

#******************************************#

#**********************************************#
# Proteins by pathway without zero bottlenecks #
#**********************************************#

data <- dataSet[!dataSet$occurrences==0,]

# Count the proteins by pathway
proteinsCount <- as.data.frame(table(data$pathway), stringsAsFactors=FALSE)

# Order the data by totalSpecies
proteinsCount <- proteinsCount[order(proteinsCount$Freq),]

# Rename the group columns
names(proteinsCount)[names(proteinsCount) == "Var1"] <- "pathway"
names(proteinsCount)[names(proteinsCount) == "Freq"] <- "proteinsCount"

# Order the dataSet by the protein count
data$pathway <- factor(data$pathway, levels = proteinsCount$pathway[order(proteinsCount$proteinsCount)])

proteinByPathway1 <- ggplot(data, aes(fill=bottleneck_classification, x=pathway)) +
    # Add the bars
    geom_bar(position="stack", stat="count", width = 0.75) +

    # Chart visual properties
    xlab("") +
    #xlab("Pathways") +
    ylab("") +
    #ylab("Proteins Count") +
    ggtitle("") + theme_bw() +
    guides(fill=guide_legend(title="Proteins Classification")) +
    scale_fill_manual(values = c("#ED553B", "#3CAEA3", "#20639B", "#173F5F")) +
    theme(axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        #axis.title.x = element_text(face="bold", size=20, margin = margin(t = 15, r = 0, b = 0, l = 0)),
        #axis.text.x = element_text(size=9, angle = 90, hjust = 1),
        axis.text.x = element_blank(),
        #axis.title.y = element_text(face="bold", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        #axis.title.y = element_text(face="bold", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", size=20, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=18),
        legend.title = element_text(face="bold", size=16),
        legend.text = element_text(size=16),
        legend.position='top') +
  annotation_custom(grobTree(textGrob("A", x=0.02,  y=0.90, hjust=0, gp=gpar(col="black", fontsize=30, fontface="bold"))))

#ggsave(paste0("./output/statistics/descriptive/proteinByPathway.png"), width = 30, height = 20, units = "cm")

proteinByPathway2 <- ggplot(data, aes(fill=bottleneck_classification, x=pathway)) +
  # Add the bars
  geom_bar(position="stack", stat="count", width = 0.75) +

  # Highlighted bars
  gghighlight(pathway == '00402' | pathway == '00232' | pathway == '00061', unhighlighted_colour = alpha("steelblue", 0.4)) +

  # Add labels to highlighted bars
  geom_label(stat='count', aes(y = ..count..+20, label = pathway), label.size = 1.5, hjust = 1, vjust = 1, fill = "#173F5F", colour = "white", alpha= 1) +

  # Chart visual properties
  xlab("Pathways") +
  ylab("") +
  #ylab("Proteins Count") +
  ggtitle("") + theme_bw() +
  guides(fill=guide_legend(title="Proteins Classification")) +
  scale_fill_manual(values = c("#ED553B", "#3CAEA3", "#20639B", "#173F5F")) +
  theme(axis.title.x = element_text(face="bold", size=20, margin = margin(t = 15, r = 0, b = 0, l = 0)),
        #axis.text.x = element_text(size=9, angle = 90, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face="bold", size=20, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size=18),
        legend.title = element_text(face="bold", size=16),
        legend.text = element_text(size=16),
        legend.position='top') +
  annotation_custom(grobTree(textGrob("B", x=0.02,  y=0.90, hjust=0, gp=gpar(col="black", fontsize=30, fontface="bold"))))

#ggsave(paste0("./output/statistics/descriptive/proteinByPathwayHighlight.png"), width = 30, height = 20, units = "cm")

figure <- ggarrange(proteinByPathway1, proteinByPathway2, heights = c(3, 3), ncol = 1, nrow = 2, align = "v", legend = "top", common.legend = TRUE)

annotate_figure(figure,
                left = text_grob("Proteins Count", color = "#000000", size=20, face = "bold", rot = 90))

ggsave(paste0("./output/statistics/descriptive/proteinByPathwayArrange.png"), width = 30, height = 20, units = "cm")

#******************************************#

#***********************************#
# Bottleneck x Betweenness X Degree #
#***********************************#

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

#**********************************************#
# Organisms by pathway without zero bottlenecks #
#**********************************************#

data <- dataSet[!dataSet$occurrences==0,]

# Aggregate pathways total species
orgByPath <- aggregate(data$totalSpecies, by=list(data$pathway), FUN=mean, stringsAsFactors=FALSE)

# Rename the group columns
names(orgByPath)[names(orgByPath) == "x"] <- "totalSpecies"
names(orgByPath)[names(orgByPath) == "Group.1"] <- "pathway"

# Order the data by totalSpecies
orgByPath <- orgByPath[order(orgByPath$totalSpecies),]

orgByPath$pathway <- factor(orgByPath$pathway, levels = orgByPath$pathway[order(orgByPath$totalSpecies)])

## Use n equally spaced breaks to assign each value to n-1 equal sized bins
colors <- c()
intervals <- 4
rangeVal <- c(seq(1, nrow(orgByPath), ceiling(nrow(orgByPath) / intervals)), nrow(orgByPath))
pallete <- c("#ED553B", "#3CAEA3", "#20639B", "#173F5F")

for (idx in 1:length(pallete)) {
  colors[rangeVal[idx]:rangeVal[idx+1]] <- pallete[idx]
}

# Classify the data
orgByPath$group <- ''
orgByPath$color <- ''

for (idx in 1:(length(rangeVal)-1) ) {
  # Set the group color
  orgByPath[rangeVal[idx]:rangeVal[idx+1],]$color <- colors[rangeVal[idx]:rangeVal[idx+1]]

  # Set the group value
  orgByPath[rangeVal[idx]:rangeVal[idx+1],]$group <- paste0('< ', orgByPath[rangeVal[idx+1],]$totalSpecies, ' organisms')
}

# Plot rganisms by pathways
ggplot(orgByPath) +
  # Add the bars
  geom_bar(aes(x=pathway, y=totalSpecies, fill=group), color='#f6f6f6', stat="identity") +

  # Add labels to bars group
  #geom_label(data = orgByPath[rangeVal,], aes(y = (totalSpecies + 200), x=pathway, label=totalSpecies, fill=group), label.size = 1.5, hjust = 1, vjust = 1, colour = "white", alpha= 1) +

  # Chart visual properties
  xlab("Pathways") +
  ylab("Organism Count") +
  ggtitle("") +
  scale_fill_manual(values = unique(colors[rangeVal])) +
  guides(fill=guide_legend(title="Pathways Group")) +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=20, margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.text.x = element_blank(),
        axis.title.y = element_text(face="bold", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(size=18),
        legend.title = element_text(face="bold", size=16),
        legend.text = element_text(size=16),
        legend.position = 'top') +
  geom_vline(xintercept = which.min(abs(orgByPath$totalSpecies - 2000)), linetype='dashed', color="red", size=0.5) +
  geom_text(aes(x=which.min(abs(orgByPath$totalSpecies - 2000)), label="\nLess than 2000 organisms", y=4500), colour="#173F5F", angle=90,
            text=element_text(size=20))

ggsave(paste0("./output/statistics/descriptive/organismsByPathway.png"), width = 30, height = 20, units = "cm")

#***************************#
# Step 3: Correlation study #
#***************************#

# Generate the correlation study
generateCorrelationStudy(dataSet, removeZeroBottlenecks_ = TRUE, verbose_ = TRUE)
