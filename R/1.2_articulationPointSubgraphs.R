#**********************************************************************************************#
# Pipeline to perform the evaluation of subgraphs related to articulation points disconnection #
#**********************************************************************************************#

# ---- IMPORT SECTION ----

# 1.2_articulationPointSubgraphs.R #

#' This is the pipeline script to perform
#' the evaluation of subgraphs related to articulation points disconnection
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library(ggplot2)
library(svglite)
library(RColorBrewer)
library(plyr)

library(DataExplorer)

#*******************************************************************************************#

# ---- SETTINGS SECTION ----

#*************************#
# Pipeline basic settings #
#*************************#

# Import the graphLoader functions
files.sources = NULL
files.sources[1] = paste0("./R/functions", "/", "graphFunctions.R")
files.sources[2] = paste0("./R/functions", "/", "kgmlFunctions.R")
files.sources[3] = paste0("./R/functions", "/", "helperFunctions.R")
files.sources[4] = paste0("./R/functions", "/", "graphPrintFunctions.R")
files.sources[5] = paste0("./R/functions", "/", "statisticsHelper.R")
sapply(files.sources, source)

# Load the pathways by organisms data
organism2pathway <- get(load(paste0("./dictionaries", "/", "organism2pathway.RData")))
pathwayList <- get(load(paste0("./dictionaries", "/", "pathwayList.RData")))
pathwayDetail <- get(load(paste0("./dictionaries", "/", "pathwayDetail.RData")))

# Defined according to the paper plot
groupThresholdMax <- 80
groupThresholdMin <- 30

#*******************************************************************************************#

# ---- FUNCTIONS SECTION ----

#*************************#
# Pipeline main functions #
#*************************#

#' Parse the KGML file, simulates the graph disconnection via its APs and calculates the subgraphs metrics
#'
#' @param removeNoise_ Remove undesirable enzyme such as: ko, map, path, cpd or gl.
#'
#' @return This function returns nothing, just export .csv files.
#'
#' @examples
#' \dontrun{
#' calculateAPSubGraphs()
#' calculateAPSubGraphs(FALSE)
#' }
#'
#' @author
#' Igor Brandão

calculateAPSubGraphs <- function(removeNoise_=TRUE) {

  # Status message
  printMessage(paste0("PERFORMING THE PATHWAYS DISCONNECTED SUB-GRAPHS EVALUATION..."))

  # Reference pathway
  reference_pathway <- 'ec'

  # Get the list of files
  folder = paste0("./output/kgml/", reference_pathway, "/")
  kgml_list <- list.files(path=folder, pattern='*.xml')
  kgml_index <- 1

  # Define the number of available pathways
  available_pathways <- length(kgml_list)

  # Load the dicionaty
  dictionary <- read.csv(file='./output/pathwaysDictionary/dictionary.csv', header=TRUE, sep=",", stringsAsFactors=FALSE)

  # Load the protein dataSet to split the groups (<30% & >=80%)
  dataSet <- read.csv(file='./output/statistics/hypergeometric/dataSet.csv', header=TRUE, sep=",", stringsAsFactors=FALSE)

  if (is.null(dictionary) | nrow(dictionary) == 0) {
    # Save the log file
    printLog(message_='The pathways nodes dictionary could not be found. Skipping it...', file_='calculateAPSubGraphs')
    return(FALSE)
  }

  if (is.null(dataSet) | nrow(dataSet) == 0) {
    # Save the log file
    printLog(message_='The nodes dataset could not be found. Skipping it...', file_='calculateAPSubGraphs')
    return(FALSE)
  }

  # Loop 01: Run through all available pathways kgml
  lapply(kgml_list, function(file) {

    # Get the pathway code
    pathway_code <- onlyNumber(file)

    # Status message
    printMessage(paste0("SIMULATING ", pathway_code, " SUB-GRAPHS DISCONNECTION [", kgml_index, " OF ", available_pathways, "]"))

    #*************************************************##
    # Load all instances (orgs) of the current pathway #
    #*************************************************##

    # Get the list of instance files
    pathway_folder = paste0("./output/", pathway_code, "/")
    file_list <- grep(list.files(path=pathway_folder), pattern='*.csv', value=T)

    # Load all csv files at once
    big.list.of.data.frames <- lapply(file_list, function(item) {
      read.csv(file=paste0(pathway_folder, item), header=TRUE, sep=",", stringsAsFactors=FALSE)
    })

    # Combine multiple organisms data frames in one
    pathwayInstancesDataSet <- do.call(rbind, big.list.of.data.frames)

    # Remove temporaly variables
    rm(big.list.of.data.frames)

    #*************************************************************************##
    # Create the final dataSet with protein frequencies and network properties #
    #*************************************************************************##

    tryCatch({
      # Load the current pathway dataframe
      current_kgml <- KGML2Dataframe(paste0(folder, file))

      # Convert the pathway data into a graph
      pathwayGraph <- KGML2GraphDictionary(paste0(folder, file), replaceOrg=TRUE, orgToReplace=reference_pathway)

      if (is.null(pathwayGraph) | isempty(pathwayGraph)) {
        # Save the log file
        printLog(message_=paste0('The network ', pathway_code, ' [', kgml_index, ' of ', available_pathways,  '] data frame is empty. Skipping it...'), file_='calculateAPSubGraphs')

        # Increment the index
        kgml_index <<- kgml_index + 1

        return(FALSE)
      }

      #*************************##
      # Prepare the pathway data #
      #*************************##

      # Remove unnecessary data from pathway data/graph
      if (removeNoise_) {
        pathwayGraph <- removeNoise(pathwayGraph)
      }

      # Convert the graph into iGraph object
      iGraph <- igraph::graph_from_data_frame(pathwayGraph, directed = TRUE)

      #*****************************************************##
      # Run the sub-graphs disconnection and get its metrics #
      #*****************************************************##

      # Apply the node bottleneck impact
      impact <- getArticulationPointSubGraphs(pathwayGraph)

      # Check if some result was returned
      if (!is.null(impact) && length(impact) > 0) {
        # Set other data
        impact['org'] <- pathwayGraph[1:nrow(impact),]$org
        impact['pathway'] <- pathwayGraph[1:nrow(impact),]$pathway
        impact['pathway_nodes'] <- length(V(iGraph))
        impact['pathway_edges'] <- length(E(iGraph))
        impact['pathway_radius'] <- min(igraph::eccentricity(iGraph))
        impact['pathway_diameter'] <- max(igraph::eccentricity(iGraph))
        colnames(impact)[1] <- 'ap_dict_id'
        impact['ap_ec'] <- NA

        if (!is.null(impact) && nrow(impact) > 0) {
          for (idx in 1:nrow(impact)) {
            # Set the AP ec according to the dictionary
            impact[idx, 'ap_ec'] <- dictionary[dictionary$id==impact[idx, 'ap_dict_id'],]$ec
          }
        }

        # Generate the exporting dataset
        result <- impact[,(ncol(impact)-6):ncol(impact)]
        result['ap_dict_id'] <- as.integer(impact$ap_dict_id)
        result['ap_group'] <- NA
        result['ap_percentage'] <- NA
        result <- merge(result, impact[,1:(ncol(impact)-7)], by='ap_dict_id', all=T)

        # Perform the AP classification according to the reference dataset
        for (idx in 1:nrow(result)) {
          result[idx, 'ap_percentage'] <- dataSet[which(dataSet$dictID==result[idx,]$ap_dict_id),]$percentage

          if (dataSet[which(dataSet$dictID==result[idx,]$ap_dict_id),]$percentage >= groupThresholdMax) {
            result[idx, 'ap_group'] <- paste0('>=', groupThresholdMax)
          } else if (dataSet[which(dataSet$dictID==result[idx,]$ap_dict_id),]$percentage < groupThresholdMax &
                     dataSet[which(dataSet$dictID==result[idx,]$ap_dict_id),]$percentage >= groupThresholdMin) {
            result[idx, 'ap_group'] <- paste0(groupThresholdMin, '<=x<', groupThresholdMax)
          } else {
            result[idx, 'ap_group'] <- paste0('<', groupThresholdMin)
          }
        }

        #****************************************#
        # Calculates the custom periphery metric #
        #****************************************#

        # TODO: Put the calculations here!

        #****************************#
        # Prepare the data to export #
        #****************************#

        # Status message
        printMessage(paste0("EXPORTING SUB-GRAPHS DATA ", pathway_code))

        # Export the pathway data
        if (!dir.exists(file.path('./output/subGraph/'))) {
          dir.create(file.path(paste0('./output/subGraph/')), showWarnings = FALSE, mode = "0775")
        }

        if (dir.exists(file.path('./output/subGraph/'))) {
          write.csv(result, file=paste0('./output/subGraph/', kgml_index, "_", pathway_code, '.csv'))
        }

        #***********************#
        # Remove temp variables #
        #***********************#

        rm(result)
      } else {
        # Save the log file
        printLog(message_=paste0('The network ', pathway_code, ' [', kgml_index, ' of ', available_pathways,  '] doesnt have any articulation point. Skipping it...'), file_='calculateAPSubGraphs')
      }

      #***********************#
      # Remove temp variables #
      #***********************#

      rm(current_kgml, pathwayGraph, iGraph, impact)

      # Increment the index
      kgml_index <<- kgml_index + 1

    }, error=function(e) {
      printMessage(e)

      # Increment the index
      kgml_index <<- kgml_index + 1

      # Save the log file
      printLog(toString(e), file_=paste0('calculateAPSubGraphs', pathway_code))

      return(FALSE)
    })
  }) # End of Loop 01
}

#' Get all data from a folder and bind it together
#'
#' @param filename_ The name of the generated file
#' @param folderName_ The folder name that contains the necessary data
#' @param filterColumns_ Flag to determine whether or not the dataSet should be filtered into default columns
#' @param verbose_ Print every status message.
#'
#' @return This function returns a data frame containing the data from all dataSets inside a folder as csv
#'
#' @examples
#' \dontrun{
#' dataSet <- generateDataSet()
#' }
#'
#' @author
#' Igor Brandão

generateConsolidatedDataSet <- function(filename_ = '', folderName_ = 'subGraph', format_ = '.csv', filterColumns_ = TRUE, verbose_ = TRUE) {
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
  dataSet <- do.call(rbind.fill, big.list.of.data.frames)

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
      dataSet$X <- NULL
    }

    # Status message
    if (verbose_) {
      printMessage(paste0(toupper(filename_), " DATASET GENERATED WITH SUCCESS!"))
      printMessage("SAVING THE GENERATED DATASET...")
    }

    # Export the subgraph data
    if (!dir.exists(file.path(paste0('./output/', folderName_ , '/')))) {
      dir.create(file.path(paste0('./output/', folderName_ , '/')), showWarnings = FALSE, mode = "0775")
    }

    write.csv(dataSet, file = paste0('./output/', folderName_, '/', filename_, format_))

    # Return the generated dataSet
    return(dataSet)
  }
}

#*******************************************************************************************#

# ---- PIPELINE SECTION ----

#***************#
# Pipeline flow #
#***************#

#****************************************#
# Step 1: Generate the network subgraphs #
#****************************************#
#calculateAPSubGraphs()

#**************************************************#
# Step 2: Export the cnsolidated subgraph datasets #
#**************************************************#
#dataSet <- generateConsolidatedDataSet(filename_='allSubGraphs')

#*********************************************************#
# Step 3: Perform plots to explore the betweenness metric #
#*********************************************************#

# Load the dataSet
dataSet <- read.csv('./output/subGraph/allSubGraphs.csv', header=TRUE, sep=",", stringsAsFactors=FALSE)

# Filter the dataSet selecting just the groups of interest
dataSetPlot <- dataSet[dataSet$ap_group==paste0('>=', groupThresholdMax) |
                         dataSet$ap_group==paste0('<', groupThresholdMin), c('ap_group', 'betweenness')]

dataSetPlot[dataSetPlot$ap_group == paste0('>=', groupThresholdMax),]$ap_group <- paste0('\u2265', groupThresholdMax, '%')
dataSetPlot[dataSetPlot$ap_group == paste0('<', groupThresholdMin),]$ap_group  <- paste0('<', groupThresholdMin, '%')

# ---- plot1 ----

# Generate the boxplot plot
plot1 <- ggplot(dataSetPlot, aes(x=ap_group, y=betweenness, na.rm = TRUE)) +
  # Add the boxplot
  stat_boxplot(geom = "errorbar", width = 0.15) +
  #geom_boxplot(aes(group=ap_group, fill=ap_group), colour=c("#831D0C", "#03080C"), fill=c("#ED553B", "#173F5F"), position=position_dodge(0.5)) +
  geom_boxplot(aes(group=ap_group, fill=ap_group), colour="black", fill="#A4A4A4", position=position_dodge(0.5)) +

  scale_y_continuous(breaks=seq(from = 0, to = 0.5, by = 0.1)) +

  # Chart visual properties
  xlab("Articulation point frequency group") +
  ylab("Betweenness") +
  ggtitle("") +
  guides(fill=guide_legend(title="")) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", color="black", size=26, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title.x = element_text(face="bold", color="black", size=20, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(face="bold", color="black", size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", color="black", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(face="bold", color="black", size=20),
        legend.title = element_text(face="bold", size=18),
        legend.text = element_text(size=16),
        legend.position = 'right') + labs(fill = "AP group")

plot1

ggsave(paste0("./output/statistics/articulationPointCentrality/apCenterPeripheryBetweenness.jpeg"), width = 30, height = 25, units = "cm")
ggsave(paste0("./output/statistics/articulationPointCentrality/apCenterPeripheryBetweenness.svg"), width = 30, height = 25, units = "cm")

# Run the Mann-Whitney Wilcoxon test
wilcox.test(dataSetPlot[dataSetPlot$ap_group==paste0('\u2265', groupThresholdMax, '%'),]$betweenness,
            dataSetPlot[dataSetPlot$ap_group==paste0('<', groupThresholdMin, '%'),]$betweenness)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ---- plot2 ----

# Generate the violin plot
plot2 <- ggplot(dataSetPlot, aes(x=ap_group, y=betweenness, na.rm = TRUE)) +
  # Add the violin
  #geom_violin(draw_quantiles = c(0.25, 0.75), linetype = "dashed", colour="blue", fill="#A4A4A4") +
  geom_violin(draw_quantiles = 0.5, colour="red", fill="#A4A4A4", size=0.5) +
  geom_violin(aes(group=ap_group, fill=ap_group), colour="black", fill="transparent", size=0.5, position=position_dodge(0.5)) +
  stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="red", fill="red") +
  stat_boxplot(geom = "errorbar", color="blue", linetype = "dashed", size=1, width = 0.35) +

  scale_y_continuous(breaks=seq(from = 0, to = 0.5, by = 0.1)) +

  # Chart visual properties
  xlab("Articulation point frequency group") +
  ylab("Betweenness") +
  ggtitle("") +
  guides(fill=guide_legend(title="")) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", color="black", size=26, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title.x = element_text(face="bold", color="black", size=20, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(face="bold", color="black", size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", color="black", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(face="bold", color="black", size=20),
        legend.title = element_text(face="bold", size=18),
        legend.text = element_text(size=16),
        legend.position = 'right') + labs(fill = "AP group")

plot2

ggsave(paste0("./output/statistics/articulationPointCentrality/apCenterPeripheryBetweennessViolin.jpeg"), width = 30, height = 25, units = "cm")
ggsave(paste0("./output/statistics/articulationPointCentrality/apCenterPeripheryBetweennessViolin.svg"), width = 30, height = 25, units = "cm")

# Run the Mann-Whitney Wilcoxon test
wilcox.test(dataSetPlot[dataSetPlot$ap_group==paste0('\u2265', groupThresholdMax, '%'),]$betweenness,
            dataSetPlot[dataSetPlot$ap_group==paste0('<', groupThresholdMin, '%'),]$betweenness)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#*******************************************************#
# Step 4: Perform plots to explore the AP impact metric #
#*******************************************************#

# ---- plot3 ----

# Select all the columns related to the graph components size
graphComponentsSize <- dataSet[,which(grepl("subgraph[0-9]{1,2}_size", names(dataSet)))]

# Calculates each AP impact score
dataSet$impact_score <- 0

# impactScore = totalNodes - maxComponentSize
for (idx in 1:nrow(dataSet)) {
  #dataSet[idx,]$impact_score <- dataSet[idx,]$pathway_nodes - max(graphComponentsSize[idx,], na.rm = T)
  dataSet[idx,]$impact_score <- dataSet[idx,]$pathway_nodes - mean(as.numeric(graphComponentsSize[idx,]), na.rm = T)
}

# Filter the dataSet selecting just the groups of interest
dataSetImpactPlot <- dataSet[dataSet$ap_group==paste0('>=', groupThresholdMax) |
                         dataSet$ap_group==paste0('<', groupThresholdMin),
                         c('ap_dict_id', 'pathway', 'ap_group', 'pathway_nodes', 'pathway_edges', 'betweenness', 'degree', 'noSubgraphs', 'impact_score')]

# Fill the pathway code with zeros
dataSetImpactPlot <- fillPathwayCodeWithZeros(dataSetImpactPlot)

dataSetImpactPlot[dataSetImpactPlot$ap_group == paste0('>=', groupThresholdMax),]$ap_group <- paste0('\u2265', groupThresholdMax, '%')
dataSetImpactPlot[dataSetImpactPlot$ap_group == paste0('<', groupThresholdMin),]$ap_group  <- paste0('<', groupThresholdMin, '%')

# Generate the boxplot plot
plot3 <- ggplot(dataSetImpactPlot, aes(x=ap_group, y=impact_score, na.rm = TRUE)) +
  # Add the boxplot
  stat_boxplot(geom = "errorbar", width = 0.15) +
  #geom_boxplot(aes(group=ap_group, fill=ap_group), colour=c("#831D0C", "#03080C"), fill=c("#ED553B", "#173F5F"), position=position_dodge(0.5)) +
  geom_boxplot(aes(group=ap_group, fill=ap_group), colour="black", fill="#A4A4A4", position=position_dodge(0.5)) +

  #scale_y_continuous(breaks=seq(from = 0, to = 0.5, by = 0.1)) +

  # Chart visual properties
  xlab("Articulation point frequency group") +
  ylab("Articulation point impact") +
  ggtitle("") +
  guides(fill=guide_legend(title="")) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", color="black", size=26, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title.x = element_text(face="bold", color="black", size=20, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(face="bold", color="black", size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", color="black", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(face="bold", color="black", size=20),
        legend.title = element_text(face="bold", size=18),
        legend.text = element_text(size=16),
        legend.position = 'right') + labs(fill = "AP group")

plot3

# Run the Mann-Whitney Wilcoxon test
wilcox.test(dataSetImpactPlot[dataSetImpactPlot$ap_group==paste0('\u2265', groupThresholdMax, '%'),]$impact_score,
            dataSetImpactPlot[dataSetImpactPlot$ap_group==paste0('<', groupThresholdMin, '%'),]$impact_score)

ggsave(paste0("./output/statistics/articulationPointCentrality/apCenterPeripheryImpact.jpeg"), width = 30, height = 25, units = "cm")
ggsave(paste0("./output/statistics/articulationPointCentrality/apCenterPeripheryImpact.svg"), width = 30, height = 25, units = "cm")

write.csv(dataSetImpactPlot, file=paste0('./output/statistics/articulationPointCentrality/apImpact', groupThresholdMax, '.csv'), row.names = F)

# Generate a exploratory report
DataExplorer::create_report(dataSetImpactPlot)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
