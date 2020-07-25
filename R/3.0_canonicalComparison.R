#*************************************************************************************#
# Pipeline to perform the comparison between canonical network and organisms network #
#*************************************************************************************#

# ---- IMPORT SECTION ----

# 3.0_canonicalComparison.R #

#' This is the pipeline script to perform
#' the comparison between canonical network and organisms network
#'
#' @author
#' Igor Brandão

# Import the necessary libraries

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
sapply(files.sources, source)

# Set the list of model organisms according to KEGG nomenclature
# Reference => ec
# Homo sapiens => hsa
# Mus musculus => mmu
# Drosophila melanogaster => dme
# Saccharomyces cerevisiae => sce
# Caenorhabditis elegans => cel
modelOrgList = c("ec", "hsa", "mmu", "dme", "sce", "cel")
pathwayList = c("00010")

#*******************************************************************************************#

# ---- FUNCTIONS SECTION ----

#*************************#
# Pipeline main functions #
#*************************#

#' Parse the organism KGML files and generate a dataset with all pathways enzymes with it's
#' frequencies and graph metrics
#'
#' @param pathwayList_ KEGG organism code
#' @param orgList_ KEGG organism code
#' @param removeNoise_ Remove undesirable enzyme such as: ko, map, path, cpd or gl.
#'
#' @return This function does not return nothing, just export .csv files.
#'
#' @examples
#' \dontrun{
#' generateOrgNetworks()
#' generateOrgNetworks(FALSE)
#' }
#'
#' @author
#' Igor Brandão

generateOrgNetworks <- function(pathwayList_, orgList_, removeNoise_=TRUE) {

  # Load the dicionaty
  dictionary <- read.csv(file='./output/pathwaysDictionary/dictionary.csv', header=TRUE, sep=",", stringsAsFactors=FALSE)

  # Loop 01: Run through all organism in list
  lapply(orgList_, function(currentOrg) {

    # Get the organism pathway list of files
    folder = paste0("./output/kgml/", currentOrg, "/")
    kgml_list <- list.files(path=folder, pattern='*.xml')
    kgml_index <- 1

    # Filter the files according to the pathway list
    kgml_list <- sapply(1:length(kgml_list), function(x) kgml_list[grepl(pathwayList_[x], kgml_list)])
    kgml_list <- unlist(kgml_list[!is.na(unlist(kgml_list, use.names=FALSE))])

    # Loop 02: Run through all KGMLs in list
    lapply(kgml_list, function(file) {

      # Define the number of available pathways
      available_pathways <- length(kgml_list)

      # Check if the folder contains files
      if (is.null(kgml_list) | length(kgml_list) == 0) {
        # Status message
        printMessage("There aren't available pathways...")
        return(FALSE)
      }

      # Load the dataframe
      current_kgml <- KGML2Dataframe(paste0(folder, file))

      # Get the pathway code
      pathway_code <- onlyNumber(file)

      # Status message
      printMessage(paste0("GATHERING ", currentOrg, " PATHWAY ", pathway_code, " DATA [", kgml_index, " OF ", available_pathways, "]"))

      tryCatch({
        # Convert the pathway data into a graph
        pathwayGraph <- KGML2Graph(paste0(folder, file), replaceOrg=TRUE, orgToReplace=currentOrg)

        #*************************##
        # Prepare the pathway data #
        #*************************##

        # Create the pathwayData dataFrame
        pathwayData <- current_kgml$nodes

        # Remove the pathwayData unnecessary columns
        pathwayData <- pathwayData[,!(names(pathwayData) %in% c('link', 'component', 'map', 'fgcolor'))]

        # Add the default columns
        pathwayData$dictID <- NA
        pathwayData$reaction_type <- NA
        pathwayData$org <- currentOrg
        pathwayData$pathway <- pathway_code
        pathwayData$is_bottleneck <- 0

        pathwayData$freq <- 0

        pathwayData$betweenness <- NA
        pathwayData$connectivity <- NA
        pathwayData$triangles <- NA
        pathwayData$clusteringCoef <- NA
        pathwayData$closenessCoef <- NA
        pathwayData$community <- NA
        pathwayData$eigenvectorScore <- NA
        pathwayData$eccentricity <- NA
        pathwayData$radius <- NA
        pathwayData$diameter <- NA
        pathwayData$degree <- NA
        pathwayData$authorityScore <- NA
        pathwayData$hubScore <- NA

        pathwayData$bottleneck_classification <- NA

        # Remove unnecessary data from pathway data/graph
        if (removeNoise_) {
          pathwayGraph <- removeNoise(pathwayGraph)

          pathwayData <- pathwayData[!grepl("^path:", pathwayData$name),]
          pathwayData <- pathwayData[!grepl("^map:", pathwayData$name),]
          pathwayData <- pathwayData[!grepl("^cpd:", pathwayData$name),]
          pathwayData <- pathwayData[!grepl("^gl:", pathwayData$name),]
          pathwayData <- pathwayData[!grepl("^ko:", pathwayData$name),]
        }

        # Assign the reaction type to each node
        for (idx in 1:nrow(pathwayData)) {
          for (idx2 in 1:length(current_kgml$reactions$name)) {
            # Check the position of the current reaction in reactions list
            if (current_kgml$reactions$name[idx2] %in% pathwayData[idx,]$reaction) {
              pathwayData[idx,]$reaction_type <- current_kgml$reactions$type[idx2]
              break()
            }
          }
        }

        # Get the graph properties
        graphProperties <- getGraphProperties(pathwayGraph)

        # Perform the graph bottleneck calculation
        iGraph <- igraph::graph_from_data_frame(pathwayGraph, directed = FALSE)
        graphBottleneck <- igraph::as_ids(getGraphBottleneck(iGraph, FALSE))

        # Assign the bottlenecks
        pathwayData$is_bottleneck[which(pathwayData$name %in% graphBottleneck)] <- 1

        # Assign the graph properties to each node
        for (idx in 1:nrow(pathwayData)) {
          # Prepare the node name to be compared
          pattern <- gsub('/', '|', pathwayData$name[idx])
          pattern <- gsub(" ","", pattern)

          # Find which rows in graphProperties match with pathwayData
          rowsToMerge <- which(grepl(pattern, graphProperties$node))[1]

          pathwayData[idx,]$betweenness <- graphProperties[rowsToMerge,]$betweenness
          pathwayData[idx,]$connectivity <- graphProperties[rowsToMerge,]$connectivity
          pathwayData[idx,]$triangles <- graphProperties[rowsToMerge,]$triangles
          pathwayData[idx,]$clusteringCoef <- graphProperties[rowsToMerge,]$clusteringCoef
          pathwayData[idx,]$closenessCoef <- graphProperties[rowsToMerge,]$closenessCoef
          pathwayData[idx,]$community <- graphProperties[rowsToMerge,]$community
          pathwayData[idx,]$eigenvectorScore <- graphProperties[rowsToMerge,]$eigenvectorScore
          pathwayData[idx,]$eccentricity <- graphProperties[rowsToMerge,]$eccentricity
          pathwayData[idx,]$radius <- graphProperties[rowsToMerge,]$radius
          pathwayData[idx,]$diameter <- graphProperties[rowsToMerge,]$diameter
          pathwayData[idx,]$degree <- graphProperties[rowsToMerge,]$degree
          pathwayData[idx,]$authorityScore <- graphProperties[rowsToMerge,]$authorityScore
          pathwayData[idx,]$hubScore <- graphProperties[rowsToMerge,]$hubScore

          # Assign the bottlenecks for kegg orthology (ko)
          if (strcmp(currentOrg, 'ko')) {
            # Check if at least one bottleneck was found
            if (length(graphBottleneck) > 0) {
              if (grepl(pattern, graphBottleneck)) {
                pathwayData[idx,]$is_bottleneck <- 1
              }
            }
          }
        }

        # Apply the node classification
        pathwayData$bottleneck_classification <- classifyBottleneck(pathwayData)$bottleneck_classification

        #***************************##
        # Count the enzyme frequency #
        #***************************##

        # Enzyme color identification
        enzyme_present_color_org <- '#BFFFBF'
        enzyme_present_color_ref <- '#BFBFFF'
        enzyme_missing_color <- '#FFFFFF'

        # Increment the enzyme frequency
        pathwayData[pathwayData$bgcolor == enzyme_present_color_org || pathwayData$bgcolor == enzyme_present_color_ref,]$freq <-
          pathwayData[pathwayData$bgcolor == enzyme_present_color_org || pathwayData$bgcolor == enzyme_present_color_ref,]$freq + 1

        #*****************************##
        # Get the enzyme dictionary ID #
        #*****************************##

        # Assign the dictionary ID to each node
        rows <- nrow(pathwayData)
        for (nodeIdx in 1:rows) {
          # Setup the current node
          currentNode <- pathwayData[nodeIdx,]

          # Find the current node into the dictionary
          if (is.na(currentNode$reaction)) {
            dictId1 <- dictionary[(dictionary$x == currentNode$x) & (dictionary$y == currentNode$y) &
                                    (dictionary$ec == currentNode$name), ]$id
          } else {
            dictId1 <- dictionary[(dictionary$x == currentNode$x) & (dictionary$y == currentNode$y) &
                                    (dictionary$reaction == currentNode$reaction) & (dictionary$ec == currentNode$name), ]$id
          }

          # Remove NAs
          dictId1 <- dictId1[!is.na(dictId1)]

          # Check if the dictionary contains the node, if the didctionary ID is empty, it means that
          # the node refers to a pathway connection (e.g:path:00020) or it is a compound
          if (is.null(dictId1) | isempty(dictId1)) {
            next()
          } else {
            pathwayData[nodeIdx, 'dictID'] <- dictId1
          }
        }

        #***********************************************************#
        # Export the organism specific data for the current pathway #
        #***********************************************************#

        # Status message
        printMessage(paste0("EXPORTING ", currentOrg, " PATHWAY ", pathway_code, " DATA [", kgml_index, " OF ", available_pathways, "]"))

        # Create the canonical analysis folder
        if (!dir.exists(file.path('./output/canonicalNetworkComparison/'))) {
          dir.create(file.path(paste0('./output/canonicalNetworkComparison/')), showWarnings = FALSE, mode = "0775")
        }

        # Create the org folder
        if (!dir.exists(file.path(paste0('./output/canonicalNetworkComparison/', currentOrg)))) {
          dir.create(file.path(paste0('./output/canonicalNetworkComparison/', currentOrg)), showWarnings = FALSE, mode = "0775")
        }

        # Save the org pathwayData into its folder
        if (dir.exists(paste0('./output/canonicalNetworkComparison/', currentOrg))) {
          write.csv(pathwayData, file=paste0('./output/canonicalNetworkComparison/', currentOrg, '/', pathway_code, '.csv'), row.names=FALSE)
        }

        #**********************************#
        # Generate the interactive network #
        #**********************************#
        graphDictionary <- KGML2GraphDictionary(paste0(folder, file), replaceOrg=TRUE, orgToReplace=currentOrg)
        graphDictionary <- removeNoise(graphDictionary)

        printInteractiveNetwork(pathway_code, currentOrg, pathwayData, graphDictionary)

      }, error=function(e) {
        printMessage(e)

        # Save the log file
        printLog(toString(e), file_=paste0('generatePathwayDataFromKGML', pathway_code))

        return(FALSE)
      })
    }) # End of Loop 02

    # Increment the index
    kgml_index <<- kgml_index + 1
  }) # End of Loop 01

  # Function finished with success
  return(TRUE)
}

#' Function to generate interactive networks
#'
#' @param removeNoise_ Remove undesirable enzyme such as: map, path, cpd or gl.
#'
#' @return This function does not return nothing, just export files.
#'
#' @examples
#' \dontrun{
#' printInteractiveNetwork()
#' printInteractiveNetwork(FALSE)
#' }
#'
#' @author
#' Igor Brandão

printInteractiveNetwork <- function(pathwayCode_, org_, pathwayData_, pathwayGraph_) {

  # Assign the parameters
  pathwayCode <- pathwayCode_
  currentOrg <- org_
  pathwayData <- pathwayData_
  pathwayGraph <- pathwayGraph_

  #--------------------------#
  # [GENERATING THE NETWORK] #
  #--------------------------#

  # Rename the AP classification
  if (nrow(pathwayData[pathwayData$bottleneck_classification=='HB',]) > 0) {
    pathwayData[pathwayData$bottleneck_classification=='HB',]$bottleneck_classification <- 'HAP'
  }

  if (nrow(pathwayData[pathwayData$bottleneck_classification=='HNB',]) > 0) {
    pathwayData[pathwayData$bottleneck_classification=='HNB',]$bottleneck_classification <- 'HUB'
  }

  if (nrow(pathwayData[pathwayData$bottleneck_classification=='NHB',]) > 0) {
    pathwayData[pathwayData$bottleneck_classification=='NHB',]$bottleneck_classification <- 'AP'
  }

  if (nrow(pathwayData[pathwayData$bottleneck_classification=='NHNB',]) > 0) {
    pathwayData[pathwayData$bottleneck_classification=='NHNB',]$bottleneck_classification <- 'Others'
  }

  # Rename the bottleneck column
  names(pathwayData)[names(pathwayData) == "bottleneck_classification"] <- "AP_classification"

  # Adjust the column percentage
  pathwayData$percentage <- 100

  # Generate the network
  generatedNetwork <- generateInteractiveNetwork(network_=pathwayGraph, networkProperties_=pathwayData, pathway_=pathwayCode, org_=currentOrg,
                                                 pathwayDetail_=NULL, dynamicNetwork_=FALSE)

  #---------------------------------#
  # [EXPORT THE NETWORK PROPERTIES] #
  #---------------------------------#

  # Generate an igraph object to extract topological features
  networkGlobalDataGraph <- igraph::graph_from_data_frame(pathwayGraph, directed = FALSE)

  networkGlobalTemp <- data.frame(pathway="", nodes="", edges="", total_species="", node_highest_impact="", disconnected_nodes="",
                                  community="", mean_degree="", mean_betweenness="", ap_number="", hap_number="", hub_number="", others_number="",
                                  stringsAsFactors=FALSE)

  # Annotate the network caracteristics
  networkGlobalTemp$pathway <- pathwayCode
  networkGlobalTemp$nodes <- length(V(networkGlobalDataGraph))
  networkGlobalTemp$edges <- length(E(networkGlobalDataGraph))
  networkGlobalTemp$total_species <- mean(pathwayData$totalSpecies)
  networkGlobalTemp$node_highest_impact <- pathwayData[which.max(pathwayData$bottleneckImpact),]$name
  networkGlobalTemp$disconnected_nodes <- max(pathwayData$bottleneckImpact)
  networkGlobalTemp$community <- max(pathwayData$community)
  networkGlobalTemp$mean_degree <- mean(pathwayData$degree)
  networkGlobalTemp$mean_betweenness <- mean(pathwayData$betweenness)
  networkGlobalTemp$ap_number <- nrow(pathwayData[pathwayData$AP_classification=="AP",])
  networkGlobalTemp$hap_number <- nrow(pathwayData[pathwayData$AP_classification=="HAP",])
  networkGlobalTemp$hub_number <- nrow(pathwayData[pathwayData$AP_classification=="HUB",])
  networkGlobalTemp$others_number <- nrow(pathwayData[pathwayData$AP_classification=="Others",])

  networkGlobalData <- rbind(networkGlobalData, networkGlobalTemp)

  #---------------------------------#
  # [EXPORT THE NETWORK PROPERTIES] #
  #---------------------------------#

  if (dir.exists(paste0('./output/canonicalNetworkComparison/', currentOrg))) {
    write.csv(networkGlobalData, file=paste0('./output/canonicalNetworkComparison/', currentOrg, '/', pathwayCode, '_properties.csv'), row.names=FALSE)
  }

  if (dir.exists(file.path(paste0('./output/canonicalNetworkComparison/', currentOrg)))) {
    filename <- paste0(pathwayCode, '.html')

    # Save the HTML file
    visSave(generatedNetwork, file = filename, selfcontained = TRUE)

    if (file.exists(filename)) {
      # Copy the file into correct directory
      file.copy(filename, paste0('./output/canonicalNetworkComparison/', currentOrg, '/', filename), overwrite = TRUE)

      # Remove the generated file
      file.remove(filename)
    } else {
      printMessage(paste0("Network file not found. Skipping it..."))
      return(FALSE)
    }
  }

  # Function finished with success
  return(TRUE)
}

#*******************************************************************************************#

# ---- PIPELINE SECTION ----

#***************#
# Pipeline flow #
#***************#
generateOrgNetworks(pathwayList, modelOrgList, removeNoise_=TRUE)
