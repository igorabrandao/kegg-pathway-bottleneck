#*************************************************************************#
# Pipeline to perform the gene essentiality evaluation of model organisms #
#*************************************************************************#

# ---- IMPORT SECTION ----

# 3.0_essentialGenesModelOrg.R #

#' This is the pipeline script to perform
#' the gene essentiality evaluation of model organisms
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

# Load the pathways by organisms data
essentialGenesModelOrg <- get(load(paste0("./dictionaries", "/", "essentialGenesModelOrg.RData")))
rm(table)

# Set the list of model organisms according to KEGG nomenclature
# Mus musculus => mmu
# Drosophila melanogaster => dme
# Saccharomyces cerevisiae => sce
# Caenorhabditis elegans => cel
modelOrgList = c("mmu", "dme", "sce", "cel")

#*******************************************************************************************#

# ---- FUNCTIONS SECTION ----

#*************************#
# Pipeline main functions #
#*************************#

#' Parse the organism KGML files and generate a dataset with all pathways enzymes with it's
#' frequencies and graph metrics
#'
#' @param org_ KEGG organism code
#' @param removeNoise_ Remove undesirable enzyme such as: ko, map, path, cpd or gl.
#'
#' @return This function does not return nothing, just export .csv files.
#'
#' @examples
#' \dontrun{
#' generateOrgDataFromKGML()
#' generateOrgDataFromKGML(FALSE)
#' }
#'
#' @author
#' Igor Brandão

generateOrgDataFromKGML <- function(org_, removeNoise_=TRUE) {

  # Check if the folder contains files
  if (is.null(org_) | length(org_) == 0) {
    # Status message
    printMessage("The organism code wasn't informed...")
    return(FALSE)
  }

  # Get the organism pathway list of files
  folder = paste0("./output/kgml/", org_, "/")
  kgml_list <- list.files(path=folder, pattern='*.xml')
  kgml_index <- 1

  # Define the number of available pathways
  available_pathways <- length(kgml_list)

  # Check if the folder contains files
  if (is.null(kgml_list) | length(kgml_list) == 0) {
    # Status message
    printMessage("There aren't available pathways...")
    return(FALSE)
  }

  # Loop 01: Run through all available pathways kgml
  lapply(kgml_list, function(file) {

    # Load the dataframe
    current_kgml <- KGML2Dataframe(paste0(folder, file))

    # Get the pathway code
    pathway_code <- onlyNumber(file)

    # Status message
    printMessage(paste0("GATHERING ", org_, " PATHWAY ", pathway_code, " DATA [", kgml_index, " OF ", available_pathways, "]"))

    tryCatch({
      # Convert the pathway data into a graph
      pathwayGraph <- KGML2Graph(paste0(folder, file), replaceOrg=TRUE, orgToReplace=org_)

      #*************************##
      # Prepare the pathway data #
      #*************************##

      # Create the pathwayData dataFrame
      pathwayData <- current_kgml$nodes[,c('entryID', 'name', 'type', 'reaction', 'bgcolor')]

      # Add the default columns
      pathwayData$reaction_type <- NA
      pathwayData$org <- org_
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

      # Remove the duplicated nodes
      pathwayData <- pathwayData[!duplicated(c("name")),]

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
        if (strcmp(org_, 'ko')) {
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
      enzyme_present_color <- '#BFFFBF'
      enzyme_missing_color <- '#FFFFFF'

      # Increment the enzyme frequency
      pathwayData[pathwayData$bgcolor == enzyme_present_color,]$freq <-
        pathwayData[pathwayData$bgcolor == enzyme_present_color,]$freq + 1

      #****************************#
      # Prepare the data to export #
      #****************************#

      # Status message
      printMessage(paste0("EXPORTING ", org_, " PATHWAY ", pathway_code, " DATA [", kgml_index, " OF ", available_pathways, "]"))

      # Export the pathway data
      if (!dir.exists(file.path('./output/essentialGenes/'))) {
        dir.create(file.path(paste0('./output/essentialGenes/')), showWarnings = FALSE, mode = "0775")
      }

      if (!dir.exists(file.path(paste0('./output/essentialGenes/', org_)))) {
        dir.create(file.path(paste0('./output/essentialGenes/', org_)), showWarnings = FALSE, mode = "0775")
      }

      if (dir.exists(paste0('./output/essentialGenes/', org_))) {
        save(pathwayData, file=paste0('./output/essentialGenes/', org_, '/', kgml_index, "_", pathway_code, '.RData'), compress = "xz")
      }
    }, error=function(e) {
      printMessage(e)

      # Save the log file
      printLog(toString(e), file_=paste0('generatePathwayDataFromKGML', pathway_code))

      return(FALSE)
    })

    # Increment the index
    kgml_index <<- kgml_index + 1
  }) # End of Loop 01

  # Function finished with success
  return(TRUE)
}

#*******************************************************************************************#

# ---- PIPELINE SECTION ----

#***************#
# Pipeline flow #
#***************#

#****************************************#
# Step 1: Generate the network subgraphs #
#****************************************#

