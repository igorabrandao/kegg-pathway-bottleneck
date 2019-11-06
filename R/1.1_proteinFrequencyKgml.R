#**********************************************************************#
# Pipeline to perform the enzyme frequencies count via KEGG kgml files #
#**********************************************************************#

# ---- IMPORT SECTION ----

# 1.1_proteinFrequencyKgml.R #

#' This is the pipeline script to perform
#' the enzymes frequencies counting via KEGG kgml files
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
files.sources[3] = paste0("./R/functions", "/", "graphPrintFunctions.R")
files.sources[4] = paste0("./R/functions", "/", "helperFunctions.R")
sapply(files.sources, source)

# Load the pathways by organisms data
organism2pathway <- get(load(paste0("./dictionaries", "/", "organism2pathway.RData")))
pathwayList <- get(load(paste0("./dictionaries", "/", "pathwayList.RData")))

#*******************************************************************************************#

# ---- FUNCTIONS SECTION ----

#*************************#
# Pipeline main functions #
#*************************#

getPathwayEnzymeKGML <- function(removeNoise_=TRUE) {

  # Reference pathway
  reference_pathway <- 'ko'

  # Get the list of files
  folder = paste0("./output/kgml/", reference_pathway, "/")
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
    current_kgml <- KGML2Dataframe(paste0(folder, kgml_list[1])) # APAGAR

    # Get the pathway code
    pathway_code <- onlyNumber(kgml_list[1])

    # Status message
    printMessage(paste0("COUNTING ", pathway_code, " ENZYMES FREQUENCIES [", kgml_index, " OF ", available_pathways, "]"))

    tryCatch({
      # Convert the pathway data into a graph
      #pathwayGraph <- KGML2Graph(paste0(folder, file), replaceOrg=TRUE, orgToReplace=reference_pathway)
      pathwayGraph <- KGML2Graph(paste0(folder, kgml_list[1]), replaceOrg=TRUE, orgToReplace=reference_pathway) # APAGAR

      #*************************##
      # Prepare the pathway data #
      #*************************##

      # Create the pathwayData dataFrame
      pathwayData <- current_kgml$nodes[,c('entryID', 'name', 'type', 'reaction')]

      # Add the default columns
      pathwayData$org <- reference_pathway
      pathwayData$pathway <- pathway_code
      pathwayData$is_bottleneck <- 0

      pathwayData$freq <- 0
      pathwayData$freq_mean <- 0
      pathwayData$total_species <- 0
      pathwayData$percentage <- 0
      pathwayData$percentage_mean <- 0

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
        pathwayData <- removeNoise(pathwayData)
        pathwayGraph <- removeNoise(pathwayGraph)
      }

      # Get the graph properties
      graphProperties <- getGraphProperties(pathwayGraph)

      # Perform the graph bottleneck calculation
      iGraph <- igraph::graph_from_data_frame(pathwayGraph, directed = FALSE)
      graphBottleneck <- igraph::as_ids(getGraphBottleneck(iGraph, FALSE))

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

        # Assign the bottlenecks
        if (grepl(pattern, graphBottleneck)) {
          pathwayData[idx,]$is_bottleneck <- 1
        }
      }

      # Apply the node classification
      pathwayData$bottleneck_classification <- classifyBottleneck(pathwayData)$bottleneck_classification

      #****************##
      # Handle org data #
      #****************##

      # Check which species have the current pathway
      org_list <- sapply(organism2pathway, function(org) {
        pathway_code %in% org
      })

      # Retrieve just the species identification
      org_list <- names(org_list[org_list==TRUE])

      # Define the number of species that have the current pathway
      available_orgs <- length(org_list)
      pathwayData$total_species <- available_orgs
      org_index <- 1

      # Enzyme color identification
      enzyme_present_color <- '#BFFFBF'
      enzyme_missing_color <- '#FFFFFF'

      # Loop 02: Run through all organisms that have the current pathway
      lapply(org_list, function(org) {
        # Status message
        printMessage(paste0("PROCESSING ", org, " SPECIE [", org_index, " OF ", available_orgs, "]"))

        tryCatch({
          # Get the organism kgml file
          org_folder = paste0("./output/kgml/", org, "/", org, pathway_code, '.xml')

          # Load the dataframe
          org_kgml <- KGML2Dataframe(org_folder)

          # Create the pathwayData dataFrame
          orgData <- org_kgml$nodes[,c('entryID', 'name', 'type', 'reaction', 'bgcolor')]

          # Remove unnecessary data from org data
          if (removeNoise_) {
            orgData <- orgData[!grepl("^path:", orgData$name),]
            orgData <- orgData[!grepl("^map:", orgData$name),]
            orgData <- orgData[!grepl("^cpd:", orgData$name),]
            orgData <- orgData[!grepl("^gl:", orgData$name),]
            orgData <- orgData[!grepl("^ko:", orgData$name),]
          }

          #***************************************##
          # Count the enzyme frequency by organism #
          #***************************************##

          # Verify which enzyme are presented for the current organism
          for (idx in 1:nrow(orgData)) {
            # Check if the enzyme is present
            if (strcmp(orgData[idx,]$bgcolor, enzyme_present_color)) {
              # Increment the enzyme frequency
              pathwayData[pathwayData$entryID==orgData[idx,]$entryID,]$freq <<-
                pathwayData[pathwayData$entryID==orgData[idx,]$entryID,]$freq + 1
            }
          }
        }, error=function(e) {
          print(paste0('The org ', org, ' could no be processed. Skipping it...'))
          return(FALSE)
        })

        # Increment the index
        org_index <<- org_index + 1
      }) # End of Loop 02

      #********************##
      # Handle other metric #
      #********************##

      # Calculate the frequency percentual
      pathwayData$freq_mean <- mean(pathwayData$freq)

      # Calculate the frequency percentage
      pathwayData$percentage <- (pathwayData$freq / pathwayData$total_species) * 100

      # Calculate the percentual mean frequency
      pathwayData$percentage_mean <- mean(pathwayData$percentage)

      #****************************#
      # Prepare the data to export #
      #****************************#

      # Status message
      printMessage(paste0("EXPORTING PATHWAY ", pathway_code))

      # Rename the nodes column
      names(pathwayData)[names(pathwayData) == "name"] <- "enzyme"

      # Export the pathway data
      if (!dir.exists(file.path('./output/totalFrequency/'))) {
        dir.create(file.path(paste0('./output/totalFrequency/')), showWarnings = FALSE, mode = "0775")
      }

      if (dir.exists(file.path('./output/totalFrequency/'))) {
        save(pathwayData, file=paste0('./output/totalFrequency/', kgml_index, "_", pathway_code, '.RData'), compress = "xz")
      }

      # Increment the index
      kgml_index <<- kgml_index + 1

    }, error=function(e) {
      print(paste0('The pathway ', pathway_code, ' could no be processed. Skipping it...'))
      return(FALSE)
    })
  }) # End of Loop 01

  # Function finished with success
  return(TRUE)
}

#*******************************************************************************************#

# ---- PIPELINE SECTION ----

#***************#
# Pipeline flow #
#***************#

#**********************************#
# Step 1: Get all pathways enzymes #
#**********************************#

# [TEST ONLY]
#lapply(1:5, getPathwayEnzymes, replaceEmptyGraph_=FALSE)
