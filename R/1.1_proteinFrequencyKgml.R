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
#files.sources[3] = paste0("./R/functions", "/", "graphPrintFunctions.R")
files.sources[3] = paste0("./R/functions", "/", "helperFunctions.R")
sapply(files.sources, source)

# Load the pathways by organisms data
organism2pathway <- get(load(paste0("./dictionaries", "/", "organism2pathway.RData")))
pathwayList <- get(load(paste0("./dictionaries", "/", "pathwayList.RData")))
pathwayDetail <- get(load(paste0("./dictionaries", "/", "pathwayDetail.RData")))

#*******************************************************************************************#

# ---- FUNCTIONS SECTION ----

#*************************#
# Pipeline main functions #
#*************************#

#' Parse the KGML file and generate the pathway dataset with enzyme frequencies and graph metrics
#'
#' @param removeNoise_ Remove undesirable enzyme such as: ko, map, path, cpd or gl.
#'
#' @return This function does not return nothing, just export .RData files.
#'
#' @examples
#' \dontrun{
#' generatePathwayDataFromKGML()
#' generatePathwayDataFromKGML(FALSE)
#' }
#'
#' @author
#' Igor Brandão

generatePathwayDataFromKGML <- function(removeNoise_=TRUE) {

  # Reference pathway
  reference_pathway <- 'ec'

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

    # Get the pathway code
    pathway_code <- onlyNumber(file)

    # Status message
    printMessage(paste0("COUNTING ", pathway_code, " ENZYMES FREQUENCIES [", kgml_index, " OF ", available_pathways, "]"))

    tryCatch({
      # Convert the pathway data into a graph
      pathwayGraph <- KGML2Graph(paste0(folder, file), replaceOrg=TRUE, orgToReplace=reference_pathway)

      #*************************##
      # Prepare the pathway data #
      #*************************##

      # Create the pathwayData dataFrame
      pathwayData <- current_kgml$nodes[,c('entryID', 'name', 'type', 'reaction')]

      # Add the default columns
      pathwayData$reaction_type <- NA
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

      # Assign the bottlenecks for enzyme code (ec)
      if (strcmp(reference_pathway, 'ec')) {
        pathwayData$is_bottleneck[which(pathwayData$name %in% graphBottleneck)] <- 1
      }

      # Remove the duplicated nodes
      pathwayData <- pathwayData[!duplicated(pathwayData$name),]

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
        if (strcmp(reference_pathway, 'ko')) {
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
        printMessage(paste0("PROCESSING ", org, " SPECIE FOR PATHWAY", pathway_code," [", org_index, " OF ", available_orgs, "]"))

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
          printMessage(e)

          # Save the log file
          printLog(toString(e), file_=paste0('generatePathwayDataFromKGML', org))

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

      # Rename the rows name
      row.names(pathwayData) <- pathwayData$name

      # Rename the nodes column
      names(pathwayData)[names(pathwayData) == "name"] <- "node"

      # Export the pathway data
      if (!dir.exists(file.path('./output/totalFrequency/'))) {
        dir.create(file.path(paste0('./output/totalFrequency/')), showWarnings = FALSE, mode = "0775")
      }

      if (dir.exists(file.path('./output/totalFrequency/'))) {
        save(pathwayData, file=paste0('./output/totalFrequency/', kgml_index, "_", pathway_code, '.RData'), compress = "xz")
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

generatePathwayGraphFromKGML <- function(removeNoise_=TRUE) {

  # Status message
  printMessage(paste0("GENERATING THE PATHWAYS GRAPHS"))

  # Reference pathway
  reference_pathway <- 'ec'

  # Get the list of files
  folder = paste0("./output/kgml/", reference_pathway, "/")
  kgml_list <- list.files(path=folder, pattern='*.xml')
  kgml_index <- 1

  # Define the number of available pathways
  available_pathways <- length(kgml_list)

  # Loop 01: Run through all available pathways kgml
  lapply(kgml_list, function(file) {

    # Get the pathway code
    pathway_code <- onlyNumber(file)

    # Status message
    printMessage(paste0("GENERATING GRAPH OF PATHWAY ", pathway_code, " [", kgml_index, " OF ", available_pathways, "]"))

    #****************************##
    # Generate the pathway graph #
    #****************************##

    tryCatch({
      # Convert the pathway data into a graph
      pathwayGraph <- KGML2GraphDictionary(paste0(folder, file), replaceOrg=TRUE, orgToReplace=reference_pathway)

      #**************************##
      # Prepare the pathway graph #
      #**************************##

      # Remove unnecessary data from pathway data/graph
      if (removeNoise_) {
        pathwayGraph <- removeNoise(pathwayGraph)
      }

      if (is.null(pathwayGraph) | isempty(pathwayGraph)) {
        # Save the log file
        printLog(message_='The pathwayGraph data frame is empty. Skipping it...', file_='generatePathwayGraphFromKGML')

        # Increment the index
        kgml_index <<- kgml_index + 1

        return(FALSE)
      }

      #****************************#
      # Prepare the data to export #
      #****************************#

      # Status message
      printMessage(paste0("EXPORTING PATHWAY ", pathway_code))

      # Export the pathway data
      if (!dir.exists(file.path('./output/allGraphs/'))) {
        dir.create(file.path(paste0('./output/allGraphs/')), showWarnings = FALSE, mode = "0775")
      }

      if (dir.exists(file.path('./output/allGraphs/'))) {
        write.csv(pathwayGraph, file=paste0('./output/allGraphs/', kgml_index, "_", pathway_code, '.csv'))
      }

      #***********************#
      # Remove temp variables #
      #***********************#

      rm(pathwayGraph)

      # Increment the index
      kgml_index <<- kgml_index + 1

    }, error=function(e) {
      printMessage(e)

      # Increment the index
      kgml_index <<- kgml_index + 1

      # Save the log file
      printLog(toString(e), file_=paste0('generatePathwayGraphFromKGML', pathway_code))

      return(FALSE)
    })
  }) # End of Loop 01

  return(TRUE)
}

#' Parse the KGML file and generate every pathway dataset of each
#' organism with enzyme frequency
#'
#' @param removeNoise_ Remove undesirable enzyme such as: ko, map, path, cpd or gl.
#'
#' @return This function does not return nothing, just export .csv files.
#'
#' @examples
#' \dontrun{
#' generateOrganismPathwayDataFromKGML()
#' generateOrganismPathwayDataFromKGML(FALSE)
#' }
#'
#' @author
#' Igor Brandão

generateOrganismPathwayDataFromKGML <- function(removeNoise_=TRUE) {

  # Reference pathway
  reference_pathway <- 'ec'

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

  # Enzyme color identification
  enzyme_present_color <- '#BFFFBF'
  enzyme_missing_color <- '#FFFFFF'

  # Loop 01: Run through all available pathways kgml
  lapply(kgml_list, function(file) {

    # Get the pathway code
    pathway_code <- onlyNumber(file)

    # Status message
    printMessage(paste0("PROCESSING ", pathway_code, " [", kgml_index, " OF ", available_pathways, "]"))

    # Set the pathway dir
    current_dir <- paste0('./output/', pathway_code)

    if (!dir.exists(file.path(current_dir))) {
      dir.create(file.path(current_dir), showWarnings = FALSE, mode = "0775")
    }

    tryCatch({

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
      org_index <- 1

      # Loop 02: Run through all organisms that have the current pathway
      lapply(org_list, function(org) {

        # Status message
        printMessage(paste0("PROCESSING ", org, " SPECIE FOR PATHWAY ", pathway_code," [", org_index, " OF ", available_orgs, "]"))

        tryCatch({
          # Get the organism kgml file
          org_folder = paste0("./output/kgml/", org, "/", org, pathway_code, '.xml')

          # Load the dataframe
          org_kgml <- KGML2Dataframe(org_folder)

          # Create the pathwayData dataFrame
          orgData <- org_kgml$nodes
          orgData$freq <- 0
          orgData$pathway <- pathway_code

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

          # Perform the enzyme frequency counting
          enzyme_present <- orgData[orgData$bgcolor==enzyme_present_color,]
          enzyme_missing <- orgData[orgData$bgcolor==enzyme_missing_color,]

          if (nrow(enzyme_present) > 0) {
            orgData[orgData$bgcolor==enzyme_present_color,]$freq <- 1
          }

          if (nrow(enzyme_missing) > 0) {
            orgData[orgData$bgcolor==enzyme_missing_color,]$freq <- 0
          }

          rm(enzyme_present, enzyme_missing)

          #**********************************#
          # Export the organism pathway data #
          #**********************************#

          # Export the pathway data
          if (dir.exists(file.path(current_dir))) {
            write.csv(orgData, file=paste0(current_dir, '/', org_index, "_", org, '.csv'))
          }
        }, error=function(e) {
          printMessage(e)

          # Save the log file
          printLog(toString(e), file_=paste0('generateOrganismPathwayDataFromKGML', org, pathway_code))

          return(FALSE)
        })

        # Increment the index
        org_index <<- org_index + 1

        #***********************#
        # Remove temp variables #
        #***********************#

        rm(org_folder, org_kgml, orgData)

      }) # End of Loop 02
    }, error=function(e) {
      printMessage(e)

      # Save the log file
      printLog(toString(e), file_=paste0('generateOrganismPathwayDataFromKGML', org, pathway_code))

      return(FALSE)
    })

    # Increment the index
    kgml_index <<- kgml_index + 1

    #***********************#
    # Remove temp variables #
    #***********************#

    rm(pathway_code, current_dir, org_list, available_orgs)

  }) # End of Loop 01

  # Function finished with success
  return(TRUE)
}

generatePathwayAllNodes <- function(removeNoise_=TRUE) {

  # Status message
  printMessage(paste0("GENERATING ALL PATHWAYS NODES LIST"))

  # Reference pathway
  reference_pathway <- 'ec'

  # Get the list of files
  folder = paste0("./output/kgml/", reference_pathway, "/")
  kgml_list <- list.files(path=folder, pattern='*.xml')
  kgml_index <- 1

  # Define the number of available pathways
  available_pathways <- length(kgml_list)

  # Set the dictionary dir
  current_dir <- './output/allNodes'

  if (!dir.exists(file.path(current_dir))) {
    dir.create(file.path(current_dir), showWarnings = FALSE, mode = "0775")
  }

  # Loop 01: Run through all available pathways kgml
  lapply(kgml_list, function(file) {

    # Get the pathway code
    pathway_code <- onlyNumber(file)

    # Status message
    printMessage(paste0("PROCESSING ", pathway_code, " [", kgml_index, " OF ", available_pathways, "]"))

    # Get the list of files
    folder = paste0("./output/", pathway_code, "/")
    file_list <- grep(list.files(path=folder), pattern='*.csv', value=T)

    # Load all csv files at once
    big.list.of.data.frames <- lapply(file_list, function(item) {
      read.csv(file=paste0(folder, item), header=TRUE, sep=",", stringsAsFactors=FALSE)
    })

    # Combine multiple data frames in one
    allNodes <- do.call(rbind, big.list.of.data.frames)

    # Export allNodes
    if (dir.exists(file.path(current_dir))) {
      write.csv(allNodes, file=paste0(current_dir, '/', pathway_code, '.csv'))
    }

    # Remove temporaly variables
    rm(big.list.of.data.frames, allNodes, pathway_code, folder, file_list)

    # Increment the index
    kgml_index <<- kgml_index + 1
  })
}

generateNodesDictionary <- function() {

  # Status message
  printMessage(paste0("GENERATING THE PATHWAYS DICTIONARY"))

  # Reference pathway
  reference_pathway <- 'ec'

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

  # Prepare the dictionary
  dictionary <- data.frame(id=NULL, pathway=NULL, x=NULL, y=NULL, reaction=NULL, ec=NULL, stringsAsFactors=FALSE)
  dictionary_index <- 1

  # Set the dictionary dir
  current_dir <- './output/pathwaysDictionary'

  if (!dir.exists(file.path(current_dir))) {
    dir.create(file.path(current_dir), showWarnings = FALSE, mode = "0775")
  }

  # Loop 01: Run through all available pathways kgml
  lapply(kgml_list, function(file) {

    # Get the pathway code
    pathway_code <- onlyNumber(file)

    # Status message
    printMessage(paste0("PROCESSING ", pathway_code, " [", kgml_index, " OF ", available_pathways, "]"))

    tryCatch({

      # Load the dataframe
      current_kgml <- KGML2Dataframe(paste0(folder, file))

      #********************************##
      # Load the pathway all nodes list #
      #********************************##

      # Get the list of files
      dataSet <- read.csv(file=paste0('./output/allNodes/', pathway_code, '.csv'), header=TRUE, sep=",", stringsAsFactors=FALSE)

      if (is.null(dataSet) | nrow(dataSet) == 0) {
        return(FALSE)
      } else {
        # Remove unnecessary data from pathway data
        dataSet <- dataSet[!grepl("^path:", dataSet$name),]
        dataSet <- dataSet[!grepl("^map:", dataSet$name),]
        dataSet <- dataSet[!grepl("^cpd:", dataSet$name),]
        dataSet <- dataSet[!grepl("^gl:", dataSet$name),]
        dataSet <- dataSet[!grepl("^ko:", dataSet$name),]
        dataSet <- dataSet[!grepl("^undefined:", dataSet$name),]
        dataSet <- dataSet[!grepl("^group:", dataSet$type),]
      }

      #************************##
      # Handle the pathway data #
      #************************##

      # Create the reference pathway dataFrame
      pathwayData <- current_kgml$nodes

      # Remove unnecessary data from pathway data
      pathwayData <- pathwayData[!grepl("^path:", pathwayData$name),]
      pathwayData <- pathwayData[!grepl("^map:", pathwayData$name),]
      pathwayData <- pathwayData[!grepl("^cpd:", pathwayData$name),]
      pathwayData <- pathwayData[!grepl("^gl:", pathwayData$name),]
      pathwayData <- pathwayData[!grepl("^ko:", pathwayData$name),]
      pathwayData <- pathwayData[!grepl("^undefined:", pathwayData$name),]

      #*********************************************##
      # Generate the dictionary from current pathway #
      #*********************************************##

      # Count the reference pathway rows
      rows <- nrow(pathwayData)

      # Loop over the reference pathway rows
      for (idx in 1:rows) {
        # Get the graphical <x,y> enzyme attributes
        x <- pathwayData[idx,]$x
        y <- pathwayData[idx,]$y

        if (!is.na(x) && !is.na(y)) {
          # Get the unique reaction for the specific <x,y> enzyme
          react <- unlist(unique(dataSet[dataSet$x == x & dataSet$y == y, ]$reaction))[1]

          # Apply the information into the dictionary
          dictionary[dictionary_index, 'id'] <<- dictionary_index
          dictionary[dictionary_index, 'pathway'] <<- pathway_code
          dictionary[dictionary_index, 'x'] <<- x
          dictionary[dictionary_index, 'y'] <<- y

          # Apply the reaction info into the dictionary
          if (is.null(react) | length(react) == 0 | is.na(react)) {
            # Apply the reference pathway reaction
            dictionary[dictionary_index, 'reaction'] <<- pathwayData[idx,]$reaction
          } else {
            # Use the organism pathway reaction
            dictionary[dictionary_index, 'reaction'] <<- react
          }

          # Apply the EC info into the dictionary
          dictionary[dictionary_index, 'ec'] <<- pathwayData[idx,]$name

          # Increment the dictionary index
          dictionary_index <<- dictionary_index + 1
        }
      }

      #***********************#
      # Remove temp variables #
      #***********************#

      rm(pathway_code, pathwayData, current_kgml)

    }, error=function(e) {
      printMessage(e)

      # Save the log file
      printLog(toString(e), file_=paste0('generateNodesDictionary', pathway_code))

      return(FALSE)
    })

    # Increment the index
    kgml_index <<- kgml_index + 1
  }) # End of Loop 01

  # Export the dictionary
  if (dir.exists(file.path(current_dir))) {
    write.csv(dictionary, file=paste0(current_dir, '/dictionary.csv'))
  }

  # Function finished with success
  return(TRUE)
}

generatePathwayFrequencyFromOrganismData <- function(removeNoise_=TRUE) {

  # Status message
  printMessage(paste0("PERFORMING THE PATHWAYS PROTEINS FREQUENCIES COUNT"))

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

  if (is.null(dictionary) | nrow(dictionary) == 0) {
    # Save the log file
    printLog(message_='The pathways nodes dictionary could not be found. Skipping it...',
             file_='generatePathwayFrequencyFromOrganismData')

    return(FALSE)
  }

  # Loop 01: Run through all available pathways kgml
  lapply(kgml_list, function(file) {

    # Get the pathway code
    pathway_code <- onlyNumber(file)

    # Status message
    printMessage(paste0("COUNTING ", pathway_code, " ENZYMES FREQUENCIES [", kgml_index, " OF ", available_pathways, "]"))

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

      #*************************##
      # Prepare the pathway data #
      #*************************##

      # Create the pathwayData dataFrame
      pathwayData <- current_kgml$nodes

      # Remove the pathwayData unnecessary columns
      pathwayData <- pathwayData[,!(names(pathwayData) %in% c('component', 'map'))]

      # Add the default columns
      pathwayData$dictID <- NA
      pathwayData$reaction_type <- NA
      pathwayData$org <- reference_pathway
      pathwayData$pathway <- pathway_code
      pathwayData$is_bottleneck <- 0
      pathwayData$bottleneckDisconnectedComponents <- 0
      pathwayData$bottleneckImpact <- 0
      pathwayData$bottleneckNormalizedImpact <- 0

      # Protein frequency columns
      pathwayData$occurrences <- 0
      pathwayData$totalSpecies <- 0
      pathwayData$percentage <- 0

      # Network metrics
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
        pathwayData <- pathwayData[!(pathwayData$type=="group"),]
        pathwayData <- pathwayData[!(pathwayData$graphicalType=="line"),]

        pathwayGraph <- removeNoise(pathwayGraph)
      }

      if (is.null(pathwayData) | nrow(pathwayData) == 0) {
        # Save the log file
        printLog(message_='The pathwayData data frame is empty. Skipping it...',
                 file_='generatePathwayFrequencyFromOrganismData')

        # Increment the index
        kgml_index <<- kgml_index + 1

        return(FALSE)
      }

      # Assign the reaction type to each node
      if (length(current_kgml$reactions) > 0) {
        for (idx in 1:nrow(pathwayData)) {
          for (idx2 in 1:length(current_kgml$reactions$name)) {
            # Check the position of the current reaction in reactions list
            if (current_kgml$reactions$name[idx2] %in% pathwayData[idx,]$reaction) {
              pathwayData[idx,]$reaction_type <- current_kgml$reactions$type[idx2]
              break()
            }
          }
        }
      }

      # Assign the dictionary ID to each node
      rows <- nrow(pathwayData)
      for (nodeIdx in 1:rows) {
        # Setup the current node
        currentNode <- pathwayData[nodeIdx,]

        # Find the current node into the dictionary
        dictId1 <- dictionary[(dictionary$x == currentNode$x) & (currentNode$y == currentNode$y) &
                                (dictionary$reaction == currentNode$reaction) & (dictionary$ec == currentNode$name), ]$id

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

      #*************************##
      # Get the graph properties #
      #*************************##

      graphProperties <- getGraphProperties(pathwayGraph)

      # Perform the graph bottleneck calculation
      iGraph <- igraph::graph_from_data_frame(pathwayGraph, directed = FALSE)
      graphBottleneck <- igraph::as_ids(getGraphBottleneck(iGraph, FALSE))

      # Assign the bottlenecks for enzyme code (ec)
      if (strcmp(reference_pathway, 'ec')) {
        for (bottleneckIdx in 1:length(graphBottleneck)) {
          pathwayData$is_bottleneck[which(pathwayData$dictID == graphBottleneck[bottleneckIdx])] <- 1
        }
      }

      #*****************************************##
      # Assign the graph properties to each node #
      #*****************************************##

      for (idx in 1:nrow(pathwayData)) {
        # Find which rows in graphProperties match with pathwayData
        rowsToMerge <- which(pathwayData[idx,]$dictID == graphProperties$node)

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
        if (strcmp(reference_pathway, 'ko')) {
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

      # Apply the node bottleneck impact
      impact <- getArticulationPointImpact(pathwayGraph)

      if (!is.null(impact) && nrow(impact) > 0) {
        for (idx in 1:nrow(impact)) {
          # Get the current impact item
          currentImpact <- impact[idx,]

          # Set the impacts
          pathwayData$bottleneckImpact[pathwayData$dictID==currentImpact$ap] <- currentImpact$impact
          pathwayData$bottleneckDisconnectedComponents[pathwayData$dictID==currentImpact$ap] <- currentImpact$noComponents
        }

        # Calculates the normalized impact
        for (idx in 1:nrow(impact)) {
          # Get the current impact item
          currentImpact <- impact[idx,]

          # Calculate the normalized impact
          pathwayMaxImpact <- max(pathwayData$bottleneckImpact)

          # Min frequency in a pathway
          pathwayMinImpact <- min(pathwayData$bottleneckImpact)

          # Normalized frequency for each protein
          pathwayData$bottleneckNormalizedImpact[pathwayData$dictID==currentImpact$ap] <-
            (pathwayData$bottleneckImpact[pathwayData$dictID==currentImpact$ap]-pathwayMinImpact)/
            (pathwayMaxImpact-pathwayMinImpact)
        }
      }

      #***************************************************##
      # Count the protein frequencies using the dictionary #
      #***************************************************##

      if (!is.null(pathwayData) & !is.null(pathwayInstancesDataSet)) {
        # Group all instances nodes by (x, y, reaction) variables
        proteinsCount <- aggregate(pathwayInstancesDataSet$freq, by=list(pathwayInstancesDataSet$x,
                                                                         pathwayInstancesDataSet$y,
                                                                         pathwayInstancesDataSet$reaction),
                                   FUN=sum, stringsAsFactors=FALSE)

        # Rename the group columns
        names(proteinsCount)[names(proteinsCount) == "x"] <- "occurrences"
        names(proteinsCount)[names(proteinsCount) == "Group.1"] <- "x"
        names(proteinsCount)[names(proteinsCount) == "Group.2"] <- "y"
        names(proteinsCount)[names(proteinsCount) == "Group.3"] <- "reaction"

        for (idx in 1:nrow(proteinsCount)) {
          # Get the group parameters
          x <- proteinsCount[idx,'x']
          y <- proteinsCount[idx,'y']
          react <- proteinsCount[idx,'reaction']

          # Find the current node into the dictionary
          current_ec <- dictionary[dictionary$x == x & dictionary$y == y & dictionary$reaction == react, ]$ec
          current_reaction <- dictionary[dictionary$x == x & dictionary$y == y & dictionary$reaction == react, ]$reaction

          # Check if the current node data was found
          if (length(current_ec) != 0 & length(current_reaction) != 0) {
            if (!is.na(current_ec) & !is.null(current_ec) & !is.na(current_reaction) & !is.null(current_reaction)) {
              # Assign the frequency into the pathwayData
              pathwayData[pathwayData$x == x & pathwayData$y == y & pathwayData$reaction == current_reaction &
                            pathwayData$name == current_ec, ]$occurrences <- proteinsCount[idx,'occurrences']
            }
          }
        }
      }

      #*******************##
      # Handle org metrics #
      #*******************##

      # Check which species have the current pathway
      org_list <- sapply(organism2pathway, function(org) {
        pathway_code %in% org
      })

      # Retrieve just the species identification
      org_list <- names(org_list[org_list==TRUE])

      # Define the number of species that have the current pathway
      available_orgs <- length(org_list)
      pathwayData$totalSpecies <- available_orgs

      #**********************##
      # Handle others metrics #
      #**********************##

      # Calculate the occurrences percentage
      pathwayData$percentage <- (pathwayData$occurrences / pathwayData$totalSpecies) * 100

      #****************************#
      # Prepare the data to export #
      #****************************#

      # Status message
      printMessage(paste0("EXPORTING PATHWAY ", pathway_code))

      # Export the pathway data
      if (!dir.exists(file.path('./output/totalFrequency/'))) {
        dir.create(file.path(paste0('./output/totalFrequency/')), showWarnings = FALSE, mode = "0775")
      }

      if (dir.exists(file.path('./output/totalFrequency/'))) {
        write.csv(pathwayData, file=paste0('./output/totalFrequency/', kgml_index, "_", pathway_code, '.csv'))
      }

      #***********************#
      # Remove temp variables #
      #***********************#

      rm(iGraph, graphProperties, graphBottleneck, pathwayData, current_kgml, pathwayInstancesDataSet,
         org_list, available_orgs, proteinsCount)

      # Increment the index
      kgml_index <<- kgml_index + 1

    }, error=function(e) {
      printMessage(e)

      # Increment the index
      kgml_index <<- kgml_index + 1

      # Save the log file
      printLog(toString(e), file_=paste0('generatePathwayFrequencyFromOrganismData', pathway_code))

      return(FALSE)
    })
  }) # End of Loop 01
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

printInteractiveNetwork <- function(removeNoise_=TRUE) {

  # Reference pathway
  reference_pathway <- 'ec'

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

    # Get the pathway code
    pathway_code <- onlyNumber(file)

    # Status message
    printMessage(paste0("GENERATING PATHWAY ", pathway_code, " INTERATIVE NETWORK [", kgml_index, " OF ", available_pathways, "]"))

    tryCatch({
      # Get the network properties
      pathwayData <- paste0('./output/totalFrequency/', kgml_index, '_', pathway_code, '.RData')

      if (file.exists(pathwayData)) {
        pathwayData <- get(load(file=pathwayData))
      } else {
        printMessage(paste0("The propertie file from pathway  ", pathway_code, " could no be found. Skipping it..."))
        return(FALSE)
      }

      # Convert the pathway data into a graph
      pathwayGraph <- KGML2Graph(paste0(folder, file), replaceOrg=TRUE, orgToReplace=reference_pathway)

      # Remove unnecessary data from pathway data/graph
      if (removeNoise_) {
        pathwayGraph <- removeNoise(pathwayGraph)
      }

      #-------------------------------#
      # [GETTING THE PATHWAY DETAILS] #
      #-------------------------------#

      # Set the pathway index
      pathway_index <- -1

      # Assign the graph properties to each node
      for (idx in 1:length(pathwayDetail)) {
        # Check the position of the current pathway in pathway detail list
        if (strcmp(pathwayDetail[[idx]]$ENTRY, paste0('map', pathway_code))) {
          pathway_index <- idx
          break()
        }
      }

      #--------------------------#
      # [GENERATING THE NETWORK] #
      #--------------------------#

      # Generate the network
      if (pathway_index == -1) {
        generatedNetwork <- generateInteractiveNetwork(pathwayGraph, pathwayData, pathway_code, NULL)
      } else {
        generatedNetwork <- generateInteractiveNetwork(pathwayGraph, pathwayData, pathway_code, pathwayDetail[[pathway_index]])
      }

      # Print the network
      print(generatedNetwork)

      # Export the network
      if (!dir.exists(file.path(paste0('./output/network/')))) {
        dir.create(file.path(paste0('./output/network/')), showWarnings = FALSE, mode = "0775")
      }

      if (dir.exists(file.path('./output/network/'))) {
        filename <- paste0(kgml_index, '_', pathway_code, '.html')

        # Save the HTML file
        visSave(generatedNetwork, file = filename, selfcontained = TRUE,
                background = "#eeefff")

        if (file.exists(filename)) {
          # Copy the file into correct directory
          file.copy(filename, paste0('./output/network/', filename), overwrite = TRUE)

          # Remove the generated file
          file.remove(filename)
        } else {
          printMessage(paste0("Network file not found. Skipping it..."))
          return(FALSE)
        }
      }
    }, error=function(e) {
      printMessage(e)

      # Save the log file
      printLog(toString(e), file_=paste0('printInteractiveNetwork', pathway_code))

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

#***************************************************************************#
# Step 1: Generate the pathways orgs intermediate files from the kgml files #
#***************************************************************************#
#generatePathwayDataFromKGML()
#generateOrganismPathwayDataFromKGML()

#************************************************************#
# Step 2: Generate the pathways all nodes list and/or graphs #
#************************************************************#
#generatePathwayAllNodes()
#generatePathwayGraphFromKGML()

#***************************************************#
# Step 3: Generate all pathways dictionary to match #
# the reference EC with each node in org pathways   #
#***************************************************#
#generateNodesDictionary()

#************************************************#
# Step 4: Perform the enzymes frequency counting #
#************************************************#
generatePathwayFrequencyFromOrganismData()

#******************************#
# Step 5: Generate the network #
#******************************#
#printInteractiveNetwork()
