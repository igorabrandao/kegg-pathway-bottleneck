#**************************************************#
# Pipeline to perform the enzyme frequencies count #
#**************************************************#

# ---- IMPORT SECTION ----

# 1.0_proteinFrequency.R #

#' This is the pipeline script to perform
#' the enzymes frequencies counting
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library(foreach)

#*******************************************************************************************#

# ---- SETTINGS SECTION ----

#*************************#
# Pipeline basic settings #
#*************************#

# Import the graphLoader functions
files.sources = NULL
files.sources[1] = paste0("./R/functions", "/", "graphFunctions.R")
#files.sources[2] = paste0("./R/functions", "/", "graphPrintFunctions.R")
files.sources[2] = paste0("./R/functions", "/", "helperFunctions.R")
sapply(files.sources, source)

# Load the pathways by organisms data
organism2pathway <- get(load(paste0("./dictionaries", "/", "organism2pathway.RData")))
pathwayList <- get(load(paste0("./dictionaries", "/", "pathwayList.RData")))

# Define in which specie the processing should begin
# default value 1 (the value should be >= 1)
start_of <- 1

#*******************************************************************************************#

# ---- FUNCTIONS SECTION ----

#*************************#
# Pipeline main functions #
#*************************#

#' Get the list of pathways from pathwayList and export
#' the enzyme list for each specie by pathway
#'
#' @param index_ Index from pathwayList representing a single pathway, e.g: 1 = 00010.
#' @param removeNoise_ Remove undesirable enzyme such as: map, path, cpd or gl.
#' @param replaceEmptyGraph_ If a reference pathway does not exist, it can be replaced by an
#' organism specific pathway when its TRUE, cc just ignore the pathway and anotate in log.
#' @param chunkSize_ During the requests to KEGG DB to avoid long time o waiting you can set
#' the size of list sent to http request.
#' @param specieRangeMin_ Represent the start index to look up the species
#' @param specieRangeMax_ Represent the end index point to look up the species
#' @param pathwayData_ Preloaded pathwayData
#'
#' @return This function does not return nothing, just export files.
#'
#' @examples
#' \dontrun{
#' getPathwayEnzymes(1, TRUE, TRUE, 50, 1, length(organism2pathway))
#' }
#'
#' @author
#' Igor Brandão

getPathwayEnzymes <- function(index_, removeNoise_=TRUE, replaceEmptyGraph_=TRUE, chunkSize_=50,
                              specieRangeMin_=1, specieRangeMax_=length(organism2pathway),
                              pathwayData_=NULL) {

  #******************************#
  # Get the pathway general info #
  #******************************#

  # Get the current pathway
  pathway <- pathwayList[index_,]

  # Format the pathway code
  pathway_code <- paste0('ec', pathway)

  # Count the total of species
  totalSpecies <- length(organism2pathway)

  #***********************************#
  # Load all enzymes from its pathway #
  #***********************************#

  # Get the enzyme list from pathway
  if (is.null(pathwayData_)) {
    # Status message
    printMessage(paste0("COUNTING ", pathway, " ENZYMES FREQUENCIES [", index_, " OF ", nrow(pathwayList), "]"))

    # Perform the web request to retrieve the data
    pathwayData <- pathwayToDataframe(pathway_code, FALSE)
  } else {
    # Use the data from the parameter
    pathwayData <- pathwayData_
  }

  # Handle empty graph
  if (is.null(pathwayData) | length(pathwayData) == 0) {
    pathwayData <- data.frame(node1 = NA, org = 'ec', pathway = pathway, is_bottleneck = 0,
                       is_presented = 0, stringsAsFactors = FALSE)

    return(pathwayData)

  } else {
    #*************************##
    # Prepare the pathway data #
    #*************************##

    # Remove unnecessary data before bottleneck calculation
    if (removeNoise_) {
      pathwayData <- pathwayData[!grepl("^path:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^path:", pathwayData$node2),]
      pathwayData <- pathwayData[!grepl("^map:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^map:", pathwayData$node2),]
      pathwayData <- pathwayData[!grepl("^cpd:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^cpd:", pathwayData$node2),]
      pathwayData <- pathwayData[!grepl("^gl:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^gl:", pathwayData$node2),]
    }

    # Get the graph properties
    graphProperties <- getGraphProperties(pathwayData)

    # Calculates the network bottleneck
    iGraph <- igraph::graph_from_data_frame(pathwayData, directed = FALSE)

    # Perform the graph bottleneck calculation
    graphBottleneck <- igraph::as_ids(getGraphBottleneck(iGraph, FALSE))

    # Convert the node2 column into node1 index_s
    aux <- unique(c(pathwayData$node1, pathwayData$node2))
    auxorg <- pathwayData$org[1]
    auxpathway <- pathwayData$pathway[1]

    # Add a new column to the enzymeFrquency dataFrame
    if (is.null(graphProperties)) {
      pathwayData <- data.frame(node1 = aux, org = auxorg, pathway = auxpathway, is_bottleneck = 0,
                                is_presented = 0, betweenness = NA, connectivity = NA, triangles = NA, clusteringCoef = NA,
                                closenessCoef = NA, community = NA, eigenvectorScore = NA, eccentricity = NA,
                                radius = NA, diameter = NA, stringsAsFactors = FALSE)
    } else {
      pathwayData <- data.frame(node1 = aux, org = auxorg, pathway = auxpathway, is_bottleneck = 0,
                                is_presented = 0, betweenness = graphProperties$betweenness, connectivity = graphProperties$connectivity,
                                triangles = graphProperties$triangles, clusteringCoef = graphProperties$clusteringCoef,
                                closenessCoef = graphProperties$closenessCoef, community = graphProperties$community,
                                eigenvectorScore = graphProperties$eigenvectorScore, eccentricity = graphProperties$eccentricity,
                                radius = graphProperties$radius, diameter = graphProperties$diameter,
                                stringsAsFactors = FALSE)
    }

    # Assign the bottlenecks
    pathwayData$is_bottleneck[which(pathwayData[,1] %in% graphBottleneck)] <- 1
  }

  #************************************#
  # Set the parameters for each specie #
  #************************************#

  # Loop over the organism list
  enzymeList <- foreach::foreach(idx = seq.int(specieRangeMin_, specieRangeMax_),
                                 .export=c('printMessage', 'pathwayToDataframe', 'as_ids', 'str_replace',
                                           'getGraphBottleneck', 'convertEntrezToECWithoutDict',
                                           'getPathwayHighlightedGenes'),
                                 .combine = "rbind") %do%
  {
    tryCatch({
      # Get its name
      specie <- names(organism2pathway[idx])

      # Status message
      printMessage(paste0("<<< Working on ", specie, pathway, " pathway... >>>"))

      # Get the pathway graph and change the column org with the current specie
      temp <- pathwayData
      temp[,c('org')] <- specie

      #**********************************************#
      # Check if the enzyme appears into the pathway #
      #**********************************************#

      # Get the highlighted enzymes list
      highlighted_enzymes <- getPathwayHighlightedGenes(paste0(specie, pathway))

      if (!is.null(highlighted_enzymes)) {
        # Concat the org string
        highlighted_enzymes <- paste(specie, highlighted_enzymes, sep=":")

        # Convert it into dataframe
        highlighted_enzymes <- as.data.frame(highlighted_enzymes, stringsAsFactors = FALSE)

        # Convert the highlighted list into EC number
        highlighted_enzymes <- convertEntrezToECWithoutDict(highlighted_enzymes[,c(1)], chunkSize_, FALSE)

        # Remove the duplicates
        highlighted_enzymes <- highlighted_enzymes[!duplicated(highlighted_enzymes),]

        # Get just the enzyme number without specie
        current_enzyme <- gsub("^[[:alpha:]]*(.*$)", "\\1", str_replace(temp$node1, ":", ""))

        # Verify if the current enzyme is highlighted and set its status
        temp$is_presented[which(current_enzyme %in% highlighted_enzymes)] <- 1

      } else {
        err <- paste0("Pathview library doesn't support the pathway: ", specie, pathway)

        # Save the log file
        if (dir.exists(file.path('./log/'))) {
          write(err, file=paste0('./log/', format(Sys.time(), "%Y%m%d_%H%M%S_"), specie, pathway, '.txt'))
        }
      }

      # Save the pathway intermediary data
      names(temp)[names(temp) == "node1"] <- "ec"

      if (!dir.exists(file.path(paste0('./output/', pathway)))) {
        dir.create(file.path(paste0('./output/', pathway)), showWarnings = FALSE, mode = "0775")
      }

      if (dir.exists(file.path('~/data3/'))) {
        dir.create(file.path(paste0('~/data3/kegg-pathway-bottleneck/output/', pathway)), showWarnings = FALSE, mode = "0775")
        save(temp, file=paste0('~/data3/kegg-pathway-bottleneck/output/', pathway, '/', idx, '_',  specie, '.RData'))
      }

      save(temp, file=paste0('./output/', pathway, '/', idx, '_',  specie, '.RData'))

      # Export the highlighted enzymes
      if (!dir.exists(file.path('./output/highlightedEnzymes/'))) {
        dir.create(file.path(paste0('./output/highlightedEnzymes/')), showWarnings = FALSE, mode = "0775")
      }

      if (!dir.exists(file.path(paste0('./output/highlightedEnzymes/', pathway)))) {
        dir.create(file.path(paste0('./output/highlightedEnzymes/', pathway)), showWarnings = FALSE, mode = "0775")
      }

      if (!is.null(highlighted_enzymes) && length(highlighted_enzymes) != 0) {
        tempHighlight <- data.frame(highlighted_enzymes = highlighted_enzymes, specie = specie,
                                    pathway = pathway, timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

        save(tempHighlight, file=paste0('./output/highlightedEnzymes/', pathway, '/', idx, '_',  specie, '.RData'))
      }

      # Return the specie [FOREACH]
      return(temp)

    }, error=function(e) {
      # Status message
      printMessage(paste0('The pathway ', specie, pathway, ' could no be processed. View the log file for more information. Skipping it...'))

      # Save the error message
      err <- conditionMessage(e)

      # Save the log file
      if (dir.exists(file.path('./log/'))) {
        write(err, file=paste0('./log/', format(Sys.time(), "%Y%m%d_%H%M%S_"), specie, pathway, '.txt'))
      }

      # Add a new column to the enzymeFrquency dataFrame
      temp <- data.frame(node1 = NA, org = specie, pathway = pathway, is_bottleneck = 0,
                         is_presented = 0, stringsAsFactors = FALSE)

      return(temp)
    })

  } # END OF FOREACH

  #****************************#
  # Prepare the data to export #
  #****************************#

  # Rename the nodes column
  names(enzymeList)[names(enzymeList) == "node1"] <- "ec"

  # Remove duplicated index_s based on EC column
  enzymeList <- enzymeList[!duplicated(enzymeList[c("ec", "org", "pathway")]),]

  # Reindex the enzymeList index_s
  rownames(enzymeList) <- 1:nrow(enzymeList)

  # Export the pathway data
  if (dir.exists(file.path('./output/'))) {
    save(enzymeList, file=paste0('./output/', index_, "_", pathway, '.RData'))
  }

  if (dir.exists(file.path('~/data3/'))) {
    save(enzymeList, file=paste0('~/data3/kegg-pathway-bottleneck/output/', index_, "_", pathway, '.RData'))
  }

  # Function finished with success
  return(TRUE)
}

#' Get the list of pathways from pathwayList and export
#' the total frequency for each pathway
#'
#' @param index_ Index from pathwayList representing a single pathway, e.g: 1 = 00010.
#' @param resumeInfo_ If TRUE export only summary data without each specie columns.
#'
#' @return This function does not return nothing, just export files.
#'
#' @examples
#' \dontrun{
#' getTotalFrequency(1)
#' }
#'
#' @author
#' Igor Brandão

getTotalFrequency <- function(index_, resumeInfo_=TRUE) {

  # Get the current pathway
  pathway <- pathwayList[index_,]

  # Status message
  printMessage(paste0("CALCULATING TOTAL FREQUENCY OF PATHWAY ", pathway))

  # Data frame to merge all data
  enzymeList <- NULL

  # Get the list of files
  folder = paste0("./output/", pathway, "/")
  file_list <- list.files(path=folder, pattern='*.RData')

  # Define the number of species
  total_species <- length(file_list)

  # Check if the folder contains files
  if (is.null(file_list) | length(file_list) == 0) {
    return(FALSE)
  }

  # Load all files at once
  big.list.of.data.frames <- lapply(file_list, function(file) {
    get(load(file=paste0(folder, file)))
  })

  # Combine multiple data frames in one
  enzymeList <- do.call(rbind, big.list.of.data.frames)

  # Remove temporaly variables
  rm(big.list.of.data.frames)

  #-------------------------------#
  # [TOTAL FREQUENCY CALCULATION] #
  #-------------------------------#

  # Get the enzymeList as a 3 way table
  enzymeTotalFrequency <- table(enzymeList$ec,enzymeList$org,enzymeList$is_presented)

  tryCatch({
    # Get just the table sinalyzing the presence of the enzyme
    enzymeTotalFrequency <- enzymeTotalFrequency[,,2]

    # Count the enzymes frequencies and transform it into a dataFrame
    enzymeTotalFrequency <- as.data.frame.matrix(enzymeTotalFrequency, stringsAsFactors = FALSE)

    # Sum each frequency
    enzymeTotalFrequency$freq <- rowSums(enzymeTotalFrequency)

    # Calculate the total frequency
    enzymeTotalFrequency$max_freq <- max(enzymeTotalFrequency$freq)

    # Count valids species (have highlights != 0)
    enzymeTotalFrequency$total_species <- (total_species - sum(sapply(enzymeTotalFrequency, function(x) all(x == 0))))

    # Calculate the frequency percentage
    enzymeTotalFrequency$percentage <- (enzymeTotalFrequency$freq / enzymeTotalFrequency$total_species) * 100

    # Calculate the mean frequency (%)
    enzymeTotalFrequency$mean <- mean(enzymeTotalFrequency$percentage)

    # Calculate the standard deviation of frequency
    enzymeTotalFrequency$std <- sd(enzymeTotalFrequency$percentage)

    #-------------------#
    # [PATHWAY METRICS] #
    #-------------------#

    # Set the is_bottleneck flag from enzymeList
    selectedEC <- enzymeList[enzymeList$ec%in%row.names(enzymeTotalFrequency),]

    # Remove duplicates EC
    selectedEC <- selectedEC[!duplicated(selectedEC[,c('ec')]),]

    # Align the column is_bottleneck
    mergeTemp <- merge(enzymeTotalFrequency, selectedEC, by.x=0, by.y="ec", all.x=T)
    enzymeTotalFrequency$is_bottleneck <- mergeTemp$is_bottleneck

    # Add the metric columns
    enzymeTotalFrequency$betweenness <- mergeTemp$betweenness
    enzymeTotalFrequency$connectivity <- mergeTemp$connectivity
    enzymeTotalFrequency$triangles <- mergeTemp$triangles
    enzymeTotalFrequency$clusteringCoef <- mergeTemp$clusteringCoef
    enzymeTotalFrequency$closenessCoef <- mergeTemp$closenessCoef
    enzymeTotalFrequency$community <- mergeTemp$community
    enzymeTotalFrequency$eigenvectorScore <- mergeTemp$eigenvectorScore
    enzymeTotalFrequency$eccentricity <- mergeTemp$eccentricity
    enzymeTotalFrequency$radius <- mergeTemp$radius
    enzymeTotalFrequency$diameter <- mergeTemp$diameter
    enzymeTotalFrequency$degree <- mergeTemp$degree
    enzymeTotalFrequency$authorityScore <- mergeTemp$authorityScore
    enzymeTotalFrequency$hubScore <- mergeTemp$hubScore
    enzymeTotalFrequency$bottleneck_classification <- mergeTemp$bottleneck_classification
    enzymeTotalFrequency$pathway <- pathway

    # Remove the quotes from column name
    colnames(enzymeTotalFrequency) <- gsub("\"", "", colnames(enzymeTotalFrequency))

    # Remove temp var
    rm(selectedEC, mergeTemp)

    # Check if its necessary remove species columns
    if (resumeInfo_) {
      enzymeTotalFrequency <- enzymeTotalFrequency[,(ncol(enzymeTotalFrequency)-20):ncol(enzymeTotalFrequency)]
    }

    # Export the pathway data
    if (!dir.exists(file.path(paste0('./output/totalFrequency/', pathway)))) {
      dir.create(file.path(paste0('./output/totalFrequency/')), showWarnings = FALSE, mode = "0775")
    }

    if (dir.exists(file.path('./output/totalFrequency/'))) {
      save(enzymeTotalFrequency, file=paste0('./output/totalFrequency/', index_, "_", pathway, '.RData'))
    }

    if (dir.exists(file.path('~/data3/'))) {
      save(enzymeTotalFrequency, file=paste0('~/data3/kegg-pathway-bottleneck/output/', index_, "_", pathway, '.RData'))
    }

    # Function finished with success
    return(TRUE)
  }, error=function(e) {
    # Status message
    printMessage(paste0("WARNING: THE PATHWAY ", pathway, " DOESNT HAVE ANY PRESENTED PROTEIN."))

    # Function finished with error
    return(NULL)
  })
}

#' Recalculates the pathway properties
#'
#' @param index_ Index from pathwayList representing a single pathway, e.g: 1 = 00010.
#' @param removeNoise_ Remove undesirable enzyme such as: map, path, cpd or gl.
#'
#' @return This function does not return nothing, just export files.
#'
#' @examples
#' \dontrun{
#' reapplyGraphProperties(1)
#' }
#'
#' @author
#' Igor Brandão

reapplyGraphProperties <- function(index_, removeNoise_=TRUE) {

  # Get the current pathway
  pathway <- pathwayList[index_,]

  # Status message
  printMessage(paste0("RECALCULATING ", pathway, " GRAPH METRICS"))

  # Get the list of files
  folder = paste0("./output/", pathway, "/")
  file_list <- list.files(path=folder, pattern='*.RData')

  # Check if the folder contains files
  if (is.null(file_list) | length(file_list) == 0) {
    return(FALSE)
  }

  # Format the pathway code
  pathway_code <- paste0('ec', pathway)

  # Get the enzyme list from pathway
  pathwayData <- pathwayToDataframe(pathway_code, FALSE)

  # Handle empty graph
  if (is.null(pathwayData) | length(pathwayData) == 0) {
    return(FALSE)
  } else {
    # Remove unnecessary data before properties calculation
    if (removeNoise_) {
      pathwayData <- pathwayData[!grepl("^path:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^path:", pathwayData$node2),]
      pathwayData <- pathwayData[!grepl("^map:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^map:", pathwayData$node2),]
      pathwayData <- pathwayData[!grepl("^cpd:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^cpd:", pathwayData$node2),]
      pathwayData <- pathwayData[!grepl("^gl:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^gl:", pathwayData$node2),]
    }

    # Get the graph properties
    graphProperties <- getGraphProperties(pathwayData)

    # Change dataFrames inside the files
    lapply(file_list, function(file) {
      # Load the dataframe
      temp <- get(load(file=paste0(folder, file)))

      # Update its columns
      temp$betweenness <- graphProperties$betweenness
      temp$connectivity <- graphProperties$connectivity
      temp$triangles <- graphProperties$triangles
      temp$clusteringCoef <- graphProperties$clusteringCoef
      temp$closenessCoef <- graphProperties$closenessCoef
      temp$community <- graphProperties$community
      temp$eigenvectorScore <- graphProperties$eigenvectorScore
      temp$eccentricity <- graphProperties$eccentricity
      temp$radius <- graphProperties$radius
      temp$diameter <- graphProperties$diameter
      temp$degree <- graphProperties$degree
      temp$authorityScore <- graphProperties$authorityScore
      temp$hubScore <- graphProperties$hubScore

      # Classify the bottleneck
      temp <- classifyBottleneck(temp)

      # Export the pathway data
      if (dir.exists(file.path('./output/'))) {
        save(temp, file=paste0(folder, file))
      }
    })
  }

  # Function finished with success
  return(TRUE)
}

#' Check each specie for each pathway and try to fill the species that doesn't
#' have any enzymes present
#'
#' @param index_ Index from pathwayList representing a single pathway, e.g: 1 = 00010.
#'
#' @return This function does not return nothing, just export files.
#'
#' @examples
#' \dontrun{
#' fillMissingEnzymesPresence(1)
#' }
#'
#' @author
#' Igor Brandão

fillMissingEnzymesPresence <- function(index_) {

  # Get the current pathway
  pathway <- pathwayList[index_,]

  # Format the pathway code
  pathway_code <- paste0('ec', pathway)

  # Status message
  printMessage(paste0("FILLING MISSING ENZYME PRESENCE IN PATHWAY ", pathway))

  # Get the list of files
  folder = paste0("./output/", pathway, "/")
  file_list <- list.files(path=folder, pattern='*.RData')

  # Check if the folder contains files
  if (is.null(file_list) | length(file_list) == 0) {
    return(FALSE)
  }

  #***********************************#
  # Load all enzymes from its pathway #
  #***********************************#

  # Try to load the pathway data 10 times
  for (i in 1:30) {
    # Get the enzyme list from pathway
    pathwayData <- pathwayToDataframe(pathway_code, FALSE)

    if (!is.null(pathwayData)) {
      break
    }
  }

  # Check if the pathwayData contains data
  if (is.null(pathwayData) | length(pathwayData) == 0) {
    return(FALSE)
  }

  # Verify each specie (by file)
  lapply(file_list, function(file) {

    tryCatch({
      # Set the specie index
      file_idx <- as.numeric(gsub("([0-9]+).*$", "\\1", file))
      print(file)

      # Load the dataframe
      temp <- get(load(file=paste0(folder, file)))

      # Check if the file is valid
      if (is.null(temp) | length(temp) == 0) {
        return(FALSE)
      } else {
        # Count the enzymes presence
        presence_count <- (nrow(temp) - sum(sapply(temp$is_presented, function(x) all(x == 0))))

        # Verify if the specie doesn't have any presence
        if (presence_count == 0) {
          # Try to get the enzymes presence info
          getPathwayEnzymes(index_, removeNoise_=TRUE, replaceEmptyGraph_=TRUE, chunkSize_=50,
                                    specieRangeMin_=file_idx, specieRangeMax_=file_idx,
                                    pathwayData_=pathwayData)
        }
      }
    }, error=function(e) {
      # If some error happened, reload the specie dataset
      getPathwayEnzymes(index_, removeNoise_=TRUE, replaceEmptyGraph_=TRUE, chunkSize_=50,
                        specieRangeMin_=file_idx, specieRangeMax_=file_idx,
                        pathwayData_=pathwayData)
    })
  })

  # Function finished with success
  return(TRUE)
}

#' Function to generate correlation matrix
#'
#' @param index_ Index from pathwayList representing a single pathway, e.g: 1 = 00010.
#'
#' @return This function does not return nothing, just export files.
#'
#' @examples
#' \dontrun{
#' generateCorrelationStudy(1)
#' }
#'
#' @author
#' Igor Brandão

generateCorrelationStudy <- function(index_) {

  # Get the current pathway
  pathway <- pathwayList[index_,]

  # Status message
  printMessage(paste0("GENERATING ", pathway, " CORRELATION STUDY"))

  # Get the network properties
  propertyFile <- paste0('./output/totalFrequency/', index_, '_', pathway, '.RData')

  if (file.exists(propertyFile)) {
    networkProperties <- get(load(file=propertyFile))
  } else {
    printMessage(paste0("The propertie file from pathway  ", pathway, " could no be found. Skipping it..."))
    return(FALSE)
  }

  #---------------------------#
  # [GENERATING CORRELATIONS] #
  #---------------------------#

  # Create the correlations
  correlations <- cor(networkProperties)

  # Export the network
  if (!dir.exists(file.path(paste0('./output/correlation/')))) {
    dir.create(file.path(paste0('./output/correlation/')), showWarnings = FALSE, mode = "0775")
  }

  if (dir.exists(file.path('./output/correlation/'))) {
    save(correlations, file=paste0('./output/correlation/', index_, "_", pathway, '.RData'))
  }

  # Function finished with success
  return(TRUE)
}

#' Function to generate interactive networks
#'
#' @param index_ Index from pathwayList representing a single pathway, e.g: 1 = 00010.
#' @param removeNoise_ Remove undesirable enzyme such as: map, path, cpd or gl.
#'
#' @return This function does not return nothing, just export files.
#'
#' @examples
#' \dontrun{
#' printInteractiveNetwork(1)
#' }
#'
#' @author
#' Igor Brandão

printInteractiveNetwork <- function(index_, removeNoise_=TRUE) {

  # Get the current pathway
  pathway <- pathwayList[index_,]

  # Status message
  printMessage(paste0("GENERATING ", pathway, " INTERATIVE NETWORK"))

  # Format the pathway code
  pathway_code <- paste0('ec', pathway)

  # Get the enzyme list from pathway
  pathwayData <- pathwayToDataframe(pathway_code, FALSE)

  # Handle empty graph
  if (is.null(pathwayData) | length(pathwayData) == 0) {
    return(FALSE)
  } else {
    # Remove unnecessary data before properties calculation
    if (removeNoise_) {
      pathwayData <- pathwayData[!grepl("^path:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^path:", pathwayData$node2),]
      pathwayData <- pathwayData[!grepl("^map:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^map:", pathwayData$node2),]
      pathwayData <- pathwayData[!grepl("^cpd:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^cpd:", pathwayData$node2),]
      pathwayData <- pathwayData[!grepl("^gl:", pathwayData$node1),]
      pathwayData <- pathwayData[!grepl("^gl:", pathwayData$node2),]
    }

    # Get the network properties
    propertyFile <- paste0('./output/totalFrequency/', index_, '_', pathway, '.RData')

    if (file.exists(propertyFile)) {
      networkProperties <- get(load(file=propertyFile))
    } else {
      printMessage(paste0("The propertie file from pathway  ", pathway, " could no be found. Skipping it..."))
      return(FALSE)
    }

    #--------------------------#
    # [GENERATING THE NETWORK] #
    #--------------------------#

    # Generate the network
    generatedNetwork <- generateInteractiveNetwork(pathwayData, networkProperties, pathway)

    # Print the network
    print(generatedNetwork)

    # Export the network
    if (!dir.exists(file.path(paste0('./output/network/')))) {
      dir.create(file.path(paste0('./output/network/')), showWarnings = FALSE, mode = "0775")
    }

    if (dir.exists(file.path('./output/network/'))) {
      filename <- paste0(index_, '_', pathway, '.html')

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
  }

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
#lapply(16:20, getPathwayEnzymes, replaceEmptyGraph_=FALSE)

# Call the function for all pathways
#lapply(start_of:nrow(pathwayList), getPathwayEnzymes, replaceEmptyGraph_=FALSE)

#*******************************************************************************************#

#********************************************#
# Step 2: Group all files into one dataframe #
# Step 3: Count the enzyme frequency         #
#********************************************#

# [TEST ONLY]
#lapply(1:1, getTotalFrequency)

# Call the function for all pathways
#lapply(start_of:nrow(pathwayList), getTotalFrequency)

#*******************************************************************************************#

#********************************#
# [OPTIONAL]                     #
# Recalculates the graph metrics #
#********************************#

# [TEST ONLY]
#lapply(1:5, reapplyGraphProperties)

# Call the function for all pathways
#lapply(start_of:nrow(pathwayList), reapplyGraphProperties)

#********************************#
# [OPTIONAL]                     #
# Fill missing enzymes presence  #
#********************************#

# [TEST ONLY]
lapply(16:20, fillMissingEnzymesPresence)

# Call the function for all pathways
#lapply(start_of:nrow(pathwayList), fillMissingEnzymesPresence)

#*******************************************************************************************#

#****************************************#
# Step 4: Generate the correlation study #
#****************************************#

# [TEST ONLY]
#lapply(1:100, generateCorrelationStudy)

# Call the function for all pathways
#lapply(start_of:nrow(pathwayList), generateCorrelationStudy)

#*******************************************************************************************#

#******************************#
# Step 5: Generate the network #
#******************************#

# [TEST ONLY]
#lapply(1:1, printInteractiveNetwork)

# Call the function for all pathways
#lapply(start_of:nrow(pathwayList), printInteractiveNetwork)
