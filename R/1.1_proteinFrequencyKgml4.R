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

getPathwayEnzymeKGML <- function(removeNoise_=TRUE) {

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
      print(paste0('The pathway ', pathway_code, ' could no be processed. Skipping it...'))
      return(FALSE)
    })

    # Increment the index
    kgml_index <<- kgml_index + 1
  }) # End of Loop 01

  # Function finished with success
  return(TRUE)
}

generateOrganismData <- function(removeNoise_=TRUE) {

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

  # APAGAR
  kgml_list <- kgml_list[41:60]
  kgml_index <- 41

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
          org_graph <- KGML2Graph(org_folder, replaceOrg=TRUE, orgToReplace=org)

          # Create the pathwayData dataFrame
          orgData <- org_kgml$nodes[,c('entryID', 'name', 'bgcolor')]
          orgData$freq <- 0
          orgData$pathway <- pathway_code

          # Remove unnecessary data from org data
          if (removeNoise_) {
            orgData <- orgData[!grepl("^path:", orgData$name),]
            orgData <- orgData[!grepl("^map:", orgData$name),]
            orgData <- orgData[!grepl("^cpd:", orgData$name),]
            orgData <- orgData[!grepl("^gl:", orgData$name),]
            orgData <- orgData[!grepl("^ko:", orgData$name),]

            org_graph <- removeNoise(org_graph)
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

          #****************************#
          # Prepare the data to export #
          #****************************#

          # Output 1
          pathwayOrg <- data.frame(pathway = orgData$pathway, org = org, entryID = orgData$entryID, freq = orgData$freq,
                              stringsAsFactors = FALSE)

          # Generate the list of unique entrez
          nodeList <- unique(c(org_graph$node1, org_graph$node2))
          entryIDList <- c()

          # Match the entrez list with the orgData$name entrez strings
          orgData$name <- gsub(" / ", "|", orgData$name)

          for (entrez in nodeList) {
            for (idx in 1:nrow(orgData)) {
              if ((grepl(orgData[idx,]$name, entrez))) {
                entryIDList <- c(entryIDList, orgData[idx,]$entryID)
                break()
              }
            }
          }

          # Remove the entrez prefix
          nodeList <- gsub(paste0(org, ":"), "", nodeList)

          # Output 2
          if (length(nodeList) == length(entryIDList)) {
            pathwayOrgWithEntrez <- data.frame(pathway = pathway_code, org = org, entryID = entryIDList, entrez = nodeList,
                                stringsAsFactors = FALSE)
          } else {
            pathwayOrgWithEntrez <- data.frame(pathway = pathway_code, org = org, entryID = NA, entrez = nodeList,
                                stringsAsFactors = FALSE)

            # Status message
            err <- paste0('Error retriving the entryID list from ', org, ' into ', pathway_code, ' pathway. Skipping it...')
            print(err)

            # Save the log file
            if (!dir.exists(file.path('./log/'))) {
              dir.create(file.path('./log/'), showWarnings = FALSE, mode = "0775")
            }

            write(err, file=paste0('./log/', format(Sys.time(), "%Y%m%d_%H%M%S_"), org, pathway_code, '.txt'))
          }

          # Export the pathway data
          if (dir.exists(file.path(current_dir))) {
            save(pathwayOrg, file=paste0(current_dir, '/', org_index, "_", org, '.RData'), compress = "xz")
            save(pathwayOrgWithEntrez, file=paste0(current_dir, '/', org_index, "_", org, '_entrez', '.RData'), compress = "xz")
          }
        }, error=function(e) {
          # Status message
          err <- paste0('The org ', org, ' could no be processed. Skipping it...')
          print(err)

          # Save the log file
          if (!dir.exists(file.path('./log/'))) {
            dir.create(file.path('./log/'), showWarnings = FALSE, mode = "0775")
          }

          write(err, file=paste0('./log/', format(Sys.time(), "%Y%m%d_%H%M%S_"), org, pathway_code, '.txt'))

          return(FALSE)
        })

        # Increment the index
        org_index <<- org_index + 1

        #***********************#
        # Remove temp variables #
        #***********************#

        rm(org_folder, org_kgml, org_graph, orgData)

      }) # End of Loop 02
    }, error=function(e) {
      # Status message
      err <- paste0('The pathway ', pathway_code, ' could no be processed. Skipping it...')
      print(err)

      # Save the log file
      if (!dir.exists(file.path('./log/'))) {
        dir.create(file.path('./log/'), showWarnings = FALSE, mode = "0775")
      }

      write(err, file=paste0('./log/', format(Sys.time(), "%Y%m%d_%H%M%S_"), org, pathway_code, '.txt'))

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
      print(paste0('The pathway ', pathway_code, ' could no be processed. Skipping it...'))
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

#***************************************************************#
# Step 1: Generate the pathways frequencies from the kgml files #
#***************************************************************#
#getPathwayEnzymeKGML()
generateOrganismData()

#******************************#
# Step 2: Generate the network #
#******************************#
#printInteractiveNetwork()
