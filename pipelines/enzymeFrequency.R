####################################################
# Pipeline to perform the enzyme frequencies count #
####################################################

# enzymeList.R #

#' This is the pipeline script to perform
#' the enzymes frequencies counting
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library(KEGGREST) # graph handler
library(KEGGgraph)
library(pathview)
library(igraph)
library(RCurl) # http connections
library(rvest) # web scraping
library(stringr) # regex manipulation
library(pracma) # string manipulation
library(foreach)

#-------------------------------------------------------------------------------------------#

###########################
# Pipeline basic settings #
###########################

# Import the graphLoader functions
files.sources = NULL
files.sources[1] = paste0("./R", "/", "graphLoader.R")
files.sources[2] = paste0("./R", "/", "graphBottleneck.R")
sapply(files.sources, source)

# Load the pathways by organisms data
organism2pathway <- get(load(paste0("./dictionnaires", "/", "organism2pathway.RData")))
pathwayList <- get(load(paste0("./dictionnaires", "/", "pathwayList.RData")))

# Define in which specie the processing should begin
# default value 1 (the value should be >= 1)
start_of <- 1

# Auxiliar function to generate messages
printMessage <- function(message_) {
  cat("\n")
  print("------------------------------------------------")
  print(message_)
  print("------------------------------------------------")
  cat("\n")
}

#-------------------------------------------------------------------------------------------#

###########################
# Pipeline main functions #
###########################

getPathwayEnzymes <- function(index_, removeNoise_=TRUE, replaceEmptyGraph_=TRUE, chunkSize_=50) {

  ################################
  # Get the pathway general info #
  ################################

  # Get the current pathway
  pathway <- pathwayList[index_,]

  # Format the pathway code
  pathway_code <- paste0('ec', pathway)

  # Count the total of species
  totalSpecies <- length(organism2pathway)
  totalSpecies <- 10 # TODO: REMOVER DEPOIS

  # Status message
  printMessage(paste0("COUNTING ", pathway, " ENZYMES FREQUENCIES [", index_, " OF ", nrow(pathwayList), "]"))

  #####################################
  # Load all enzymes from its pathway #
  #####################################

  # Get the enzyme list from pathway
  pathwayData <- pathwayToDataframe(pathway_code, FALSE)

  # Handle empty graph
  if (is.null(pathwayData) | length(pathwayData) == 0) {
    pathwayData <- data.frame(node1 = NA, org = specie, pathway = pathway, is_bottleneck = 0,
                       is_presented = 0, stringsAsFactors = FALSE)

    return(temp)
  }

  ######################################
  # Set the parameters for each specie #
  ######################################

  # Loop over the organism list
  enzymeList <- foreach::foreach(idx = seq.int(1, totalSpecies),
                                 .export=c('printMessage', 'pathwayToDataframe', 'as_ids', 'str_replace',
                                           'getGraphBottleneck', 'convertEntrezToECWithoutDict',
                                           'getPathwayHighlightedGenes'),
                                 .combine = "rbind") %do%
  {
    tryCatch({
      # Get its name
      specie <- names(organism2pathway[idx])

      # Get the pathway graph and change the column org with the current specie
      temp <- pathwayData
      temp[,c('org')] <- specie

      # Status message
      printMessage(paste0("<<< Working on ", specie, pathway, " pathway... >>>"))

      # Remove unnecessary data before bottleneck calculation
      if (removeNoise_) {
        temp <- temp[!grepl("^path:", temp$node1),]
        temp <- temp[!grepl("^path:", temp$node2),]
        temp <- temp[!grepl("^map:", temp$node1),]
        temp <- temp[!grepl("^map:", temp$node2),]
        temp <- temp[!grepl("^cpd:", temp$node1),]
        temp <- temp[!grepl("^cpd:", temp$node2),]
        temp <- temp[!grepl("^gl:", temp$node1),]
        temp <- temp[!grepl("^gl:", temp$node2),]
      }

      # Calculates the network bottleneck
      iGraph <- igraph::graph_from_data_frame(temp, directed = FALSE)

      # Perform the graph bottleneck calculation
      graphBottleneck <- igraph::as_ids(getGraphBottleneck(iGraph, FALSE))

      # Convert the node2 column into node1 index_s
      aux <- unique(c(temp$node1, temp$node2))
      auxorg <- temp$org[1]
      auxpathway <- temp$pathway[1]

      # Add a new column to the enzymeFrquency dataFrame
      temp <- data.frame(node1 = aux, org = auxorg, pathway = auxpathway, is_bottleneck = 0,
                         is_presented = 0, stringsAsFactors = FALSE)

      # Assign the bottlenecks
      temp$is_bottleneck[which(temp[,1] %in% graphBottleneck)] <- 1

      ################################################
      # Check if the enzyme appears into the pathway #
      ################################################

      # Get the highlighted enzymes list
      highlighted_enzymes <- getPathwayHighlightedGenes(paste0(specie, pathway), allMapped_=TRUE)

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

      # Save the intermediary data
      if (!dir.exists(file.path(paste0('./output/', pathway)))) {
        dir.create(file.path(paste0('./output/', pathway)), showWarnings = FALSE)
      }

      if (dir.exists(file.path('~/data3/'))) {
        dir.create(file.path(paste0('~/data3/kegg-pathway-bottleneck/output/', pathway)), showWarnings = FALSE)
        save(temp, file=paste0('~/data3/kegg-pathway-bottleneck/output/', pathway, '/', idx, '_',  specie, '.RData'))
      }

      save(temp, file=paste0('./output/', pathway, '/', idx, '_',  specie, '.RData'))

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

  ##############################
  # Prepare the data to export #
  ##############################

  # Rename the nodes column
  names(enzymeList)[names(enzymeList) == "node1"] <- "ec"

  # Remove duplicated index_s based on EC column
  enzymeList <- enzymeList[!duplicated(enzymeList[c("ec", "org", "pathway")]),]

  # Reindex the enzymeList index_s
  rownames(enzymeList) <- 1:nrow(enzymeList)

  # Export the specie data
  if (dir.exists(file.path('./output/'))) {
    save(enzymeList, file=paste0('./output/', pathway, '.RData'))
  }

  if (dir.exists(file.path('~/data3/'))) {
    save(enzymeList, file=paste0('~/data3/kegg-pathway-bottleneck/output/', pathway, '.RData'))
  }

  # Function finished with success
  return(TRUE)
}

#-------------------------------------------------------------------------------------------#

####################################
# Step 1: Get all pathways enzymes #
####################################

#lapply(start_of:1, getPathwayEnzymes, replaceEmptyGraph_=FALSE)
lapply(start_of:start_of, getPathwayEnzymes, replaceEmptyGraph_=FALSE)

#enzymeList <- do.call("rbind", enzymeList)

lapply(start_of:nrow(pathwayList), getPathwayEnzymes, replaceEmptyGraph_=FALSE)

#-------------------------------------------------------------------------------------------#

##############################################
# Step 2: Group all files into one dataframe #
##############################################

# Data frame to merge all data
enzymeList <- NULL

# Get the list of files
folder = "./output/"
file_list <- list.files(path=folder, pattern='*.RData')

# Load all files at once
big.list.of.data.frames <- lapply(file_list, function(file) {
  get(load(file=paste0(folder, file)))
})

# Combine multiple data frames in one
enzymeList <- do.call(rbind, big.list.of.data.frames)

# Remove temporaly variables
rm(big.list.of.data.frames)

#-------------------------------------------------------------------------------------------#

######################################
# Step 3: Count the enzyme frequency #
######################################

#-------------------#
# [TOTAL FREQUENCY] #
#-------------------#

# Get the enzymeList as a 3 way table
enzymeTotalFrequency <- table(enzymeList$ec,enzymeList$org,enzymeList$is_presented)

# Get just the table sinalyzing the presence of the enzyme
enzymeTotalFrequency <- enzymeTotalFrequency[,,2]

# Count the enzymes frequencies and transform it into a dataFrame
enzymeTotalFrequency <- as.data.frame.matrix(enzymeTotalFrequency, stringsAsFactors = FALSE)

# Add columns according to the specie list
enzymeTotalFrequency <- cbind(enzymeTotalFrequency, as.list(c('is_bottleneck', 'freq', 'total_freq')))

# Remove the quotes from column name
colnames(enzymeTotalFrequency) <- gsub("\"", "", colnames(enzymeTotalFrequency))

# Set the is_bottleneck flag from enzymeList (not working)
sapply(enzymeTotalFrequency, function(idx)
  enzymeTotalFrequency[idx, ]$is_bottleneck <- enzymeList[enzymeList$ec == row.names(enzymeTotalFrequency[idx, ]),]$is_bottleneck[1]
)


#------------------------#
# [FREQUENCY BY PATHWAY] #
#------------------------#

# Filter just the enzymes with some frequency
enzymeFrequencyByPathway <- enzymeList[enzymeList$is_presented == 1,]

# Count the enzyme frequency by pathway code
enzymeFrequencyByPathway <- aggregate(x = enzymeFrequencyByPathway,
                                      by = list(enzymeFrequencyByPathway$ec, enzymeFrequencyByPathway$pathway),
                                      FUN = length)

# Remove unnecessary columns
enzymeFrequencyByPathway[3:6] <- list(NULL)

# Rename the columns
names(enzymeFrequencyByPathway)[names(enzymeFrequencyByPathway) == "Group.1"] <- "ec"
names(enzymeFrequencyByPathway)[names(enzymeFrequencyByPathway) == "Group.2"] <- "pathway"
names(enzymeFrequencyByPathway)[names(enzymeFrequencyByPathway) == "is_presented"] <- "frequency"

#---------------------------#
# [FREQUENCY BY BOTTLENECK] #
#---------------------------#

#-----------------------#
# [FREQUENCY BY SPECIE] #
#-----------------------#

