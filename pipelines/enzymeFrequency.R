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
library(igraph)
library(RCurl) # http connections
library(rvest) # web scraping
library(stringr) # regex manipulation
library(pracma) # string manipulation


library(progress)
library(parallel)
library(snow)
library(doSNOW)
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

# Define in which specie the processing should begin
# default value 1 (the value should be >= 1)
start_of <- 1

# Empty enzyme frequency dataFrame
enzymeList <- NULL
pathwaysNotExtracted <- data.frame(org = character(0), pathway = character(0))

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

getPathwayEnzymes <- function(row, removeNoise=TRUE, replaceEmptyGraph=TRUE) {

  #################################
  # Get the organism general info #
  #################################

  # Get the current organism
  org <- organism2pathway[row]

  # Get its name
  specie <- names(org)

  ###########################
  # Set parallelism options #
  ###########################

  nCores <- parallel::detectCores()
  cl <- snow::makeSOCKcluster(nCores)
  on.exit(snow::stopCluster(cl))
  doSNOW::registerDoSNOW(cl)

  progress <- function() {
    pb$tick()
  }

  opts <- list(progress = progress)
  ntasks <- length(unlist(org))

  pb <- progress::progress_bar$new(format = "running [:bar] :percent elapsed time :elapsed",
                                   total = ntasks, clear = FALSE,
                                   width = 60)

  # Status message
  printMessage(paste0("COUNTING ", specie, " ENZYMES FREQUENCIES [", row, " OF ", length(organism2pathway), "]"))

  # Loop over the current organism pathways code
  result <- foreach::foreach(idx = seq.int(1, ntasks), .export=c('printMessage', 'pathwayToDataframe', 'as_ids',
                                                                     'getGraphBottleneck', 'convertEntrezToECWithoutDict',
                                                                     'pathwaysNotExtracted'),
                                 .combine = "rbind", .options.snow = opts) %dopar%
  {
    #####################################
    # Load all enzymes from its pathway #
    #####################################

    # Set the flag extraction by Entrez
    extracted_by_entrez <- FALSE

    # Get the pathway object
    pathway <- unlist(org)[idx]

    # Format the pathway code
    pathway_code <- paste0('ec', pathway)

    # Status message
    printMessage(paste0("<<< Requesting ", specie, ": ", pathway_code, "... >>>"))

    # Get the enzyme list from pathway
    temp <- pathwayToDataframe(pathway_code, TRUE, specie)

    # Handle empty graph
    if (replaceEmptyGraph) {
      if (is.null(temp)) {
        # Try to find the specific pathway for an organism
        pathway_code_tmp <- paste0(specie, pathway)

        # Status message
        printMessage(paste0("<<< Requesting specific pathway for ", pathway_code_tmp, "... >>>"))

        # Receive the specie data as Entrez
        temp <- pathwayToDataframe(pathway_code_tmp, TRUE, specie)

        # Set the specific flag
        extracted_by_entrez <- TRUE
      }
    } else {
      # Save the not extracted pathway
      not_extracted_tmp<-data.frame(specie, pathway)
      names(not_extracted_tmp)<-c("org", "pathway")

      # Add each pathways not extracted
      pathwaysNotExtracted <<- rbind(pathwaysNotExtracted, not_extracted_tmp)
    }

    # Check if the organism has the current pathway
    if (!is.null(temp)) {
      # Remove unnecessary data before bottleneck calculation
      if (removeNoise) {
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

      # Convert the node2 column into node1 rows
      aux <- unique(c(temp$node1, temp$node2))
      auxorg <- temp$org[1]
      auxpathway <- temp$pathway[1]

      # Add a new column to the enzymeFrquency dataFrame
      temp <- data.frame(node1 = aux, org = auxorg, pathway = auxpathway, is_bottleneck = 0,
                         is_presented = 0, stringsAsFactors = FALSE)

      # Assign the bottlenecks
      temp$is_bottleneck[which(temp[,1] %in% graphBottleneck)] <- 1

      # If necessary convert Entrez to EC
      if (extracted_by_entrez) {
        # Status message
        printMessage(paste0("<<< Converting Entrez to EC for pathway: ", pathway_code_tmp, "... >>>"))

        # Remove duplicated rows based on entrez column
        temp <- temp[!(duplicated(temp$node1) & duplicated(temp$pathway)),]

        # Convert Entrez to EC
        temp$node1 <- convertEntrezToECWithoutDict(temp$node1, 50, TRUE)
      }

      # Add each pathways enzymes into enzymeList dataFrame
      # enzymeList <<- rbind(enzymeList, temp)
    }

    return(temp)
  }

  return(result)
}

checkGenePresence <- function(row) {
  # Check if the current pathway is the same of the previous one
  if ( current_pathway != paste0(enzymeList$org[row], enzymeList$pathway[row]) ) {
    # Update the current pathway
    current_pathway <<- paste0(enzymeList$org[row], enzymeList$pathway[row])

    # Status message
    printMessage(paste0("<<< Working on ", current_pathway, " pathway... >>>"))

    # Get the highlighted enzymes list
    highlighted_enzymes <<- getPathwayHighlightedGenes(current_pathway, allMapped_=TRUE)

    # Concat the org string
    highlighted_enzymes <<- paste(enzymeList$org[row], highlighted_enzymes, sep=":")

    # Convert it into dataframe
    highlighted_enzymes <<- as.data.frame(highlighted_enzymes, stringsAsFactors = FALSE)

    # Convert the highlighted list into EC number
    highlighted_enzymes <<- convertEntrezToECWithoutDict(highlighted_enzymes[,c(1)], 50, TRUE)

    # Remove the duplicates
    highlighted_enzymes <<- highlighted_enzymes[!duplicated(highlighted_enzymes),]
  }

  # Get just the enzyme number without specie
  current_enzyme <- gsub("^[[:alpha:]]*(.*$)", "\\1", str_replace(enzymeList$ec[row], ":", ""))

  # Verify if the current enzyme is highlighted and set its status
  if (current_enzyme %in% highlighted_enzymes) {
    enzymeList$is_presented[row] <<- 1
  }
}

#-------------------------------------------------------------------------------------------#

####################################
# Step 1: Get all pathways enzymes #
####################################

enzymeList <- lapply(start_of:2, getPathwayEnzymes, replaceEmptyGraph=FALSE)
enzymeList <- do.call("rbind", enzymeList)
#enzymeList <- sapply(start_of:length(organism2pathway), function(idx) getPathwayEnzymes(idx, replaceEmptyGraph=FALSE))

####################################
# Step 2: Clear duplicates enzymes #
####################################

# Rename the nodes column
names(enzymeList)[names(enzymeList) == "node1"] <- "ec"

# Remove duplicated rows based on entrez column
enzymeList <- enzymeList[!duplicated(enzymeList[c("ec","pathway")]),]

# Reindex the enzymeList rows
rownames(enzymeList) <- 1:nrow(enzymeList)

########################################################
# Step 3: Check if the enzyme appears into the pathway #
########################################################

current_pathway <- ""
highlighted_enzymes <- NULL

# Add a new column to the enzymeFrquency dataFrame
enzymeList$is_presented <- 0

# Order the dataFrame
enzymeList <- enzymeList[order(enzymeList$org, enzymeList$pathway),]

sapply(start_of:length(unlist(enzymeList$ec)), function(idx) checkGenePresence(idx))

# Remove intermediary variables
rm(current_pathway, highlighted_enzymes)

############################################
# Step 4: Count the enzyme total frequency #
############################################

# Filter just the enzymes with some frequency
enzymeTotalFrequency <- enzymeList[enzymeList$is_presented == 1,]

# Count the enzymes frequencies and transform it into a dataFrame
enzymeTotalFrequency <- as.data.frame(table(enzymeTotalFrequency$ec), stringsAsFactors = FALSE)
names(enzymeTotalFrequency)[names(enzymeTotalFrequency) == "Var1"] <- "ec"

# TODO: Calcular a frequência dos genes presentes por rota

# Count the enzyme frequency by pathway code
enzymeFrequencyByPathway

# APAGAR DEPOIS!
# Count the proteins by its orthologous_group and specie
result <- setNames(aggregate(enzymeList[,c('ec')],
                             by=list(enzymeList$pathway), length),
                   c('ec', 'pathway', 'frequency'))

aggregate(x = enzymeList,
          by = list(enzymeList$ec, enzymeList$pathway),
          FUN = length)


#enzymeList[which(enzymeList$ec %in% enzymeTotalFrequency$ec),]
