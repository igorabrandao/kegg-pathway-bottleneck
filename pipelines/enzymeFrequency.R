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

getPathwayEnzymes <- function(row, removeNoise_=TRUE, replaceEmptyGraph_=TRUE) {

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
  result <- foreach::foreach(idx = seq.int(1, ntasks), .export=c('printMessage', 'pathwayToDataframe', 'as_ids', 'str_replace',
                                                                     'getGraphBottleneck', 'convertEntrezToECWithoutDict',
                                                                     'getPathwayHighlightedGenes'),
                                 .combine = "rbind", .options.snow = opts) %dopar%
  {
    #sink("log.txt", append=TRUE)

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
    if (is.null(temp) | length(temp) == 0) {
      if (replaceEmptyGraph_) {
        # Try to find the specific pathway for an organism
        pathway_code_tmp <- paste0(specie, pathway)

        # Status message
        printMessage(paste0("<<< Requesting specific pathway for ", pathway_code_tmp, "... >>>"))

        # Receive the specie data as Entrez
        temp <- pathwayToDataframe(pathway_code_tmp, TRUE, specie)

        # Set the specific flag
        extracted_by_entrez <- TRUE
      }
      else {
        temp <- data.frame(node1 = NA, org = specie, pathway = pathway, is_bottleneck = 0,
                           is_presented = 0, stringsAsFactors = FALSE)

        return(temp)
      }
    }

    # Check if the organism has the current pathway
    if (!is.null(temp) & length(temp) != 0) {
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

      # Convert the node2 column into node1 rows
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

      # If necessary convert Entrez to EC
      if (extracted_by_entrez) {
        # Remove the duplicates
        highlighted_enzymes <- highlighted_enzymes[!duplicated(highlighted_enzymes),]

        # Get just the enzyme number without specie
        current_enzyme <- gsub("^[[:alpha:]]*(.*$)", "\\1", str_replace(temp$node1, ":", ""))

        # Verify if the current enzyme is highlighted and set its status
        if (is.null(highlighted_enzymes) | length(highlighted_enzymes) == 0) {
          temp$is_presented[which(current_enzyme %in% highlighted_enzymes)] <- 1
        }

        # Status message
        printMessage(paste0("<<< Converting Entrez to EC for pathway: ", pathway_code_tmp, "... >>>"))

        # Remove duplicated rows based on entrez column
        temp <- temp[!(duplicated(temp$node1) & duplicated(temp$pathway)),]

        # Convert Entrez to EC
        temp$node1 <- convertEntrezToECWithoutDict(temp$node1, 50, TRUE)
      } else {
        # Convert the highlighted list into EC number
        highlighted_enzymes <- convertEntrezToECWithoutDict(highlighted_enzymes[,c(1)], 50, TRUE)

        # Remove the duplicates
        highlighted_enzymes <- highlighted_enzymes[!duplicated(highlighted_enzymes),]

        # Get just the enzyme number without specie
        current_enzyme <- gsub("^[[:alpha:]]*(.*$)", "\\1", str_replace(temp$node1, ":", ""))

        # Verify if the current enzyme is highlighted and set its status
        if (is.null(highlighted_enzymes) | length(highlighted_enzymes) == 0) {
          temp$is_presented[which(current_enzyme %in% highlighted_enzymes)] <- 1
        }
      }


      printMessage(paste0("is_present indexes: ", str(which(current_enzyme %in% highlighted_enzymes))))
    }

    # Return the specie [FOREACH]
    return(temp)
  }

  doSNOW::stopCluster(c1)

  # Return the entire dataSet [FUNCTION]
  return(result)
}

#-------------------------------------------------------------------------------------------#

####################################
# Step 1: Get all pathways enzymes #
####################################

enzymeList <- lapply(start_of:2, getPathwayEnzymes, replaceEmptyGraph_=FALSE)
enzymeList <- do.call("rbind", enzymeList)
#enzymeList <- sapply(start_of:length(organism2pathway), function(idx) getPathwayEnzymes(idx, replaceEmptyGraph_=FALSE))

####################################
# Step 2: Clear duplicates enzymes #
####################################

# Rename the nodes column
names(enzymeList)[names(enzymeList) == "node1"] <- "ec"

# Remove duplicated rows based on entrez column
enzymeList <- enzymeList[!duplicated(enzymeList[c("ec", "org", "pathway")]),]

# Reindex the enzymeList rows
rownames(enzymeList) <- 1:nrow(enzymeList)

############################################
# Step 3: Count the enzyme total frequency #
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


enzymeList[which(enzymeList$ec %in% enzymeTotalFrequency$ec),]
