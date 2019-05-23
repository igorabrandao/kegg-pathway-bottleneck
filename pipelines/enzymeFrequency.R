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
enzymeList <- data.frame()
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

####################################
# Step 1: Get all pathways enzymes #
####################################

getPathwayEnzymes <- function(row, removeNoise=TRUE, replaceEmptyGraph=TRUE) {
  #################################
  # Get the organism general info #
  #################################

  # Get the current organism
  org <- organism2pathway[row]

  # Get its name
  specie <- names(org)

  # Status message
  printMessage(paste0("COUNTING ", specie, " ENZYMES FREQUENCIES [", row, " OF ", length(organism2pathway), "]"))

  # Loop over the current organism pathways code
  for(idx in start_of:length(unlist(org))) {
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
    if (is.null(temp)) {
      if (replaceEmptyGraph) {
        # Try to find the specific pathway for an organism
        pathway_code_tmp <- paste0(specie, pathway)

        # Status message
        printMessage(paste0("<<< Requesting specific pathway for ", pathway_code_tmp, "... >>>"))

        # Receive the specie data as Entrez
        temp <- pathwayToDataframe(pathway_code_tmp, TRUE, specie)

        # Set the specific flag
        extracted_by_entrez <- TRUE
      } else {
        # Save the not extracted pathway
        not_extracted_tmp<-data.frame(specie, pathway)
        names(not_extracted_tmp)<-c("org", "pathway")

        # Add each pathways not extracted
        pathwaysNotExtracted <<- rbind(pathwaysNotExtracted, not_extracted_tmp)
      }
    }

    # Check if the organism has the current pathway
    if (!is.null(temp)) {
      # Remove unnecessary data before bottleneck calculation
      if (removeNoise) {
        temp <- temp[!grepl("^path:", temp$node1),]
        temp <- temp[!grepl("^path:", temp$node2),]
        temp <- temp[!grepl("^map:", temp$node1),]
        temp <- temp[!grepl("^map:", temp$node2),]
      }

      # Calculates the network bottleneck
      iGraph <- igraph::graph_from_data_frame(temp, directed = FALSE)

      # Perform the graph bottleneck calculation
      graphBottleneck <- as_ids(getGraphBottleneck(iGraph, FALSE))

      # Convert the node2 column into node1 rows
      temp_df <- as.data.frame(temp[,c(2, 3, 4)], stringsAsFactors = FALSE)
      temp$node2 <- NULL
      names(temp_df)[names(temp_df) == "node2"] <- "node1"
      temp <- rbind(temp, temp_df)

      # Add a new column to the enzymeFrquency dataFrame
      temp$is_bottleneck <- 0

      # Assign the bottlenecks
      temp$is_bottleneck[which(temp[,1] %in% graphBottleneck)] <- 1

      # If necessary convert Entrez to EC
      if (extracted_by_entrez) {
        # Status message
        printMessage(paste0("<<< Converting Entrez to EC for pathway: ", pathway_code_tmp, "... >>>"))

        # Convert Entrez to EC
        temp$node1 <- convertEntrezToECWithoutDict(temp$node1, 50, TRUE)
      }

      # Add each pathways enzymes into enzymeList dataFrame
      enzymeList <<- rbind(enzymeList, temp)
    }
  }
}

sapply(start_of:length(organism2pathway), function(idx) getPathwayEnzymes(idx))

#-------------------------------------------------------------------------------------------#

####################################
# Step 2: Clear duplicates enzymes #
####################################

# Rename the nodes column
names(enzymeList)[names(enzymeList) == "node1"] <- "ec"

# Remove duplicated rows based on entrez column
enzymeList <- enzymeList[!duplicated(enzymeList[c("ec","pathway")]),]

# Reindex the enzymeList rows
rownames(enzymeList) <- 1:nrow(enzymeList)

#-------------------------------------------------------------------------------------------#

########################################################
# Step 3: Check if the enzyme appears into the pathway #
########################################################

current_pathway <- ""
highlighted_enzymes <- NULL

# Add a new column to the enzymeFrquency dataFrame
enzymeList$is_presented <- 0

# Order the dataFrame
enzymeList <- enzymeList[order(enzymeList$org, enzymeList$pathway),]

checkGenePresence <- function(row) {
  # Check if the current pathway is the same of the previous one
  if ( current_pathway != paste0(enzymeList$org[row], enzymeList$pathway[row]) ) {
    # Update the current pathway
    current_pathway <<- paste0(enzymeList$org[row], enzymeList$pathway[row])

    # Status message
    cat("\n")
    print(paste0("<<< Working on ", current_pathway, " pathway... >>>"))
    cat("\n")

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

sapply(start_of:length(unlist(enzymeList$ec)), function(idx) checkGenePresence(idx))

# Remove intermediary variables
rm(current_pathway, highlighted_enzymes)

#-------------------------------------------------------------------------------------------#

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

