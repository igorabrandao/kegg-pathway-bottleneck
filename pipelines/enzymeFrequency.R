####################################################
# Pipeline to perform the enzyme frequencies count #
####################################################

# enzymeFrequency.R #

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

# Function to detect null values
is.not.null <- function(x) !is.null(x)

# Empty enzyme frequency dataFrame
enzymeFrequency <- data.frame()

#-------------------------------------------------------------------------------------------#

####################################
# Step 1: Get all pathways enzymes #
####################################

# Loop over species matrix
for(row in start_of:length(organism2pathway)) {
# for(row in start_of:1) {
  #################################
  # Get the organism general info #
  #################################

  # Get the current organism
  org <- organism2pathway[row]

  # Get its name
  specie <- names(org)

  # Status message
  cat("\n")
  print("------------------------------------------------")
  print(paste0("COUNTING ", specie, " ENZYMES FREQUENCIES [", row, " OF ", length(organism2pathway), "]"))
  print("------------------------------------------------")
  cat("\n")

  # Loop over the current organism pathways code
  for(idx in start_of:length(unlist(org))) {
    #####################################
    # Load all enzymes from its pathway #
    #####################################

    # Get the pathway object
    pathway <- unlist(org)[idx]

    # Format the pathway code
    pathway_code <- paste0('ec', pathway)

    # Status message
    cat("\n")
    print(paste0("<<< Requesting ", pathway_code, "... >>>"))
    cat("\n")

    # Get the enzyme list from pathway
    temp <- pathwayToDataframe(pathway_code, TRUE, specie)

    # Check if the organism has the current pathway
    if (is.not.null(temp)) {
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

      # Add each pathways enzymes into enzymeFrequency dataFrame
      enzymeFrequency <- rbind(enzymeFrequency, temp)
    }
  }
}

# Remove intermediary variables
rm(temp, temp_df, idx)

#-------------------------------------------------------------------------------------------#

####################################
# Step 2: Clear duplicates enzymes #
####################################

# Rename the nodes column
names(enzymeFrequency)[names(enzymeFrequency) == "node1"] <- "ec"

# Remove duplicated rows based on entrez column
enzymeFrequency <- enzymeFrequency[!duplicated(enzymeFrequency[c("ec","pathway")]),]

# Reindex the enzymeFrequency rows
rownames(enzymeFrequency) <- 1:nrow(enzymeFrequency)

#-------------------------------------------------------------------------------------------#

########################################################
# Step 3: Check if the enzyme appears into the pathway #
########################################################

current_pathway <- ""
highlighted_enzymes <- NULL

# Add a new column to the enzymeFrquency dataFrame
enzymeFrequency$belong_to_pathway <- 0

# Order the dataFrame
enzymeFrequency <- enzymeFrequency[order(enzymeFrequency$org, enzymeFrequency$pathway),]

# Loop over the enzymeFrequency dataFrame
for(row in start_of:length(unlist(enzymeFrequency$ec))) {

  # Check if the current pathway is the same of the previous one
  if ( current_pathway != paste0(enzymeFrequency$org[row], enzymeFrequency$pathway[row]) ) {
    # Update the current pathway
    current_pathway <- paste0(enzymeFrequency$org[row], enzymeFrequency$pathway[row])

    # Get the highlighted enzymes list
    highlighted_enzymes <- getPathwayHighlightedGenes(current_pathway, allMapped_=TRUE)

    # Concat the org string
    highlighted_enzymes <- paste(enzymeFrequency$org[row], highlighted_enzymes, sep=":")

    # Convert it into dataframe
    highlighted_enzymes <- as.data.frame(highlighted_enzymes, stringsAsFactors = FALSE)

    # Convert the highlighted list into EC number
    highlighted_enzymes <- convertEntrezToECWithoutDict(highlighted_enzymes[,c(1)], 50, TRUE)

    # Remove the duplicates
    highlighted_enzymes <- highlighted_enzymes[!duplicated(highlighted_enzymes),]
  }

  # Get just the enzyme number without specie
  current_enzyme <- gsub("^[[:alpha:]]*(.*$)", "\\1", str_replace(enzymeFrequency$ec[row], ":", ""))

  # Verify if the current enzyme is highlighted and set its status
  if (current_enzyme %in% highlighted_enzymes) {
    enzymeFrequency$belong_to_pathway[row] <- 1
  }
}

# Remove intermediary variables
rm(current_pathway, highlighted_enzymes, row)

#-------------------------------------------------------------------------------------------#

############################################
# Step 4: Count the enzyme total frequency #
############################################

# Count the enzymes frequencies and transform it into a dataFrame
enzymeTotalFrequency <- as.data.frame(table(enzymeFrequency[,c(5)]), stringsAsFactors = FALSE)
