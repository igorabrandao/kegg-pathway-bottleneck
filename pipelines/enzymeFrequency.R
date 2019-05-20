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
    pathway_code <- paste0(specie, pathway)

    # Status message
    cat("\n")
    print(paste0("<<< Requesting ", pathway_code, "... >>>"))
    cat("\n")

    # Get the enzyme list from pathway
    temp <- pathwayToDataframe(pathway_code)

    # Check if the organism has the current pathway
    if (is.not.null(temp)) {
      # Add each pathways enzymes into enzymeFrequency dataFrame
      enzymeFrequency <- rbind(enzymeFrequency, temp)
    }
  }
}

# Convert the node2 column into node1 rows
temp_df <- as.data.frame(enzymeFrequency[,c(2, 3, 4)], stringsAsFactors = FALSE)
enzymeFrequency$node2 <- NULL
names(temp_df)[names(temp_df) == "node2"] <- "node1"
enzymeFrequency <- rbind(enzymeFrequency, temp_df)

# Remove intermediary variables
rm(temp, temp_df, idx)

#-------------------------------------------------------------------------------------------#

####################################
# Step 2: Clear duplicates enzymes #
####################################

# Rename [node1 -> entrez] column name
names(enzymeFrequency)[names(enzymeFrequency) == "node1"] <- "entrez"

# Remove duplicated rows based on entrez column
enzymeFrequency <- enzymeFrequency[!duplicated(enzymeFrequency[c("entrez","pathway")]),]

# Reindex the enzymeFrequency rows
rownames(enzymeFrequency) <- 1:nrow(enzymeFrequency)

#-------------------------------------------------------------------------------------------#

######################################
# Step 3: Convert the entrez into EC #
######################################

# Empty EC list dataFrame
ec_list_df <- data.frame("ec" = character(0), stringsAsFactors = FALSE)

# Check if enzyme frequency dataFrame is not null
if (is.not.null(enzymeFrequency)) {
  # Call the function to convert the entrez code into EC
  ec_list_df <- convertEntrezToECWithoutDict(enzymeFrequency[,c(1)], 50, TRUE)
}

# Add the EC column into ezyme frequency dataFrame
enzymeFrequency = cbind(enzymeFrequency, ec_list_df)

# Remove intermediary variables
rm(ec_list_df)

#-------------------------------------------------------------------------------------------#

########################################################
# Step 4: Check if the enzyme appears into the pathway #
########################################################

current_pathway <- ""
highlighted_enzymes <- NULL

# Add a new column to the enzymeFrquency dataFrame
enzymeFrequency$belong_to_pathway <- 0

# Order the dataFrame
enzymeFrequency <- enzymeFrequency[order(enzymeFrequency$org, enzymeFrequency$pathway),]

# Loop over the enzymeFrequency dataFrame
for(row in start_of:length(unlist(enzymeFrequency$entrez))) {
  # Check if the current pathway is the same of the previous one
  if ( current_pathway != paste0(enzymeFrequency$org[row], enzymeFrequency$pathway[row]) ) {
    # Update the current pathway
    current_pathway <- paste0(enzymeFrequency$org[row], enzymeFrequency$pathway[row])

    # Get the highlighted enzymes list
    highlighted_enzymes <- getPathwayHighlightedGenes(current_pathway, allMapped_=TRUE)
  }

  # Get just the enzyme number without specie
  current_enzyme <- gsub("^[[:alpha:]]*(.*$)", "\\1", str_replace(enzymeFrequency$entrez[row], ":", ""))

  # Verify if the current enzyme is highlighted and set its status
  if (current_enzyme %in% highlighted_enzymes) {
    enzymeFrequency$belong_to_pathway[row] <- 1
  }
}

# Remove intermediary variables
rm(current_pathway, highlighted_enzymes, row)

#-------------------------------------------------------------------------------------------#

###############################################
# Step 5: Check if the enzyme is a bottleneck #
###############################################

#Para calcular o bottleneck, pegar o objeto iGraph
#e aplicar a função de calculo de bottleneckss diretamente
#sem necessidade de filtrar os nós pintados para cada organismo

# TODO: Test code (remove later)
xpto <- getPathwayHighlightedGenes("hsa00010", allMapped_=TRUE)
xpto2 <- getPathwayHighlightedGenes("hsa00010", allMapped_=FALSE)

iGraph <- pathwayToDataframe("hsa00010")
iGraph <- igraph::graph_from_data_frame(iGraph, directed = FALSE)

# Vertex communites
iGraph <- setGraphCommunity(iGraph)

# Vertex betweenness
iGraph <- setGraphBetwenness(iGraph)

# Vertex closeness
iGraph <- setGraphCloseness(iGraph)

# Vertex clustering
iGraph <- setGraphClustering(iGraph)

# Perform the graph bottleneck calculation
graphBottleneck <- getGraphBottleneck(iGraph, TRUE)




current_pathway <- ""
pathway_bottlenecks <- NULL

# Add a new column to the enzymeFrquency dataFrame
enzymeFrequency$is_bottleneck <- 0

# Loop over the enzymeFrequency dataFrame
for(row in start_of:length(unlist(enzymeFrequency$entrez))) {
  # Check if the current pathway is the same of the previous one
  if ( current_pathway != paste0(enzymeFrequency$org[row], enzymeFrequency$pathway[row]) ) {
    # Update the current pathway
    current_pathway <- paste0(enzymeFrequency$org[row], enzymeFrequency$pathway[row])

    # Get the pathway bottlenecks list
    #pathway_bottlenecks <- getPathwayHighlightedGenes(current_pathway)
  }
}

# Remove intermediary variables
rm(current_pathway, highlighted_enzymes, row)

#-------------------------------------------------------------------------------------------#

############################################
# Step 6: Count the enzyme total frequency #
############################################

# Count the enzymes frequencies and transform it into a dataFrame
enzymeTotalFrequency <- as.data.frame(table(enzymeFrequency[,c(5)]), stringsAsFactors = FALSE)
