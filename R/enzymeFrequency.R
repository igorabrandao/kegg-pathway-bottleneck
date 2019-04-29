##################################################
# Script to perform the enzyme frequencies count #
##################################################

# enzymeFrequency ####

#' This is the pipeline script to perform
#' the enzymes frequencies counting
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library(KEGGREST)
library(igraph)

##############################################

# Import the graphLoader functions
files.sources = paste0("./R", "/", "graphLoader.R")
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

##############################################

# Loop over species matrix
# for(row in start_of:length(organism2pathway)) {
for(row in start_of:1) {
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
    print(paste0("<<< ", pathway_code, " >>>"))
    cat("\n")

    # Add each pathways enzymes into enzymeFrequency dataFrame
    enzymeFrequency <- rbind(enzymeFrequency, pathwayToDataframe(pathway_code))
  }
}
