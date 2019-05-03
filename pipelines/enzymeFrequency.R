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
library(KEGGREST)
library(igraph)
library(RCurl)
library(rvest)
library(stringr)

#-------------------------------------------------------------------------------------------#

###########################
# Pipeline basic settings #
###########################

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

# Remove the temp variable
rm(temp)

#-------------------------------------------------------------------------------------------#

######################################
# Step 2: Convert the entrez into EC #
######################################

# Empty EC list dataFrame
ec_list_df <- data.frame("EC" = character(0), stringsAsFactors = FALSE)

# Check if enzyme frequency dataFrame is not null
if (is.not.null(enzymeFrequency)) {
  # Call the function to convert the entrez code into EC
  ec_list_df <- convertEntrezToECWithoutDict(enzymeFrequency[,c(1)], 100, TRUE)
}

# Add the EC column into ezyme frequency dataFrame
ec_list_df[(nrow(ec_list_df) + 1):nrow(enzymeFrequency),] = NA # it's temporaly until fix the conversion function
enzymeFrequency = cbind(enzymeFrequency, ec_list_df)

#-------------------------------------------------------------------------------------------#

############################################
# Step 3: Count the enzyme total frequency #
############################################

# Count the enzymes frequencies and transform it into a dataFrame
enzymeTotalFrequency <- as.data.frame(table(enzymeFrequency[,c(5)]), stringsAsFactors = FALSE)
