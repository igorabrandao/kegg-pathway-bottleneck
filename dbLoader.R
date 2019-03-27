###################################################
# Script to insert the pathway data automatically #
# into DB                                         #
###################################################

# dbLoader ####

#' This is an automation script to insert KEGG data into mysql DB
#'
#' @author
#' Igor Brandão

# Import the necessary libraries

# Import all graph bottleneck library
files.sources = paste0("./R", "/", list.files(path = "./R"))
sapply(files.sources, source)

##############################################
# Define which pathway will be extracted
code <- "00010"

##############################################

# First of all, load all species
data(korg, package = "pathview", verbose = FALSE)

# Loop over species matrix
for(row in 1:nrow(korg)) {
  # TODO: Run the script here
}

###########################
# Insert the pathway data #
###########################

# Used for test purpose (needs to run over the loop)
specie <- korg[1, "kegg.code"]
pathway <- paste0(specie, code)

# Get the gene list
genesList <- KEGGREST::keggLink(specie, pathway)
genesList <- gsub("^[[:alpha:]]*:(.*$)", "\\1", genesList)

# Insert the pathway data into DB
pathway_data <- list(specie=specie, code=code, gene_count=as.numeric(length(genesList)))
last_insert_id <- insertPathway(pathway_data)

########################
# Insert the gene data #
########################
