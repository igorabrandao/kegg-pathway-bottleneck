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

# Load the KEGG pathway and convert it into iGraph object
iGraph <- graph_from_data_frame(pathwayToDataframe(pathway))

# Get the pathway bottleneck
graphBottleneck <- names(getGraphBottleneck(iGraph, TRUE))

# Get the gene list
genesList <- KEGGREST::keggLink(specie, pathway)
genesList <- gsub("^[[:alpha:]]*:(.*$)", "\\1", genesList)

# Get the highlighted gene list
highlightedGenesList <- getPathwayHighlightedGenes(pathway)

# Insert the pathway data into DB
pathway_data <- list(specie=specie, code=code, gene_count=as.numeric(length(genesList)))
last_insert_id <- insertPathway(pathway_data)

########################
# Insert the gene data #
########################

# Loop over each gene
for(gene in genesList) {
  # Default flags
  is_bottleneck = 0
  belong_to_specie = 0

  # Verify if the current gene is a bottleneck
  if (paste0(specie, ":", gene) %in% graphBottleneck) {
    is_bottleneck = 1
  }

  # Verify if the current gene belong to the current specie
  if (gene %in% highlightedGenesList) {
    belong_to_specie = 1
  }

  # Aggregate the data into a list
  gene_data <- list(pathway_id=last_insert_id, entrez=paste0(specie, ":", gene),
                    is_bottleneck=is_bottleneck, belong_to_specie=belong_to_specie)

  # Insert the gene data into DB
  insertGene(gene_data)
}
