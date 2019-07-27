##################################################
# Example how to use the graph bottlenck package #
##################################################

# bottleneck ####

#' This is an example how to use the graph bottleneck package
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library(KEGGREST)
library(igraph)
library(RedeR)
library(RColorBrewer)

# Import all graph bottleneck library
files.sources = paste0("./R", "/", list.files(path = "./R"))
sapply(files.sources, source)

##############################################
# Define which pathway will be analysed
prefix <- "ec"
code <- "00260"
pathway <- paste0(prefix, code)

# Load EC from KO dictionnaire
ec_dictionnaire <-  get(load(paste0("./dictionnaires", "/", "KO2EC.RData")))

# Load the KO dictionnaire data
ko_dictionnaire <- get(load(paste0("./dictionnaires", "/", "KO", code, ".RData")))

# Load the pathways by organisms data
organism2pathway <- get(load(paste0("./dictionnaires", "/", "organism2pathway.RData")))

##############################################

# Load the KEGG reference pathway
# referencePathway <- getReferencePathway(code, ec_dictionnaire)
# iGraph <- graph_from_data_frame(referencePathway, directed=FALSE)

# Load the KEGG pathway and convert it into iGraph object
iGraph <- graph_from_data_frame(pathwayToDataframe(pathway))

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

# Print the graph in RedPort
printBottleneckInRedPort(iGraph, graphBottleneck)

# Export KEGG  pathway image
printBottleneckPathwayImage(pathway, ECToEntrez(graphBottleneck, ec_dictionnaire, ko_dictionnaire), TRUE)
