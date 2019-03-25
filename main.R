##################################################
# Example how to use the graph bottlenck package #
##################################################

# main ####

#' This is an example how to use the graph bottleneck package
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library(igraph)
library(RedeR)
library(RColorBrewer)

# Import all graph bottleneck library
files.sources = paste0("./R", "/", list.files(path = "./R"))
sapply(files.sources, source)

##############################################
# Define which pathway will be analysed
prefix <- "hsa"
code <- "00010"
pathway <- paste0(prefix, code)

##############################################

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
printBottleneckPathwayImage(pathway, graphBottleneck, TRUE)
