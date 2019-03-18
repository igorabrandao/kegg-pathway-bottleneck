# Author:
# Kegg Pathway Bottlenecks
# 14 March 2019

# LOADING PACKAGES
library(KEGGgraph)
library(KEGG.db)
library(Rgraphviz)
library(igraph)
library(RedeR)
library(RColorBrewer)

##########################
# Load Kegg Pathway Data #
##########################

tmp <- tempfile()

# Note: retrieveKGML uses a try-download mechanism (since the KEGGgraph version
# 1.1.2 ) to retrieve the KGML file from remote KEGG FTP server
pName <- "p53 signaling pathway"
pId <- mget(pName, KEGGPATHNAME2ID)[[1]]
retrieveKGML(pId, organism="cel", destfile=tmp, method="wget", quiet=TRUE)

# First we read in KGML file for human MAPK signaling pathway (with KEGG ID hsa04010)
mapkKGML <- system.file("extdata/hsa04010.xml", package="KEGGgraph")

####################################################
# Conversion from Kegg Pathway into a graph object #
####################################################

# Once the file is ready, we parse the KGML file into an object of KEGGpathway, which can
# be later converted into the graph object
mapkpathway <- parseKGML(mapkKGML)
mapkpathway

# Convert the the pathway into a graph object
mapkG <- KEGGpathway2Graph(mapkpathway, expandGenes=TRUE)
mapkG

########################################
# Convert the graph object into iGraph #
########################################

# Get the data attribute from graph
aux <- names(mapkG@edgeData@data)

# Converts it into a dataFrame
aux <- as.data.frame(aux, stringsAsFactors = FALSE)

# Split the node1->node2 into 2 columns
aux$node2 <- gsub("^.*\\|(.*)$", "\\1", aux$aux)
colnames(aux)[1] <- "node1"
aux$node1 <- gsub("^(.*)\\|.*$", "\\1", aux$node1)

# iGraph object transformation
iGraph <- graph_from_data_frame(aux, directed = TRUE)

#######################
# Community detection #
#######################

# greedy method (hiearchical, fast method)
c3 = cluster_edge_betweenness(iGraph)

# define the cluster attribute
V(iGraph)$group <- membership(c3)

# modularity measure [optional]
modularity(c3)

###############################
# Edge betwenness calculation #
###############################

# calculates the betweenness
V(iGraph)$betweenness <- betweenness(iGraph, normalized = TRUE)

# normalizes betweenness
V(iGraph)$betweenness <- V(iGraph)$betweenness/max(V(iGraph)$betweenness)

################################
# Display the graph in RedPort #
################################

# Set the color palette according to the communities count
myColors <- rainbow(length(c3))

# Set the graph attributes
V(iGraph)$nodeLineColor <- myColors[V(iGraph)$group]
E(iGraph)$edgeColor <- "grey80"

pal <- brewer.pal(9, "YlOrRd")
color_col <- colorRampPalette(pal)(10)

iGraph <- att.setv(g = iGraph, from = "betweenness", to = "nodeColor",
                   cols = color_col, na.col = "grey80", breaks = seq(0, 1, 0.1))

# Set the betweenness into label name [warning: bad visualization]
V(iGraph)$nodeAlias <- paste0(names(V(iGraph)), " | ", V(iGraph)$betweenness)

# Set the graph direction
E(iGraph)$arrowDirection <- 1
V(iGraph)$nodeLineWidth <- 5

# Create RedPort object
rdp <- RedPort()

# Open the connection
calld(rdp)

# Add the graph into the RedPort
addGraph(rdp, iGraph)

# Add the legend with color scale
addLegend.color(rdp, colvec=iGraph$legNodeColor$scale, size=15, labvec=iGraph$legNodeColor$legend, 
                title="Betweenness Centrality Scale (BCS)")

# Relax the graph visualization
relax(rdp)

# Clear the RedPort [optional]
resetd(rdp)

##################################
# Graph minimum cut [bottleneck] #
##################################

min_cut(iGraph, value.only=FALSE)

###########################
# General info [optional] #
###########################

# Top 10 most betweenness vertices
toprBetweenness <- sort(V(iGraph)$betweenness, decreasing=TRUE)[1:10]
toprBetweenness

###############
# Data export #
###############

# Save the iGraph object
save(iGraph, file = "iGraph.RData", compress = "xz")

# Load the iGraph object
load("iGraph.RData")