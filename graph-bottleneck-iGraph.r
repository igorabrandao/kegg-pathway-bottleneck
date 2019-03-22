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

# First we read in KGML file for human MAPK signaling pathway (with KEGG ID hsa00010)
mapkKGML <- system.file("extdata/hsa00010.xml", package="KEGGgraph")

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

#################################
# Vertex closeness #
#################################

# closeness refers to how connected a node is to its neighbors
V(iGraph)$closeness <- closeness(iGraph, vids=V(iGraph))

#################################
# Vertex clustering coefficient #
#################################

# calculates the local clustering coefficient for each vertex
V(iGraph)$clustering <- transitivity(iGraph, type = c("local"), vids = NULL,
                                      weights = NULL, isolates = c("NaN", "zero"))

# calculates the global clustering coefficient
transitivity(iGraph, type = c("global"), vids = NULL,
                              weights = NULL, isolates = c("NaN", "zero"))

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
V(iGraph)$nodeAlias <- paste0(names(V(iGraph)), " | B: ", V(iGraph)$betweenness
                              , " | Clu: ", V(iGraph)$clustering, " | Clo: ", V(iGraph)$closeness)

V(iGraph)$nodeAlias <- names(V(iGraph))

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

####################
# Grapg Bottleneck #
####################

# *** APPROACH 01 ***

# Graph minimum cut [bottleneck]
min_cut(iGraph, value.only=FALSE)

# A subgraph of a connected graph is a minimum spanning tree if it is tree, 
# and the sum of its edge weights are the minimal among all tree subgraphs of the graph. 
# A minimum spanning forest of a graph is the graph consisting of the minimum spanning 
# trees of its components.
mst(iGraph)

# --------------------------------------------------------------------------------------

# *** APPROACH 02 ***
# Definition of hubs and bottlenecks.
# We defined hubs as all proteins that are in the top 20% of the degree distribution 
# (i.e., proteins that have the 20% highest number of neighbors). 

# Accordingly, we defined bottlenecks as the proteins that are in the top 20% in terms of 
# betweenness. Varying this cutoff from 10% to 40% had no significant impact on results
betweenness_percentual_rate <- 0.2

topBetweenness <- sort(V(iGraph)$betweenness, decreasing=TRUE)
top20Betweenness <- topBetweenness[1:as.integer(length(topBetweenness) * betweenness_percentual_rate)]

# Print the top 20%
V(iGraph)[match(top20Betweenness, V(iGraph)$betweenness)]
top20Betweenness

# --------------------------------------------------------------------------------------

# *** APPROACH 03 ***

# Articulation points [cut vertices]

# They are vertices whose removal increases the number of connected 
# components in a graph
articulation_points(iGraph)

###############
# Data export #
###############

# Save the iGraph object
save(iGraph, file = "iGraph.RData", compress = "xz")

# Load the iGraph object
load("iGraph.RData")

# --------------------------------------------------------------------------------------

# XLSX export
library("xlsx")

# Get th data frame
iGraph_df <- as.data.frame(list(Vertex=V(iGraph)$name, Betweenness=V(iGraph)$betweenness, 
              Clustering=V(iGraph)$clustering, Closeness=V(iGraph)$closeness, 
              Group=V(iGraph)$group), stringsAsFactors=FALSE)

# Write the data set in a workbook
write.xlsx(iGraph_df, file = "igraph.xlsx",
           sheetName = "IGRAPH-DATA", append = FALSE)
