########################################################
### LIBRARY IMPORTING
########################################################

# Import KEGGgraph library
library(KEGGgraph)

# Import KEGG.db library
library(KEGG.db)

# Import Rgraphviz library
library(Rgraphviz)

# Import RBGL library
library(RBGL)

########################################################
### DATA LOADING
########################################################

tmp <- tempfile()

# Note: retrieveKGML uses a try-download mechanism (since the KEGGgraph version
# 1.1.2 ) to retrieve the KGML file from remote KEGG FTP server
pName <- "p53 signaling pathway"
pId <- mget(pName, KEGGPATHNAME2ID)[[1]]
retrieveKGML(pId, organism="cel", destfile=tmp, method="wget", quiet=TRUE)

# First we read in KGML file for human MAPK signaling pathway (with KEGG ID hsa04010)
mapkKGML <- system.file("extdata/hsa04010.xml", package="KEGGgraph")

########################################################
### CONVERSION FROM KEGG PATHWAY INTO A GRAPH
########################################################

# Once the file is ready, we parse the KGML file into an object of KEGGpathway, which can
# be later converted into the graph object
mapkpathway <- parseKGML(mapkKGML)
mapkpathway

# Convert the the pathway into a graph object
mapkG <- KEGGpathway2Graph(mapkpathway, expandGenes=TRUE)
mapkG

########################################################
### DEGREE ATTRIBUTES QUERYING
########################################################

# Calculates the edges in and outs
mapkGoutdegrees <- sapply(edges(mapkG), length)
mapkGindegrees <- sapply(inEdges(mapkG), length)
topouts <- sort(mapkGoutdegrees, decreasing=T)
topins <- sort(mapkGindegrees, decreasing=T)

# Display the top ins & outs
topouts[1:10]
topins[1:10]

########################################################
### GRAPH PLOTTING
########################################################

# Get the in out from vertices
outs <- sapply(edges(mapkG), length) > 0
ins <- sapply(inEdges(mapkG), length) > 0
ios <- outs | ins

# translate the KEGG IDs into Gene Symbol
if(require(org.Hs.eg.db)) {
    ioGeneID <- translateKEGGID2GeneID(names(ios))
    nodesNames <- sapply(mget(ioGeneID, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
} else {
    nodesNames <- names(ios)
}

# Set the vertices name
names(nodesNames) <- names(ios)

# Set the names as labels
nAttrs <- list();
nAttrs$fillcolor <- makeAttr(mapkG, "lightgrey", list(orange=names(ios)[ios]))
nAttrs$label <- nodesNames

# Plot the graph
plot(
    mapkG, "neato", nodeAttrs=nAttrs,
    attrs = list(node=list(fillcolor="lightgreen",
    width = "0.75", shape="ellipse"),
    edge = list(arrowsize="0.7")),
    main = "HSA00010 - Glycolysis / Gluconeogenesis"
)

########################################################
### GRAPH PLOTTING WITH BETWEENNESS
########################################################

# Calculates the betweeness
bcc <- brandes.betweenness.centrality(mapkG)
rbccs <- bcc$relative.betweenness.centrality.vertices[1L,]
toprbccs <- sort(rbccs,decreasing=TRUE)[1:10]

toprbccName <- names(toprbccs)
toprin <- sapply(toprbccName, function(x) inEdges(mapkG)[x])
toprout <- sapply(toprbccName, function(x) edges(mapkG)[x])
toprSubnodes <- unique(unname(c(unlist(toprin), unlist(toprout), toprbccName)))
toprSub <- subGraph(toprSubnodes, mapkG)

nAttrs <- list()
tops <- toprbccs
topLabels <- lapply(toprbccName, function(x) x); names(topLabels) <- tops
nAttrs$label <- makeAttr(toprSub, "", topLabels)
nAttrs$fillcolor <- makeAttr(toprSub, "lightblue", list(orange=toprbccName))
nAttrs$width <- makeAttr(toprSub,"",list("0.8"=toprbccName))

plot(toprSub, "twopi", nodeAttrs=nAttrs, attrs=list(graph=list(start=2)))
