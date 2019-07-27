#*****************************************#
# Functions to handle graph visualization #
#*****************************************#

# graphPrintFunctions.R #

# ---- IMPORT SECTION ----

#' This is set of functions to handle graph visualization
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library(igraph)
library(RedeR)
library(RColorBrewer)
library(visNetwork)
library(scales)

#*******************************************************************************************#

# ---- GRAPH PRINT SECTION ----

#' Function to generate interactive networks based on redPort library
#'
#' @param iGraph_ iGraph object.
#' @param bottleneck_ List containing the bottlenecks to be colored.
#' @param verbose_ Flag to display or not the processing messages.
#'
#' @return This function does not return nothing, just open a redPort
#' window instance..
#'
#' @author
#' Igor Brandão

printBottleneckInRedPort <- function(iGraph_, bottleneck_, verbose_=FALSE) {
  # Set the color palette according to the communities count
  myColors <- rainbow(length(unique(V(iGraph_)$group)))

  # Set the graph attributes
  V(iGraph_)$nodeLineColor <- myColors[V(iGraph_)$group]
  E(iGraph_)$edgeColor <- "grey80"

  pal <- brewer.pal(9, "YlOrRd")
  color_col <- colorRampPalette(pal)(10)

  iGraph_ <- att.setv(g = iGraph_, from = "betweenness", to = "nodeColor",
                      cols = color_col, na.col = "grey80", breaks = seq(0, 1, 0.1))

  # Set the betweenness into label name [warning: bad visualization]
  V(iGraph_)$nodeAlias <- paste0(names(V(iGraph_)), " | B: ", V(iGraph_)$betweenness
                                 , " | Clu: ", V(iGraph_)$clustering, " | Clo: ", V(iGraph_)$closeness)

  V(iGraph_)$nodeAlias <- names(V(iGraph_))

  # Set the graph direction
  E(iGraph_)$arrowDirection <- 1
  V(iGraph_)$nodeLineWidth <- 5

  # Create RedPort object
  rdp <- RedPort()

  # Open the connection
  calld(rdp)

  # Add the graph into the RedPort
  addGraph(rdp, iGraph_)

  # Add the legend with color scale
  addLegend.color(rdp, colvec=iGraph_$legNodeColor$scale, size=15, labvec=iGraph_$legNodeColor$legend,
                  title="Betweenness Centrality Scale (BCS)")

  # Relax the graph visualization
  relax(rdp)

  # Select the bottlenecks
  selectNodes(rdp, names(bottleneck_))

  # print log
  if (verbose_) {
    print("Bottleneck visualization in RedPort generated successfully!")
  }
}

#' Get the image from a given KEGG pathway
#'
#' Given a KEGG pathway ID, this function saves its image in the current
#' working directory.
#'
#' @param pathway A KEGG pathway ID.
#'
#' @param IDs Character vector containing ENTREZ or KO identifiers.
#'
#' @return This function returns an image (PNG).
#'
#' @examples
#' \dontrun{
#' getPathwayImage("hsa00010")
#' getPathwayImage("ko00010")
#' }
#'
#' @importFrom KEGGREST keggLink
#' @importFrom pathview pathview
#'
#' @author
#' Diego Morais / Igor Brandão

printBottleneckPathwayImage <- function(pathway_, bottleneck_, verbose_=FALSE) {

  # Basic info
  species <- gsub("^([[:alpha:]]*).*$", "\\1", pathway)

  pathway <- gsub("^[[:alpha:]]*(.*$)", "\\1", pathway_)
  #bottleneck <- gsub("^[[:alpha:]]*:(.*$)", "\\1", names(bottleneck_))
  bottleneck <- bottleneck_
  data(bods, package = "pathview", verbose = FALSE)

  # Generate the pathway with bottlenecks
  img <- pathview::pathview(gene.data = bottleneck, pathway.id = pathway,
        species = species,
        out.suffix = "_bottleneck",
        high = list(gene = "#FF6961"),
        kegg.native = TRUE,
        same.layer = FALSE,
        new.signature = FALSE,
        plot.col.key = FALSE,
        map.symbol = TRUE, # bug
        gene.annotpkg = NA, # bug
        map.null = FALSE) # cpd size

  # Remove the generated xml file
  invisible(suppressWarnings(file.remove(paste0(species, pathway, ".xml"))))

  # print log
  if (verbose_) {
    print(paste0("Pathway ", pathway_, " visualization generated successfully!"))
  }
}

#' Function to generate interactive networks based on visNetwork library
#'
#' @param network_ Network data frame with all properties.
#' @param networkProperties_ Contains main information about the network nodes.
#' @param pathway_ Network name.
#'
#' @return This function does not return nothing, just export files.
#'
#' @examples
#' \dontrun{
#' printInteractiveNetwork(pathwayData, properties)
#' }
#'
#' @author
#' Igor Brandão

generateInteractiveNetwork <- function(network_, networkProperties_, pathway_="") {

  # Color pallet
  pal <- brewer.pal(9, "YlOrRd")
  pal2 <- brewer.pal(8, "Dark2")

  # Calculates the network bottleneck
  iGraph <- igraph::graph_from_data_frame(network_, directed = FALSE)

  # Convert the iGraph object toVisNetworkData
  data <- toVisNetworkData(iGraph)

  # Create vis object
  vis.nodes <- networkProperties_
  vis.links <- data$edges

  vis.nodes$color.border <- "white"
  vis.nodes$borderWidth <- 0

  # Apply the border color by bottleneck status
  vis.nodes$color.border[which(vis.nodes$is_bottleneck == 0)] <- "white"
  vis.nodes$color.border[which(vis.nodes$is_bottleneck == 1)] <- "blue"
  vis.nodes$borderWidth[which(vis.nodes$is_bottleneck == 0)] <- 0
  vis.nodes$borderWidth[which(vis.nodes$is_bottleneck == 1)] <- 4

  vis.nodes$shadow <- TRUE # Nodes will drop shadow
  vis.nodes$id     <- row.names(vis.nodes) # Node ID
  vis.nodes$label  <- paste0(row.names(vis.nodes), "\n(", vis.nodes$bottleneck_classification, ")") # Node label
  vis.nodes$title  <- paste0(vis.nodes$bottleneck_classification, " - ", row.names(vis.nodes), " degree: ", vis.nodes$degree, " betweenness: ", vis.nodes$betweenness) # Text on click
  vis.nodes$borderWidth <- 2 # Node border width

  # Properties when node highlighted
  vis.nodes$color.highlight.background <- "orange"
  vis.nodes$color.highlight.border <- "darkred"

  betweennessScaleValues <- 1

  tryCatch({
    # Generates the background color scale
    betweennessScaleValues <- cut(vis.nodes$betweenness, breaks = seq(min(vis.nodes$betweenness),
                                  max(vis.nodes$betweenness), len = 100),
                                  include.lowest = TRUE)

  }, error=function(e) {})

  # Apply the background color scale
  vis.nodes$color.background <- colorRampPalette(pal)(99)[betweennessScaleValues]

  # Apply node size according to its frequency
  vis.nodes$size <- scales::rescale(vis.nodes$freq, to=c(10, 30))

  # Set network links properties
  vis.links$width <- 1 # line width
  vis.links$color <- "gray"    # line color
  vis.links$arrows <- "middle" # arrows: 'from', 'to', or 'middle'
  vis.links$smooth <- FALSE    # should the edges be curved?
  vis.links$shadow <- FALSE    # edge shadow

  # Generate the visNetwor object
  visNetworkObj <- visNetwork(nodes = vis.nodes, edges = vis.links,
             background="#eeefff", width = '1200px', height = '800px',
             main=paste0("Pathway ", pathway_),
             submain=paste0("Nodes: ", length(V(iGraph)), " Edges: ", length(E(iGraph))),
             footer= "Note: nodes sizes are related to its frequencies")

  # Define the legend groups
  visNetworkObj <- visGroups(visNetworkObj, groupname = "Bottleneck", shape = "star",
                       color = list(background = "gray", border="black"))

  visNetworkObj <- visGroups(visNetworkObj, groupname = "Non-bottleneck", shape = "dot",
                       color = list(background = "tomato", border="black"))

  # Add a legend
  visNetworkObj <- visLegend(visNetworkObj, enabled = TRUE, useGroups = TRUE,
            main="Legend", position="left", ncol=1)

  # Add custom options
  visNetworkObj <- visOptions(visNetworkObj, highlightNearest = TRUE, selectedBy = "bottleneck_classification")

  # Generate the network
  return(visNetworkObj)
}
