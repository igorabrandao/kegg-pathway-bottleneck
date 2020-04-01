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
library(RColorBrewer)
library(visNetwork)
library(scales)

#*******************************************************************************************#

# ---- GRAPH PRINT SECTION ----

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
#' @param pathwayDetail_ The details related to the pathway.
#' @param dynamicNetwork_ Determins whether or not the network can be manipulated.
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

generateInteractiveNetwork <- function(network_, networkProperties_, pathway_="", pathwayDetail_=NULL, dynamicNetwork_=FALSE) {

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
  vis.nodes$borderWidth[which(vis.nodes$is_bottleneck == 0)] <- 2 # Node border width
  vis.nodes$borderWidth[which(vis.nodes$is_bottleneck == 1)] <- 4 # AP Node border width

  vis.nodes$id     <- vis.nodes$dictID # Node ID
  vis.nodes$label  <- paste0(vis.nodes$name, "\n(", vis.nodes$AP_classification, ")") # Node label
  vis.nodes$title  <- paste0("Node: ", vis.nodes$name, "<hr>",
                             "Classification: ", vis.nodes$AP_classification, "<br>",
                             "Is AP: ", ifelse(vis.nodes$is_bottleneck==1, 'Yes', 'No') , "<br>",
                             "AP impact: ", vis.nodes$bottleneckImpact, "<br>",
                             "Disconnected components: ", vis.nodes$bottleneckDisconnectedComponents, "<hr>",
                             "Community: ", vis.nodes$community, "<br>",
                             "Degree: ", vis.nodes$degree, "<br>",
                             "Betweenness: ", format(round(vis.nodes$betweenness, 4), nsmall = 4), "<br>",
                             "Clustering coefficient: ", format(round(vis.nodes$clusteringCoef, 4), nsmall = 4), "<br>",
                             "Closeness coefficient: ", format(round(vis.nodes$closenessCoef, 4), nsmall = 4), "<br>",
                             "Authority score: ", format(round(vis.nodes$authorityScore, 4), nsmall = 4), "<br>",
                             "Hub score: ", format(round(vis.nodes$hubScore, 4), nsmall = 4), "<hr>",
                             "Frequency: ", format(round(vis.nodes$percentage, 2), nsmall = 2), "%") # Text on click
  vis.nodes$shadow <- TRUE # Nodes will drop shadow

  # Properties when node highlighted
  vis.nodes$color.highlight.background <- "orange"
  vis.nodes$color.highlight.border <- "darkred"

  vis.nodes$color.highlight.background[which(vis.nodes$is_bottleneck == 1)] <- "#20639B"
  vis.nodes$color.highlight.border[which(vis.nodes$is_bottleneck == 1)] <- "#173F5F"

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
  vis.nodes$size <- scales::rescale(vis.nodes$occurrences, to=c(10, 30))

  # Set network links properties
  vis.links$width <- 1 # line width
  vis.links$color <- "gray"    # line color
  vis.links$arrows <- "middle" # arrows: 'from', 'to', or 'middle'
  vis.links$smooth <- FALSE    # should the edges be curved?
  vis.links$shadow <- FALSE    # edge shadow

  # Paint the edges bridges in red
  G <- igraph::graph_from_data_frame(vis.links[,1:2], directed = FALSE)
  num_comp <- length(decompose.graph(G))

  for (i in 1:length(E(G))) {
    G_sub <- delete.edges(G, i)
    if ( length( decompose.graph(G_sub) ) > num_comp ) vis.links$color[i] <- "red"
  }

  # Generate the visNetwor object
  if (is.null(pathwayDetail_) | length(pathwayDetail_) == 0) {
    visNetworkObj <- visNetwork(nodes = vis.nodes, edges = vis.links,
                                background="#ffffff", width = '100%', height = '100vh',
                                main=paste0("Pathway ", pathway_),
                                submain=paste0("<b>Nodes:</b> ", length(V(iGraph)), " <b>Edges:</b> ", length(E(iGraph))))
  } else {
    visNetworkObj <- visNetwork(nodes = vis.nodes, edges = vis.links,
                                background="#ffffff", width = '100%', height = '100vh',
                                main=paste0("Pathway ", pathway_, " - ", pathwayDetail_$NAME),
                                submain=paste0(
                                  "<br> <b>Description:</b> ", pathwayDetail_$DESCRIPTION,
                                  "<br><br> <b>Class:</b> ", pathwayDetail_$CLASS,
                                  "<br><br> <b>Nodes:</b> ", length(V(iGraph)), " <b>Edges:</b> ", length(E(iGraph))
                                ))
  }

  visNetworkObj <- visNetwork(nodes = vis.nodes, edges = vis.links, background="#ffffff", width = '100%', height = '85vh')

  # Define the legend groups
  #visNetworkObj <- visGroups(visNetworkObj, groupname = "Bottleneck", shape = "star", color = list(background = "gray", border="black"))
  #visNetworkObj <- visGroups(visNetworkObj, groupname = "Non-bottleneck", shape = "dot", color = list(background = "tomato", border="black"))

  # Add a legend
  #visNetworkObj <- visLegend(visNetworkObj, enabled = TRUE, useGroups = TRUE, main="Legend", position="left", ncol=1)

  # Generate a dynamic network
  if (dynamicNetwork_) {
    visNetworkObj <- visExport(visNetworkObj, type = "png", name = "network", label = paste0("Export as png"), background = "#fff",
              float = "right", style = NULL,loadDependencies = TRUE)

    visNetworkObj <- visExport(visNetworkObj, type = "jpeg", name = "network", label = paste0("Export as jpeg"), background = "#fff",
              float = "left", style = NULL,loadDependencies = TRUE)

    visNetworkObj <- visExport(visNetworkObj, type = "pdf", name = "network", label = paste0("Export as pdf"), background = "#fff",
              float = "right", style = NULL,loadDependencies = TRUE)

    # Add custom physics
    visNetworkObj <- visPhysics(visNetworkObj, stabilization = TRUE, solver = 'forceAtlas2Based',
                                forceAtlas2Based = list(gravitationalConstant = -75, avoidOverlap = 0.3))

    # Add custom options
    visNetworkObj <- visOptions(visNetworkObj, autoResize = TRUE, manipulation = TRUE, selectedBy = 'AP_classification',
                                highlightNearest = list(enabled = T, degree = 2, hover = T))

    # Add interaction
    visNetworkObj <- visInteraction(visNetworkObj, navigationButtons = TRUE, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE,
                                    keyboard = TRUE, hideEdgesOnDrag = TRUE, tooltipDelay = 0)
  } else {
    # Static network
    visNetworkObj <- visPhysics(visNetworkObj, stabilization = TRUE, solver = 'forceAtlas2Based',
                                forceAtlas2Based = list(gravitationalConstant = -75, avoidOverlap = 0.3))
  }

  # Generate the network
  return(visNetworkObj)
}
