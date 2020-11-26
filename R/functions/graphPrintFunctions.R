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
library(dplyr)

# Graph handling
library(igraph)
library(visNetwork)

# Graph plot
library(ggraph)

# Graph layouts
library(graphlayouts)
library(oaqc)

# Color pallete and scale
library(RColorBrewer)
library(viridis)
library(scales)

# Image export
library(svglite)

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

generateInteractiveNetwork <- function(network_, networkProperties_, pathway_="", org_="", pathwayDetail_=NULL, dynamicNetwork_=FALSE) {

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

  # Define the nodes group
  vis.nodes$group = vis.nodes$is_bottleneck

  if (sum(vis.nodes$group == 1, na.rm = T) > 0) {
    vis.nodes[vis.nodes$group == 1,]$group <- 'AP'
  }

  if (sum(vis.nodes$group == 0, na.rm = T) > 0) {
    vis.nodes[vis.nodes$group == 0,]$group <- 'Non-AP'
  }

  # Set the initial nodes attributes
  vis.nodes$color.border <- "white"
  vis.nodes$borderWidth <- 0

  # Apply the border color by bottleneck status
  vis.nodes$color.border[which(vis.nodes$is_bottleneck == 0)] <- "white"
  vis.nodes$color.border[which(vis.nodes$is_bottleneck == 1)] <- "blue"
  vis.nodes$borderWidth[which(vis.nodes$is_bottleneck == 0)] <- 2 # Node border width
  vis.nodes$borderWidth[which(vis.nodes$is_bottleneck == 1)] <- 4 # AP Node border width

  vis.nodes$id     <- vis.nodes$dictID # Node ID
  vis.nodes$label  <- paste0(vis.nodes$name, "\n(", vis.nodes$AP_classification, ")") # Node label
  vis.nodes$title  <- paste0("EC: ", vis.nodes$name, "<br>",
                             "Entrez: ", vis.nodes$entrez, "<hr>",
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
                             "Frequency: ", format(round(vis.nodes$percentage, 2), nsmall = 2), "% <hr>",
                             "More info: ", vis.nodes$link) # Text on click
  vis.nodes$shadow <- TRUE # Nodes will drop shadow

  # Properties when node highlighted
  vis.nodes$color.highlight.background <- "orange"
  vis.nodes$color.highlight.border <- "darkred"

  vis.nodes$color.highlight.background[which(vis.nodes$is_bottleneck == 1)] <- "#20639B"
  vis.nodes$color.highlight.border[which(vis.nodes$is_bottleneck == 1)] <- "#173F5F"

  betweennessScaleValues <- 1

  # Replace NA values
  if (sum(is.na(vis.nodes$betweenness), na.rm = T) > 0) {
    vis.nodes[is.na(vis.nodes$betweenness),]$betweenness = 0
  }

  tryCatch({
    # Generates the background color scale
    betweennessScaleValues <- cut(vis.nodes$betweenness, breaks = seq(min(vis.nodes$betweenness),
                                  max(vis.nodes$betweenness), len = 100),
                                  include.lowest = TRUE)

  }, error=function(e) {})

  # Apply the background color scale
  vis.nodes$color.background <- colorRampPalette(pal)(99)[betweennessScaleValues]

  # Apply node size according to its frequency
  vis.nodes$size <- scales::rescale(vis.nodes$percentage, to=c(10, 30))

  # Set network links properties
  vis.links$width <- 1 # line width
  vis.links$arrows <- "middle" # arrows: 'from', 'to', or 'middle'
  vis.links$smooth <- TRUE    # should the edges be curved?
  vis.links$shadow <- FALSE    # edge shadow

  # line color
  vis.links$color <- NA

  if (sum(is.na(vis.links$reaction1Status), na.rm = T) > 0) {
    vis.links[is.na(vis.links$reaction1Status),]$reaction1Status = 'reversible'
  }

  if (sum(vis.links$reaction1Status == 'reversible', na.rm = T) > 0) {
    vis.links[vis.links$reaction1Status == 'reversible',]$edge_color <- "gray"
  }

  if (sum(vis.links$reaction1Status == 'irreversible', na.rm = T) > 0) {
    vis.links[vis.links$reaction1Status == 'irreversible',]$edge_color <- "darkred"
  }

  # Line title
  vis.links$title <- paste0("Reaction: ", vis.links$reaction1, "<br>",
                            "Status: ", vis.links$reaction1Status, "<br>",
                            "Node1: ", vis.links$ec1, "<br>",
                            "Node2: ", vis.links$ec2, "<br>") # Text on click

  # Paint the edges bridges in red
  #G <- igraph::graph_from_data_frame(vis.links[,1:2], directed = FALSE)
  #num_comp <- length(decompose.graph(G))

  #for (i in 1:length(E(G))) {
  #  G_sub <- delete.edges(G, i)
  #  if ( length( decompose.graph(G_sub) ) > num_comp ) vis.links$color[i] <- "red"
  #}

  # Generate the visNetwor object
  if (is.null(pathwayDetail_) | length(pathwayDetail_) == 0) {
    visNetworkObj <- visNetwork(nodes = vis.nodes, edges = vis.links,
                                background="#ffffff", width = '100%', height = '100vh',
                                main=paste0("Pathway ", org_, pathway_),
                                submain=paste0("<b>Nodes:</b> ", length(V(iGraph)), " <b>Edges:</b> ", length(E(iGraph))))
  } else {
    visNetworkObj <- visNetwork(nodes = vis.nodes, edges = vis.links,
                                background="#ffffff", width = '100%', height = '100vh',
                                main=paste0("Pathway ", org_, pathway_, " - ", pathwayDetail_$NAME),
                                submain=paste0(
                                  "<br> <b>Description:</b> ", pathwayDetail_$DESCRIPTION,
                                  "<br><br> <b>Class:</b> ", pathwayDetail_$CLASS,
                                  "<br><br> <b>Nodes:</b> ", length(V(iGraph)), " <b>Edges:</b> ", length(E(iGraph))
                                ))
  }

  #visNetworkObj <- visNetwork(nodes = vis.nodes, edges = vis.links, background="#ffffff", width = '100%', height = '85vh')

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

#' Function to generate a static network based on ggraph library
#'
#' @param network_ Network data frame with all properties.
#' @param networkProperties_ Contains main information about the network nodes.
#' @param pathway_ Network name.
#' @param org_ Organism name.
#'
#' @return This function does not return nothing, just export files.
#'
#' @examples
#' \dontrun{
#' generateStaticNetwork(pathwayData, properties)
#' }
#'
#' @author
#' Igor Brandão
#'
generateStaticNetwork <- function(network_, networkProperties_, pathway_="", org_="", customLayout_="sparse_stress") {

  #----------------------------#
  # [ORGANIZZE THE GRAPH DATA] #
  #----------------------------#
  vertices <- networkProperties_
  relations <- network_

  #--------------------------#
  # [PREPARE THE GRAPH DATA] #
  #--------------------------#

  # Remove unnecessary data
  vertices$X <- NULL
  vertices$x <- NULL
  vertices$y <- NULL

  # Adjust the columns names
  names(relations)[names(relations) == "node1"] <- "from"
  names(relations)[names(relations) == "node2"] <- "to"

  # Change the relation from -> to names
  relations$from <- relations$entryID1
  relations$to <- relations$entryID2

  #---------------------#
  # [VERTEX AESTHETICS] #
  #---------------------#

  # Set the node label
  vertices$label  <- vertices$name

  # Set the initial nodes aesthetic attributes
  vertices$color.border <- "white"
  vertices$borderWidth <- 0

  # Apply the border color by bottleneck status
  if (sum(vertices$is_bottleneck == 0, na.rm = T) > 0) {
    vertices$color.border[which(vertices$is_bottleneck == 0)] <- "#ffffff"
    vertices$borderWidth[which(vertices$is_bottleneck == 0)] <- 1 # Node border width
  }

  if (sum(vertices$is_bottleneck == 1, na.rm = T) > 0) {
    vertices$color.border[which(vertices$is_bottleneck == 1)] <- "#005b96"
    vertices$borderWidth[which(vertices$is_bottleneck == 1)] <- 2 # AP Node border width
  }

  # Vertex drop shadow
  vertices$shadow <- TRUE # Nodes will drop shadow

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # ********************************************** #
  # Vertex aesthetics binded to scalar attributes
  # ********************************************** #

  # Vertex background color scale according to the betweenness
  betweennessScaleValues <- 1

  # Replace NA values
  if (sum(is.na(vertices$betweenness), na.rm = T) > 0) {
    vertices[is.na(vertices$betweenness),]$betweenness = 0
  }

  tryCatch({
    # Generates the background color scale
    betweennessScaleValues <- cut(vertices$betweenness, breaks = seq(min(vertices$betweenness),
        max(vertices$betweenness), len = 100), include.lowest = TRUE)

  }, error=function(e) {})

  # Color pallet
  pal <- brewer.pal(9, "YlOrBr")

  # Apply the background color scale
  vertices$color.background <- colorRampPalette(pal)(99)[betweennessScaleValues]

  # Apply node size according to its frequency
  vertices$vertex_size <- scales::rescale(vertices$percentage, to=c(0, 100))

  #--------------------#
  # [EDGES AESTHETICS] #
  #--------------------#

  # Set network links properties
  relations$edge_width <- 0.1 # line width
  relations$edge_arrows <- "middle" # arrows: 'from', 'to', or 'middle'
  relations$edge_smooth <- TRUE    # should the edges be curved?
  relations$edge_shadow <- FALSE    # edge shadow

  # line color
  relations$edge_color <- NA

  if (sum(is.na(relations$reaction1Status), na.rm = T) > 0) {
    relations[is.na(relations$reaction1Status),]$reaction1Status = 'reversible'
  }

  if (sum(relations$reaction1Status == 'reversible', na.rm = T) > 0) {
    relations[relations$reaction1Status == 'reversible',]$edge_color <- "gray"
  }

  if (sum(relations$reaction1Status == 'irreversible', na.rm = T) > 0) {
    relations[relations$reaction1Status == 'irreversible',]$edge_color <- "darkred"
  }

  #-------------------------#
  # [SET THE FACTORS ORDER] #
  #-------------------------#

  # Order the vertices by the entryID
  vertices$entryID <- factor(vertices$entryID, levels = vertices$entryID[order(vertices$entryID)])

  #-----------------------------#
  # [GENERATE THE GRAPH OBJECT] #
  #-----------------------------#
  iGraph <- igraph::graph_from_data_frame(relations, directed = FALSE, vertices = vertices)

  #------------------#
  # [PLOT THE GRAPH] #
  #------------------#
  # Validate the used layout since the sparse stress require additional parameters
  if (customLayout_ == "sparse_stress") {
    staticGraph <- ggraph(iGraph, layout = customLayout_, pivots = nrow(vertices), weights = NA)
  } else {
    staticGraph <- ggraph(iGraph, layout = customLayout_)
  }

  staticGraph <- staticGraph +
    # Edges
    geom_edge_fan(aes(colour = reaction1Status),
                  edge_alpha = 0.5,
                  angle_calc = 'along',
                  label_dodge = unit(2.5, 'mm'),
                  arrow = arrow(length = unit(2, 'mm'), type = 'closed'),
                  end_cap = circle(3, 'mm')) +

    # Nodes
    geom_node_point(aes(fill = betweenness, size = vertex_size, stroke = borderWidth, colour = AP_classification),
                    shape=21, alpha = 1) +

    # Nodes label
    geom_node_text(aes(filter = vertex_size >= 0, label = label),
                   size = 3, family="serif", repel = TRUE, check_overlap = TRUE,
                   nudge_y = -0.19) +

    # Nodes customizations
    scale_fill_gradientn("Betweenness", colours = brewer.pal(9, "YlOrBr"), limits=c(min(vertices$betweenness), max(vertices$betweenness))) +
    scale_color_manual("AP classification", values = c('blue', 'black')) +

    # Edges customizations
    scale_edge_color_manual("Reaction status", values = c('#ff7b7b', 'grey66')) +
    scale_edge_width_continuous(range = c(0.2,3)) +

    # Theming
    theme_graph() +
    theme(legend.position = "right") +

    # Rename the legends
    labs(size = "Enzymes frequency (%)") +

    # Set the legends order
    guides(size = guide_legend(order = 1),
           colour = guide_legend(order = 2))

  staticGraph

  #---------------------------#
  # [EXPORT THE GRAPH FIGURE] #
  #---------------------------#
  filename <- paste0(org_, pathway_)

  # Export the network
  if (!dir.exists(file.path(paste0('./output/network/static')))) {
    dir.create(file.path(paste0('./output/network/static')), showWarnings = FALSE, mode = "0775")
  }

  if (dir.exists(file.path(paste0('./output/network/static')))) {
    ggsave(paste0('./output/network/static/', filename, '.png'), width = 35, height = 20, units = "cm")

    # Export the svg
    if (!dir.exists(file.path(paste0('./output/network/static/svg')))) {
      dir.create(file.path(paste0('./output/network/static/svg')), showWarnings = FALSE, mode = "0775")
    }

    ggsave(paste0('./output/network/static/svg/', filename, '.svg'), width = 35, height = 20, units = "cm")
  }
}
