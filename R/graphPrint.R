########################################
# Functions to print the iGraph_ object #
########################################

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
}

# getPathwayImage ####

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
#' Diego Morais

printBottleneckPathwayImage <- function(pathway, IDs = NULL) {
  species <- gsub("^([[:alpha:]]*).*$", "\\1", pathway)
  if(is.null(IDs)){
    IDs <- KEGGREST::keggLink(species, pathway)
    IDs <- gsub("^[[:alpha:]]*:(.*$)", "\\1", IDs)
  }
  pathway <- gsub("^[[:alpha:]]*(.*$)", "\\1", pathway)
  data(bods, package = "pathview", verbose = FALSE)
  now <- format(Sys.time(), "%Y%m%dT%H%M%S")
  img <- suppressWarnings(suppressMessages(pathview::pathview(gene.data = IDs,
                                                              pathway.id = pathway,
                                                              out.suffix = now,
                                                              species = species,
                                                              high = list(gene = "darkseagreen1"),
                                                              kegg.native = TRUE,
                                                              same.layer = FALSE,
                                                              new.signature = FALSE,
                                                              plot.col.key = FALSE,
                                                              map.symbol = FALSE, # bug
                                                              gene.annotpkg = NA, # bug
                                                              map.null = FALSE))) # cpd size
  invisible(suppressWarnings(file.remove(paste0(species, pathway, ".xml"))))
  return(NULL)
}
