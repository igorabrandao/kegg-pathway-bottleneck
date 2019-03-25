#########################################
# Functions to print the iGraph_ object #
#########################################

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
#' Diego Morais / Igor Brandão

printBottleneckPathwayImage <- function(pathway_, bottleneck_, verbose_=FALSE) {

  # Basic info
  species <- gsub("^([[:alpha:]]*).*$", "\\1", pathway)

  pathway <- gsub("^[[:alpha:]]*(.*$)", "\\1", pathway_)
  bottleneck <- gsub("^[[:alpha:]]*:(.*$)", "\\1", names(bottleneck_))
  data(bods, package = "pathview", verbose = FALSE)

  # Generate the pathway with bottlenecks
  img <- pathview::pathview(gene.data = bottleneck, pathway.id = pathway,
        species = species,
        out.suffix = "_bottleneck",
        high = list(gene = "#FF6961"),
        kegg.native = TRUE,
        same.layer = FALSE,
        new.signature = FALSE,
        plot.col.key = TRUE,
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
