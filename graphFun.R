# pathway2dataframe ####

#' Get the edges from a given KEGG pathway
#'
#' Given a KEGG pathway ID, this function returns a data.frame ready to create
#' an igraph object.
#'
#' @param pathway A KEGG pathway ID.
#' 
#' @return This function returns a data.frame containing the edges from a
#' KEGG pathway.
#'
#' @examples
#' \dontrun{
#' df <- pathway2dataframe("hsa00010")
#' df2 <- pathway2dataframe("ko00010")
#' df3 <- pathway2dataframe("mmu00010")
#' }
#'
#' @importFrom KEGGREST keggGet
#' @importFrom KEGGgraph parseKGML
#' @importFrom KEGGgraph KEGGpathway2Graph
#'
#' @author
#' Diego Morais

pathway2dataframe <- function(pathway){
  genesOnly <- !grepl("^ko", pathway)
  kgml <- suppressMessages(KEGGREST::keggGet(pathway, "kgml"))
  mapkpathway <- KEGGgraph::parseKGML(kgml)
  mapkG <- KEGGgraph::KEGGpathway2Graph(mapkpathway, genesOnly)
  rm(pathway, mapkpathway, kgml)
  aux <- names(mapkG@edgeData@data)
  aux <- as.data.frame(aux, stringsAsFactors = FALSE)
  aux$node2 <- gsub("^.*\\|(.*)$", "\\1", aux$aux)
  colnames(aux)[1] <- "node1"
  aux$node1 <- gsub("^(.*)\\|.*$", "\\1", aux$node1)
  return(aux)
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

getPathwayImage <- function(pathway, IDs = NULL){
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

# graphProperties ####

#' Get some graph properties from a directed graph
#'
#' Given a data.frame of directed edges, this function computes the connectivity,
#' clustering coefficient, and betweenness.
#'
#' @param edges A data.frame containing directed edges from a KEGG graph.
#' 
#' @return This function returns a data.frame containing three columns: connectivity,
#' clustering coefficient, and betweenness.
#'
#' @examples
#' \dontrun{
#' df <- graphProperties(pathway2dataframe("hsa00010"))
#' }
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph betweenness
#' @importFrom igraph count_triangles
#' @importFrom igraph transitivity
#'
#' @author
#' Diego Morais

graphProperties <- function(edges){
  g <- igraph::graph_from_data_frame(edges, directed = TRUE)
  betweenness <- igraph::betweenness(g, normalized = TRUE)
  result <- data.frame(node = names(betweenness), betweenness = betweenness,
                       stringsAsFactors = FALSE)
  rownames(result) <- NULL
  k <- as.data.frame(table(edges$node1))
  result$connectivity <- 0
  result$connectivity <- k[match(result$node, k$Var1), 2]
  result$connectivity[is.na(result$connectivity)] <- 0
  result$triangles <- vapply(result$node, function(x){
    as.integer(igraph::count_triangles(g, vids = x))
  }, integer(1))
  result$clusteringCoef <- igraph::transitivity(g, vids = result$node,
                                                isolates = "zero",
                                                type = "local")
  return(result)
}
