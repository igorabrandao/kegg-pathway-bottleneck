#############################################
# Functions to handle iGraph object loading #
#############################################

# pathwayToDataframe ####

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
#' df <- pathwayToDataframe("hsa00010")
#' df2 <- pathwayToDataframe("ko00010")
#' df3 <- pathwayToDataframe("mmu00010")
#' }
#'
#' @importFrom KEGGREST keggGet
#' @importFrom KEGGgraph parseKGML
#' @importFrom KEGGgraph KEGGpathway2Graph
#'
#' @author
#' Diego Morais

pathwayToDataframe <- function(pathway) {
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

# getGraphProperties ####

#' Get some graph properties from a directed graph
#'
#' Given a data.frame of directed edges, this function computes the connectivity,
#' clustering coefficient, and betweenness.
#'
#' @param iGraph_ A data.frame containing directed edges from a KEGG graph.
#' 
#' @return This function returns a data.frame containing three columns: connectivity,
#' clustering coefficient, and betweenness.
#'
#' @examples
#' \dontrun{
#' df <- getGraphProperties(pathwayToDataframe("hsa00010"))
#' }
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph betweenness
#' @importFrom igraph count_triangles
#' @importFrom igraph transitivity
#' @importFrom igraph closeness
#'
#' @author
#' Diego Morais

getGraphProperties <- function(iGraph_) {
  g <- igraph::graph_from_data_frame(iGraph_, directed = TRUE)
  betweenness <- igraph::betweenness(g, normalized = TRUE)
  result <- data.frame(node = names(betweenness), betweenness = betweenness,
                       stringsAsFactors = FALSE)
  rownames(result) <- NULL
  k <- as.data.frame(table(iGraph_$node1))
  result$connectivity <- 0
  result$connectivity <- k[match(result$node, k$Var1), 2]
  result$connectivity[is.na(result$connectivity)] <- 0
  result$triangles <- vapply(result$node, function(x){
    as.integer(igraph::count_triangles(g, vids = x))
  }, integer(1))
  result$clusteringCoef <- igraph::transitivity(g, vids = result$node,
                                                isolates = "zero",
                                                type = "local")
  result$closenessCoef <- igraph::closeness(g, vids=result$node)
  return(result)
}
