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
#' Diego Morais & Igor Brandão

pathwayToDataframe <- function(pathway_) {
  genesOnly <- !grepl("^ko", pathway_)
  kgml <- suppressMessages(KEGGREST::keggGet(pathway_, "kgml"))
  mapkpathway <- KEGGgraph::parseKGML(kgml)
  mapkG <- KEGGgraph::KEGGpathway2Graph(mapkpathway, genesOnly)
  aux <- names(mapkG@edgeData@data)
  aux <- as.data.frame(aux, stringsAsFactors = FALSE)
  aux$node2 <- gsub("^.*\\|(.*)$", "\\1", aux$aux)
  colnames(aux)[1] <- "node1"
  aux$node1 <- gsub("^(.*)\\|.*$", "\\1", aux$node1)
  aux$org <- gsub("^([[:alpha:]]*).*$", "\\1", pathway_)
  aux$pathway <- gsub("^[[:alpha:]]*(.*$)", "\\1", pathway_)
  rm(pathway_, mapkpathway, kgml)
  return(aux)
}

# getReferencePathway ####

#' Get the reference KEGG pathway
#'
#' Given a KEGG pathway ID and a list of KOs/ECs, this function returns a
#' data.frame ready to create an igraph object with ECs.
#'
#' @param pathway A KEGG pathway ID.
#' @param ko_ec_dictionnaire_ Data.frame containing the KOs/ECs equivalence
#'
#' @return This function returns a data.frame containing the edges from a
#' KEGG pathway by its ECs.
#'
#' @examples
#' \dontrun{
#' df <- pathwayToDataframe("00010", KO2EC)
#' }
#'
#' @importFrom KEGGREST keggGet
#' @importFrom KEGGgraph parseKGML
#' @importFrom KEGGgraph KEGGpathway2Graph
#'
#' @author
#' Diego Morais / Igor Brandão
#'

getReferencePathway <- function(pathway_, ko_ec_dictionnaire_) {
  # Load the ECs list from KEGG according to a pathway CODE
  ecs <- keggLink("ec", paste0("map", pathway_))
  ecs <- unname(ecs)

  # Load the KEGG XML according to the KO reference pathway
  kgml <- suppressMessages(KEGGREST::keggGet(paste0("ko", pathway_), "kgml"))

  # Parse the XML data
  mapkpathway <- KEGGgraph::parseKGML(kgml)

  # Get the graph object
  mapkG <- KEGGgraph::KEGGpathway2Graph(mapkpathway, FALSE)

  # Remove intermediary variables
  rm(mapkpathway, kgml)

  # Get the edge data and convert it into a dataFrame
  graphData <- names(mapkG@edgeData@data)
  graphData <- as.data.frame(graphData, stringsAsFactors = FALSE)

  # Change the data.frame columns name
  graphData$node2 <- gsub("^.*\\|(.*)$", "\\1", graphData$graphData)
  colnames(graphData)[1] <- "node1"
  graphData$node1 <- gsub("^(.*)\\|.*$", "\\1", graphData$node1)

  # Remove ko: prefix
  graphData$node1 <- gsub("ko:", "", graphData$node1)
  graphData$node2 <- gsub("ko:", "", graphData$node2)

  # Filter: get just the row in which its KOs belong to KO/EC dictionary
  graphData <- graphData[graphData$node1 %in% names(ko_ec_dictionnaire_) &
                           graphData$node2 %in% names(ko_ec_dictionnaire_),]

  # Run through all data in first column
  for(ko in unique(graphData$node1)) {
    # Count how many ECs are related to the current KO
    flag <- ifelse(length(KO2EC[[ko]])>1, TRUE, FALSE)

    if (flag) {
      # Which rows should be replaced
      toReplace <- which(graphData$node1==ko)

      # Receive all column 2 entries according toReplace
      node2 <- graphData[toReplace, 2]
      node2 <- sort(node2)

      # Get the clean data from KO/EC dictionary
      info <- unname(unlist(ko_ec_dictionnaire_[ko]))

      # Create a temporary data.frame to "replace" the original data
      temp <- data.frame(node1 = rep(info, length(toReplace)),
                         node2 = sort(rep(node2, length(info))), stringsAsFactors = FALSE)

      # Remove the rows in original data.frame that should be translated to ECs
      graphData <- graphData[-toReplace,]

      # Add the temporary data.frame to the original one
      graphData <- rbind(graphData, temp)
    } else {
      graphData$node1[which(graphData$node1==ko)] <- ko_ec_dictionnaire_[ko]
    }
  }

  # Run through all data in first column
  for(ko in unique(graphData$node2)) {
    # Count how many ECs are related to the current KO
    flag <- ifelse(length(KO2EC[[ko]])>1, TRUE, FALSE)

    if (flag) {
      # Which rows should be replaced
      toReplace <- which(graphData$node2==ko)

      # Receive all column 1 entries according toReplace
      node1 <- unlist(graphData[toReplace, 1])
      node1 <- sort(node1)

      # Get the clean data from KO/EC dictionary
      info <- unname(unlist(ko_ec_dictionnaire_[ko]))

      # Create a temporary data.frame to "replace" the original data
      temp <- data.frame(node2 = rep(info, length(toReplace)),
                         node1 = sort(rep(node1, length(info))), stringsAsFactors = FALSE)

      # Invert the data frame column order
      temp <- temp[, c(2, 1)]

      # Remove the rows in original data.frame that should be translated to ECs
      graphData <- graphData[-toReplace,]

      # Add the temporary data.frame to the original one
      graphData <- rbind(graphData, temp)
    } else {
      graphData$node2[which(graphData$node2==ko)] <- ko_ec_dictionnaire_[ko]
    }
  }

  # Remove duplicates
  graphData <- unique(graphData)
  temp <- data.frame(node1 = unlist(graphData$node2),
                     node2 = unlist(graphData$node1), stringsAsFactors = FALSE)
  graphData <- rbind(graphData, temp)
  # Re-index data.frame
  rownames(graphData) <- NULL
  return(graphData)
}

entrezToEC <- function(entrez_, ko_dictionnaire_, ec_dictionnaire_) {
  load("dictionnaires/KO00010.RData")
  load("dictionnaires/KO2EC.RData")
  ko_dictionnaire_ <- ENTREZ2KO
  ec_dictionnaire_ <- KO2EC
  return(ec_dictionnaire_[[ko_dictionnaire_[[as.character(entrez_)]]]])
}

# unlist(entrezToECMulti(c(2821, 669)))
entrezToECMulti <- Vectorize(entrezToEC, vectorize.args = "entrez_")

ECToEntrez <- function(ec_list_) {
  # Unlist the
  unlistEC <- unlist(ec_dictionnaire)
  unlistKO <- unlist(ko_dictionnaire)

  # Get the KO from EC
  ko_list <- names(which(unlistEC == as_ids(ec_list_)))

  # Retrieve the entrez list
  entrez_list <- list()

  for (ko in ko_list) {
    entrez_list <- append(entrez_list, names(which(unlistKO == ko)))
  }

  return(unlist(entrez_list))
}

# getPathwayHighlightedGenes ####

#' Get the image from a given KEGG pathway
#'
#' Given a KEGG pathway ID, this function saves its image in the current
#' working directory.
#'
#' @param pathway A KEGG pathway ID.
#'
#' @param IDs Character vector containing ENTREZ or KO identifiers.
#'
#' @return This function saves an image (PNG) in the current working directory
#' and returns the identifiers of the highlighted nodes.
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

getPathwayHighlightedGenes <- function(pathway_, IDs = NULL) {
  img <- NULL
  species <- gsub("^([[:alpha:]]*).*$", "\\1", pathway_)

  if(is.null(IDs)) {
    IDs <- KEGGREST::keggLink(species, pathway_)
    IDs <- gsub("^[[:alpha:]]*:(.*$)", "\\1", IDs)
  }

  pathway <- gsub("^[[:alpha:]]*(.*$)", "\\1", pathway_)
  now <- format(Sys.time(), "%Y-%m-%d--%H-%M-%S")

  tryCatch({
    img <- suppressWarnings(suppressMessages(
      pathview::pathview(gene.data = IDs,
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
                         map.null = FALSE)
    )) # cpd size
  }, warning = function(w) {
    # warning-handler-code
  }, error = function(e) {
    # error-handler-code
  }, finally = {
    # cleanup-code
  })

  # Remove generated files
  rm(gene.idtype.bods, korg, cpd.simtypes, gene.idtype.list)
  invisible(suppressWarnings(file.remove(paste0(species, pathway, ".xml"))))
  invisible(suppressWarnings(file.remove(paste0(species, pathway, ".png"))))
  invisible(suppressWarnings(file.remove(paste0(species, pathway, ".", now, ".png"))))

  # Get the highlighted genes
  highlightedGenes <- list(unique(img$plot.data.gene$kegg.names))

  # Convert the list into a vector
  highlightedGenes <- unlist(highlightedGenes, use.names=FALSE)

  # Return the highlighted genes
  return(highlightedGenes)
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
