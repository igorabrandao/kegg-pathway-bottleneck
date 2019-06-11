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

pathwayToDataframe <- function(pathway_, replaceOrg=FALSE, orgToReplace='') {
  tryCatch({
    # Determine the genesOnly parameter (Default true)
    genesOnly <- !grepl("^ko|^ec", pathway_)

    # Request the graph data from KEGG
    kgml <- suppressMessages(KEGGREST::keggGet(pathway_, "kgml"))

    # Convert it into graph object
    mapkpathway <- KEGGgraph::parseKGML(kgml)
    mapkG <- KEGGgraph::KEGGpathway2Graph(mapkpathway, genesOnly)

    # Get the node data
    aux <- names(mapkG@edgeData@data)

    # If node data is empty, use the edge data
    if (is.null(aux) | length(aux) == 0) {
      return(NULL)
    }

    # Adjust the columns
    aux <- as.data.frame(aux, stringsAsFactors = FALSE)
    aux$node2 <- gsub("^.*\\|(.*)$", "\\1", aux$aux)
    colnames(aux)[1] <- "node1"
    aux$node1 <- gsub("^(.*)\\|.*$", "\\1", aux$node1)

    if (length(unlist(aux)) > 0) {
      # It means the organism has the pathway
      if (!replaceOrg) {
        aux$org <- gsub("^([[:alpha:]]*).*$", "\\1", pathway_)
      } else {
        aux$org <-orgToReplace
      }

      aux$pathway <- gsub("^[[:alpha:]]*(.*$)", "\\1", pathway_)
      rm(pathway_, mapkpathway, kgml)
      return(aux)
    } else {
      # It means the organism doesn't have the pathway
      rm(pathway_, mapkpathway, kgml)
      return(NULL)
    }
  }, error=function(e) {
    print(paste0('The pathway ', pathway_, ' could no be found. Skipping it...'))
    return(NULL)
  })
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

# getPathwayHighlightedGenes ####

#' Get the image from a given KEGG pathway
#'
#' Given a KEGG pathway ID, this function saves its image in the current
#' working directory.
#'
#' @param pathway_ A KEGG pathway ID.
#'
#' @return This function returns a list of highlighted enzymes according
#' to the pathway
#'
#' @examples
#' \dontrun{
#' getPathwayImage("hsa00010")
#' }
#'
#' @importFrom KEGGREST keggLink
#' @importFrom pathview pathview
#'
#' @author
#' Igor Brandão
#'
getPathwayHighlightedGenes <- function(pathway_, genesOnly_=TRUE) {
  tryCatch({
    # Request the graph data from KEGG
    kgml <- suppressMessages(KEGGREST::keggGet(pathway_, "kgml"))

    # Convert it into graph object
    mapkpathway <- KEGGgraph::parseKGML(kgml)
    mapkG <- KEGGgraph::KEGGpathway2Graph(mapkpathway)

    # Get the node data
    aux <- names(mapkG@edgeData@data)

    # If node data is empty, use the edge data
    if (is.null(aux) | length(aux) == 0) {
      return(NULL)
    }

    # Adjust the columns
    aux <- as.data.frame(aux, stringsAsFactors = FALSE)
    aux$node2 <- gsub("^.*\\|(.*)$", "\\1", aux$aux)
    colnames(aux)[1] <- "node1"
    aux$node1 <- gsub("^(.*)\\|.*$", "\\1", aux$node1)

    # Convert the node2 column into node1 index_s
    aux <- unique(c(aux$node1, aux$node2))

    aux <- gsub("^[[:alpha:]]*(.*$)", "\\1", str_replace(aux, ":", ""))

    return(unique(aux))

  }, error=function(e) {
    print(paste0('It wasnt possible to retrieve the highlights from pathway ', pathway_, ' . Skipping it...'))
    return(NULL)
  })
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
#' Diego Morais / Igor Brandão

getGraphProperties <- function(iGraph_) {

  tryCatch({
    # Convert the dataFrame to iGraph object
    g <- igraph::graph_from_data_frame(iGraph_, directed = TRUE)

    # Calculates betweenness
    # Measure of centrality in a graph based on shortest paths
    betweenness <- igraph::betweenness(g, normalized = TRUE)

    # Define the result dataFrame
    result <- data.frame(node = names(betweenness), betweenness = betweenness,
                         stringsAsFactors = FALSE)

    rownames(result) <- NULL

    k <- as.data.frame(table(iGraph_$node1))

    # Calculates the connectivity
    # # minimum number of elements (nodes or edges) that need to be removed to separate
    # the remaining nodes into isolated subgraphs.
    result$connectivity <- 0
    result$connectivity <- k[match(result$node, k$Var1), 2]
    result$connectivity[is.na(result$connectivity)] <- 0

    # Calculates the triangles
    # How many triangles a vertex is part of, in a graph, or just list the triangles of a graph.
    result$triangles <- vapply(result$node, function(x){
      as.integer(igraph::count_triangles(g, vids = x))
    }, integer(1))

    # Calculates the clustering coefficient
    # Transitivity measures the probability that the adjacent vertices of a vertex are connected
    result$clusteringCoef <- igraph::transitivity(g, vids = result$node,
                                                  isolates = "zero",
                                                  type = "local")

    # Calculates the closeness coefficient
    # Measures how many steps is required to access every other vertex from a given vertex.
    result$closenessCoef <- suppressWarnings(igraph::closeness(g, vids=result$node))

    # Calculates the vertex communities
    # greedy method (hiearchical, fast method)
    c3 = cluster_edge_betweenness(g)
    result$community <- as.integer(membership(c3))

    # Calculates the Eigenvector Centrality Scores of Network Positions
    # It is a measure of the influence of a node in a network.
    eigen_centrality <- igraph::eigen_centrality(g, directed = TRUE)
    result$eigenvectorScore <- unlist(eigen_centrality[1]) # Just the scores

    # Calculates the eccentricity of a vertex
    # It is defined as the maximum distance of one vertex from other vertex
    result$eccentricity <- igraph::eccentricity(g, vids=result$node)

    # Calculates the radius of a vertex (entire graph)
    # The smallest eccentricity in a graph is called its radius
    result$radius <- igraph::radius(g)

    # Calculates the diameter of a vertex (entire graph)
    # The diameter of a graph is the length of the longest geodesic
    result$diameter <- igraph::diameter(g, directed = TRUE)

    # Return the result data frame
    return(result)

  }, error=function(e) {
    print(paste0('It wasnt possible to retrieve properties from the graph. Skipping it...'))
    return(NULL)
  })
}

##################################
# Identifiers conversion section #
##################################

entrezToEC <- function(entrez_, ko_dictionnaire_, ec_dictionnaire_) {
  return(ec_dictionnaire_[[ko_dictionnaire_[[as.character(entrez_)]]]])
}

# unlist(entrezToECMulti(c(2821, 669)))
entrezToECMulti <- Vectorize(entrezToEC, vectorize.args = "entrez_")

ECToEntrez <- function(ec_list_, ec_dictionnaire_, ko_dictionnaire_) {
  # Unlist the dictionaries
  unlistEC <- unlist(ec_dictionnaire_)
  unlistKO <- unlist(ko_dictionnaire_)

  # Get the KO from EC
  ko_list <- names(which(unlistEC == as_ids(ec_list_)))

  # Retrieve the entrez list
  entrez_list <- list()

  for (ko in ko_list) {
    entrez_list <- append(entrez_list, names(which(unlistKO == ko)))
  }

  return(unlist(entrez_list))
}

# convertEntrezToECWithoutDict ####

#' Receive an entrez list and return an EC dataFrame
#'
#' Given a list of entrez, this function computes converts it via web scraping
#' into EC dataFrame.
#'
#' @param entrez_list_ A list containing entrez codes.
#' @param chunk_size_ How many entrez should be processed at a time.
#' @param verbose_ Print every status message.
#'
#' @return This function returns a data.frame containing one column: EC code.
#'
#' @examples
#' \dontrun{
#' ec_df <- convertEntrezToECWithoutDict(enzymeTotalFrequency)
#' }
#'
#' @importFrom RCurl
#' @importFrom rvest
#' @importFrom stringr
#'
#' @author
#' Igor Brandão

convertEntrezToECWithoutDict <- function(entrez_list_, chunk_size_=50, verbose_=FALSE, pathway_name_='') {

  # Perform an auxiliar scraping into mygene plataform to convert Entrez to EC
  #' @param entrez_number_ Entrez number withou specie
  #' @examples
  #' \dontrun{
  #' ec_result <- auxiliarScrap(entrez_number_)
  #' }
  #
  auxiliarScrap <- function (entrez_number_) {
    # Perform an auxiliar scraping
    scraping_aux <- getURL(paste0("http://mygene.info/v3/gene/", entrez_number_, "?fields=ec%2Cname"))

    # Get the raw EC list
    ec_value <- unlist(str_extract_all(toString(scraping_aux), "(\\d+)(?:\\.(\\d+)\\.(\\d+)\\.(\\d+))"))

    # Return the result
    return(paste(ec_value, collapse = ' / '));
  }

  tryCatch({
    if (verbose_) {
      cat("\n")
      print("------------------------------------------------")
      print("INITIALIZING THE ENTREZ -> EC CONVERSION")
      print("------------------------------------------------")
      cat("\n")
    }

    # Empty EC list dataFrame
    ec_list_df <- data.frame("ec" = character(0), stringsAsFactors = FALSE)

    # Remove the names from list to get just the values
    entrez_list_ <- unname(entrez_list_)

    # Break the list into chunks
    chunk_size <- chunk_size_
    chunked_entrez_list <- split(entrez_list_, ceiling(seq_along(entrez_list_)/chunk_size))

    # Log variable
    log <- ""

    # Loop over each chunk
    for(idx in 1:length(chunked_entrez_list)) {
      # Store just the number code of each entrez
      current_entrez_list <- str_extract(unlist(chunked_entrez_list[idx]), "\\-*\\d+\\.*\\d*")

      # Format the entrez list to be requested
      request_param <- paste0("https://www.kegg.jp/dbget-bin/www_bget?",
                              paste(unlist(chunked_entrez_list[idx]), collapse= "+", sep=""))

      if (verbose_) {
        cat("\n")
        print("------------------------------------------------")
        print(paste0("RUNNING CHUNK [", idx, " OF ", length(chunked_entrez_list), "] WITH SIZE: ", chunk_size_))
        cat("\n")
        print(paste0("REQUEST: ", request_param))
        print("------------------------------------------------")
        cat("\n")

        log <- paste0(log, "\n\n")
        log <- paste0(log, "------------------------------------------------\n")
        log <- paste0(log, (paste0("RUNNING CHUNK [", idx, " OF ", length(chunked_entrez_list), "] WITH SIZE: ", chunk_size_)))
        log <- paste0(log, "\n")
        log <- paste0(log, paste0("REQUEST: ", request_param))
        log <- paste0(log, "\n------------------------------------------------")
        log <- paste0(log, "\n\n")
      }

      # Get the entire KEGG webpage
      scraping <- getURL(request_param)

      # Get the raw EC list
      ec_list <- unlist(str_extract_all(toString(scraping), "\\[EC:(.*?)\\]"))

      # Get the raw Entrez list
      entrez_list <- unlist(str_extract_all(toString(scraping), "<code><nobr>+([0-9])+&nbsp;&nbsp"))
      entrez_list <- str_extract(unlist(entrez_list), "\\-*\\d+\\.*\\d*")

      # Run each entrez to define how it gonne be converted
      for (item in 1:length(current_entrez_list)) {
        # Check if the KEGG html contains the current entrez
        entrez_position <- which(current_entrez_list[item]==entrez_list)

        # Use KEGG data
        if (length(entrez_position) > 0) {
          # Apply the regex to remove trash from EC string
          ec_list[entrez_position] <- str_replace_all(ec_list[entrez_position], "\\[EC:(.*?)\\>", "")
          ec_list[entrez_position] <- str_replace_all(ec_list[entrez_position], "</a>\\]", "")
          ec_list[entrez_position] <- str_replace_all(ec_list[entrez_position], "</a>(.*?)\\>", " / ")
          ec_list[entrez_position] <- str_replace_all(ec_list[entrez_position], "</a>", " / ")
          ec_list[entrez_position] <- str_replace_all(ec_list[entrez_position], "]", "")

          if (verbose_) {
            log <- paste0(log, paste0((item), ") ", unlist(chunked_entrez_list[idx])[entrez_position], " -> ", ec_list[entrez_position]))
            log <- paste0(log, "\n")
          }

          # Add the EC item to the dataFrame
          if (!is.null(ec_list[entrez_position]) & !is.na(ec_list[entrez_position])) {
            ec_list_df[nrow(ec_list_df) + 1,] = paste(ec_list[entrez_position], collapse = ' / ')
          } else {
            ec_list_df[nrow(ec_list_df) + 1,] = auxiliarScrap(current_entrez_list[entrez_position])
          }
        } else {
          # Alternative method
          ec_list_df[nrow(ec_list_df) + 1,] = auxiliarScrap(current_entrez_list[item])
        }
      }
    }

    if (verbose_) {
      # Write the log into file
      write(log, file = paste0("./log/entrez_ec_conversion_log_", pathway_name_, ".txt"))

      cat("\n")
      print("------------------------------------------------")
      print("ENTREZ -> EC CONVERSION FINISHED WITH SUCCESS!")
      print("------------------------------------------------")
      cat("\n")
    }

    # Return the dataFrame containing the EC list
    return(ec_list_df)

  }, error=function(e){
    print(paste0('Some error happened during the [Entrez -> EC] convertion. Skipping it...'))

    return(NULL)
  })
}
