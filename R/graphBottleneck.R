############################################
# Functions to manipulate graph properties #
# related to  bottlencks                   #
############################################

###########
# Setters #
###########

# setGraphCommunity ####

#' Given an iGraph object, this function set each vertex communities
#'
#' @param iGraph_ An iGraph object.
#'
#' @return This function returns the same iGraph object with an extra 'group' column.
#'
#' @examples
#' \dontrun{
#' iGraphObj <- setGraphCommunity(iGraphObj)
#' }
#'
#' @importFrom igraph cluster_edge_betweenness
#' @importFrom igraph membership
#' @importFrom igraph modularity
#'
#' @author
#' Igor Brandão

setGraphCommunity <- function(iGraph_, verbose_=FALSE) {
  # greedy method (hiearchical, fast method)
  c3 = cluster_edge_betweenness(iGraph_)

  # define the cluster attribute
  V(iGraph_)$group <- membership(c3)

  # modularity measure
  if (verbose_) {
    print(modularity(c3))
    print("Graph community setted successfully!")
  }

  # return the iGraph object
  return(iGraph_)
}

setGraphBetwenness <- function(iGraph_, normalized_=TRUE, verbose_=FALSE) {
  # calculates the betweenness
  V(iGraph_)$betweenness <- betweenness(iGraph_, normalized = normalized_)

  # normalizes betweenness
  if (normalized_) {
    V(iGraph_)$betweenness <- V(iGraph_)$betweenness/max(V(iGraph_)$betweenness)
  }

  # print the betweenness
  if (verbose_) {
    print(V(iGraph_)$betweenness)
    print("Graph betweenness setted successfully!")
  }

  # return the iGraph object
  return(iGraph_)
}

setGraphCloseness <- function(iGraph_, verbose_=FALSE) {
  # closeness refers to how connected a node is to its neighbors
  V(iGraph_)$closeness <- closeness(iGraph_, vids=V(iGraph_))

  # print the closeness
  if (verbose_) {
    print(V(iGraph_)$closeness)
    print("Graph closeness setted successfully!")
  }

  # return the iGraph object
  return(iGraph_)
}

setGraphClustering <- function(iGraph_, verbose_=FALSE) {
  # calculates the local clustering coefficient for each vertex
  V(iGraph_)$clustering <- transitivity(iGraph_, type = c("local"), vids = NULL,
                                       weights = NULL, isolates = c("NaN", "zero"))

  # print the clustering
  if (verbose_) {
    print(V(iGraph_)$clustering)
    print("Graph clustering setted successfully!")
  }

  # return the iGraph object
  return(iGraph_)
}

# --------------------------------------------------------------------------------------

###########
# Getters #
###########

getTopBetweenness <- function(iGraph_, betweenness_percentual_rate_=0.2, verbose_=FALSE) {
  topBetweenness <- sort(V(iGraph)$betweenness, decreasing=TRUE)
  topBetweennessWithRate <- topBetweenness[1:as.integer(length(topBetweenness) * betweenness_percentual_rate_)]

  # print the clustering
  if (verbose_) {
    print(V(iGraph)[match(topBetweennessWithRate, V(iGraph)$betweenness)])
    print(topBetweennessWithRate)
  }

  return(topBetweennessWithRate)
}

getGraphBottleneck <- function(iGraph_, verbose_=FALSE) {
  articulation_points <- articulation_points(iGraph_)

  # print the articulation points
  if (verbose_) {
    print(articulation_points)
    print("Graph bottlenecks calculated successfully!")
  }

  return(articulation_points)
}
