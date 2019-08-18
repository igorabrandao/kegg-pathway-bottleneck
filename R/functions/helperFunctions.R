#********************************#
# Package with general functions #
#********************************#

# ---- IMPORT SECTION ----

# helperFunctions.R #

# Import the necessary libraries


#*******************************************************************************************#

# ---- MESSAGE HANDLER ----

#' Function to generate dynamic messages
#'
#' @param message_ String representing the content of the message.
#'
#' @return This function does not return nothing, just export files.
#'
#' @examples
#' \dontrun{
#' printMessage('Hello World!')
#' }
#'
#' @author
#' Igor Brandão

printMessage <- function(message_) {
  cat("\n")
  print("------------------------------------------------")
  print(message_)
  print("------------------------------------------------")
  cat("\n")
}

#*******************************************************************************************#

# ---- GRAPH HANDLER ----

#' Function to apply attributes into iGraph objects
#'
#' @param graph_ iGraph object.
#' @param default_ String representing the default value of the attribute.
#' @param valNodeList_ List of values to apply into the graph nodes.
#'
#' @return This function does not return nothing, just export files.
#'
#' @examples
#' \dontrun{
#' nAttrs$fillcolor <- makeAttr(graph, "lightblue", list(orange=toprbccName))
#' }
#'
#' @author
#' Igor Brandão

makeAttr <- function(graph_, default_, valNodeList_) {
  tmp <- nodes(graph_)
  x <- rep(default_, length(tmp)); names(x) <- tmp

  if(!missing(valNodeList_)) {
    stopifnot(is.list(valNodeList_))
    allnodes <- unlist(valNodeList_)
    stopifnot(all(allnodes %in% tmp))
    for(i in seq(valNodeList)) {
      x[valNodeList_[[i]]] <- names(valNodeList_)[i]
    }
  }
  return(x)
}

#*******************************************************************************************#

# ---- STRING HANDLER ----

#' Function to perform trim operation
#'
#' @param str_ The string to be trimmed.
#'
#' @return Returns string w/o leading or trailing whitespace
#'
#' @examples
#' \dontrun{
#'  myDummy$country <- trim(myDummy$country)
#' }
#'
#' @author
#' Igor Brandão

trim <- function (str_) {
  gsub("^\\s+|\\s+$", "", str_)
}
