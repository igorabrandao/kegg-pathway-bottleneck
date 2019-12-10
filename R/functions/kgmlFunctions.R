#*************************************#
# Functions to handle KEGG kgml files #
#*************************************#

# kgmlFunctions.R #

# ---- IMPORT SECTION ----

#' This is the set of functions to handle the KEGG kgml files
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library("XML")

#*******************************************************************************************#

# ---- KGML FUNCTIONS ----

.xmlChildrenWarningFree <- function(xmlNode)
{
  # Verify if the XML node has a child node
  if (is.null(xmlNode$children))
    return(NULL)
  return(XML::xmlChildren(xmlNode))
}

parsePathwayInfo <- function(root)
{
  # Declare the pathwayInfo dataframe
  pathwayInfo <- data.frame(name = NA, org = NA, number = NA, title = NA, image = NA,
                            link = NA, stringsAsFactors = FALSE)

  # Access the entry attributes
  attrs <- xmlAttrs(root)

  # Retrieve the 1-dimensional attributes
  name <- attrs[["name"]]
  org <- attrs[["org"]]
  number <- attrs[["number"]]
  title <- getNamedElement(attrs, "title")
  image <- getNamedElement(attrs, "image")
  link <- getNamedElement(attrs, "link")

  # Set the pathwayInfo attributes
  pathwayInfo$name <- name
  pathwayInfo$org <- org
  pathwayInfo$number <- number
  pathwayInfo$title <- title
  pathwayInfo$image <- image
  pathwayInfo$link <- link

  # Return the pathwayInfo dataFrame
  return(pathwayInfo)
}

parseGraphics <- function(graphics)
{
  # Declare the graphic dataframe
  g <- data.frame(x = NA, y = NA, graphicalType = NA, width = NA, height = NA, fgcolor = NA,
                  bgcolor = NA, stringsAsFactors = FALSE)

  # Check if the graphics contains data
  if (is.null(graphics))
    return(g)

  # Access the entry attributes
  attrs <- xmlAttrs(graphics)

  # Set the graphical attributes
  g$x <- as.integer(getNamedElement(attrs, "x"))
  g$y <- as.integer(getNamedElement(attrs, "y"))
  g$graphicalType <- getNamedElement(attrs, "type")
  g$width <-as.integer(getNamedElement(attrs, "width"))
  g$height <- as.integer(getNamedElement(attrs, "height"))
  g$fgcolor <- getNamedElement(attrs, "fgcolor")
  g$bgcolor <- getNamedElement(attrs, "bgcolor")

  # Return the graphics dataFrame
  return(g)
}

parseEntry <- function(entry)
{
  # Declare the node dataframe
  newNode <- data.frame(entryID = NA, name = NA, type = NA, link = NA, component = NA,
                            reaction = NA, map = NA, x = NA, y = NA, graphicalType = NA,
                            width = NA, height = NA, fgcolor = NA, bgcolor = NA, stringsAsFactors = FALSE)

  # Access the entry attributes
  attrs <- xmlAttrs(entry)

  # Retrieve the 1-dimensional attributes
  entryID <- attrs[["id"]]
  name <- unname(unlist(strsplit(attrs["name"], " ")))
  name <- paste(name, collapse = " / ")
  type <- attrs[["type"]]
  link <- getNamedElement(attrs, "link")
  reaction <- getNamedElement(attrs, "reaction")
  map <- getNamedElement(attrs, "map")
  graphics <- xmlChildren(entry)$graphics

  # Retrieve the N-dimensional attributes
  g <- parseGraphics(graphics)

  # Verify the node type
  if (type != "group") {
    # Set the node attributes
    newNode$entryID <- entryID
    newNode$name <- name
    newNode$type <- type
    newNode$link <- link
    newNode$reaction <- reaction
    newNode$map <- map

    # Set the graphical attributes
    newNode$x <- g$x
    newNode$y <- g$y
    newNode$graphicalType <- g$graphicalType
    newNode$width <- g$width
    newNode$height <- g$height
    newNode$fgcolor <- g$fgcolor
    newNode$bgcolor <- g$bgcolor
  }
  else if (type == "group") {
    children <- xmlChildren(entry)
    children <- children[names(children) == "component"]

    if (length(children) == 0) {
      component <- as.character(NA)
    }
    else {
      component <- sapply(children, function(x) {
        if (xmlName(x) == "component") {
          return(xmlAttrs(x)["id"])
        }
        else {
          return(as.character(NA))
        }
      })
    }

    component <- unname(unlist(component))
    component <- paste(component, collapse = " / ")

    # Set the node attributes
    newNode$entryID <- entryID
    newNode$name <- name
    newNode$type <- type
    newNode$link <- link
    newNode$component <- component
    newNode$reaction <- reaction
    newNode$map <- map

    # Set the graphical attributes
    newNode$x <- g$x
    newNode$y <- g$y
    newNode$graphicalType <- g$graphicalType
    newNode$width <- g$width
    newNode$height <- g$height
    newNode$fgcolor <- g$fgcolor
    newNode$bgcolor <- g$bgcolor
  }

  # Return the node as a dataframe
  return(newNode)
}

parseSubType <- function(subtype)
{
  # Declare the subtype dataframe
  subtypeDf <- data.frame(name = NA, value = NA, stringsAsFactors = FALSE)

  # Access the subtype attributes
  attrs <- xmlAttrs(subtype)

  # Retrieve the 1-dimensional attributes
  subtypeDf$name <- attrs[["name"]]
  subtypeDf$value <- attrs[["value"]]

  # Return the subtype dataFrame
  return(subtypeDf)
}

parseRelation <- function(relation)
{
  # Declare the newEdge dataframe
  newEdge <- data.frame(entry1 = NA, entry2 = NA, type = NA, subtype = NA, stringsAsFactors = FALSE)

  # Access the relation attributes
  attrs <- xmlAttrs(relation)

  # Retrieve the 1-dimensional attributes
  newEdge$entry1 <- attrs[["entry1"]]
  newEdge$entry2 <- attrs[["entry2"]]
  newEdge$type <- attrs[["type"]]

  # Retrieve the N-dimensional attributes
  subtypeNodes <- xmlChildren(relation)
  subtypes <- sapply(subtypeNodes, parseSubType)
  subtypes <- paste(subtypes, collapse = " / ")

  # Apply the subtypes into the dataFrame
  newEdge$subtype <- subtypes

  # Return the edge dataFrame
  return(newEdge)
}

parseReaction <- function(reaction)
{
  # Declare the newReaction dataframe
  newReaction <- data.frame(name = NA, type = NA, children = NA, childrenNames = NA, substrateIndices = NA,
                            productIndices = NA, substrateName = NA, productName = NA, stringsAsFactors = FALSE)

  # Access the reaction attributes
  attrs <- xmlAttrs(reaction)

  # Retrieve the 1-dimensional attributes
  newReaction$name <- attrs[["name"]]
  newReaction$type <- attrs[["type"]]

  reactionChildren <- xmlChildren(reaction)
  newReaction$children <- paste(as.character(reactionChildren), collapse = " / ")

  childrenNames <- paste(as.character(names(xmlChildren(reaction))), collapse = " / ")
  newReaction$childrenNames <- childrenNames

  substrateIndices <- grep("^substrate$", childrenNames)
  productIndices <- grep("^product$", childrenNames)

  newReaction$substrateIndices <- paste(as.character(grep("^substrate$", childrenNames)), collapse = " / ")
  newReaction$productIndices <- paste(as.character(grep("^product$", childrenNames)), collapse = " / ")

  # Retrieve the N-dimensional attributes
  substrateAltName <- paste(as.character(vector("character", length(substrateIndices))), collapse = " / ")
  productAltName <- paste(as.character(vector("character", length(productIndices))), collapse = " / ")

  # Adjust the names
  substrateAltName <- paste(substrateAltName, collapse = " / ")
  productAltName <- paste(productAltName, collapse = " / ")

  newReaction$substrateName <- substrateAltName
  newReaction$productName <- productAltName

  for (i in seq(along = substrateIndices)) {
    ind <- substrateIndices[i]
    substrate <- children[[ind]]
    substrateName[i] <- xmlAttrs(substrate)[["name"]]
    substrateAltName[i] <- as.character(NA)
    substrateChildren <- .xmlChildrenWarningFree(substrate)
    if (!is.null(substrateChildren)) {
      substrateAlt <- substrateChildren$alt
      substrateAltName[i] <- xmlAttrs(substrateAlt)[["name"]]
    }
  }
  for (i in seq(along = productIndices)) {
    ind <- productIndices[i]
    product <- children[[ind]]
    productName[i] <- xmlAttrs(product)[["name"]]
    productChildren <- .xmlChildrenWarningFree(product)
    productAltName[i] <- as.character(NA)
    if (!is.null(productChildren)) {
      productAlt <- productChildren$alt
      productAltName[i] <- xmlAttrs(productAlt)[["name"]]
    }
  }

  # Adjust the names
  substrateAltName <- paste(substrateAltName, collapse = " / ")
  productAltName <- paste(productAltName, collapse = " / ")

  newReaction$substrateAltName <- substrateAltName
  newReaction$productAltName <- productAltName

  # Return the reaction dataFrame
  return(newReaction)
}

#*******************************************************************************************#

# ---- MAIN FLOW FUNCTION ----

#' Get the KEGG kgml file and parse its information
#'
#' Given a KEGG kgml file, this function parse this file and returns a data.frame contaning
#' all relevante information related to this pathway
#'
#' @param kgml_ Name of KGML file.
#'
#' @return This function returns a list of data.frames containing all relevante information
#' related to a pathway
#'
#' @examples
#' \dontrun{
#' lst <- kgmlToDataframe('ko00010.xml')
#' }
#'
#' @importFrom 'XML'
#'
#' @author
#' Igor Brandão

KGML2Dataframe <- function(kgml_) {
  # Verify and process the KGML file
  tryCatch(doc <- xmlTreeParse(kgml_, getDTD = FALSE, error = xmlErrorCumulator(immediate = FALSE)),
           error = function(e) {
             fileSize <- file.info(kgml_)$size[1]
             bytes <- sprintf("%d byte%s", fileSize, ifelse(fileSize >
                                                              1, "s", ""))
             msg <- paste("The file", kgml_, "seems not to be a valid KGML file\n")
             if (fileSize < 100L)
               msg <- paste(msg, "[WARNING] File size (", bytes,
                            ") is unsually small; please check.\n", sep = "")
             msg <- paste(msg, "\nDetailed error messages from",
                          "XML::xmlTreeParse:\n", sep = "")
             cat(msg)
             stop(e)
           })

  # Retrieve the XML root
  r <- xmlRoot(doc)

  # Retrieve the XML root children
  childnames <- sapply(xmlChildren(r), xmlName)

  # Define which children are entries
  isEntry <- childnames == "entry"

  # Define which children are relations
  isRelation <- childnames == "relation"

  # Define which children are reactions
  isReaction <- childnames == "reaction"

  # Retrieve the pathway info as dataFrame
  kegg.pathwayinfo <- parsePathwayInfo(r)

  # Retrieve the pathway nodes as dataFrame
  big.list.of.data.frames <- lapply(r[isEntry], parseEntry)
  kegg.nodes <- do.call(rbind, big.list.of.data.frames)
  rm(big.list.of.data.frames)

  # Retrieve the pathway edges as dataFrame
  big.list.of.data.frames <- lapply(r[isRelation], parseRelation)
  kegg.edges <- do.call(rbind, big.list.of.data.frames)
  rm(big.list.of.data.frames)

  # Retrieve the pathway reactions as dataFrame
  big.list.of.data.frames <- lapply(r[isReaction], parseReaction)
  kegg.reactions <- do.call(rbind, big.list.of.data.frames)
  rm(big.list.of.data.frames)

  # Create a list with pathway dataFrames info
  pathway <- list(pathwayinfo = kegg.pathwayinfo, nodes = kegg.nodes, edges = kegg.edges,
                  reactions = kegg.reactions)

  # Return the list with pathway dataFrames info
  return(pathway)
}

#' Get the KEGG kgml file and convert it into graph
#'
#' Given a KEGG kgml file, this function parse this file and returns the graph structure
#'
#' @param kgml_ Name of KGML file.
#' @param replaceOrg_ Flag to determine with wheter or not the organism name will be changed.
#' @param orgToReplace_ The organism identification.
#'
#' @return This function returns a list of data.frames containing all relevante information
#' related to a pathway
#'
#' @examples
#' \dontrun{
#' graph <- kgmlToDataframe('ko00010.xml')
#' graph <- kgmlToDataframe('hsa00010.xml', TRUE, 'hsa')
#' }
#'
#' @importFrom KEGGgraph
#' @importFrom 'XML'
#'
#' @author
#' Igor Brandão

KGML2Graph <- function(kgml_, replaceOrg_=FALSE, orgToReplace_='') {
  # Get the pathway name
  pathway_name <- onlyNumber(kgml_)

  tryCatch({
    # Convert it into graph object
    mapkG <- KEGGgraph::parseKGML2Graph(kgml_, genesOnly=FALSE)

    # Get the node data
    aux <- names(mapkG@edgeData@data)

    # If node data is empty, use the edge data
    if (is.null(aux) | length(aux) == 0) {
      tryCatch({
        # Try another way to get the node data
        aux <- mapkG@nodes
      }, error=function(e) {
        return(NULL)
      })

      if (is.null(aux) | length(aux) == 0) {
        return(NULL)
      }
    }

    # Adjust the columns
    aux <- as.data.frame(aux, stringsAsFactors = FALSE)
    aux$node2 <- gsub("^.*\\|(.*)$", "\\1", aux$aux)
    colnames(aux)[1] <- "node1"
    aux$node1 <- gsub("^(.*)\\|.*$", "\\1", aux$node1)

    if (length(unlist(aux)) > 0) {
      # It means the organism has the pathway
      if (!replaceOrg_) {
        aux$org <- gsub("^([[:alpha:]]*).*$", "\\1", pathway_name)
      } else {
        aux$org <- orgToReplace_
      }

      aux$pathway <- gsub("^[[:alpha:]]*(.*$)", "\\1", pathway_name)
      return(aux)
    } else {
      # It means the organism doesn't have the pathway
      return(NULL)
    }
  }, error=function(e) {
    print(paste0('The pathway ', pathway_name, ' could no be converted. Skipping it...'))
    return(NULL)
  })
}

KGML2GraphDictionary <- function(kgml_, replaceOrg_=FALSE, orgToReplace_='') {

  # Load the dicionaty
  dictionary <- read.csv(file='./output/pathwaysDictionary/dictionary.csv', header=TRUE, sep=",", stringsAsFactors=FALSE)

  if (is.null(dictionary) | nrow(dictionary) == 0) {
    # Save the log file
    printLog(message_='The pathways nodes dictionary could not be found. Skipping it...',
             file_='KGML2GraphDictionary')

    return(FALSE)
  }

  # Create the KGML dataframe
  pathway <- KGML2Dataframe(kgml_)

  if (is.null(pathway) | isempty(pathway)) {
    # Save the log file
    printLog(message_='The pathway data could not be found. Skipping it...', file_='KGML2GraphDictionary')

    return(FALSE)
  }

  # Generate the hand-made graph with dictionary identifiers
  pathwayGraph <- data.frame(node1='', node2='', ec1='', entryID1='', reaction1='', x1='', y1='',
                             ec2='', entryID2='', reaction2='', x2='',  y2='', org='', pathway='', stringsAsFactors = FALSE)

  #**************************************##
  # Fill the pathway graph with entry IDs #
  #**************************************##

  rows <- nrow(pathway$edges)

  if (is.null(rows) | isempty(rows)) {
    return(NULL)
  } else {
    for (edgeIdx in 1:rows) {
      pathwayGraph[edgeIdx, 'node1'] <- pathway$edges[edgeIdx,]$entry1
      pathwayGraph[edgeIdx, 'node2'] <- pathway$edges[edgeIdx,]$entry2
    }
  }

  # Set basic info
  if (!replaceOrg_) {
    pathwayGraph['org'] <- 'ec'
  } else {
    pathwayGraph['org'] <- orgToReplace_
  }

  pathwayGraph['pathway'] <- pathway$pathwayinfo$number

  #**************************************************##
  # Convert the entry IDs according to the dictionary #
  #**************************************************##

  rows <- nrow(pathwayGraph)
  for (nodeIdx in 1:rows) {
    # Setup the current node
    currentNode1 <- c()
    currentNode2 <- c()

    # Get the info related to the current node
    currentNode1$entryID <- pathwayGraph[nodeIdx, 'node1'] # node1: entry ID
    currentNode2$entryID <- pathwayGraph[nodeIdx, 'node2'] # node2: entry ID
    currentNode1$ec <- pathway$nodes[pathway$nodes$entryID==currentNode1$entryID,]$name # node1: EC
    currentNode2$ec <- pathway$nodes[pathway$nodes$entryID==currentNode2$entryID,]$name # node2: EC
    currentNode1$reaction <- pathway$nodes[pathway$nodes$entryID==currentNode1$entryID,]$reaction # node1: reaction
    currentNode2$reaction <- pathway$nodes[pathway$nodes$entryID==currentNode2$entryID,]$reaction # node2: reaction
    currentNode1$x <- pathway$nodes[pathway$nodes$entryID==currentNode1$entryID,]$x # node1: graph x
    currentNode2$x <- pathway$nodes[pathway$nodes$entryID==currentNode2$entryID,]$x # node2: graph x
    currentNode1$y <- pathway$nodes[pathway$nodes$entryID==currentNode1$entryID,]$y # node1: graph y
    currentNode2$y <- pathway$nodes[pathway$nodes$entryID==currentNode2$entryID,]$y # node2: graph y

    #*********************************************##
    # Convert the entry IDs into the dictionary ID #
    #*********************************************##

    # Find the current node into the dictionary
    if (is.na(currentNode1$reaction)) {
      dictId1 <- dictionary[dictionary$x == currentNode1$x & dictionary$y == currentNode1$y, ]$id
    } else {
      dictId1 <- dictionary[dictionary$x == currentNode1$x & dictionary$y == currentNode1$y & dictionary$reaction == currentNode1$reaction, ]$id
    }

    if (is.na(currentNode2$reaction)) {
      dictId2 <- dictionary[dictionary$x == currentNode2$x & dictionary$y == currentNode2$y, ]$id
    } else {
      dictId2 <- dictionary[dictionary$x == currentNode2$x & dictionary$y == currentNode2$y & dictionary$reaction == currentNode2$reaction, ]$id
    }

    # Check if the dictionary contains the node, if the didctionary ID is empty, it means that
    # the node refers to a pathway connection (e.g:path:00020) or it is a compound
    if (is.null(dictId1) | isempty(dictId1)) {
      # Remove the current row
      pathwayGraph[nodeIdx,] <- NA
      next()
    } else {
      pathwayGraph[nodeIdx, 'node1'] <- dictId1
    }

    if (is.null(dictId2) | isempty(dictId2)) {
      # Remove the current row
      pathwayGraph[nodeIdx,] <- NA
      next()
    } else {
      pathwayGraph[nodeIdx, 'node2'] <- dictId2
    }

    #***************************##
    # Set graph other attributes #
    #***************************##

    pathwayGraph[nodeIdx, 'ec1'] <- currentNode1$ec
    pathwayGraph[nodeIdx, 'ec2'] <- currentNode2$ec
    pathwayGraph[nodeIdx, 'entryID1'] <- currentNode1$entryID
    pathwayGraph[nodeIdx, 'entryID2'] <- currentNode2$entryID
    pathwayGraph[nodeIdx, 'reaction1'] <- currentNode1$reaction
    pathwayGraph[nodeIdx, 'reaction2'] <- currentNode2$reaction
    pathwayGraph[nodeIdx, 'x1'] <- currentNode1$x
    pathwayGraph[nodeIdx, 'x2'] <- currentNode2$x
    pathwayGraph[nodeIdx, 'y1'] <- currentNode1$y
    pathwayGraph[nodeIdx, 'y2'] <- currentNode2$y
  }

  # Remove the NA rows
  if ( (any(is.na(pathwayGraph$reaction1)) | any(is.na(pathwayGraph$reaction2))) &
       (any(is.na(pathwayGraph$node1)) | any(is.na(pathwayGraph$node2))) ) {
    # Just the complete cases
    pathwayGraph <- pathwayGraph[complete.cases(pathwayGraph), ]

    # Reindex the dataframe
    if (nrow(pathwayGraph) > 0) {
      rownames(pathwayGraph) <- 1:nrow(pathwayGraph)
    } else {
      pathwayGraph <- NULL
    }
  }

  # Return the pathway graph
  return(pathwayGraph)
}
