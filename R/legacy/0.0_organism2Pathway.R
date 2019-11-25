#*******************************************************************#
# Functions to load metabolic pathways and organisms data from KEGG #
#*******************************************************************#

# ---- PATHWAYS RETRIEVING SECTION ----

# 0.0_organism2Pathway.R #

#' This is an automation script to retrieve metabolic pathways and organism from KEGG
#'
#' @author
#' Diego Morais & Igor Brandão

# Import the necessary libraries
library(KEGGREST)

# ---- SETTINGS SECTION ----

#*************************#
# Pipeline basic settings #
#*************************#

# Import the graphLoader functions
files.sources = NULL
files.sources[1] = paste0("./R/functions", "/", "helperFunctions.R")
sapply(files.sources, source)

#*******************************************************************************************#

#********************#
# Metabolic pathways #
#********************#

# Get all the metabolic pathways code
pathwayList <- names(keggList("pathway"))
pathwayList <- gsub("path:", "", pathwayList)

# Get just the numeric part
pathwayList <- gsub("[^0-9.]", "", pathwayList)

# Convert into dataFrame object
pathwayList <- data.frame(pathway = pathwayList, stringsAsFactors = FALSE)

# Save the pathwayList
save(pathwayList, file = "./dictionaries/pathwayList.RData", compress = "xz")

#***********#
# Organisms #
#***********#

# Get all organisms code
organismList <- as.data.frame(keggList("organism"), stringsAsFactors = FALSE)
save(organismList, file = "./dictionaries/organismList.RData", compress = "xz")

# Get just the organism identification
org <- unique(organismList$organism)

#****************************#
# Metabolic pathways details #
#****************************#

# Load all pathway detail at once
pathwayDetail <- lapply(unlist(pathwayList), function(pathwayCode) {
  # Status message
  printMessage(paste0("GETTING PATHWAY ", pathwayCode, " DETAILS"))

  tryCatch({
    # Retrieve the pathway detail from KEEG API
    keggGet(paste0('map', pathwayCode))
  }, error=function(e) {
    printMessage(paste0("Pathway ", pathwayCode, " could no be found. Skipping it..."))
  })
})

# Count the number of pathways that have information (valid pathways)
length(which(lengths(pathwayDetail, use.names = TRUE)==1))

# Remove one dimension from the list
pathwayDetail <- lapply(pathwayDetail, function(item) {
  item[[1]]
})

# Save the pathwayDetail (list)
save(pathwayDetail, file = "./dictionaries/pathwayDetail.RData", compress = "xz")

#********************************#
# Metabolic pathways x Organisms #
#********************************#

# Get all metabolic pathways by organism
organism2pathway <- lapply(org, function(i) {
  message(i)
  flagError <- FALSE
  tryCatch({
    res <- unique(keggLink("pathway", i))
  }, warning = function(w) {
    # warning-handler-code
  }, error = function(e) {
    flagError <<- TRUE
  }, finally = {
    # cleanup-code
  })
  if (flagError) {
    return(NULL)
  }
  res <- gsub(paste0("path:", i), "", res)
  res <- list(unlist(res))
  names(res) <- i
  return(res)
})

# Rename the listitems
names(organism2pathway) <- org

# Remove one dimension from the list
organism2pathway <- lapply(organism2pathway, function(item) {
  item[[1]]
})

# Remove the NULL cases
#organism2pathway[sapply(organism2pathway, is.null)] <- NULL

# Save the organism2pathway (list)
save(organism2pathway, file = "./dictionaries/organism2pathway.RData", compress = "xz")
