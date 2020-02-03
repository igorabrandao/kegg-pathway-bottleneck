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

# ---- Pathways ----

# Get all the metabolic pathways code
pathwayList <- names(keggList("pathway"))
pathwayList <- gsub("path:", "", pathwayList)

# Get just the numeric part
pathwayList <- gsub("[^0-9.]", "", pathwayList)

# Convert into dataFrame object
pathwayList <- data.frame(pathway = pathwayList, stringsAsFactors = FALSE)

# Save the pathwayList
save(pathwayList, file = "./dictionaries/pathwayList.RData", compress = "xz")

#*******************************************************************************************#

# ---- Organisms ----

# Get all organisms code
organismList <- as.data.frame(keggList("organism"), stringsAsFactors = FALSE)
save(organismList, file = "./dictionaries/organismList.RData", compress = "xz")

# Get just the organism identification
org <- unique(organismList$organism)

#*******************************************************************************************#

# ---- Pathways details ----

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

#*******************************************************************************************#

# ---- Pathways x Organisms ----

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

# Rename the list items
names(organism2pathway) <- org

# Remove one dimension from the list
organism2pathway <- lapply(organism2pathway, function(item) {
  item[[1]]
})

# Remove the NULL cases
#organism2pathway[sapply(organism2pathway, is.null)] <- NULL

# Save the organism2pathway (list)
save(organism2pathway, file = "./dictionaries/organism2pathway.RData", compress = "xz")

#*******************************************************************************************#

# ---- Pathways x Phylogeny ----

# Load the pathway and org dictionary (if exists)
organismList <- get(load(paste0("./dictionaries", "/", "organismList.RData")))
phylogenyPerPathway <- get(load(paste0("./dictionaries", "/", "pathwayList.RData")))
organism2pathway <- get(load(paste0("./dictionaries", "/", "organism2pathway.RData")))

# Split the phylogeny from the dictionary
uniquePhylogeny <- data.frame(do.call('rbind', strsplit(as.character(organismList$phylogeny),';',fixed=TRUE)))
phylogeny1 <- data.frame(lapply(unique(uniquePhylogeny$X1), as.character), stringsAsFactors=FALSE)
phylogeny2 <- data.frame(lapply(unique(uniquePhylogeny$X2), as.character), stringsAsFactors=FALSE)
phylogeny <- cbind(phylogeny1, phylogeny2)

# Add phylogeny columns to the pathway list based on orgs classification
for (idx in 1:ncol(phylogeny)) {
  phylogenyPerPathway[phylogeny[1,idx]] <- 0
}

# Count the organism phylogeny per pathway
sapply(organismList$organism, function(org_) {
  # Get the pathway list from the specif org
  orgIdx <- which(names(organism2pathway) == org_)
  orgPathwayList <- unlist(organism2pathway[orgIdx])

  # Get the philogeny of a specific org
  orgPhylogeny <-
    organismList[organismList$organism == org_, ]$phylogeny

  # Finally perform the phylogeny count
  phylogenies <-
    unlist(strsplit(as.character(orgPhylogeny), ';', fixed = TRUE))

  for (orgPathwayListIdx in 1:length(orgPathwayList)) {
    # Add-up the count
    phylogenyPerPathway[phylogenyPerPathway$pathway == orgPathwayList[orgPathwayListIdx],
                        names(phylogenyPerPathway) %in% phylogenies] <<-
      phylogenyPerPathway[phylogenyPerPathway$pathway == orgPathwayList[orgPathwayListIdx],
                          names(phylogenyPerPathway) %in% phylogenies] + 1
  }
})
# end loop

# Save the phylogenyPerPathway dataFrame as CSV
write.csv(phylogenyPerPathway, file=paste0('./dictionaries/phylogenyPerPathway.csv'))

#*******************************************************************************************#
