###########################################
# Functions to load all organism pathways #
###########################################

# dbLoader ####

#' This is an automation script to insert KEGG data into mysql DB
#'
#' @author
#' Diego Morais & Igor Brandão

# Import the necessary libraries
library(KEGGREST)

######################
# Metabolic pathways #
######################

# Get all the metabolic pathways code
allPathways <- names(keggList("pathway"))
allPathways <- gsub("path:", "", allPathways)

# save(allPathways, file = "allPathways.RData", compress = "xz")

#############
# Organisms #
#############

# Get all organisms code
org <- as.data.frame(keggList("organism"), stringsAsFactors = FALSE)
org <- unique(org$organism)

# save(org, file = "allOrganisms.RData", compress = "xz")

##################################
# Metabolic pathways x Organisms #
##################################

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

organism2pathway <- lapply(organism2pathway, function(i) {
  i[[1]]
})
names(organism2pathway) <- org
organism2pathway[sapply(organism2pathway, is.null)] <- NULL
# save(organism2pathway, file = "organism2pathway.RData", compress = "xz")
