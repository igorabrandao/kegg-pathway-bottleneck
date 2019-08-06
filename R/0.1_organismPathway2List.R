#*************************************************#
# Adapt the organism2pathway data to pathway list #
#*************************************************#

# ---- PATHWAYS CONVERSION TO LIST SECTION ----

# 0.1_organismPathway2List.R #

#' This is the pipeline script to perform
#' the conversion from organism2pathway to pathway list
#'
#' @author
#' Igor Brandão

#-------------------------------------------------------------------------------------------#

# Load the pathways by organisms data
organism2pathway <- get(load(paste0("./dictionaries", "/", "organism2pathway.RData")))

# Convert the list to dataFrame
pathwayList <- data.frame(pathway=unlist(organism2pathway), stringsAsFactors=FALSE)

# Remove the duplicates
pathwayList <-as.data.frame(pathwayList[!duplicated(pathwayList),], stringsAsFactors=FALSE)
names(pathwayList) <- "pathway"

# Export into a new file
save(pathwayList, file='./dictionaries/pathwayList.RData')
