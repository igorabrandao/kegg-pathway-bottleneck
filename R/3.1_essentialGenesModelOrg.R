#*************************************************************************#
# Pipeline to perform the gene essentiality evaluation of model organisms #
#*************************************************************************#

# ---- IMPORT SECTION ----

# 3.0_essentialGenesModelOrg.R #

#' This is the pipeline script to perform
#' the gene essentiality evaluation of model organisms
#'
#' @author
#' Igor Brandão

# Import the necessary libraries
library(biomaRt)

#*******************************************************************************************#

# ---- SETTINGS SECTION ----

#*************************#
# Pipeline basic settings #
#*************************#

# Import the graphLoader functions
files.sources = NULL
files.sources[1] = paste0("./R/functions", "/", "graphFunctions.R")
files.sources[2] = paste0("./R/functions", "/", "kgmlFunctions.R")
files.sources[3] = paste0("./R/functions", "/", "helperFunctions.R")
files.sources[4] = paste0("./R/functions", "/", "graphPrintFunctions.R")
sapply(files.sources, source)

# Load the pathways by organisms data
essentialGenesModelOrg <- get(load(paste0("./dictionaries", "/", "essentialGenesModelOrg.RData")))

# Set the list of model organisms according to KEGG nomenclature
# Mus musculus => mmu
# Drosophila melanogaster => dme
# Saccharomyces cerevisiae => sce
# Caenorhabditis elegans => cel
modelOrgList = c("mmu", "dme", "sce", "cel")
biomaRtOrgs = c("mmusculus_gene_ensembl", "dmelanogaster_gene_ensembl", "scerevisiae_gene_ensembl", "celegans_gene_ensembl")

#*******************************************************************************************#

# ---- FUNCTIONS SECTION ----

#*************************#
# Pipeline main functions #
#*************************#

#' Parse the organism KGML files and generate the list of essential geness {by its lethality status}
#'
#' @param orgList_ List of KEGG organism code
#'
#' @return This function does not return nothing, just export .csv files.
#'
#' @examples
#' \dontrun{
#' generateOrgDataFromKGML(c("mmu", "dme", "sce", "cel"))
#' }
#'
#' @author
#' Igor Brandão
#'
generateOrgGeneList <- function(orgList_) {

  # Load the dicionaty
  dictionary <- read.csv(file='./output/pathwaysDictionary/dictionary.csv', header=TRUE, sep=",", stringsAsFactors=FALSE)

  # Loop 01: Run through all organism in list
  lapply(orgList_, function(currentOrg) {

    # Get the organism pathway list of files
    folder = paste0("./output/kgml/", currentOrg, "/")
    kgml_list <- list.files(path=folder, pattern='*.xml')

    # Define the number of available pathways
    available_pathways <- length(kgml_list)
    kgml_index <- 1

    # Check if the folder contains files
    if (is.null(kgml_list) | length(kgml_list) == 0) {
      # Status message
      printMessage(paste0("There aren't available pathways for ", currentOrg, "..."))
      return(FALSE)
    }

    # Loop 02: Run through all KGMLs in list
    lapply(kgml_list, function(file) {

      # Load the dataframe
      current_kgml <- KGML2Dataframe(paste0(folder, file))

      # Get the pathway code
      pathway_code <- onlyNumber(file)

      # Status message
      printMessage(paste0("GATHERING ", currentOrg, " PATHWAY ", pathway_code, " DATA [", kgml_index, " OF ", available_pathways, "]"))

      # Skip the global metabolic map
      if (pathway_code %in% c('01100')) {
        return(FALSE)
      }

      tryCatch({
        # Convert the pathway data into a graph
        pathwayGraph <- KGML2Graph(paste0(folder, file), replaceOrg=TRUE, orgToReplace=currentOrg)

        #*************************##
        # Prepare the pathway data #
        #*************************##

        # Create the pathwayData dataFrame
        pathwayData <- current_kgml$nodes

        # Remove the pathwayData unnecessary columns
        pathwayData <- pathwayData[,!(names(pathwayData) %in% c('link', 'component', 'map', 'fgcolor'))]

        # Add the default columns
        pathwayData$entrez <- NA
        pathwayData$dictID <- NA
        pathwayData$reaction_type <- NA
        pathwayData$org <- currentOrg
        pathwayData$pathway <- pathway_code
        pathwayData$is_bottleneck <- 0

        pathwayData$freq <- 0

        pathwayData$betweenness <- NA
        pathwayData$connectivity <- NA
        pathwayData$triangles <- NA
        pathwayData$clusteringCoef <- NA
        pathwayData$closenessCoef <- NA
        pathwayData$community <- NA
        pathwayData$eigenvectorScore <- NA
        pathwayData$eccentricity <- NA
        pathwayData$radius <- NA
        pathwayData$diameter <- NA
        pathwayData$degree <- NA
        pathwayData$authorityScore <- NA
        pathwayData$hubScore <- NA

        pathwayData$bottleneck_classification <- NA

        # Remove unnecessary data from pathway data/graph
        pathwayGraph <- removeNoise(pathwayGraph)

        pathwayData <- pathwayData[!grepl("^path:", pathwayData$name),]
        pathwayData <- pathwayData[!grepl("^map:", pathwayData$name),]
        pathwayData <- pathwayData[!grepl("^cpd:", pathwayData$name),]
        pathwayData <- pathwayData[!grepl("^gl:", pathwayData$name),]
        pathwayData <- pathwayData[!grepl("^ko:", pathwayData$name),]

        # Assign the reaction type to each node
        for (idx in 1:nrow(pathwayData)) {
          for (idx2 in 1:length(current_kgml$reactions$name)) {
            if (is.null(current_kgml$reactions$name[idx2]) | is.na(pathwayData[idx,]$reaction)) {
              next()
            }

            # Check the position of the current reaction in reactions list
            if (current_kgml$reactions$name[idx2] %in% pathwayData[idx,]$reaction) {
              pathwayData[idx,]$reaction_type <- current_kgml$reactions$type[idx2]
              break()
            }
          }
        }

        # Get the graph properties
        graphProperties <- getGraphProperties(pathwayGraph)

        # Perform the graph bottleneck calculation
        iGraph <- igraph::graph_from_data_frame(pathwayGraph, directed = FALSE)
        graphBottleneck <- igraph::as_ids(getGraphBottleneck(iGraph, FALSE))

        # Assign the bottlenecks
        pathwayData$is_bottleneck[which(pathwayData$name %in% graphBottleneck)] <- 1

        # Assign the graph properties to each node
        for (idx in 1:nrow(pathwayData)) {
          # Prepare the node name to be compared
          pattern <- gsub('/', '|', pathwayData$name[idx])
          pattern <- gsub(" ","", pattern)

          # Find which rows in graphProperties match with pathwayData
          rowsToMerge <- which(grepl(pattern, graphProperties$node))[1]

          pathwayData[idx,]$betweenness <- graphProperties[rowsToMerge,]$betweenness
          pathwayData[idx,]$connectivity <- graphProperties[rowsToMerge,]$connectivity
          pathwayData[idx,]$triangles <- graphProperties[rowsToMerge,]$triangles
          pathwayData[idx,]$clusteringCoef <- graphProperties[rowsToMerge,]$clusteringCoef
          pathwayData[idx,]$closenessCoef <- graphProperties[rowsToMerge,]$closenessCoef
          pathwayData[idx,]$community <- graphProperties[rowsToMerge,]$community
          pathwayData[idx,]$eigenvectorScore <- graphProperties[rowsToMerge,]$eigenvectorScore
          pathwayData[idx,]$eccentricity <- graphProperties[rowsToMerge,]$eccentricity
          pathwayData[idx,]$radius <- graphProperties[rowsToMerge,]$radius
          pathwayData[idx,]$diameter <- graphProperties[rowsToMerge,]$diameter
          pathwayData[idx,]$degree <- graphProperties[rowsToMerge,]$degree
          pathwayData[idx,]$authorityScore <- graphProperties[rowsToMerge,]$authorityScore
          pathwayData[idx,]$hubScore <- graphProperties[rowsToMerge,]$hubScore

          # Assign the bottlenecks for kegg orthology (ko)
          if (strcmp(currentOrg, 'ko')) {
            # Check if at least one bottleneck was found
            if (length(graphBottleneck) > 0) {
              if (grepl(pattern, graphBottleneck)) {
                pathwayData[idx,]$is_bottleneck <- 1
              }
            }
          }
        }

        # Apply the node classification
        pathwayData$bottleneck_classification <- classifyBottleneck(pathwayData)$bottleneck_classification

        #***************************##
        # Count the enzyme frequency #
        #***************************##

        # Enzyme color identification
        enzyme_present_color_org <- '#BFFFBF'
        enzyme_present_color_ref <- '#BFBFFF'
        enzyme_missing_color <- '#FFFFFF'

        # Increment the enzyme frequency
        pathwayData[pathwayData$bgcolor == enzyme_present_color_org || pathwayData$bgcolor == enzyme_present_color_ref,]$freq <-
          pathwayData[pathwayData$bgcolor == enzyme_present_color_org || pathwayData$bgcolor == enzyme_present_color_ref,]$freq + 1

        #*****************************##
        # Get the enzyme dictionary ID #
        #*****************************##

        if (!strcmp(currentOrg, 'ec')) {
          pathwayData$entrez <- pathwayData$name
        }

        # Assign the dictionary ID to each node
        rows <- nrow(pathwayData)
        for (nodeIdx in 1:rows) {
          # Setup the current node
          currentNode <- pathwayData[nodeIdx,]

          # Find the current node into the dictionary
          if (is.na(currentNode$reaction)) {
            dictId <- dictionary[(dictionary$x == currentNode$x) & (dictionary$y == currentNode$y), ]
          } else {
            dictId <- dictionary[(dictionary$x == currentNode$x) & (dictionary$y == currentNode$y) &
                                   (dictionary$reaction == currentNode$reaction), ]
          }

          # Remove NAs
          dictId <- dictId[!is.na(dictId$ec),]

          # Check if the dictionary contains the node, if the didctionary ID is empty, it means that
          # the node refers to a pathway connection (e.g:path:00020) or it is a compound
          if (is.null(dictId) | isempty(dictId) | nrow(dictId) == 0) {
            pathwayData[nodeIdx, 'name'] <- NA
            next()
          } else {
            pathwayData[nodeIdx, 'dictID'] <- dictId$id

            if (!strcmp(currentOrg, 'ec')) {
              pathwayData[nodeIdx, 'name'] <- dictId$ec
            }
          }
        }

        #***********************************************************#
        # Export the organism specific data for the current pathway #
        #***********************************************************#

        # Status message
        printMessage(paste0("EXPORTING ", currentOrg, " PATHWAY ", pathway_code, " DATA [", kgml_index, " OF ", available_pathways, "]"))

        # Create the canonical analysis folder
        if (!dir.exists(file.path('./output/essentialGenes/'))) {
          dir.create(file.path(paste0('./output/essentialGenes/')), showWarnings = FALSE, mode = "0775")
        }

        # Create the org folder
        if (!dir.exists(file.path(paste0('./output/essentialGenes/', currentOrg)))) {
          dir.create(file.path(paste0('./output/essentialGenes/', currentOrg)), showWarnings = FALSE, mode = "0775")
        }

        # Save the org pathwayData into its folder
        if (dir.exists(paste0('./output/essentialGenes/', currentOrg))) {
          write.csv(pathwayData, file=paste0('./output/essentialGenes/', currentOrg, '/', pathway_code, '.csv'), row.names=FALSE)
        }

        #*********************#
        # Increment the index #
        #*********************#
        kgml_index <<- kgml_index + 1

      }, error=function(e) {
        printMessage(e)

        # Save the log file
        printLog(toString(e), file_=paste0('generateOrgGeneList', pathway_code))

        return(FALSE)
      })
    }) # End of Loop 02

  }) # End of Loop 01

  # Function finished with success
  return(TRUE)
}

#' Convert the ensembl peptide ID to ensembl ID
#'
#' @param biomaRtOrgs_ List of biomart organism datasets
#'
#' @return This function does not return nothing, just export .RData file.
#'
#' @examples
#' \dontrun{
#' generateOrgDataFromKGML(c("mmusculus_gene_ensembl", "dmelanogaster_gene_ensembl", "scerevisiae_gene_ensembl", "celegans_gene_ensembl"))
#' }
#'
#' @author
#' Igor Brandão
#'
convertPeptideToEnsembl <- function(biomaRtOrgs_) {

  # Load the pathways by organisms data
  essentialGenesModelOrg <- get(load(paste0("./dictionaries", "/", "essentialGenesModelOrg.RData")))

  # Add the ensembl_gene_id column
  essentialGenesModelOrg$ensembl_gene_id <- NA
  essentialGenesModelOrg$entrezgene_id <- NA
  essentialGenesModelOrg$entrezgene_accession <- NA

  # Loop counters
  index <- 1
  available_orgs <- length(biomaRtOrgs_)

  # Loop 01: Run through all available pathways kgml
  lapply(biomaRtOrgs_, function(biomaRtOrg) {
    # Status message
    printMessage(paste0("CONVERTING ", biomaRtOrg, " DATASET [", index, " OF ", available_orgs, "]"))

    # Retrieve the ensembl info
    ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = biomaRtOrg)

    # Filter the ensembl ID
    ids <- getBM(filters = "ensembl_peptide_id",
             attributes = c("ensembl_peptide_id", "ensembl_gene_id", 'entrezgene_id', 'entrezgene_accession'),
             values = essentialGenesModelOrg$ensembl_peptide_id, mart = ensembl)

    # Remove NA cases
    ids[ids == ""] <- NA
    ids <- na.omit(ids)

    # Apply the ensembl IDs and the entrez_id by the peptide ID
    for (idx in 1:nrow(ids)) {
      # Ensembl ID
      essentialGenesModelOrg[essentialGenesModelOrg$ensembl_peptide_id == ids[idx,]$ensembl_peptide_id,]$ensembl_gene_id <<-
        ids[idx,]$ensembl_gene_id

      # Entrez ID
      essentialGenesModelOrg[essentialGenesModelOrg$ensembl_peptide_id == ids[idx,]$ensembl_peptide_id,]$entrezgene_id <<-
        ids[idx,]$entrezgene_id

      # Entrez accession
      essentialGenesModelOrg[essentialGenesModelOrg$ensembl_peptide_id == ids[idx,]$ensembl_peptide_id,]$entrezgene_accession <<-
        ids[idx,]$entrezgene_accession
    }

    # Increment the index
    index <<- index + 1

  }) # End of Loop 01

  # Update the dataSet
  save(essentialGenesModelOrg, file = "./dictionaries/essentialGenesModelOrg2.RData", compress = "xz")
}

#' Function to aggregate all org datasets into one
#'
#' @param orgList_ List of KEGG organism code
#'
#' @return This function does not return nothing, just export .csv files.
#'
#' @examples
#' \dontrun{
#' joinOrgDatasets(c("mmu", "dme", "sce", "cel"))
#' }
#'
#' @author
#' Igor Brandão
#'
joinOrgDatasets <- function(orgList_) {

  # Loop counters
  index <- 1
  available_orgs <- length(orgList_)

  # Loop 01: Run through all organism in list
  lapply(orgList_, function(currentOrg) {

    # Status message
    printMessage(paste0("MERGING ", currentOrg, " datasets [", index, " OF ", available_orgs, "]"))

    # Get the list of files
    folder = paste0("./output/essentialGenes/", currentOrg, "/")
    file_list <- grep(list.files(path=folder), pattern='*.csv', value=T)

    # Load all csv files at once
    big.list.of.data.frames <- lapply(file_list, function(item) {
      read.csv(file=paste0(folder, item), header=TRUE, sep=",", stringsAsFactors=FALSE)
    })

    # Combine multiple data frames in one
    allNodes <- do.call(rbind, big.list.of.data.frames)

    # Export allNodes
    write.csv(allNodes, file=paste0('./output/essentialGenes/', currentOrg, '.csv'))

    # Remove temporaly variables
    rm(big.list.of.data.frames, allNodes, folder, file_list)

    # Increment the index
    index <<- index + 1
  })
}

#' Function to match the organisms genes by its entrez ID with the ensembl ID
#' from gene essential list
#'
#' @param orgList_ List of KEGG organism code
#'
#' @return This function does not return nothing, just export .csv files.
#'
#' @examples
#' \dontrun{
#' matchOrgEntrezWithEnsembl(c("mmu", "dme", "sce", "cel"))
#' }
#'
#' @author
#' Igor Brandão
#'
matchOrgEntrezWithEnsembl <- function(orgList_) {

  # Loop counters
  index <- 1
  available_orgs <- length(orgList_)

  # Loop 01: Run through all organism in list
  lapply(orgList_, function(currentOrg) {

    # Status message
    printMessage(paste0("MATCHING THE ", currentOrg, " GENE LIST WITH THE LETHALITY LIST [", index, " OF ", available_orgs, "]"))

    # Get the list of files
    orgGenes <- read.csv(paste0("./output/essentialGenes/", currentOrg, ".csv"), header=TRUE, sep=",", stringsAsFactors=FALSE)

    # Extra fields to orgGenes df
    orgGenes$lethal_nonlethal <- NA
    orgGenes$entrezgene_accession <- NA

    # Loop 02: Run through all organism gene list
    idx = 1
    apply(orgGenes, 1, function(currentGene) {
      # Current lethality status
      lethalityStatus <- c()
      accessionIDs <- c()

      # Generate a temporary entrez list for the current gene
      currentEntrezList <- orgGenes[idx,]$entrez

      # Remove the org code from entrez
      currentEntrezList <- gsub(paste0(currentOrg, ':'), '', currentEntrezList)

      # Split the entrez list
      currentEntrezList <- unlist(str_split(currentEntrezList, " / "))

      # Loop 03: Run through all entrez in the current gene
      for (idx2 in 1:length(currentEntrezList)) {
        # Retrieve the lethality status for the current entrez
        currentEntrez <- essentialGenesModelOrg[which(essentialGenesModelOrg$entrezgene_id == currentEntrezList[idx2]),]

        # Append the status to the current gene status
        lethalityStatus <- c(lethalityStatus, currentEntrez$lethal_nonlethal)
        accessionIDs <- c(accessionIDs, currentEntrez$entrezgene_accession)
      }

      #*******************************************************************************************************************************#
      # Verify the lethalityStatus, if at least one entrez is classified as lethal the whole group is lethal                          #
      #                                                                                                                               #
      # This situation reflects the case of an enzyme (EC) that has several subunits (several entrez id or mgids). The hypothesis     #
      # is that the loss of a subunit probably affects the structure of the enzyme as a whole and, therefore, if you find a subunit   #
      # that generates lethality, this alteration is also likely to be lethal as a whole                                              #
      #*******************************************************************************************************************************#
      orgGenes[idx,]$lethal_nonlethal <<- ifelse('lethal' %in% lethalityStatus, 'lethal', 'nonlethal')
      orgGenes[idx,]$entrezgene_accession <<- toString(accessionIDs)

      # Increment the index
      idx <<- idx + 1
    })

    # Export the organisms genes with the lethality status
    write.csv(orgGenes, file=paste0('./output/essentialGenes/', currentOrg, '.csv'))

    # Remove temporaly variables
    rm(orgGenes)

    # Increment the index
    index <<- index + 1
  })
}

#*******************************************************************************************#

# ---- PIPELINE SECTION ----

#***************#
# Pipeline flow #
#***************#

#******************************************************#
# Step 1: Convert the genes from peptide ID to ensembl #
# Warning: to heavy, run just once                     #
#******************************************************#
convertPeptideToEnsembl(biomaRtOrgs)

#*************************************************************************#
# Step 2: Generate the organisms list of genes for all available pathways #
#*************************************************************************#
generateOrgGeneList(modelOrgList)

#*****************************************************************#
# Step 3: Join organisms list of genes for all available pathways #
#*****************************************************************#
joinOrgDatasets(modelOrgList)

#******************************************************************************************#
# Step 4: Match the organisms genes with the list of lethal genes (essentialGenesModelOrg) #
#******************************************************************************************#
matchOrgEntrezWithEnsembl(modelOrgList)

#*******************************************************#
# Step 5: Perform plots to explore the lethality metric #
#*******************************************************#
