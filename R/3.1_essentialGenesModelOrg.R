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
library(biomartr)

library(ggplot2)
library(svglite)
library(RColorBrewer)
library(plyr)

library(ggpubr)
library(gghighlight)
library(grid)
library(GGally)

library(DataExplorer)

#*******************************************************************************************#

# ---- SETTINGS SECTION ----

#*************************#
# Pipeline basic settings #
#*************************#

# Import the graphLoader functions
files.sources = NULL
files.sources[1] = paste0("./R/functions", "/", "graphFunctions.R")
files.sources[2] = paste0("./R/functions", "/", "kgmlFunctions.R")
files.sources[3] = paste0("./R/functions", "/", "statisticsHelper.R")
files.sources[4] = paste0("./R/functions", "/", "helperFunctions.R")
sapply(files.sources, source)

# Load the pathways by organisms data
#essentialGenesModelOrg <- get(load(paste0("./dictionaries", "/", "essentialGenesModelOrg.csv")))
essentialGenesModelOrg <- read.csv(file='./dictionaries/essentialGenesModelOrg.csv', header=TRUE, sep=",", stringsAsFactors=FALSE)

# Set the list of model organisms according to KEGG nomenclature
modelOrgList = unique(essentialGenesModelOrg$org)
modelOrgList = c("hsa", "spo", "dme", "pau", "eco", "mtv", "sey", "sao", "hin", "mpn") # Lista prioridade Clóvis
biomaRtOrgs = c("hsapiens_gene_ensembl", "", "dmelanogaster_gene_ensembl", "", "", "", "", "", "", "")

# List the biomart orgs
#ensembl=useMart("ensembl")
#biomaRt::listDatasets(ensembl)

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

#' Function to add the organism mart dataset name into the dictionary
#'
#' @param orgList_ List of KEGG organism name
#'
#' @return This function does not return nothing, just export .csv file.
#'
#' @examples
#' \dontrun{
#' generateOrgDataFromKGML(c("Homo sapiens", "Mus musculus"))
#' }
#'
#' @author
#' Igor Brandão
#'
addOrgMartDatasetName <- function(orgList_) {

  # Load the pathways by organisms data
  essentialGenesModelOrg <- read.csv(file='./dictionaries/essentialGenesModelOrg.csv', header=TRUE, sep=",", stringsAsFactors=FALSE)

  # Add the ensembl_gene_id column
  essentialGenesModelOrg$biomartDataset <- NA

  # Loop counters
  index <- 1
  available_orgs <- length(orgList_)

  # Loop 01: Run through all available pathways kgml
  lapply(orgList_, function(org) {
    # Status message
    printMessage(paste0("GETTING ", org, " BIOMART DATASET NAME [", index, " OF ", available_orgs, "]"))

    orgAttributes <- NULL

    tryCatch({
      # Retrieve the org attributes
      orgAttributes <- biomartr::organismBM(organism = org)
    }, error=function(e) {
      # Save the log file
      printLog(toString(e), file_=paste0('addOrgMartDatasetName', org))
    })

    # Increment the index
    index <<- index + 1

    if (is.null(orgAttributes) | length(orgAttributes) == 0) {
      return(FALSE)
    } else {
      datasetName <- grep(orgAttributes$dataset, pattern='*_gene_ensembl', value=T)

      if (is.null(datasetName) | length(datasetName) == 0) {
        return(FALSE)
      } else {
        printMessage(datasetName)

        # Add the biomart dataset name to the essential genes list
        essentialGenesModelOrg[essentialGenesModelOrg$sciName == org,]$biomartDataset <<- datasetName
      }
    }
  }) # End of Loop 01

  # Update the dataSet
  write.csv(essentialGenesModelOrg, file=paste0('./dictionaries/essentialGenesModelOrg2.csv'), row.names = F)
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

    # Clear the unnecessary fields
    #Precisamos do Nome do organismo, codigo kegg, pathway, EnsGenID e is_ap
    allNodes <- allNodes[,c('org', 'pathway', 'name', 'entrez', 'is_bottleneck', 'reaction_type')]

    # Export allNodes
    write.csv(allNodes, file=paste0('./output/essentialGenes/', currentOrg, '.csv'), row.names = F)

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
    orgGenes$ensembl_gene_id <- NA

    # Loop 02: Run through all organism gene list
    idx = 1
    apply(orgGenes, 1, function(currentGene) {
      # Current lethality status
      lethalityStatus <- c()
      accessionIDs <- c()
      ensemblIDs <- c()

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
        ensemblIDs <- c(ensemblIDs, currentEntrez$ensembl_gene_id)
      }

      #*******************************************************************************************************************************#
      # Verify the lethalityStatus, if at least one entrez is classified as lethal the whole group is lethal                          #
      #                                                                                                                               #
      # This situation reflects the case of an enzyme (EC) that has several subunits (several entrez id or mgids). The hypothesis     #
      # is that the loss of a subunit probably affects the structure of the enzyme as a whole and, therefore, if you find a subunit   #
      # that generates lethality, this alteration is also likely to be lethal as a whole                                              #
      #*******************************************************************************************************************************#
      #  E, DE (?), ES (?), D (?) - essential
      #*******************************************************************************************************************************#
      #  NE - nonessential
      #*******************************************************************************************************************************#
      # Condicionais:
      # growth defective (GD)
      # growth advantaged (GA)
      # F: fitness in vitro
      # E-infection: required for single infection to lung of mice
      # F-infection: Potential intermediate attenuation during single infection
      # E-co-infection: required for co-infection to lung of mice
      # F-co-infection: Potential intermediate attenuation during co-infection
      #*******************************************************************************************************************************#
      # Desconhecido:
      # ND, U
      #*******************************************************************************************************************************#
      # Incerto
      # S
      #*******************************************************************************************************************************#

      orgGenes[idx,]$lethal_nonlethal <<- ifelse('lethal' %in% lethalityStatus, 'lethal', 'nonlethal')
      orgGenes[idx,]$entrezgene_accession <<- toString(accessionIDs)
      orgGenes[idx,]$ensembl_gene_id <<- toString(ensemblIDs)

      # Increment the index
      idx <<- idx + 1
    })

    # Export the organisms genes with the lethality status
    write.csv(orgGenes, file=paste0('./output/essentialGenes/', currentOrg, '.csv'), row.names = F)

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
addOrgMartDatasetName(unique(essentialGenesModelOrg$sciName))

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

# ---- PLOT SECTION ----

# Get the list of files
folder = paste0("./output/essentialGenes/")
file_list <- grep(list.files(path=folder), pattern='*.csv', value=T)

# Load all csv files at once
big.list.of.data.frames <- lapply(file_list, function(item) {
  read.csv(file=paste0(folder, item), header=TRUE, sep=",", stringsAsFactors=FALSE)
})

# Combine multiple data frames in one
orgGenes <- do.call(rbind, big.list.of.data.frames)
rm(big.list.of.data.frames)

# For now we'll use just the MMU data since its better annotated
# TODO: Remove the next line
orgGenes <- orgGenes[orgGenes$org == 'mmu',]

write.csv(orgGenes, file=paste0('./output/essentialGenes/mmuGenesList.csv'))

# OR

# Load the gene list dataSet with the reaction and lethality classification
orgGenes <- read.csv(file='./output/essentialGenes/mmuGenesWithClassification.csv', header=TRUE, sep=",", stringsAsFactors=FALSE)

#**************#
# Data summary #
#**************#

# A more advanced and complete way to see the structure of our dataset
str(orgGenes)

# Summarize the data
print(summary(orgGenes))
print(summary(as.factor(orgGenes$reaction_type)))
print(summary(as.factor(orgGenes$is_bottleneck)))
print(summary(as.factor(orgGenes$lethal_nonlethal)))

#**************#
# Summary plot #
#**************#

# Define which dataSet column will be displayed into the plot
columns = c("reaction_type", "is_bottleneck", "betweenness", "degree", "lethal_nonlethal")
columnLabels = c("Reaction type", "AP status", "Betweenness", "Degree", "Lethality status")
plotTitle = paste0("Data summary")

# ---- plot1 ----

# Drawing a scatterplot matrix of freq, totalSpecies, percentage, and is_bottleneck using the pairs function
plot1 <- ggpairs(orgGenes, columns = columns, columnLabels = columnLabels, title = plotTitle,
                 mapping = aes(color = lethal_nonlethal),
                 lower = list(
                   continuous = "smooth",
                   combo = "facetdensity"
                 ), cardinality_threshold = 1000) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_bw()

# Recolor the matrix
for (i in 1:plot1$nrow) {
  for (j in 1:plot1$ncol) {
    plot1[i, j] <- plot1[i, j] +
      scale_fill_manual(values = c("#ED553B", "#173F5F")) +
      scale_color_manual(values = c("#ED553B", "#173F5F"))
  }
}

plot1

ggsave(paste0("./output/essentialGenes/plots/dataSummary.jpeg"), width = 30, height = 25, units = "cm")
ggsave(paste0("./output/essentialGenes/plots/dataSummary.svg"), width = 30, height = 25, units = "cm")

#*****************#
# Reaction status #
#*****************#

# Add the gene classification
orgGenes$gene_reaction_group <- NA

# Loop 01: Run through all organism gene list
# Warning: it's heavy!
idx = 1
apply(orgGenes, 1, function(currentGene) {
  # Apply the gene classification
  if (is.na(orgGenes[idx,]$reaction_type)) {
    orgGenes[idx,]$gene_reaction_group <<- 0
  }
  else if (orgGenes[idx,]$reaction_type == 'irreversible') {
    orgGenes[idx,]$gene_reaction_group <<- 1
  } else if (orgGenes[idx,]$reaction_type == 'reversible') {
    orgGenes[idx,]$gene_reaction_group <<- 2
  }

  # Increment the index
  idx <<- idx + 1
})

#******************#
# Lethality status #
#******************#

# Add the gene classification
orgGenes$gene_classification <- NA
orgGenes$gene_classification_group <- NA

# Loop 01: Run through all organism gene list
# Warning: it's heavy!
idx = 1
apply(orgGenes, 1, function(currentGene) {
  # Apply the gene classification
  if (orgGenes[idx,]$lethal_nonlethal == 'lethal' & orgGenes[idx,]$is_bottleneck == 1) {
    orgGenes[idx,]$gene_classification <<- 'AP lethal'
    orgGenes[idx,]$gene_classification_group <<- 1
  } else if (orgGenes[idx,]$lethal_nonlethal == 'lethal' & orgGenes[idx,]$is_bottleneck == 0) {
    orgGenes[idx,]$gene_classification <<- 'Non AP lethal'
    orgGenes[idx,]$gene_classification_group <<- 2
  } else if (orgGenes[idx,]$lethal_nonlethal == 'nonlethal' & orgGenes[idx,]$is_bottleneck == 1) {
    orgGenes[idx,]$gene_classification <<- 'AP non lethal'
    orgGenes[idx,]$gene_classification_group <<- 3
  } else {
    orgGenes[idx,]$gene_classification <<- 'Non AP non lethal'
    orgGenes[idx,]$gene_classification_group <<- 4
  }

  # Increment the index
  idx <<- idx + 1
})

# Rewrite the gene list dataSet with the reaction and lethality classification
write.csv(orgGenes, file=paste0('./output/essentialGenes/orgGenesWithClassification.csv'), row.names = F)

# Load the gene list dataSet with the reaction and lethality classification
orgGenes <- read.csv(file='./output/essentialGenes/orgGenesWithClassification.csv', header=TRUE, sep=",", stringsAsFactors=FALSE)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ---- plot2 ----

# Status message
printMessage("Lethality status plot")

# Generate the boxplot plot
plot2 <- ggplot(orgGenes, aes(fill=gene_classification, x=org), ymin = -Inf, ymax = Inf) +
  # Add the bars
  geom_bar(position="dodge", stat="count") +

  scale_y_continuous(breaks=seq(from = 0, to = nrow(orgGenes), by = 2500)) +
  scale_fill_manual(values = c("#ED553B", "#3CAEA3", "#a98600", "#173F5F")) +

  # Chart visual properties
  xlab("Genes classification") +
  ylab("Genes count") +
  ggtitle("") +
  guides(fill=guide_legend(title="")) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", color="black", size=26, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title.x = element_text(face="bold", color="black", size=20, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(face="bold", color="black", size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", color="black", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(face="bold", color="black", size=20),
        legend.title = element_text(face="bold", size=18),
        legend.text = element_text(size=16),
        legend.position = 'right') + labs(fill = "Classification")

plot2

ggsave(paste0("./output/essentialGenes/plots/genesClassificationOrg.jpeg"), width = 30, height = 25, units = "cm")
ggsave(paste0("./output/essentialGenes/plots/genesClassificationOrg.svg"), width = 30, height = 25, units = "cm")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#******************************#
# AP status x Lethality status #
#******************************#

# ---- plot3 ----

# Status message
printMessage("AP status x Lethality status")

# Run the Mann-Whitney Wilcoxon test
apStatus_lethalityStatus_test <- wilcox.test(orgGenes[orgGenes$gene_classification=='AP lethal',]$gene_classification_group,
            orgGenes[orgGenes$gene_classification=='AP non lethal',]$gene_classification_group)

# Generate the boxplot plot
plot3 <- ggplot(orgGenes[orgGenes$gene_classification=='AP lethal' | orgGenes$gene_classification=='AP non lethal',],
                         aes(fill=gene_classification, x=org), ymin = -Inf, ymax = Inf) +
  # Add the bars
  geom_bar(position="dodge", stat="count") +

  scale_y_continuous(breaks=seq(from = 0, to = nrow(orgGenes), by = 500)) +
  scale_fill_manual(values = c("#ED553B", "#3CAEA3")) +

  # Chart visual properties
  xlab("Genes classification") +
  ylab("Genes count") +
  ggtitle("") +
  guides(fill=guide_legend(title="")) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", color="black", size=26, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title.x = element_text(face="bold", color="black", size=20, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(face="bold", color="black", size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", color="black", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(face="bold", color="black", size=20),
        legend.title = element_text(face="bold", size=18),
        legend.text = element_text(size=16),
        legend.position = 'right') + labs(fill = "Classification")

annotate_figure(plot3, text_grob("p-value < 2.2e-16", x=0.2,  y=-2, hjust=0, color = "#000000", size=12))

ggsave(paste0("./output/essentialGenes/plots/genesApClassificationOrg.jpeg"), width = 30, height = 25, units = "cm")
ggsave(paste0("./output/essentialGenes/plots/genesApClassificationOrg.svg"), width = 30, height = 25, units = "cm")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ---- plot4 ----

# Status message
printMessage("Reaction type plot")

reactionDf <- orgGenes[!is.na(orgGenes$reaction_type),]

# Generate the boxplot plot
plot4 <- ggplot(reactionDf, aes(fill=reaction_type, x=gene_classification), ymin = -Inf, ymax = Inf) +
  # Add the bars
  geom_bar(position="dodge", stat="count") +

  scale_y_continuous(breaks=seq(from = 0, to = nrow(orgGenes), by = 250)) +
  scale_fill_manual(values = c("#ED553B", "#3CAEA3", "#a98600", "#173F5F")) +

  # Chart visual properties
  xlab("Reaction classification per gene type") +
  ylab("Genes count") +
  ggtitle("") +
  guides(fill=guide_legend(title="")) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", color="black", size=26, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title.x = element_text(face="bold", color="black", size=20, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(color="black", size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", color="black", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(face="bold", color="black", size=20),
        legend.title = element_text(face="bold", size=18),
        legend.text = element_text(size=16),
        legend.position = 'right') + labs(fill = "Classification")

plot4

ggsave(paste0("./output/essentialGenes/plots/genesReactionsClassification.jpeg"), width = 30, height = 25, units = "cm")
ggsave(paste0("./output/essentialGenes/plots/genesReactionsClassification.svg"), width = 30, height = 25, units = "cm")

#::::::::::::::::::::::::::::::::

reactionDf <- reactionDf[reactionDf$is_bottleneck==1,]

print(nrow(reactionDf[reactionDf$reaction_type == 'irreversible',]))
print(nrow(reactionDf[reactionDf$reaction_type == 'reversible',]))

# Generate the boxplot plot
plot4_2 <- ggplot(reactionDf, aes(fill=reaction_type, x=reaction_type), ymin = -Inf, ymax = Inf) +
  # Add the bars
  geom_bar(position="dodge", stat="count") +

  scale_y_continuous(breaks=seq(from = 0, to = nrow(orgGenes), by = 250)) +
  scale_fill_manual(values = c("#ED553B", "#3CAEA3", "#a98600", "#173F5F")) +

  # Chart visual properties
  xlab("Reaction classification") +
  ylab("APs count") +
  ggtitle("") +
  guides(fill=guide_legend(title="")) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", color="black", size=26, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title.x = element_text(face="bold", color="black", size=20, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(color="black", size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", color="black", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(face="bold", color="black", size=20),
        legend.title = element_text(face="bold", size=18),
        legend.text = element_text(size=16),
        legend.position = 'right') + labs(fill = "Classification")

plot4_2

ggsave(paste0("./output/essentialGenes/plots/apsReactionsClassification.jpeg"), width = 30, height = 25, units = "cm")
ggsave(paste0("./output/essentialGenes/plots/apsReactionsClassification.svg"), width = 30, height = 25, units = "cm")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ---- plot5 ----

# Status message
printMessage("Betweenness by gene classification")

# Run the Kruskal-Wallis test
gene_classification_betweenness_test <- kruskal.test(betweenness ~ gene_classification, data = orgGenes)
pairwise.wilcox.test(orgGenes$betweenness, orgGenes$gene_classification, p.adjust.method = "BH")

# Remove NA cases and non sense betweenness
betweennessDataset <- orgGenes[!is.na(orgGenes$betweenness),]
betweennessDataset <- betweennessDataset[betweennessDataset$betweenness >= 0 & betweennessDataset$betweenness <= 1,]

# Keep just the AP genes data
betweennessDataset <- betweennessDataset[betweennessDataset$is_bottleneck == 1,]

# Remove the outliers
betweennessDataset <- betweennessDataset[betweennessDataset$betweenness >= 0 & betweennessDataset$betweenness <= 0.45,]

# Print some metrics related to AP lethal and AP non lethal
print(mean(betweennessDataset[betweennessDataset$gene_classification_group==1,]$betweenness))
print(mean(betweennessDataset[betweennessDataset$gene_classification_group==3,]$betweenness))

# Generate the boxplot plot
plot5 <- ggplot(betweennessDataset, aes(x=gene_classification, y=betweenness, na.rm = TRUE)) +
  # Add the boxplot
  stat_boxplot(geom = "errorbar", width = 0.15) +
  geom_boxplot(aes(group=gene_classification, fill=gene_classification),
               colour=c("#C92D12", "#1E5953"),
               fill=c("#ED553B", "#3CAEA3"), position=position_dodge(0.5)) +

  # Chart visual properties
  xlab("Genes classification") +
  ylab("Betweenness") +
  ggtitle("") +
  guides(fill=guide_legend(title="")) +
  theme_bw() +
  theme(plot.title = element_text(face="bold", color="black", size=26, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title.x = element_text(face="bold", color="black", size=20, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(color="black", size=18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", color="black", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(face="bold", color="black", size=20),
        legend.title = element_text(face="bold", size=18),
        legend.text = element_text(size=16),
        legend.position = 'top') + labs(fill = "Classification")

annotate_figure(plot5, text_grob('Kruskal-Wallis, p-value < 2.2e-16', x=0.2,  y=-2, hjust=0, color = "#000000", size=12))

ggsave(paste0("./output/essentialGenes/plots/genesClassificationBetweenness.jpeg"), width = 30, height = 25, units = "cm")
ggsave(paste0("./output/essentialGenes/plots/genesClassificationBetweenness.svg"), width = 30, height = 25, units = "cm")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ---- plot6 ----

# Status message
printMessage("Genes AP, lethal and irreversible distribution")

genesApLethalIrreversibleList <- orgGenes[orgGenes$reaction_type == 'irreversible' &
                                            orgGenes$is_bottleneck == 1 & orgGenes$lethal_nonlethal == 'lethal',]

genesApLethalIrreversibleList <- genesApLethalIrreversibleList[!is.na(genesApLethalIrreversibleList$reaction_type) &
                                                               !is.na(genesApLethalIrreversibleList$is_bottleneck) &
                                                                 !is.na(genesApLethalIrreversibleList$lethal_nonlethal),]

genesApLethalIrreversibleList <- genesApLethalIrreversibleList[,c('pathway', 'org', 'dictID', 'name', 'entrez', 'ensembl_gene_id', 'reaction', 'reaction_type', 'is_bottleneck', 'betweenness', 'degree',
                                                                  'community', 'lethal_nonlethal', 'gene_classification')]

# Fill the pathway code with zeros
genesApLethalIrreversibleList <- fillPathwayCodeWithZeros(genesApLethalIrreversibleList)

View(genesApLethalIrreversibleList)
DataExplorer::create_report(genesApLethalIrreversibleList)

# Write the genesApLethalIrreversibleList
write.csv(genesApLethalIrreversibleList, file=paste0('./output/essentialGenes/genesApLethalIrreversibleList.csv'), row.names = F)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

dfReport <- orgGenes[,c('reaction_type', 'betweenness', 'is_bottleneck', 'lethal_nonlethal', 'gene_classification')]
DataExplorer::create_report(dfReport)
