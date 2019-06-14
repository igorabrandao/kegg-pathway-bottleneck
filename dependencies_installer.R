##########################
# Dependencies installer #
##########################

# dependencies_installer.R #

#' Install the necessary dependencies to run the main
#' pipeline
#'
#' @author
#' Igor Brandão

# BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# KEGGREST
BiocManager::install("KEGGREST")

# KEGGREST
BiocManager::install("KEGGgraph")

# igraph
install.packages("igraph")

# RCurl
install.packages("RCurl")

# rvest
install.packages("rvest")

# stringr
install.packages("stringr")

# pracma
install.packages("pracma")

# foreach
install.packages("foreach")
