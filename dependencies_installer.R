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
install.packages("igraph", repos = "http://cran.us.r-project.org")

# RCurl
install.packages("RCurl", repos = "http://cran.us.r-project.org")

# rvest
install.packages("rvest", repos = "http://cran.us.r-project.org")

# stringr
install.packages("stringr", repos = "http://cran.us.r-project.org")

# pracma
install.packages("pracma", repos = "http://cran.us.r-project.org")

# foreach
install.packages("foreach", repos = "http://cran.us.r-project.org")
