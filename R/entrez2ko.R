##########################################
# Script to convert entrez genes into ko #
##########################################

# entre2ko ####

#' This is the pipeline script to perform
#' the enzymes frequencies counting
#'
#' @author
#' Diego Morais

# Import the necessary libraries
library(KEGGREST)
library(RCurl)
library(rvest)

##############################################

entrez2KoDict <- function() {

  kos <- keggLink("ko", "ko00010")

  # puxar entrez usando ko
  # https://www.kegg.jp/kegg-bin/view_ortholog_table?orthology=K01810

  kos <- gsub("ko:", "", kos)
  kos <- unname(kos)

  df <- getURL(paste0("https://www.kegg.jp/kegg-bin/view_ortholog_table?orthology=",
                      paste(kos, collapse = "+")))
  df <- gsub("<br />", "-|-", df, fixed = TRUE)

  df <- read_html(df)
  df <- html_node(df, "table")
  df <- html_table(df, header = TRUE, fill = TRUE)

  # save(df, file = "RAWKO2ENTREZ.RData", compress = "xz")

  df[1, ] <- gsub("^([[:alnum:]]*).*$", "\\1",df[1, ])
  colnames(df) <- df[1, ]
  df <- df[-1,]
  table <- df

  # foram feitas apenas 100 consultas das 103

  kos <- kos[101:103]
  df <- getURL(paste0("https://www.kegg.jp/kegg-bin/view_ortholog_table?orthology=",
                      paste(kos, collapse = "+")))
  df <- gsub("<br />", "-|-", df, fixed = TRUE)
  df <- read_html(df)
  df <- html_node(df, "table")
  df <- html_table(df, header = TRUE, fill = TRUE)
  df[1, ] <- gsub("^([[:alnum:]]*).*$", "\\1",df[1, ])
  colnames(df) <- df[1, ]
  df <- df[-1,]
  table2 <- df

  # organizando dados

  table <- table[, -c(1:3)]
  table2 <- table2[, -c(1:3)]
  table[table == ""] <- NA
  table2[table2 == ""] <- NA

  # save(table, table2, file = "RAWKO2ENTREZ.RData", compress = "xz")

  rm(kos)

  table2 <- t(table2)
  df <- NULL

  for(i in 1:nrow(table2)){
    aux <- na.omit(unique(unlist(strsplit(table2[i,], split = "-|-", fixed = T))))
    aux <- list(aux)
    id <- rownames(table2)[i]
    names(aux) <- id
    entrez <- unname(unlist(aux))
    dfaux <- data.frame(KO = rep(id, length(entrez)), ENTREZ = entrez, stringsAsFactors = F)
    df <- rbind(df, dfaux)
  }
  table <- t(table)
  for(i in 1:nrow(table)){
    aux <- na.omit(unique(unlist(strsplit(table[i,], split = "-|-", fixed = T))))
    aux <- list(aux)
    id <- rownames(table)[i]
    names(aux) <- id
    entrez <- unname(unlist(aux))
    dfaux <- data.frame(KO = rep(id, length(entrez)), ENTREZ = entrez, stringsAsFactors = F)
    df <- rbind(df, dfaux)
  }
  rm(entrez, i, id, dfaux, aux, table, table2)

  aux <- split(df, df$ENTREZ)
  aux <- lapply(aux, function(i){
    i$KO
  })
  ENTREZ2KO <- aux
  save(ENTREZ2KO, file = "ENTREZ2KO.RData", compress = "xz")
}
