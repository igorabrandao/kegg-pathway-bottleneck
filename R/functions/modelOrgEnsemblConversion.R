library(biomaRt)

# Load the pathways by organisms data
essentialGenesModelOrg <- get(load(paste0("./dictionaries", "/", "essentialGenesModelOrg.RData")))
rm(table)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")

#listDatasets(ensembl)

mouseIDs <- getBM(filters = "ensembl_peptide_id",
             attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
             values = essentialGenesModelOrg$ensembl_peptide_id, mart = ensembl)

mouseIDs[mouseIDs == ""] <- NA
mouseIDs <- na.omit(mouseIDs)

essentialGenesModelOrg = merge(essentialGenesModelOrg, mouseIDs, by='ensembl_peptide_id', all=T)









save(essentialGenesModelOrg, file = "./dictionaries/essentialGenesModelOrg.RData", compress = "xz")
