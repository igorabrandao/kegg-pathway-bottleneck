#*********#
# CRÍTICA #
#*********#

# The authors posit that articulation points may be essential enzymes in pathways, but there is not significant
# analysis of this hypothesis. It would be informative to more broadly compare AP enzymes to the DEG. Perhaps a
# figure comparing the percentage of APs and non-APs that are essential would provide additional insights that
# increase the utility of this AP analysis.

#*************************#
# IDEIAS DE COMO RESOLVER #
#*************************#

# 1) AMPLIAR A DISCUSSÃO DE ESSENCIALIDADE. ACRESCENTAR NOS DADOS DA IARA E TAMBÉM AMPLIAR A DISCUSSÃO DO DEG
# 2) COMPOR UMA FIGURA COM A AVALIAÇÃO DE ESSENCIALIDADE DE M MUSCULUS E E. COLI.
# 3) DISCUTIR O CONCEITO DE ESSENCIALIDADE DE ROTA. OU SEJA, SE UMA ROTA NÃO FOR ESSENCIAL, TALVEZ O AP DESSA ROTA TAMBÉM NÃO SEJA.

#*********************************************#
# CLASSIFICAÇÕES DE ESSENCIALIDADE DO DATASET #
#*********************************************#

# ----------------------------------------------------------------------#
#   E, DE (?), ES (?), D (?) - essential
# ----------------------------------------------------------------------#
#  NE - nonessential
# ----------------------------------------------------------------------#
#  Condicionais:
#  growth defective (GD)
#  growth advantaged (GA)
#  F: fitness in vitro
#  E-infection: required for single infection to lung of mice
#  F-infection: Potential intermediate attenuation during single infection
#  E-co-infection: required for co-infection to lung of mice
#  F-co-infection: Potential intermediate attenuation during co-infection
# ----------------------------------------------------------------------#
#  Desconhecido:
#  ND, U
# ----------------------------------------------------------------------#
#  Incerto
#  S
# ----------------------------------------------------------------------#

#*********#
# SOLUÇÃO #
#*********#

# ---- PIPELINE ----
library(ggplot2)
library(svglite)

#***************************************************************************#
# ----  Passo 1: Carregar os dados para a análise ----
#                                                                           #
# Nota: lembre de ajustar o caminho do arquivo de acordo com o seu ambiente #
#***************************************************************************#

# Carrega o dataSet com os genes do OGEE
dataSetOGEE <- read.csv("./revisao_jbc/revisor1/major1/t1_OGEE.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# Carrega o dataSet com os genes do KEGG + métricas de rede
dataSetOrgGenes <- read.csv("./revisao_jbc/revisor1/major1/t2_orgGenes.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# Carrega o dataSet com o dicionário (ensemvlID -> Entrez ID)
ensemblEntrezDictionary <- read.csv("./revisao_jbc/revisor1/major1/t3_ensemblEntrezDictionary.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

#***************************************************************************#
# ----  Passo 2: Match dos dataSets através do dicionário ----
#***************************************************************************#

# Faz o match entre os dataSets
dataSetMatch <- dataSetOrgGenes
dataSetMatch <- cbind(dataSetMatch, sciName = NA)
dataSetMatch <- cbind(dataSetMatch, kingdom = NA)
dataSetMatch <- cbind(dataSetMatch, ensemblID = NA)
dataSetMatch <- cbind(dataSetMatch, symbols = NA)
dataSetMatch <- cbind(dataSetMatch, datasets = NA)
dataSetMatch <- cbind(dataSetMatch, datasetIDs = NA)
dataSetMatch <- cbind(dataSetMatch, essentiality.status = NA)
dataSetMatch <- cbind(dataSetMatch, essentiality.consensus = NA)
dataSetMatch <- cbind(dataSetMatch, biomartDataset = NA)

# Percorre todos os genes do KEGG e faz um match com os dados do OGEE
# Aviso: esse trecho é muito pesado!
idx = 1
apply(dataSetOrgGenes, 1, function(x) {
  # Busca o item do dicionário correspondente ao gene atual
  currentDictionaryEntry <- ensemblEntrezDictionary[ensemblEntrezDictionary$ensemblID == dataSetOGEE[idx,]$ensemblID,]

  # Se não houver correspondências com o dicionário pula o gene atual
  if (nrow(currentDictionaryEntry) > 0) {
    # Find which rows in graphProperties match with pathwayData
    rowsToMerge <- which(grepl(paste0('\\', dataSetOGEE[idx,]$org, ':', currentDictionaryEntry$entrezGeneID, '\\b'), dataSetOrgGenes$entrez))

    print(paste0('Index: ', idx))
    print(paste0('Rows to merge: ', rowsToMerge))

    # Se houver alguma correspondência, faz o match com o gene correspondente
    if (length(rowsToMerge) != 0) {
      for (idx2 in 1:length(rowsToMerge)) {
        print(paste0('Current row: ', rowsToMerge[idx2]))

        dataSetMatch[rowsToMerge[idx2],]$sciName <<- dataSetOGEE[idx,]$sciName
        dataSetMatch[rowsToMerge[idx2],]$kingdom <<- dataSetOGEE[idx,]$kingdom
        dataSetMatch[rowsToMerge[idx2],]$ensemblID <<- dataSetOGEE[idx,]$ensemblID
        dataSetMatch[rowsToMerge[idx2],]$symbols <<- dataSetOGEE[idx,]$symbols
        dataSetMatch[rowsToMerge[idx2],]$datasets <<- dataSetOGEE[idx,]$datasets
        dataSetMatch[rowsToMerge[idx2],]$datasetIDs <<- dataSetOGEE[idx,]$datasetIDs
        dataSetMatch[rowsToMerge[idx2],]$essentiality.status <<- dataSetOGEE[idx,]$essentiality.status
        dataSetMatch[rowsToMerge[idx2],]$essentiality.consensus <<- dataSetOGEE[idx,]$essentiality.consensus
        dataSetMatch[rowsToMerge[idx2],]$biomartDataset <<- dataSetOGEE[idx,]$biomartDataset
      }
    }
  }

  # Increment the index
  idx <<- idx + 1
})

# Exporta as correspondências entre os dados do KEGG e OGEE
write.csv(dataSetMatch, file=paste0('./revisao_jbc/revisor1/major1/t4_dataSetMatch.csv'), row.names = F)

#**********************************************************************************************#
# ----  Passo 3: Separar os genes essenciais e não essenciais nos datasets de APs e non-APs ----
#**********************************************************************************************#

# Carrega o dataSet de correspondência entre os dados do KEGG e OGEE
dataSetMatch <- read.csv("./revisao_jbc/revisor1/major1/t4_dataSetMatch.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# Remove casos inválidos
dataSetMatch <- dataSetMatch[!is.na(dataSetMatch$essentiality.consensus),]

# Remove os genes duplicados pelo KEGG
dataSetMatch <- unique(dataSetMatch)

# Separa os dataSets the APs e Non APs
dataSetAp <- dataSetMatch[dataSetMatch$is_bottleneck == 1,]
dataSetNonAp <- dataSetMatch[dataSetMatch$is_bottleneck == 0,]

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Filtra os APs essenciais
apEssential <- dataSetMatch[dataSetMatch$is_bottleneck == 1 & dataSetMatch$essentiality.consensus == 'Essential',]

# Filtra os APs não essenciais
apNonEssential <- dataSetMatch[dataSetMatch$is_bottleneck == 1 & dataSetMatch$essentiality.consensus == 'Nonessential',]

# Filtra os APs com essencialidade condicional
apEssentialityConditional <- dataSetMatch[dataSetMatch$is_bottleneck == 1 & dataSetMatch$essentiality.consensus == 'Conditional',]

# Filtra os não APs essenciais
nonApEssential <- dataSetMatch[dataSetMatch$is_bottleneck == 0 & dataSetMatch$essentiality.consensus == 'Essential',]

# Filtra os não APs não essenciais
nonApNonEssential <- dataSetMatch[dataSetMatch$is_bottleneck == 0 & dataSetMatch$essentiality.consensus == 'Nonessential',]

# Filtra os não APs com essencialidade condicional
nonApEssentialityConditional <- dataSetMatch[dataSetMatch$is_bottleneck == 0 & dataSetMatch$essentiality.consensus == 'Conditional',]

#************************************************************************#
# ----  Passo 4: Verificar as proporações de APs e non APs essenciais ----
#************************************************************************#

# Calcula as proporções dos APs e non APs essenciais ou não comparando com todos os genes do dataSet
genesProportionAllGenes <- data.frame(apClassification = NA, essentialityClassification = NA, proportion = NA, stringsAsFactors = F)

genesProportionAllGenes[1,] <- c('AP', 'Conditional', (nrow(apEssentialityConditional) / nrow(dataSetMatch))) # APs com essencialidade condicional
genesProportionAllGenes[2,] <- c('AP', 'Essential', (nrow(apEssential) / nrow(dataSetMatch))) # APs essenciais
genesProportionAllGenes[3,] <- c('AP', 'Non-Essential', (nrow(apNonEssential) / nrow(dataSetMatch))) # APs não essenciais

genesProportionAllGenes[4,] <- c('Non-AP', 'Conditional', (nrow(nonApEssentialityConditional) / nrow(dataSetMatch))) # Non APs com essencialidade condicional
genesProportionAllGenes[5,] <- c('Non-AP', 'Essential', (nrow(nonApEssential) / nrow(dataSetMatch))) # Non APs essenciais
genesProportionAllGenes[6,] <- c('Non-AP', 'Non-Essential', (nrow(nonApNonEssential) / nrow(dataSetMatch))) # Non APs non essenciais

# Converte as proporções para o tipo numérico
genesProportionAllGenes$proportion = as.numeric(genesProportionAllGenes$proportion)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Calcula as proporções dos APs e non APs essenciais ou não somente no mesmo grupo (AP ou non-AP)
genesProportion <- data.frame(apClassification = NA, essentialityClassification = NA, proportion = NA, stringsAsFactors = F)

genesProportion[1,] <- c('AP', 'Conditional', (nrow(apEssentialityConditional) / nrow(dataSetAp))) # APs com essencialidade condicional
genesProportion[2,] <- c('AP', 'Essential', (nrow(apEssential) / nrow(dataSetAp))) # APs essenciais
genesProportion[3,] <- c('AP', 'Non-Essential', (nrow(apNonEssential) / nrow(dataSetAp))) # APs não essenciais

genesProportion[4,] <- c('Non-AP', 'Conditional', (nrow(nonApEssentialityConditional) / nrow(dataSetNonAp))) # Non APs com essencialidade condicional
genesProportion[5,] <- c('Non-AP', 'Essential', (nrow(nonApEssential) / nrow(dataSetNonAp))) # Non APs essenciais
genesProportion[6,] <- c('Non-AP', 'Non-Essential', (nrow(nonApNonEssential) / nrow(dataSetNonAp))) # Non APs non essenciais

# Converte as proporções para o tipo numérico
genesProportion$proportion = as.numeric(genesProportion$proportion)

#**************************************#
# ----  Passo 5: Criar alguns plots ----
#**************************************#

# ----  plot 1 ----

# Plot comparando os grupos dos genes pela essencialidade comparando com todos os genes do dataSet
plot1 <- ggplot(genesProportionAllGenes) +
  # Add the bars
  geom_bar(aes(fill=essentialityClassification, x=apClassification, y=proportion), position="dodge", stat="identity") +

  # Add labels to bars group
  geom_text(aes(y=(proportion + 0.06), x=c(0.7, 1, 1.3, 1.7, 2, 2.3), label=(paste0(format(round(proportion * 100, 2), nsmall = 2), '%'))),
            size=6) +

  scale_y_continuous(limits=c(0, 1), labels = scales::percent) +
  scale_fill_manual(values = c("#a98600", "#173F5F", "#ED553B")) +

  # Chart visual properties
  xlab("Classification") +
  ylab("Proportion") +
  ggtitle("Proportion of APs and non-APs by essentiality in all genes") +
  theme_bw() +
  theme(plot.title = element_text(face="bold", color="black", size=26, margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title.x = element_text(face="bold", color="black", size=20, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(face="bold", color="black", size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", color="black", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.text.y = element_text(face="bold", color="black", size=20),
        legend.title = element_text(face="bold", size=18),
        legend.text = element_text(size=16),
        legend.position = 'right') + labs(fill = "Classification")

plot1

# Exporta a imagem para edição posterior
ggsave(paste0("./revisao_jbc/revisor1/major1/figuras/genesClassificationAllGenes.jpeg"), width = 30, height = 20, units = "cm")
ggsave(paste0("./revisao_jbc/revisor1/major1/figuras/genesClassificationAllGenes.svg"), width = 30, height = 20, units = "cm")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ----  plot 2 ----

# Plot comparando os grupos dos genes pela essencialidade comparando somente os genes do mesmo grupo
plot2 <- ggplot(genesProportion) +
  # Add the bars
  geom_bar(aes(fill=essentialityClassification, x=apClassification, y=proportion), position="dodge", stat="identity") +

  # Add labels to bars group
  geom_text(aes(y=(proportion + 0.06), x=c(0.7, 1, 1.3, 1.7, 2, 2.3), label=(paste0(format(round(proportion * 100, 2), nsmall = 2), '%'))),
            size=6) +

  scale_y_continuous(limits=c(0, 1), labels = scales::percent) +
  scale_fill_manual(values = c("#a98600", "#173F5F", "#ED553B")) +

  # Chart visual properties
  xlab("Classification") +
  ylab("Proportion") +
  ggtitle("Proportion of APs and non-APs by essentiality and group") +
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

print(paste0("Quantidade de Genes APs: ", nrow(dataSetAp)))
print(paste0("Quantidade de Genes não APs: ", nrow(dataSetNonAp)))

# Exporta a imagem para edição posterior
ggsave(paste0("./revisao_jbc/revisor1/major1/figuras/genesClassification.jpeg"), width = 30, height = 20, units = "cm")
ggsave(paste0("./revisao_jbc/revisor1/major1/figuras/genesClassification.svg"), width = 30, height = 20, units = "cm")
