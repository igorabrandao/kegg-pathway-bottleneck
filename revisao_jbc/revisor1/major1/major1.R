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

# Carrega o dataSet referente a análise deste major
dataSet <- read.csv("./revisao_jbc/revisor1/major1/dataset.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# Limpa o dataSet e remove os casos dos genes sem classificação de AP
dataSet <- dataSet[!is.na(dataSet$is_ap),]

# Separa os dataSets the APs e Non APs
dataSetAp <-  dataSet[dataSet$is_ap == 1,]
dataSetNonAp <-  dataSet[dataSet$is_ap == 0,]

# Filtra somente os genes essenciais
dataSetEssential <- dataSet[dataSet$essential %in% c('E', 'DE', 'ES', 'D'),]

#**********************************************************************************************#
# ----  Passo 2: Separar os genes essenciais e não essenciais nos datasets de APs e non-APs ----
#**********************************************************************************************#

# Nota: veja a tabela de classificação no início deste arquivo

#*****#
# APs #
#*****#

# Filtra os APs essenciais
apEssential <- dataSetAp[dataSetAp$essential %in% c('E', 'DE', 'ES', 'D'),]

# Filtra os APs não essenciais
apNonEssential <- dataSetAp[dataSetAp$essential %in% c('NE'),]

#*********#
# Non APs #
#*********#

# Filtra os não APs essenciais
nonApEssential <- dataSetNonAp[dataSetNonAp$essential %in% c('E', 'DE', 'ES', 'D'),]

# Filtra os não APs não essenciais
nonApNonEssential <- dataSetNonAp[dataSetNonAp$essential %in% c('NE'),]

#************************************************************************#
# ----  Passo 3: Verificar as proporações de APs e non APs essenciais ----
#************************************************************************#

# Calcula as proporções dos APs e non APs essenciais ou não comparando com todos os genes do dataSet
genesProportionAllGenes <- data.frame(apClassification = NA, essentialityClassification = NA, proportion = NA, stringsAsFactors = F)

# APs essenciais
genesProportionAllGenes[1,] <- c('AP', 'Essential', (nrow(apEssential) / nrow(dataSet)))

# APs não essenciais
genesProportionAllGenes[2,] <- c('AP', 'Non-Essential', (nrow(apNonEssential) / nrow(dataSet)))

# Non APs essenciais
genesProportionAllGenes[3,] <- c('Non-AP', 'Essential', (nrow(nonApEssential) / nrow(dataSet)))

# Non APs non essenciais
genesProportionAllGenes[4,] <- c('Non-AP', 'Non-Essential', (nrow(nonApNonEssential) / nrow(dataSet)))

# Converte as proporções para o tipo numérico
genesProportionAllGenes$proportion = as.numeric(genesProportionAllGenes$proportion)

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Calcula as proporções dos APs e non APs essenciais ou não somente no mesmo grupo (AP ou non-AP)
genesProportion <- data.frame(apClassification = NA, essentialityClassification = NA, proportion = NA, stringsAsFactors = F)

# APs essenciais
genesProportion[1,] <- c('AP', 'Essential', (nrow(apEssential) / nrow(dataSetAp)))

# APs não essenciais
genesProportion[2,] <- c('AP', 'Non-Essential', (nrow(apNonEssential) / nrow(dataSetAp)))

# Non APs essenciais
genesProportion[3,] <- c('Non-AP', 'Essential', (nrow(nonApEssential) / nrow(dataSetNonAp)))

# Non APs non essenciais
genesProportion[4,] <- c('Non-AP', 'Non-Essential', (nrow(nonApNonEssential) / nrow(dataSetNonAp)))

# Converte as proporções para o tipo numérico
genesProportion$proportion = as.numeric(genesProportion$proportion)

#**************************************#
# ----  Passo 4: Criar alguns plots ----
#**************************************#

# ----  plot 1 ----

# Plot comparando os grupos dos genes pela essencialidade comparando com todos os genes do dataSet
plot1 <- ggplot(genesProportionAllGenes) +
  # Add the bars
  geom_bar(aes(fill=essentialityClassification, x=apClassification, y=proportion), position="dodge", stat="identity") +

  # Add labels to bars group
  geom_text(aes(y=(proportion + 0.06), x=c(0.75, 1.20, 1.75, 2.20), label=(paste0(format(round(proportion * 100, 2), nsmall = 2), '%'))),
            size=6) +

  scale_y_continuous(limits=c(0, 1), labels = scales::percent) +
  scale_fill_manual(values = c("#ED553B", "#173F5F")) +

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
ggsave(paste0("./revisao_jbc/revisor1/major1/genesClassificationAllGenes.jpeg"), width = 30, height = 20, units = "cm")
ggsave(paste0("./revisao_jbc/revisor1/major1/genesClassificationAllGenes.svg"), width = 30, height = 20, units = "cm")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ----  plot 2 ----

# Plot comparando os grupos dos genes pela essencialidade comparando somente os genes do mesmo grupo
plot2 <- ggplot(genesProportion) +
  # Add the bars
  geom_bar(aes(fill=essentialityClassification, x=apClassification, y=proportion), position="dodge", stat="identity") +

  # Add labels to bars group
  geom_text(aes(y=(proportion + 0.06), x=c(0.75, 1.20, 1.75, 2.20), label=(paste0(format(round(proportion * 100, 2), nsmall = 2), '%'))),
            size=6) +

  scale_y_continuous(limits=c(0, 1), labels = scales::percent) +
  scale_fill_manual(values = c("#ED553B", "#173F5F")) +

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

# Exporta a imagem para edição posterior
ggsave(paste0("./revisao_jbc/revisor1/major1/genesClassification.jpeg"), width = 30, height = 20, units = "cm")
ggsave(paste0("./revisao_jbc/revisor1/major1/genesClassification.svg"), width = 30, height = 20, units = "cm")
