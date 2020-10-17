#*********#
# CRÍTICA #
#*********#

# The threshold values of <30% and >=80% used in Figures 4 and 5 appear arbitrary,
# and the selection of these values is not elaborated in the text. Is there any meaning or
# biological relevance behind these thresholds?

#*************************#
# IDEIAS DE COMO RESOLVER #
#*************************#

# 1) ANALIZAR OS QUARTIS SEPARANDO PELA FREQUÊNCIA DAS ENZIMAS

#*********#
# SOLUÇÃO #
#*********#

# ---- PIPELINE ----

#***************************************************************************#
# ----  Passo 1: Carregar os dados para a análise ----
#                                                                           #
# Nota: lembre de ajustar o caminho do arquivo de acordo com o seu ambiente #
#***************************************************************************#

# Carrega o dataSet referente a análise deste major
dataSet <- read.csv("./revisao_jbc/revisor1/dados/major2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# Importante: em todas as análises do artigo removemos as enzimas com frequência = 0%
dataSet <- dataSet[dataSet$frequency != 0,]

# Separa os grupos de enzimas APs e non-APs
dataSetAP <- dataSet[dataSet$isAp == 1,]
dataSetNonAP <- dataSet[dataSet$isAp == 0,]

# Verificação se a contagem das enzimas bate
if ((nrow(dataSetAP) + nrow(dataSetNonAP)) == nrow(dataSet)) {
  print('Contagem de enzimas OK!')
} else {
  print('Erro!')
}

#***********************************************************************#
# ----  Passo 2: Separar os datasets de APs e non-APs em quartiles ----
#***********************************************************************#

#*****#
# APs #
#*****#

# Classifica os quartiles acordo com a frequência das enzimas
dataSetAP$quartile <- NA
dataSetAP[dataSetAP$frequency < 25,]$quartile <- 'Q1'
dataSetAP[dataSetAP$frequency >= 25 & dataSetAP$frequency < 50,]$quartile <- 'Q2'
dataSetAP[dataSetAP$frequency >= 50 & dataSetAP$frequency < 75,]$quartile <- 'Q3'
dataSetAP[dataSetAP$frequency >= 75 & dataSetAP$frequency <= 100,]$quartile <- 'Q4'

# Classifica os quartiles de acordo com o tamanho do dataSet (divisões equitativas)
dataSetAP$expectedQuartile <- NA
dataSetAP$expectedQuartile <- with(dataSetAP, factor(findInterval(dataSetAP$frequency, c(-Inf,
                        quantile(dataSetAP$frequency, probs=c(0.25, .5, .75), na.rm=TRUE), Inf)),
                        labels=c("Q1","Q2","Q3","Q4")))

# Conta os quartiles do grupo de enzimas APs
count(dataSetAP$quartile)
count(dataSetAP$expectedQuartile)

#*********#
# Non APs #
#*********#

# Classifica os quartiles acordo com a frequência das enzimas
dataSetNonAP$quartile <- NA
dataSetNonAP[dataSetNonAP$frequency < 25,]$quartile <- 'Q1'
dataSetNonAP[dataSetNonAP$frequency >= 25 & dataSetNonAP$frequency < 50,]$quartile <- 'Q2'
dataSetNonAP[dataSetNonAP$frequency >= 50 & dataSetNonAP$frequency < 75,]$quartile <- 'Q3'
dataSetNonAP[dataSetNonAP$frequency >= 75 & dataSetNonAP$frequency <= 100,]$quartile <- 'Q4'

# Classifica os quartiles de acordo com o tamanho do dataSet (divisões equitativas)
dataSetNonAP$expectedQuartile <- NA
dataSetNonAP$expectedQuartile <- with(dataSetNonAP, factor(findInterval(dataSetNonAP$frequency, c(-Inf,
                        quantile(dataSetNonAP$frequency, probs=c(0.25, .5, .75), na.rm=TRUE), Inf)),
                        labels=c("Q1","Q2","Q3","Q4")))

# Conta os quartiles do grupo de enzimas não APs
count(dataSetNonAP$quartile)
count(dataSetNonAP$expectedQuartile)

#***********************************************************************#
# ----  Passo 3: Análises estatísticas ----
#***********************************************************************#

# H0: A distribuiçào das enzimas é equitativa
# H1: A distribuiçào das enzimas não é equitativa


#*******************************************************************************************#
