#*******************************************#
# Pipeline to perform the binomial analysis #
#*******************************************#

# ---- IMPORT SECTION ----

# 2.2_binomial #

#' This is the pipeline script to perform
#' the binomial analysis
#'
#' @author
#' Clóvis F. Reis / Igor Brandão

# Import the necessary libraries
library(ggplot2)
library(svglite)
#options(scipen = 999, digits = 2) # sig digits

# Import the basic functions
files.sources = NULL
files.sources[1] = paste0("./R/functions", "/", "statisticsHelper.R")
files.sources[2] = paste0("./R/functions", "/", "helperFunctions.R")
sapply(files.sources, source)

#*******************************************************************************************#

# ---- SETTINGS SECTION ----

p_value_  <- 0.01
intervalo <- 10

#*******************************************************************************************#

# ---- PIPELINE SECTION ----

#***************#
# Pipeline flow #
#***************#

# Load the dataSet
dataSet_ <- generateDataSetCSV(testName_ = 'hypergeometric')

#**************************************************************************##
# Apply the 2nd normalization:                                              #
# 100% of occurrence in a pathway can be compared with 50% of other pathway #
#**************************************************************************##

# Get the unique pathways
uniquePathways <- unique(dataSet_$pathway)

# Calculates the normalized frequency
for (item in uniquePathways) {
# Max frequency in a pathway
pathwayMaxFrequency <- max(dataSet_$percentage[dataSet_$pathway==item])

# Min frequency in a pathway
pathwayMinFrequency <- min(dataSet_$percentage[dataSet_$pathway==item])

# Normalized frequency for each protein
dataSet_$normalizedPercentage[dataSet_$pathway==item] <-
    (dataSet_$percentage[dataSet_$pathway==item]-pathwayMinFrequency)/
    (pathwayMaxFrequency-pathwayMinFrequency)
}

# Fix for NAN cases (when min and max frequency have the same values)
dataSet_$normalizedPercentage[is.nan(dataSet_$normalizedPercentage)] <- dataSet_$percentage[is.nan(dataSet_$normalizedPercentage)]/100

#********************##
# Prepare the dataset #
#********************##

# Filter dataSet from proteins with ZERO frequency
dataSet_ <- dataSet_[!dataSet_$occurrences==0, c('name', 'pathway', 'is_bottleneck', 'bottleneck_classification', 'percentage', 'normalizedPercentage')]

# Order the dataSet
dataSet_ <- dataSet_[order(dataSet_$percentage, decreasing = T),]

# Define the filename
exportFile <- "binomialDistribution"

#*****************************************##
# Generate the hypergeometric distribution #
#*****************************************##

# Create a dataFrame for the result
distribution <- data.frame(ini=0,
                           range=0,
                           bottleneck=0,
                           non_bottleneck=0,
                           binom=NA,
                           stringsAsFactors = F)

# Count the bottlenecks and non-bottlenecks
countsBase <- c(bottleneck=nrow(dataSet_[dataSet_$is_bottleneck ==1,]), non_bottleneck=nrow(dataSet_[dataSet_$is_bottleneck !=1,]))

proporcao<- countsBase[1]/(countsBase[1]+countsBase[2])

# Verification
if ((countsBase[1] +countsBase[2]) != nrow(dataSet_)) {
# Save the log file
printLog(message_='The sum of bottlenecks with non-bottlenecks is different from the dataSet total rows. Skipping it...',
            file_='hypergeometricDistribution')

return(FALSE)
}


ranges<-seq(1,101,intervalo)
#ranges<-c(1,477,1248,2014,5498)
# Loop over all dataSet
idx=1
for (idx in 1:(length(ranges)-1)) {
  # Set the cumulative range [initVal:range]
  initVal <- ranges[idx]
  range<-ranges[idx+1]

  # Temporaly dataFrame for indexed results
  countsTop <- data.frame(ini=numeric(),
                          range=numeric(),
                          bottleneck=numeric(),
                          non_bottleneck=numeric(),
                          binom=numeric(),
                          stringsAsFactors = F)

  # Retrieve the cumulative proteins
  top <- dataSet_[dataSet_$percentage>ranges[idx]&
                    dataSet_$percentage<ranges[idx+1], ]

  # Number of draws
  drawn <- nrow(top)

  # The number of articulation points in the accumulated group
  btn <- nrow(top[top$is_bottleneck == 1,])
  nbtn <- nrow(top[top$is_bottleneck == 0,])

  binom<-binom.test(btn,(btn+nbtn),proporcao,alternative = "g")
  binom<-binom$p.value

  countsTop[1,] <- t(c(initVal, range-1, btn,nbtn,binom))


  # Bind each result
  distribution <- rbind(distribution, countsTop)
  initVal <- range+1
}

# Adjust the p-value in order to compensate the accumulated error with
# Benjamini-Hochberg method
distribution$pCor <- p.adjust(distribution$binom, method = "BH")
distribution<-na.exclude(distribution)

# ---- plo1 ----

# Classify the articulation points
distribution$group <-''
distribution[distribution$pCor<=p_value_,]$group <- 'Significative'
distribution[distribution$pCor>p_value_,]$group <- 'Non-Significative'

# Plot the graph
plot1 <- ggplot() + ggtitle("") + # for the main title
  xlab("Range (%)") +
  ylab("Ratio (Articulation point/Total)") +

  geom_col(data = distribution, aes(x=range-5, y=(bottleneck/(bottleneck+non_bottleneck)), fill=group), width=9) +

  geom_hline(yintercept = proporcao, lty=2, col="#E75480") +

  scale_fill_manual(values = c("#ED553B", "#173F5F")) +

  theme_bw() +
  theme(plot.title = element_text(face="bold", size=20, hjust = 0),
       axis.title.x = element_text(face="bold", size=20, margin = margin(t = 15, r = 0, b = 0, l = 0)),
       axis.text.x = element_text(size=16),
       axis.title.y = element_text(face="bold", size=20, margin = margin(t = 0, r = 15, b = 0, l = 0)),
       axis.text.y = element_text(size=16),
       legend.title = element_text(face="bold", size=18),
       legend.text = element_text(size=16),
       legend.position = 'none') + labs(fill = "Range group")

plot1

ggsave(paste0("./output/statistics/hypergeometric/binomial.png"), width = 30, height = 20, units = "cm")
ggsave(paste0("./output/statistics/hypergeometric/binomial.svg"), width = 30, height = 20, units = "cm")

write.csv(distribution, file=paste0('./output/statistics/hypergeometric/binomial.csv'))
