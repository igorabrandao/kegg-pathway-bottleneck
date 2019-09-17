rm(list = ls())
setwd("/home/clovis/Doutorado/Projetos/GrupoDalmolin/Igor/")

load("1_00010.RData")
enzymeList<-enzymeTotalFrequency

load("2_00020.RData")
enzymeList<-rbind(enzymeList,enzymeTotalFrequency)

load("3_00030.RData")
enzymeList<-rbind(enzymeList,enzymeTotalFrequency)

hist(enzymeList$percentage[enzymeList$freq>=0], breaks = 20)
mean<-mean(enzymeList$freq[enzymeList$freq>1])
sd<-sd(enzymeList$freq[enzymeList$freq>1])

#enzymeList<-enzymeList[enzymeList$freq>=(mean-sd),]
enzymeList<-enzymeList[order(enzymeList$percentage,decreasing = T),]

percRange<- unique(round(enzymeList$percentage,1))


countsBase<-c(white=nrow(enzymeList[enzymeList$is_bottleneck ==1,]),
                       black=nrow(enzymeList[enzymeList$is_bottleneck !=1,]))
countsBase[1]+countsBase[2]

result<-data.frame(perc=numeric(),
                      white=numeric(),
                      black=numeric(),
                      drawn=numeric(),
                      freq=numeric(),
                      hyp=numeric())

top<-enzymeList[1:13,]
drawn<-nrow(top)
freq<-nrow(top[top$is_bottleneck == 1,])
phyper(q=freq,
       m=countsBase["white"],
       n=countsBase["black"],
       k=drawn,
       lower.tail = F)


perc=94
for( perc in percRange){
  countsTop<-data.frame(perc=numeric(),
                        white=numeric(),
                        black=numeric(),
                        drawn=numeric(),
                        freq=numeric(),
                        hyp=numeric())
  
  top<-enzymeList[enzymeList$percentage>=perc,]
  drawn<-nrow(top)
  freq<-nrow(top[top$is_bottleneck == 1,])
  countsTop[1,]<-t(c(perc,c(countsBase),drawn,freq,NA))
  countsTop[1,"hyp"]<-phyper(countsTop[1,"freq"],
                        countsTop[1,"white"],
                        countsTop[1,"black"],
                        countsTop[1,"drawn"],
                        lower.tail = F)
  result<-rbind(result,countsTop)
  #exclude non significant zeros
  #occured when drawn all balls avaliable on upper tail test
#  countsTop$hyp[countsTop$white+countsTop$black == countsTop$drawn &
                  #countsTop$hyp == 0]<- NA
#  countsTop<-countsTop[, c("faixa","hyp")]
  
                        
}

nrow(result[result$hyp<=0.05&round(result$hyp,5)!=0,])
(result[result$hyp<=0.01&round(result$hyp,5)!=0,])
