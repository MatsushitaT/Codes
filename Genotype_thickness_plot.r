rm(list=ls())
setwd("/Users/taqmatsu/Analysis/Cortex_thickness/GWAS")

library("ggplot2")

Top24map <- read.table("24topSNP.map",header=F, sep="")
col <- rep("character",54)
Top24ped <- read.table("24topSNP.ped",header=F, colClasses=col, sep="")

thickness <- read.table("phenotype_not3sd.txt",header=T,sep="")

for(i in 1:24)
  if(i==1){Top24table <- paste(Top24ped[,2*i+5],Top24ped[,2*i+6])} else {
  Top24table <- cbind(Top24table, paste(Top24ped[,2*i+5],Top24ped[,2*i+6]))}

colnames(Top24table) <- as.character(Top24map$V2)

Top24table <- cbind(FID=Top24ped[,1],IID=Top24ped[,2],Top24table)

Top24table <- merge(Top24table, thickness, all.x=T, all.y=F)
str(plotdata)
for(i in 3:26){
  for(j in 1:68+26){
    plotdata <- Top24table[,c(i,j)]; title <- paste(colnames(plotdata)[1],colnames(plotdata)[2],sep="_");
    plotresult <- ggplot(plotdata, aes(plotdata[,1],plotdata[,2]))+
      ylab(colnames(plotdata[2])) +
      xlab(colnames(plotdata)[1]) +
      opts(legend.position = "none") +
      geom_boxplot(aes(fill=plotdata[,1]),outlier.size=0)+
      geom_jitter(alpha=0.4)
      ggsave(paste("Genotype_thickness/",title,".pdf",sep=""), width = 5, height = 5)}} 
                       