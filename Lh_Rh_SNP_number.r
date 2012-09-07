rm(list=ls())
setwd("~/Analysis/Cortex_thickness/GWAS/")
library("ggplot2")
library("maptools")

load("08212012Big_table.rdata")

#Make table containing SNP number being over -Log10P 5
plotdata <- cbind(
  Lh=apply(data[,280:313], 2, function(x)
  sum(x>5 & !is.na(x))
        ),
  Rh=apply(data[,280:313+34], 2, function(x)
  sum(x>5 & !is.na(x))
        )
)
rownames(plotdata) <- gsub("lh_","", rownames(plotdata))
rownames(plotdata) <- gsub("_adjusted_by_age_no3sd_REC","",rownames(plotdata))
plotdata$Lh_Rh<-round(plotdata$Lh/plotdata$Rh,2)
write.csv(plotdata, file="plotdata.csv")
plotdata <-as.data.frame(plotdata)

plotdata2 <- cbind(
  Lh=apply(data[,212:245], 2, function(x)
    sum(x>5 & !is.na(x))
  ),
  Rh=apply(data[,212:245+34], 2, function(x)
    sum(x>5 & !is.na(x))
  )
)

rownames(plotdata2) <- gsub("lh_","", rownames(plotdata))
rownames(plotdata2) <- gsub("_adjusted_by_age_no3sd_DOM","",rownames(plotdata))
plotdata2$Lh_Rh<-round(plotdata2$Lh/plotdata2$Rh,2)
write.csv(plotdata2, file="plotdata2.csv")
plotdata2 <-as.data.frame(plotdata2)

plot <- ggplot(plotdata, aes(Lh, Rh, label=rownames(plotdata))) +
  geom_point() + xlim(0,600) + ylim(0,600) +
  geom_abline(intercept=0, slope=1) +
  geom_text(hjust=0.2, vjust=-0.8, size=6, angle=15)
ggsave("Lh_Rh_SNP(-Log10P>5,rec,no3sd).pdf",width=15, height=15)

plotdiv <- ggplot(plotdata, aes(Lh, Rh, label=rownames(plotdata))) +
  geom_point() + xlim(0,150) + ylim(0,150) +
  geom_abline(intercept=0, slope=1) +
  geom_text(hjust=0.2, vjust=-0.8, size=6, angle=15)
ggsave("Lh_Rh_SNP_div(-Log10P>5,rec,no3sd).pdf",width=15, height=15)

plotdiv2 <- ggplot(plotdata, aes(Lh, Rh, label=rownames(plotdata))) +
  geom_point() + xlim(0,25) + ylim(0,25) +
  geom_abline(intercept=0, slope=1) +
  geom_text(hjust=0.2, vjust=-0.8, size=6, angle=15)
ggsave("Lh_Rh_SNP_div2(-Log10P>5,rec,no3sd).pdf",width=15, height=15)

plot2 <- ggplot(plotdata2, aes(Lh, Rh, label=rownames(plotdata2)))+
  geom_point() + xlim(0,20) + ylim(0,20) +
  geom_abline(intercept=0, slope=1) +
  geom_text(hjust=0.2, vjust=-0.8, size=6, angle=15)
ggsave("Lh_Rh_SNP(-Log10P>5,dom_no3sd).pdf",width=15, height=15)