rm(list=ls())
setwd("/Users/taqmatsu/Analysis/Cortex_thickness/GWAS")
load("08212012Big_table.rdata")

data.collapse <- data[,grep("by_sex_and_age", colnames(data),fixed=T)]

for(i in 1:68){
  if(i == 1){data.collapse2 <- pmax(data.collapse[,i],
                                    data.collapse[,i+68],
                                    data.collapse[,i+68*2],na.rm=T)} else {
  data.collapse2 <- cbind(data.collapse2, 
                          pmax(data.collapse[,i],data.collapse[,i+68], data.collapse[,i+68*2],na.rm=T))}
}


for(j in 1:2){
  data.collapse2 <- cbind(data.collapse2, 
                          pmax(data.collapse[,204+j],data.collapse[,206+j],data.collapse[,208+j],na.rm=T)) }

colnames(data.collapse2) <- c(gsub("_adjusted_by_sex_and_age_DOM","",colnames(data.collapse)[1:68]),
                              "median_lh", "median_rh")

data.collapse2 <- cbind(data[,1:7],data.collapse2)



data.small2 <- data
data.small2[data.small2$MAF < 0.15 , grep("REC",names(data.small2))] <- 0
data.small2[data.small2$MAF < 0.005 , -1:-7] <- 0
cbind(cumsum(rev(table(cut(apply(data.small2[,-1:-7],1,max),9:33/2)))))
data.small2 <- data.small2[which(apply(data.small2[,-1:-7],1,max) > 7),]
data.small2 <- data.small2[order(data.small2$CHR,data.small2$BP),]
data.small2 <- cbind(data.small2[1:7],
                     data.small2[7 + which(apply(data.small2[,-1:-7],2,max) > 7)])

data.collapse3 <- data.collapse2[data.collapse2$SNP %in% data.small2$SNP,]
colnames(data.collapse3)
for(i in c(1:34+7)){
  if(i == 8){plot.data <- cbind(data.collapse3[,c(1,i,i+34)],
                                rep(gsub("lh_","",colnames(data.collapse3)[i]),nrow(data.collapse3)))
                                colnames(plot.data)[2:3] <- c("Lh","Rh")  } else {
                                    temp <- cbind(data.collapse3[,c(1,i,i+34)],
                                    rep(gsub("lh_","",colnames(data.collapse3)[i]),nrow(data.collapse3)))
                                    colnames(temp)[2:3] <- c("Lh", "Rh")
                                    plot.data <- rbind(plot.data,temp) 
                                                     }
}

colnames(plot.data)[4] <- "Region"
Lh_Rh <- ggplot(plot.data, aes(Lh, Rh, label=Region))+
          geom_point()+
          geom_text()+
          facet_wrap(~SNP)
ggsave("Pvalue_plot_Lh_Rh.pdf", width=15, height=15)