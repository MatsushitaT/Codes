rm(list=ls())
library("ggplot2")
#set workdesk
setwd("~/Analysis/Cortex_thickness/Final_dataset/")

## Read Connection
conn <- read.table("Cytoscape_data/Lateral_View.sif", header=F, sep="", as.is=T)
colnames(conn) <- c("Region1", "Connection","Region2")

conn[,"Connection"] <- gsub("Directed","Strong",conn[,"Connection"])

## Read node attribution
attri <- read.table("Cytoscape_data/Fullname.na", header=F, skip=1, sep="", as.is=T)
attri <- attri[,-2]
attri[,2] <- paste("lh_",tolower(gsub("\\+", "",attri[,2])),sep="")
attri[,2] <- gsub("\\-","",attri[,2])
colnames(attri) <- c("Abb","Full")

##Load thickness data. Including individual with age < 1 in controls
## All data is truncated to mean+- 3sd range
load("cortical.thickness.rdata")
load("control_thickness.rdata")

##Confirm name corresponding
setdiff(colnames(cortical.thickness)[1:34],attri[,"Full"])
setdiff(colnames(thick.ctl)[1:34],attri[,"Full"])

##Make heatmaps
cormat <- cor(cortical.thickness[,1:34], method="spearman")
par(oma = c(7, 0, 0, 7))
hc <- hclust(as.dist(1-cormat),method = 'ward')
heatmap(1-cormat,Rowv=as.dendrogram(hc),Colv=as.dendrogram(hc), scale='none')
dev.copy(pdf, file="cases_correlation_lh.pdf")
dev.off()

cormat.control <- cor(thick.ctl[,1:34], method="spearman")
par(oma = c(7, 0, 0, 7))
hc <- hclust(as.dist(1-cormat.control),method = 'ward')
heatmap(1-cormat.control,Rowv=as.dendrogram(hc),Colv=as.dendrogram(hc), scale='none')
dev.copy(pdf, file="controls_correlation_lh.pdf")
dev.off()

##All combination
allcmb <- t(combn(attri[,1],2))

##Excluding combination in "conn"
allcmb <- allcmb[!(paste(allcmb[,1],allcmb[,2]) %in% paste(conn[,1],conn[,3]))
                 & !(paste(allcmb[,1],allcmb[,2]) %in% paste(conn[,3], conn[,1])),]

##Make total combination table
allcmb <- cbind(as.data.frame(allcmb),"NonEdge")[,c(1,3,2)]
colnames(allcmb) <- c("Region1","Connection","Region2")
total.cmb <- rbind(conn,allcmb)

##Cor.test in each combination
correlation <- apply(total.cmb,1,function(x){
            cor.test(cortical.thickness[,colnames(cortical.thickness)==attri[attri$Abb==x[1],"Full"]],
            cortical.thickness[,colnames(cortical.thickness)==attri[attri$Abb==x[3],"Full"]],
            method="spearman")
})

correlation.control <- apply(total.cmb,1,function(x){
  cor.test(thick.ctl[,colnames(thick.ctl)==attri[attri$Abb==x[1],"Full"]],
           thick.ctl[,colnames(thick.ctl)==attri[attri$Abb==x[3],"Full"]],
           method="spearman")
})

##Linear test (case vs. control) in each combination
linear.case_ctl <- apply(total.cmb,1,function(x){tmptable <- rbind(
  data.frame(Region1=as.numeric(cortical.thickness[,colnames(cortical.thickness)==attri[attri$Abb==x[1],"Full"]]),
        Region2=as.numeric(cortical.thickness[,colnames(cortical.thickness)==attri[attri$Abb==x[3],"Full"]]),
        phenotype="case"), 
  data.frame(Region1=as.numeric(thick.ctl[,colnames(thick.ctl)==attri[attri$Abb==x[1],"Full"]]), 
        Region2=as.numeric(thick.ctl[,colnames(thick.ctl)==attri[attri$Abb==x[3],"Full"]]),
        phenotype="control"));
                                              summary(lm(Region1~Region2*phenotype, tmptable))})


##Make table including rho and p.value, and slope difference between
##cases and controls with p value

total.cmb <- data.frame(total.cmb,
                   Cases=t(sapply(correlation, "[", c("p.value","estimate"))),
                   Controls=t(sapply(correlation.control, "[", c("p.value","estimate"))),
                        t(sapply(linear.case_ctl, function(lm)
                          unlist(list("Case_control"=lm$coefficients["Region2:phenotypecontrol",c(1,4)])))))

for(i in 4:7){
  total.cmb[,i] <- as.numeric(total.cmb[,i])
}


pairwise.wilcox.test(total.cmb[,5], total.cmb[,2], p.adj = "bonf")

write.csv(total.cmb, file="Correlation_region_rho.csv")

plotresult <- ggplot(total.cmb, aes(total.cmb[,2],total.cmb[,5]))+
  ylab("Rho(cases)") +
  xlab(colnames(total.cmb)[2]) +
  opts(legend.position = "none") +
  geom_boxplot(aes(fill=total.cmb[,2]),outlier.size=0)+
  geom_jitter(alpha=0.4)
ggsave("Rho_plot(cases).pdf",width=10, height=10)

plotresult <- ggplot(total.cmb, aes(total.cmb[,2],total.cmb[,7]))+
  ylab("Rho(controls)") +
  xlab(colnames(total.cmb)[2]) +
  opts(legend.position = "none") +
  geom_boxplot(aes(fill=total.cmb[,2]),outlier.size=0)+
  geom_jitter(alpha=0.4)
ggsave("Rho_plot(controls).pdf",width=10, height=10)

plotresult <- ggplot(total.cmb, aes(-log10(total.cmb[,4]),total.cmb[,5],colour=total.cmb[,2]))+
  ylab("Rho(cases)") +
  xlab("-log10P") +
  geom_point()