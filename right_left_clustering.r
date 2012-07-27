rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")

source("~/Analysis/Cortex_thickness/Codes/cor.test.2.r")

##### READ DATA
## read table
left <- read.table("cortical-thickness-left-hemi.txt",
                   sep="\t",header=T,quote="", comment.char="",as.is=T)

## name first four columns
names(left)[1:4] <- c("ID","date","msID1","msID2")

## make EPIC ID
left$ID <- paste("EPIC",left$ID, sep="")

## clean up date
left$date <- gsub("'","",left$date)
left$date <- gsub("=","/",left$date)
## clean up msID1
left$msID1  <- gsub("'","",left$msID1)
left$msID1  <- gsub(" ","",left$msID1)
left$msID1  <- gsub("/","_",left$msID1)
## msID1 and msID2 should match
summary(paste(left$msID1,"t1v",sep="_") == left$msID2) ## must be all true


## read table
right <- read.table("cortical-thickness-right-hemi.txt",
                    sep="\t",header=T,quote="", comment.char="",as.is=T)

## name first four columns
names(right)[1:4] <- c("ID","date","msID1","msID2")

## make EPIC ID
right$ID <- paste("EPIC",right$ID, sep="")

## clean up date
right$date <- gsub("'","",right$date)
right$date <- gsub("=","/",right$date)
## clean up msID1
right$msID1  <- gsub("'","",right$msID1)
right$msID1  <- gsub(" ","",right$msID1)
right$msID1  <- gsub("/","_",right$msID1)
## msID1 and msID2 should match
summary(paste(right$msID1,"t1v",sep="_") == right$msID2) ## must be all true

## join tables
summary(left$msID1 == right$msID1) ## must be all true
summary(left$msID2 == right$msID2) ## must be all true
summary(left$ID == right$ID) ## must be all true


left.cortices <- left[-1:-4]
rownames(left.cortices) <- left$ID
colnames(left.cortices) <- gsub("_thickness", "", colnames(left.cortices))

right.cortices <- right[-1:-4]
rownames(right.cortices) <- right$ID
colnames(right.cortices) <- gsub("_thickness", "", colnames(right.cortices))

## remove non-joined data
rm(left,right)

left.cortices[left.cortices==0] <- NA
right.cortices[right.cortices==0] <- NA

## missing rate in each individual
(ind.left <- apply(left.cortices, 1, function(x)sum(is.na(x))/length(x)))
(ind.right <- apply(right.cortices, 1, function(x)sum(is.na(x))/length(x)))

## missing rate in each cortex
(cereb.left <- apply(left.cortices, 2, function(x)sum(is.na(x))/length(x)))
(cereb.right <- apply(right.cortices, 2, function(x)sum(is.na(x))/length(x)))
left.cortices <- left.cortices[ind.left<0.05,cereb.left<0.05]
right.cortices <- right.cortices[ind.right<0.05,cereb.right<0.05]

##make distance matrix by correlation (pearson)
corpearson.left <- as.dist(1-cor(left.cortices,use='pair'))
corpearson.right <- as.dist(1-cor(right.cortices,use='pair'))
##make distance matrix by correlation (spearman)
corspearman.left <- as.dist(1-cor(left.cortices,use='pair', method="spearman"))
corspearman.right <- as.dist(1-cor(right.cortices,use='pair', method="spearman"))


methods <- c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")

## plot dendrogram by each method using distance information based on the correlation coefficient (pearson)
pdf("Lt_clustering_by_pearson.pdf")
  par(ps=8)
  for(m in methods)plot(hclust(corpearson.left,method=m),main=paste(m,"Pearson", sep=":"))
dev.off()

pdf("Rt_clustering_by_pearson.pdf")
par(ps=8)
for(m in methods)plot(hclust(corpearson.right,method=m),main=paste(m,"Pearson", sep=":"))
dev.off()

## plot dendrogram by each method using distance information based on the correlation coefficient (spearman)
pdf("Lt_clustering_by_spearman.pdf")
par(ps=8)
for(m in methods)plot(hclust(corspearman.left,method=m),main=paste(m,"Spearman", sep=":"))
dev.off()

pdf("Rt_clustering_by_spearman.pdf")
par(ps=8)
for(m in methods)plot(hclust(corspearman.right,method=m),main=paste(m,"Spearman", sep=":"))
dev.off()


## Brain regions included in each cluster (3-34) by Ward method
## based on the information of the correlation coefficient (pearson), and
## frequency of combination of regions showing r>0.5 and p<0.05 in each cluster.

for(i in 3:34){
  Cor.group <- cutree(hclust(corspearman.left, method="ward"),i)
  cat("\n","Divided by ", i, "groups","\n\n",append=T,file="Lt_Grouped_by_pearson.txt")
  rpfull <- 0; total<-0
  for(numg in 1:i){
    cat("Group",numg, "\n", append=T, file="Lt_Grouped_by_pearson.txt")
    sink("Lt_Grouped_by_pearson.txt", append=T)
    print(names(Cor.group[Cor.group==numg]),quote=F)
    sink()
    if(sum(Cor.group==numg)<2){next}
    cts <- cor.test.2(left.cortices[names(which(Cor.group==numg))],"pearson")
    rplist <- rbind(r=sapply(cts,  "[[" , "estimate"), p.value=sapply(cts,  "[[" , "p.value"))
    rpfull <- rpfull+ncol(as.data.frame(rplist[,rplist[1,]>0.5&rplist[2,]<0.05]))
    total <- total+ncol(rplist)
    cat("\n", "Combinations with r > 0.5 and p <0.05: ", ncol(as.data.frame(rplist[,rplist[1,]>0.5&rplist[2,]<0.05])),
        "in ", ncol(rplist), "\n\n", append = T, file="Lt_Grouped_by_pearson.txt")}
  cat("\n","Total: ",rpfull, "in", total, "\n",append = T, file = "Lt_Grouped_by_pearson.txt")}

for(i in 3:34){
  Cor.group <- cutree(hclust(corspearman.right, method="ward"),i)
  cat("\n","Divided by ", i, "groups","\n\n",append=T,file="Rt_Grouped_by_pearson.txt")
  rpfull <- 0; total<-0
  for(numg in 1:i){
    cat("Group",numg, "\n", append=T, file="Rt_Grouped_by_pearson.txt")
    sink("Rt_Grouped_by_pearson.txt", append=T)
    print(names(Cor.group[Cor.group==numg]),quote=F)
    sink()
    if(sum(Cor.group==numg)<2){next}
    cts <- cor.test.2(right.cortices[names(which(Cor.group==numg))],"pearson")
    rplist <- rbind(r=sapply(cts,  "[[" , "estimate"), p.value=sapply(cts,  "[[" , "p.value"))
    rpfull <- rpfull+ncol(as.data.frame(rplist[,rplist[1,]>0.5&rplist[2,]<0.05]))
    total <- total+ncol(rplist)
    cat("\n", "Combinations with r > 0.5 and p <0.05: ", ncol(as.data.frame(rplist[,rplist[1,]>0.5&rplist[2,]<0.05])),
        "in ", ncol(rplist), "\n\n", append = T, file="Rt_Grouped_by_pearson.txt")}
  cat("\n","Total: ",rpfull, "in", total, "\n",append = T, file = "Rt_Grouped_by_pearson.txt")}

## Brain regions included in each cluster (3-34) by Ward method
## based on the information of the correlation coefficient (spearman), and
## frequency of combination of regions showing r>0.5 and p<0.05 in each cluster.
for(i in 3:34){
  Corspear.group <- cutree(hclust(corspearman.left, method="ward"),i)
  cat("\n","Divided by ", i, "groups","\n\n",append=T,file="Lt_Grouped_by_spearman.txt")
  rpfull <- 0; total<-0
  for(numg in 1:i){
    cat("Group",numg, "\n", append=T, file="Lt_Grouped_by_spearman.txt")
    sink("Lt_Grouped_by_spearman.txt", append=T)
    print(names(Corspear.group[Corspear.group==numg]),quote=F)
    sink()
    if(sum(Corspear.group==numg)<2){next}
    cts <- cor.test.2(left.cortices[names(which(Corspear.group==numg))],"spearman")
    rplist <- rbind(r=sapply(cts,  "[[" , "estimate"), p.value=sapply(cts,  "[[" , "p.value"))
    rpfull <- rpfull+ncol(as.data.frame(rplist[,rplist[1,]>0.5&rplist[2,]<0.05]))
    total <- total+ncol(rplist)
    rplist <- rbind(r=sapply(cts,  "[[" , "estimate"), p.value=sapply(cts,  "[[" , "p.value"))
    cat("\n", "Combinations with r > 0.5 and p <0.05: ", ncol(as.data.frame(rplist[,rplist[1,]>0.5&rplist[2,]<0.05])),
        "in ", ncol(rplist), "\n\n", append = T, file="Lt_Grouped_by_spearman.txt")}
  cat("\n","Total: ",rpfull, "in", total, "\n",append = T, file = "Lt_Grouped_by_spearman.txt")}

for(i in 3:34){
  Corspear.group <- cutree(hclust(corspearman.right, method="ward"),i)
  cat("\n","Divided by ", i, "groups","\n\n",append=T,file="Rt_Grouped_by_spearman.txt")
  rpfull <- 0; total<-0
  for(numg in 1:i){
    cat("Group",numg, "\n", append=T, file="Rt_Grouped_by_spearman.txt")
    sink("Rt_Grouped_by_spearman.txt", append=T)
    print(names(Corspear.group[Corspear.group==numg]),quote=F)
    sink()
    if(sum(Corspear.group==numg)<2){next}
    cts <- cor.test.2(right.cortices[names(which(Corspear.group==numg))],"spearman")
    rplist <- rbind(r=sapply(cts,  "[[" , "estimate"), p.value=sapply(cts,  "[[" , "p.value"))
    rpfull <- rpfull+ncol(as.data.frame(rplist[,rplist[1,]>0.5&rplist[2,]<0.05]))
    total <- total+ncol(rplist)
    rplist <- rbind(r=sapply(cts,  "[[" , "estimate"), p.value=sapply(cts,  "[[" , "p.value"))
    cat("\n", "Combinations with r > 0.5 and p <0.05: ", ncol(as.data.frame(rplist[,rplist[1,]>0.5&rplist[2,]<0.05])),
        "in ", ncol(rplist), "\n\n", append = T, file="Rt_Grouped_by_spearman.txt")}
  cat("\n","Total: ",rpfull, "in", total, "\n",append = T, file = "Rt_Grouped_by_spearman.txt")}
