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

cortical.thickness <- cbind(left[,-1:-4],right[-1:-4])
rownames(cortical.thickness) <- left$ID
colnames(cortical.thickness) <- gsub("_thickness","",colnames(cortical.thickness))

## remove non-joined data
rm(left,right)

cortical.thickness[cortical.thickness==0] <- NA

cort3sd <- cortical.thickness
cort3sd <- apply(cort3sd,2,function(x){
  m <- mean(x,na.rm=T)
  sd <- sd(x,na.rm=T)
  pmin(pmax(x,m-3*sd),m+3*sd)})

## missing rate in each individual
(ind <- apply(cortical.thickness, 1, function(x)sum(is.na(x))/length(x)))
(ind2 <- apply(cort3sd, 1, function(x)sum(is.na(x))/length(x)))

## missing rate in each cortex
(cereb <- apply(cortical.thickness, 2, function(x)sum(is.na(x))/length(x)))
(cereb2 <- apply(cort3sd, 2, function(x)sum(is.na(x))/length(x)))
cortical.thickness <- cortical.thickness[ind<0.05,cereb<0.05]
cort3sd <- cort3sd[ind2<0.05,cereb2<0.05]

tmp <- cortical.thickness[, c(rbind(1:34, 1:34+34))]
##make distance matrix by correlation (pearson)
cormat <- as.dist(1-cor(tmp,use='pair'))
##make distance matrix by correlation (spearman)
corspear <- as.dist(1-cor(tmp,use='pair', method="spearman"))
##make distance matrix by Euclidean method
eucmat <- dist(t(cortical.thickness))
##make distance matrix by Manhattan method
manmat <- dist(t(cortical.thickness), method="manhattan")

methods <- c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")

## plot dendrogram by each method using distance information based on the correlation coefficiency (pearson)
pdf("clustering_by_pearson.pdf")
  par(ps=8)
  for(m in methods)plot(hclust(cormat,method=m),main=paste(m,"Pearson", sep=":"))
dev.off()

## plot dendrogram by each method using distance information based on the correlation coefficiency (spearman)
pdf("clustering_by_spearman.pdf")
par(ps=8)
for(m in methods)plot(hclust(corspear,method=m),main=paste(m,"Spearman", sep=":"))
dev.off()

## plot dendrogram by each method using distance information by Euclidean method
pdf("clustering_by_Euc.pdf")
par(ps=8)
for(m in methods)plot(hclust(eucmat,method=m),main=paste(m,"Euclidean", sep=":"))
dev.off()

## plot dendrogram by each method using distance information by Manhattan method
pdf("clustering_by_Manhattan.pdf")
par(ps=8)
for(m in methods)plot(hclust(manmat, method=m),main=paste(m,"Manhattan", sep=":"))
dev.off()

## Brain regions included in each cluster (3-34) by Ward method
## based on the information of the correlation coefficiency (pearson), and
## frequency of combination of regions showing r>0.5 and p<0.05 in each cluster.

for(i in 3:34){
  Cor.group <- cutree(hclust(cormat, method="ward"),i)
  cat("\n","Divided by ", i, "groups","\n\n",append=T,file="Grouped_by_pearson.txt")
  for(numg in 1:i){
    cat("Group",numg, "\n", append=T, file="Grouped_by_pearson.txt")
    sink("Grouped_by_pearson.txt", append=T)
    print(names(Cor.group[Cor.group==numg]),quote=F)
    sink()
    if(sum(Cor.group==numg)<2){next}
    cts <- cor.test.2(cortical.thickness[names(which(Cor.group==numg))],"pearson")
    rplist <- rbind(r=sapply(cts,  "[[" , "estimate"), p.value=sapply(cts,  "[[" , "p.value"))
    cat("\n", "Combination with r > 0.5 and p <0.05: ", ncol(as.data.frame(rplist[,rplist[1,]>0.5&rplist[2,]<0.05])),
        "in ", ncol(rplist), "\n\n", append = T, file="Grouped_by_pearson.txt")}}

## Brain regions included in each cluster (3-34) by Ward method
## based on the information of the correlation coefficiency (spearman), and
## frequency of combination of regions showing r>0.5 and p<0.05 in each cluster.
for(i in 3:34){
  Corspear.group <- cutree(hclust(corspear, method="ward"),i)
  cat("\n","Divided by ", i, "groups","\n\n",append=T,file="Grouped_by_spearman.txt")
  for(numg in 1:i){
    cat("Group",numg, "\n", append=T, file="Grouped_by_spearman.txt")
    sink("Grouped_by_spearman.txt", append=T)
    print(names(Corspear.group[Corspear.group==numg]),quote=F)
    sink()
    if(sum(Corspear.group==numg)<2){next}
    cts <- cor.test.2(cortical.thickness[names(which(Corspear.group==numg))],"spearman")
    rplist <- rbind(r=sapply(cts,  "[[" , "estimate"), p.value=sapply(cts,  "[[" , "p.value"))
    cat("\n", "Combination with r > 0.5 and p <0.05: ", ncol(as.data.frame(rplist[,rplist[1,]>0.5&rplist[2,]<0.05])),
        "in ", ncol(rplist), "\n\n", append = T, file="Grouped_by_spearman.txt")}}


## Brain regions included in each cluster (3-9) by Ward method
## based on the Euclidean distance, and
## frequency of samples showing p<0.05 diffrence between each cluster
for(i in 3:9){
  Euc.group <- cutree(hclust(eucmat, method="ward"),i)
  oneway <- lapply(as.data.frame(t(cortical.thickness)),function(x)oneway.test(x~Euc.group))
  nump <- sum(sapply(oneway, function(x)x$p.value<0.05))
  cat("\n","Divided by ", i, "groups","\n\n",append=T,file="Grouped_by_Euc.txt")
  for(numg in 1:i){
    cat("Group",numg, "\n", append=T, file="Grouped_by_Euc.txt")
    sink("Grouped_by_Euc.txt", append=T)
    print(names(Euc.group[Euc.group==numg]),quote=F)
    sink()}
  cat("\n","p<0.05: ", nump, "in", nrow(cortical.thickness), "\n\n", append=T, file="Grouped_by_Euc.txt")}

## Brain regions included in each cluster (3-9) by Kmeans method
## and frequency of samples showing p<0.05 diffrence between each cluster
for(i in 3:9){
  ans <- kmeans(t(cortical.thickness),i, nstart=10000)
  oneway <- lapply(as.data.frame(t(cortical.thickness)),function(x)oneway.test(x~ans[[1]]))
  nump <- sum(sapply(oneway, function(x)x$p.value<0.05))
  cat("\n","Divided by ", i, "groups","\n\n",append=T,file="Grouped_by_Kmeans.txt")
  for(numg in 1:i){
    cat("Group",numg, "\n", append=T, file="Grouped_by_Kmeans.txt")
    sink("Grouped_by_Kmeans.txt", append=T)
    print(names(ans[[1]][ans[[1]]==numg]),quote=F)
    sink()}
  cat("\n","p<0.05: ", nump, "in", nrow(cortical.thickness), "\n\n", append=T, file="Grouped_by_Kmeans.txt")}



#cor.test.2(cortical.thickness[c("lh_caudalanteriorcingulate",
#                                "lh_caudalmiddlefrontal",
#                                "rh_bankssts",
#                                "lh_parstriangularis")])
#which(Euc.group==8)
#names(which(Euc.group==8))
#cts <- cor.test.2(cortical.thickness[names(which(Euc.group==8))])
#sapply(cts,  "[[" , "estimate")
#sapply(cts,  "[[" , "p.value")