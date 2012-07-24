rm(list=ls())
setwd("~/Analysis/Gene_Cortex")

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
cormat <- as.dist(1-cor(tmp,use='pair'))
eucmat <- dist(t(cortical.thickness))
manmat <- dist(t(cortical.thickness), method="manhattan")

methods <- c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")

pdf("clustering_by_cor.pdf")
  par(ps=8)
  for(m in methods)plot(hclust(cormat,method=m),main=paste(m,"Cor", sep=":"))
dev.off()

pdf("clustering_by_Euc.pdf")
par(ps=8)
for(m in methods)plot(hclust(eucmat,method=m),main=paste(m,"Euclidean", sep=":"))
dev.off()

pdf("clustering_by_Manhattan.pdf")
par(ps=8)
for(m in methods)plot(hclust(manmat, method=m),main=paste(m,"Manhattan", sep=":"))
dev.off()
table(Cor.group)

for(i in 3:9){
  Cor.group <- cutree(hclust(cormat, method="ward"),i)
  oneway <- lapply(as.data.frame(t(cortical.thickness)),function(x)oneway.test(x~Cor.group))
  nump <- sum(sapply(oneway, function(x)x$p.value<0.05))
  cat("\n","Divided by ", i, "groups","\n\n",append=T,file="Grouped_by_cor.txt")
  for(numg in 1:i){
  cat("Group",numg, "\n", append=T, file="Grouped_by_cor.txt")
  sink("Grouped_by_cor.txt", append=T)
  print(names(Cor.group[Cor.group==numg]),quote=F)
  sink()}
  cat("\n","p<0.05: ", nump, "in", nrow(cortical.thickness), "\n\n", append=T, file="Grouped_by_cor.txt")}

for(i in 3:9){
  Euc.group <- cutree(hclust(eucmat, method="ward"),i)
  oneway <- lapply(as.data.frame(t(cortical.thickness)),function(x)oneway.test(x~Cor.group))
  nump <- sum(sapply(oneway, function(x)x$p.value<0.05))
  cat("\n","Divided by ", i, "groups","\n\n",append=T,file="Grouped_by_Euc.txt")
  for(numg in 1:i){
    cat("Group",numg, "\n", append=T, file="Grouped_by_Euc.txt")
    sink("Grouped_by_Euc.txt", append=T)
    print(names(Euc.group[Euc.group==numg]),quote=F)
    sink()}
  cat("\n","p<0.05: ", nump, "in", nrow(cortical.thickness), "\n\n", append=T, file="Grouped_by_Euc.txt")}
