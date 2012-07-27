rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")

library("maptools")

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

## missing rate in each individual
(ind <- apply(cortical.thickness, 1, function(x)sum(is.na(x))/length(x)))

## missing rate in each cortex
(cereb <- apply(cortical.thickness, 2, function(x)sum(is.na(x))/length(x)))
cortical.thickness <- cortical.thickness[ind<0.05,cereb<0.05]


##Principal component analysis
Region.pca <- prcomp(cortical.thickness,scale=T)

##Caluculate component loadings
loadings <- t(t(Region.pca$rotation)*Region.pca$sdev)

##Plot
pdf(file="PC.pdf")
screeplot(Region.pca)
dev.off()

plot(Region.pca$x[,c(1:2)])
pointLabel(Region.pca$x[,1:2],labels=rownames(Region.pca$x),cex=.5)


for(i in 1:3){
  pdf(file=paste("PCA_Loading_",i,".pdf",sep=""), width=10, height=10)
  plot(loadings[,combn(c(1:3),2)[,i]])
  pointLabel(loadings[,combn(c(1:3),2)[,i]],labels=rownames(loadings),cex=.5)
  abline(h=0,lty=2,lwd=2,col='red')
  dev.off()
   }