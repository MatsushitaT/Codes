rm(list=ls())
setwd("e:/documents/cortical thickness")

##### READ DATA
## read table
left <- read.table("data/cortical-thickness-left-hemi.txt",
                   sep="\t",header=T,quote="", comment.char="",as.is=T)
## drop first column
left <- left[-1]
## name next three columns
names(left)[1:3] <- c("date","msID1","msID2")
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
right <- read.table("data/cortical-thickness-right-hemi.txt",
                   sep="\t",header=T,quote="", comment.char="",as.is=T)
## drop first column
right <- right[-1]
## name next three columns
names(right)[1:3] <- c("date","msID1","msID2")
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
cortical.thickness <- cbind(left[-1:-3],right[-1:-3])
rownames(cortical.thickness) <- left$msID1
colnames(cortical.thickness) <- gsub("_thickness","",colnames(cortical.thickness))

## remove non-joined data
rm(left,right)

## create a list of brain areas
brain.areas <- unique(sapply(strsplit(names(cortical.thickness),"_"),'[',2))
summary(colnames(cortical.thickness) ==
        c(t(outer(c("lh","rh"),brain.areas,paste, sep = "_"))) ) ## should be all true

##### PLOT HISTOGRAMS
source("e:/r/plot.linear.model.r")
dir.create("plots",F)
for(brain.area in brain.areas){
  par(mfcol = c(2,3))
  left.thickness  <- cortical.thickness[[paste("lh",sep="_",brain.area)]]
  right.thickness <- cortical.thickness[[paste("rh",sep="_",brain.area)]]
  nonMissing <- left.thickness > 0.01 & right.thickness > 0.01
  
  ## left distribution
  hist(left.thickness,
       breaks = 100, main = "left hemisphere",
       xlab = paste(brain.area,'thickness'))
  ## right distribution
  hist(right.thickness,
       breaks = 100, main = "right hemisphere",
       xlab = paste(brain.area,'thickness'))
  ## left versus right
  plot.linear.model(pch = 20,
       left.thickness,
       right.thickness,
       xlab="left",ylab="right")
  abline(0,1,col = 'red',lty = 2)
  ## left versus right (missing data removed)
  plot.linear.model(pch = 20,
       left.thickness[nonMissing],
       right.thickness[nonMissing],
       xlab="left",ylab="right")
  abline(0,1,col = 'red',lty = 2)
  ## left versus right (Smooth)
  smoothScatter(pch = 20,
       left.thickness[nonMissing],
       right.thickness[nonMissing],
       xlab="left",ylab="right")
  abline(0,1,col = 'red',lty = 2)
  ## plot histogram of differences
  difference <-  left.thickness[nonMissing]-
                right.thickness[nonMissing]
  hist(difference,
       xlab="left - right",breaks = 50,
       main = ("green = mean ,  red = median"))
  abline(v=0,lty=1,lwd=2,col='black')
  abline(v=mean  (difference),lty=2,lwd=2,col='green')
  abline(v=median(difference),lty=2,lwd=2,col='red')
  ## save plot
  savePlot(paste(sep="","plots/","distribution and left-right correlation"," ",brain.area),'jpeg')
  rm(brain.area,left.thickness,right.thickness, nonMissing,difference)
  }

##### calculate pairwise correlation
tmp <- cortical.thickness[, c(rbind(1:34, 1:34+34))]
tmp[tmp==0] <- NA
cormat <- cor(tmp,use='pair')
image(1-cormat)
hc <- hclust(as.dist(1-cormat),method = 'ward')
heatmap(1-cormat,Rowv=NA,Colv=NA, scale='none')
heatmap(1-cormat, scale='none')
heatmap(1-cormat,Rowv=as.dendrogram(hc),Colv=as.dendrogram(hc), scale='none')
savePlot("plots/correlations between thicknesses",'eps')
hist(c(cormat),breaks = 100)