rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")
library(ggplot2)

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

## make table for ploting
tmp <- cortical.thickness[, c(rbind(1:34, 1:34+34))]
tmp2 <- data.frame(tmp, median=apply(tmp,1,median))
tmp2 <- as.data.frame(t(tmp2[order(tmp2$median),]))[1:68,]##Sort ID by the median

plotdata <- NULL
for(i in 1:ncol(tmp)){
  plotdata<-rbind(plotdata, data.frame(Region = rep(names(tmp)[i],
                                       nrow(tmp)),Thickness = tmp[,i],
                                       ID = rownames(tmp)))}

plotdata2 <- NULL
for(i in 1:ncol(tmp2)){
  plotdata2<-rbind(plotdata2, data.frame(ID = rep(names(tmp2)[i],
                                       nrow(tmp2)),Thickness = tmp2[,i],
                                       Region = rownames(tmp2)))}

## make boxplot
ploting <- ggplot(plotdata, aes(Region, Thickness))
ploting + geom_boxplot(aes(fill =Region)) +
opts(legend.position = "none") +
opts(axis.text.x=theme_text(angle=90, hjust=1))
ggsave("boxplot_by_region.png", width = 10, height = 8, dpi = 300)

for(i in 1:12)
{plotdata3 <- plotdata2[(1+68*50*(i-1)):(68*50*i),]
ploting2 <- ggplot(plotdata3, aes(ID, Thickness))
ploting2 + geom_boxplot() +
  opts(legend.position = "none") +
  opts(axis.text.x=theme_text(angle=90, hjust=1))
ggsave(paste("boxplot_by_individual",i,".png",sep=""), width = 15, height = 8, dpi = 400)}

##make median histogram in each individual
med <- apply(cortical.thickness,1,median)
hist(med, main="Distribution of median (green:mean, red:median)",breaks=50)
abline(v=mean  (med),lty=2,lwd=2,col='green')
abline(v=median(med),lty=2,lwd=2,col='red')
dev.copy(pdf, file="median_histogram.pdf")
dev.off()
