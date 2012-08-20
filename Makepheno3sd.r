rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")
"%w/o%" <- function(x, y) x[!x %in% y]

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
cortical.thickness <- cortical.thickness[rownames(cortical.thickness) %w/o% c("EPIC416","EPIC293"),]

cortical.thickness <- apply(cortical.thickness,2,function(x){
  m <- mean(x,na.rm=T)
  sd <- sd(x,na.rm=T)
  pmin(pmax(x,round(m-3*sd,2)),round(m+3*sd,2))})

## make FID and IID
cortical.thickness <- data.frame(FID = rownames(cortical.thickness),IID = rownames(cortical.thickness),cortical.thickness)
save(cortical.thickness, file="cortical.thickness.rdata")

## phenotype file for PLINK
write.table(cortical.thickness, "GWAS/phenotype.txt",col.names=T, row.names=F, quote=F)

## Covariate File
clinical.data <- read.table("20120807_CoriticalThickness.csv", sep=";", header=T)
clinical.data <- clinical.data[,c(2,79:81)]
clinical.data[,"GSKID"] <- paste("EPIC",clinical.data[,"GSKID"], sep="")
covariate <- merge(cortical.thickness,clinical.data,by.x="FID",by.y="GSKID", all.x=T)
covariate <- covariate[,c(1,2,71:73)]
write.table(covariate, "GWAS/covariates.txt",col.names=T, row.names=F, quote=F)

## File for Extract individuals
write.table(cortical.thickness[,1:2], "GWAS/extract_individual.txt",col.names=F, row.names=F, quote=F)