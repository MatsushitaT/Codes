rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")

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

## make FID and IID
cortical.thickness <- cbind(FID = rownames(cortical.thickness),IID = rownames(cortical.thickness),cortical.thickness)

cort3sd <- cbind(FID = rownames(cort3sd),IID = rownames(cort3sd),cort3sd)

## phenotype file for PLINK
write.table(cortical.thickness, "phenotype.txt",col.names=T, row.names=F, quote=F)
write.table(cort3sd, "pheno3sd.txt",col.names=T, row.names=F, quote=F)
