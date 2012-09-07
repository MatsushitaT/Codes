rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")


clinical.data <- read.table("20120807_CoriticalThickness.csv", sep=";", header=T)
clinical.data <- clinical.data[,c("GSKID","gender","ageatexam","diseasecourse","diseaseduration")]
clinical.data[,"GSKID"] <- paste("EPIC",clinical.data[,"GSKID"], sep="")

##### For cases
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

##Exclude cases without data
cortical.thickness <- cortical.thickness[ind<0.05,cereb<0.05]

##Cases without clinical data
excl <- setdiff(rownames(cortical.thickness),clinical.data$GSKID)
cortical.thickness <- cortical.thickness[!(rownames(cortical.thickness) %in% excl),]

##Check outliners
hist(apply(cortical.thickness, 1, median), breaks=100, 
     main="Distribution of median (case)",xlab="median" )
dev.copy(pdf, file="histogram(cases).pdf")
dev.off()

## Exclude outliners
cortical.thickness <- cortical.thickness[c(-which.max(apply(cortical.thickness, 1, median)),
                                         -which.min(apply(cortical.thickness, 1, median))),]

## Truncated to 3sd range
cortical.thickness <- apply(cortical.thickness,2,function(x){
	m <- mean(x,na.rm=T)
	sd <- sd(x,na.rm=T)
	pmin(pmax(x,round(m-3*sd,2)),round(m+3*sd,2))})


save(cortical.thickness, file="Final_dataset/cortical.thickness.rdata")

## make FID and IID
cortical.thickness <- cbind(FID = rownames(cortical.thickness),IID = rownames(cortical.thickness),cortical.thickness)

## phenotype file for PLINK
write.table(cortical.thickness, "Final_dataset/phenotype.txt",col.names=T, row.names=F, quote=F)

## Covariate File
covariate <- merge(cortical.thickness,clinical.data,by.x="FID",by.y="GSKID", all.x=T)
covariate <- covariate[,c(1,2,71:74)]
write.table(covariate, "Final_dataset/covariates.txt",col.names=T, row.names=F, quote=F)

## File for Extract individuals
write.table(cortical.thickness[,1:2], "Final_dataset/extract_individual.txt",col.names=F, row.names=F, quote=F)

#######For controls
thick.ctl <- read.csv("ctrls_4.5_thk.csv")
thick.ctl$fsid. <- gsub("'","",thick.ctl$fsid.)
thick.ctl$TYPE. <- gsub("'","",thick.ctl$TYPE.)
thick.ctl$Gender. <- gsub("'","",thick.ctl$Gender.)

#thick.ctl$"lh_mean(unknown)" <- apply(thick.ctl[,5:39],1,mean)
#thick.ctl$"rh_mean(unknown)" <- apply(thick.ctl[,40:74],1,mean)
#thick.ctl$lh_mean <- apply(thick.ctl[,6:39],1,mean)
#thick.ctl$rh_mean <- apply(thick.ctl[,41:74],1,mean)
#thick.ctl$"lh_median(unknown)" <- apply(thick.ctl[,5:39],1,median)
#thick.ctl$"rh_median(unknown)" <- apply(thick.ctl[,40:74],1,median)
#thick.ctl$lh_median <- apply(thick.ctl[,6:39],1,median)
#thick.ctl$rh_median <- apply(thick.ctl[,41:74],1,median)

ctl.clinical <- thick.ctl[,1:4]
row.names(thick.ctl) <- thick.ctl[,1]
thick.ctl <- thick.ctl[,c(-1:-5, -40,-75:-76)]
colnames(thick.ctl)[1:34] <- gsub("\\.","",colnames(thick.ctl)[1:34])
colnames(thick.ctl)[1:34] <- paste("lh_",colnames(thick.ctl)[1:34],sep="")
colnames(thick.ctl)[35:68] <- gsub("\\.\\.1","",colnames(thick.ctl)[35:68])
colnames(thick.ctl)[35:68] <- paste("rh_",colnames(thick.ctl)[35:68],sep="")


##Check outliners
hist(apply(thick.ctl, 1, median), breaks=30,
     main="Distribution of median (controls)",
     xlab="median",xlim=c(1.884,2.95))
dev.copy(pdf,file="histogram(controls).pdf")
dev.off()

## Truncated to 3sd range
thick.ctl <- apply(thick.ctl,2,function(x){
  m <- mean(x,na.rm=T)
  sd <- sd(x,na.rm=T)
  pmin(pmax(x,round(m-3*sd,2)),round(m+3*sd,2))})

save(thick.ctl, file="Final_dataset/control_thickness.rdata")

