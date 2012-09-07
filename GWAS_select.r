rm(list=ls())
setwd("~/Analysis/Cortex_thickness/GWAS/")
regionlist <- readLines("Regionlist.txt")
model <- data.frame(naming=c("DOM","REC","ADD"), dir=c("Linear_dominant","Linear_recessive","Linear_trend"))
frq <- read.table(file="plink.frq",header=T)

for(file in regionlist){cons <- strsplit(file, "_")
                        region <- paste(cons[[1]][1],cons[[1]][2],sep="_")
                        file.name <- paste("plink.",region,".assoc.linear",sep="")
                        selmod <- cons[[1]][length(cons[[1]])]
                        mod <- as.character(model[model$naming==selmod,"dir"])
                        adjust <- gsub(paste(region,"_",sep=""),"",file);adjust <- gsub(paste("_",selmod,sep=""),"",adjust)
                        d <- read.table(paste(adjust,"/",mod,"/",file.name,sep=""),header=T)
                        d <- d[d$TEST==selmod & !is.na(d$P),c("SNP","P")]
                        d <- merge(d,frq[,c("SNP","MAF")],all.x=T, all.y=F)
                        if(selmod=="REC"){d <- d[d$MAF >= 0.15,]} else {d <- d[d$MAF >= 0.005,]}
                        write.table(d[,c("SNP","P")],file=paste(file,".txt",sep=""),col.names=F,row.names=F,quote=F)
  
}