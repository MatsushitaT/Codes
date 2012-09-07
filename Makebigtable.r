rm(list=ls())

dir <- "~/Analysis/Cortex_thickness/GWAS/"

wdlist <- c("adjusted_by_age", "adjusted_by_age_no3sd", "adjusted_by_sex_and_age", "median_adjusted_by_age", "median_adjusted_by_sex_and_age")
flist <- data.frame(folder=c("Linear_dominant", "Linear_recessive", "Linear_trend"), naming=c("DOM","REC", "ADD"))
data <- data.frame(CHR=NA, SNP=NA, BP=NA)

for(i in wdlist){
  setwd(paste(dir,i,sep=""))
  for(j in flist$folder){
    files <- list.files(as.character(j))
    for(filename in files){
      if(regexpr("^plink.*assoc.linear$",filename) <0){next}
      d <- read.table(paste(j,"/",filename, sep=""), header=T); d <- d[d$TEST==as.character(flist[flist$folder==j,"naming"]),]
      d$LogP <- -log10(d$P); d <- d[d$LogP > 4.5 & !is.na(d$LogP), c("CHR","SNP", "BP")]
      data <- merge(d,data[!is.na(data$CHR),],all=T)
      }
    }
  }


for(i in wdlist){
  setwd(paste(dir,i,sep=""))
  for(j in flist$folder){
    files <- list.files(as.character(j))
    for(filename in files){
      if(regexpr("^plink.*assoc.linear$",filename) <0){next}
      title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title)
      d <- read.table(paste(j,"/",filename, sep=""), header=T); d <- d[d$TEST==as.character(flist[flist$folder==j,"naming"]),]
      d$LogP <- -log10(d$P);
      colnames(d)[10]<-paste(title,"_",i,"_",flist[flist$folder==j,"naming"], sep="")
      data <- merge(data,d[,c(2,10)], all.x=T, all.y=F)
    }
  }
}

setwd("~/Analysis/Cortex_thickness/GWAS/")
frq <- read.table(file="plink.frq",header=T)
data <- merge(data, frq[,c(2,5)],all.x=T, all.y=F)
hwt <- read.table(file="plink.hwe", header=T);hwt <- hwt[hwt$TEST=="ALL",]
data <- merge(data, hwt[,c(2,6,9)], all.x=T, all.y=F)
colnames(data)[630] <- "HWT-P"

load("~/Analysis/snp129 x hg18 (condensed).rdata")
SNP129xHG18 <- SNP129xHG18[,c(3,7)]
data <- merge(data, SNP129xHG18,by.x="SNP", by.y="name", all.x=T, all.y=F)
data <- data[,c(1:3,628:631,4:627)]

save(data, file="Big_table_final.rdata")
write.csv(data,file="Big_table.csv")