#set workdesk
setwd("~/Analysis/Cortex_thickness/GWAS/median/adjusted_by_age/")
rm(list=ls())
load("~/Analysis/snp129 x hg18 (condensed).rdata")
folders <- list.files()

#.assoc.linear file correspondent to each phenotype without -o option (plink.***.assoc.linear)

for(workfolder in folders){
  if(workfolder == "Linear_dominant"){
    files <- list.files(workfolder)
    numSNP <- 0
    for(filename in files){
      if(regexpr("^plink.*assoc.linear$",filename) <0){next}
      title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title)
      d <- read.table(paste(workfolder,"/",filename,sep=""), header=T); d <- d[d$TEST=="DOM",]
      d$LogP <- -log10(d$P); n <- nrow(d)
      thres <- c(-log10(0.05/n),-log10(0.05/n/13),-log10(0.05/n/34))
      d <- d[!is.na(d$LogP) & d$LogP > thres[1],];
      d <- merge(d, SNP129xHG18, by.x="SNP", by.y="name", all.x=T, all.y=F);d <- d[,c(1:10,14:16)];d <- d[order(d$CHR),]
      numSNP <- numSNP + nrow(d)
      cat("\n",title,"\n",file = "SNP_above_threshold(dom).txt", append=T)
      cat("SNP number: ", n, ", ", "-Log P value thrseshold=", thres, "\n", file = "SNP_above_threshold(dom).txt", append=T)
      sink("SNP_above_threshold(dom).txt",append=T)
      for(i in 3:1){
        if(i==3){meetthres <- d[d$LogP > thres[i],];print(meetthres)
                 cat("------------------------------------------------------------------------------------------------------- -Log10P=",thres[i], "SNP number =", nrow(meetthres),"\n")} else {
                   meetthres <- d[d$LogP > thres[i] & d$LogP <= thres[i+1],]
                   print(meetthres);cat("------------------------------------------------------------------------------------------------------- -Log10P=",thres[i], "SNP number =", nrow(meetthres),"\n")}}
    cat("Total SNPs number=",nrow(d[d$LogP > thres[1],]),"\n", file = "SNP_above_threshold(dom).txt", append=T)
    sink()}
    cat("Total number of SNPs = ", numSNP, file="SNP_above_threshold(dom).txt", append=T)} else if(workfolder == "Linear_recessive"){
      files <- list.files(workfolder)
      numSNP <- 0
      for(filename in files){
        if(regexpr("^plink.*assoc.linear$",filename) <0){next}
        title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title)
        d <- read.table(paste(workfolder,"/",filename,sep=""), header=T); d <- d[d$TEST=="REC",]
        d$LogP <- -log10(d$P); n <- nrow(d)
        thres <- c(-log10(0.05/n),-log10(0.05/n/13),-log10(0.05/n/34))
        d <- d[!is.na(d$LogP) & d$LogP > thres[1],];
        d <- merge(d, SNP129xHG18, by.x="SNP", by.y="name", all.x=T, all.y=F); d <- d[,c(1:10,14:16)];d <- d[order(d$CHR),]
        numSNP <- numSNP + nrow(d)
        cat("\n",title,"\n",file = "SNP_above_threshold(rec).txt", append=T)
        cat("SNP number: ", n, ", ", "-Log P value thrseshold=", thres, "\n", file = "SNP_above_threshold(rec).txt", append=T)
        sink("SNP_above_threshold(rec).txt",append=T)
        for(i in 3:1){
          if(i==3){meetthres <- d[d$LogP > thres[i],];print(meetthres)
                   cat("------------------------------------------------------------------------------------------------------- -Log10P=",thres[i], "SNP number =", nrow(meetthres),"\n")} else {
                     meetthres <- d[d$LogP > thres[i] & d$LogP <= thres[i+1],]
                     print(meetthres);cat("------------------------------------------------------------------------------------------------------- -Log10P=",thres[i], "SNP number =", nrow(meetthres),"\n")}}
        cat("Total SNPs number=",nrow(d[d$LogP > thres[1],]),"\n", file = "SNP_above_threshold(rec).txt", append=T)
        sink()}
      cat("Total number of SNPs = ", numSNP, file="SNP_above_threshold(rec).txt", append=T)} else {
        files <- list.files(workfolder)
        numSNP <- 0
        for(filename in files){
          if(regexpr("^plink.*assoc.linear$",filename) <0){next}
          title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title)
          d <- read.table(paste(workfolder,"/",filename,sep=""), header=T); d <- d[d$TEST=="ADD",]
          d$LogP <- -log10(d$P); n <- nrow(d)
          thres <- c(-log10(0.05/n),-log10(0.05/n/13),-log10(0.05/n/34))
          d <- d[!is.na(d$LogP) & d$LogP > thres[1],];
          d <- merge(d, SNP129xHG18, by.x="SNP", by.y="name", all.x=T, all.y=F); d <- d[,c(1:10,14:16)];d <- d[order(d$CHR),]
          numSNP <- numSNP + nrow(d)
          cat("\n",title,"\n",file = "SNP_above_threshold(add).txt", append=T)
          cat("SNP number: ", n, ", ", "-Log P value thrseshold=", thres, "\n", file = "SNP_above_threshold(add).txt", append=T)
          sink("SNP_above_threshold(add).txt",append=T)
          for(i in 3:1){
            if(i==3){meetthres <- d[d$LogP > thres[i],];print(meetthres)
                     cat("------------------------------------------------------------------------------------------------------- -Log10P=",thres[i], "SNP number =", nrow(meetthres),"\n")} else {
                       meetthres <- d[d$LogP > thres[i] & d$LogP <= thres[i+1],]
                       print(meetthres);cat("------------------------------------------------------------------------------------------------------- -Log10P=",thres[i], "SNP number =", nrow(meetthres),"\n")}}
          cat("Total SNPs number=",nrow(d[d$LogP > thres[1],]),"\n", file = "SNP_above_threshold(add).txt", append=T)
          sink()}
        cat("Total number of SNPs = ", numSNP, file="SNP_above_threshold(add).txt", append=T)}
  }