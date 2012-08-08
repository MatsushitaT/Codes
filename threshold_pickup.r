#set workdesk
setwd("~/Analysis/Cortex_thickness/GWAS/adjusted_by_age/")
rm(list=ls())

folders <- list.files()

#.assoc.linear file correspondent to each phenotype without -o option (plink.***.assoc.linear)

for(workfolder in folders){
  if(workfolder == "Linear_dominant"){
    files <- list.files(workfolder)
    numSNP <- 0
    for(filename in files){
      if(regexpr("^plink.*assoc.linear$",filename) <0){next}
      title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title);
      dom <- read.table(paste(workfolder,"/",filename,sep=""), header=T);dom <- dom[dom$TEST=="DOM",];
  dom$LogP <- -log10(dom$P)
  domthres <- c(-log10(0.05/nrow(dom)),-log10(0.05/nrow(dom)/13),-log10(0.05/nrow(dom)/34))
  numSNP <- numSNP + nrow(dom[!is.na(dom$LogP) & dom$LogP > domthres[1],])
  cat("\n",title,"\n",file = "SNP_above_threshold(dom).txt", append=T)
  cat("SNP number: ", nrow(dom), ", ", "-Log P value thrseshold=", domthres, "\n", file = "SNP_above_threshold(dom).txt", append=T)
  sink("SNP_above_threshold(dom).txt",append=T)
  for(i in 3:1){if(i==3){meetthres <- dom[!is.na(dom$LogP) & dom$LogP > domthres[i],];print(meetthres);
                         cat("----------------------------------------------------------------------------------- -Log10P=",domthres[i], "SNP number =", nrow(meetthres),"\n")} else {
                           meetthres <- dom[!is.na(dom$LogP) & dom$LogP > domthres[i] & dom$LogP <= domthres[i+1],];
                           print(meetthres);cat("----------------------------------------------------------------------------------- -Log10P=",domthres[i], "SNP number =", nrow(meetthres),"\n")}}
  cat("Total SNPs number=",nrow(dom[!is.na(dom$LogP) & dom$LogP > domthres[1],]),"\n", file = "SNP_above_threshold(dom).txt", append=T)
  sink()}
cat("Total number of SNPs = ", numSNP, file="SNP_above_threshold(dom).txt", append=T)} else if(workfolder == "Linear_recessive"){
  files <- list.files(workfolder)
  numSNP <- 0
  for(filename in files){
    if(regexpr("^plink.*assoc.linear$",filename) <0){next}
    title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title);
    dom <- read.table(paste(workfolder,"/",filename,sep=""), header=T);dom <- dom[dom$TEST=="REC",];
    dom$LogP <- -log10(dom$P)
    domthres <- c(-log10(0.05/nrow(dom)),-log10(0.05/nrow(dom)/13),-log10(0.05/nrow(dom)/34))
    numSNP <- numSNP + nrow(dom[!is.na(dom$LogP) & dom$LogP > domthres[1],])
    cat("\n",title,"\n",file = "SNP_above_threshold(rec).txt", append=T)
    cat("SNP number: ", nrow(dom), ", ", "-Log P value thrseshold=", domthres, "\n", file = "SNP_above_threshold(rec).txt", append=T)
    sink("SNP_above_threshold(rec).txt",append=T)
    for(i in 3:1){if(i==3){meetthres <- dom[!is.na(dom$LogP) & dom$LogP > domthres[i],];print(meetthres);
                           cat("----------------------------------------------------------------------------------- -Log10P=",domthres[i], "SNP number =", nrow(meetthres),"\n")} else {
                             meetthres <- dom[!is.na(dom$LogP) & dom$LogP > domthres[i] & dom$LogP <= domthres[i+1],];
                             print(meetthres);cat("----------------------------------------------------------------------------------- -Log10P=",domthres[i], "SNP number =", nrow(meetthres),"\n")}}
    cat("Total SNPs number=",nrow(dom[!is.na(dom$LogP) & dom$LogP > domthres[1],]),"\n", file = "SNP_above_threshold(rec).txt", append=T)
    sink()}
  cat("Total number of SNPs = ", numSNP, file="SNP_above_threshold(rec).txt", append=T)} else {
    files <- list.files(workfolder)
    numSNP <- 0
    for(filename in files){
      if(regexpr("^plink.*assoc.linear$",filename) <0){next}
      title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title);
      dom <- read.table(paste(workfolder,"/",filename,sep=""), header=T);dom <- dom[dom$TEST=="ADD",];
      dom$LogP <- -log10(dom$P)
      domthres <- c(-log10(0.05/nrow(dom)),-log10(0.05/nrow(dom)/13),-log10(0.05/nrow(dom)/34))
      numSNP <- numSNP + nrow(dom[!is.na(dom$LogP) & dom$LogP > domthres[1],])
      cat("\n",title,"\n",file = "SNP_above_threshold(add).txt", append=T)
      cat("SNP number: ", nrow(dom), ", ", "-Log P value thrseshold=", domthres, "\n", file = "SNP_above_threshold(add).txt", append=T)
      sink("SNP_above_threshold(add).txt",append=T)
      for(i in 3:1){if(i==3){meetthres <- dom[!is.na(dom$LogP) & dom$LogP > domthres[i],];print(meetthres);
                             cat("----------------------------------------------------------------------------------- -Log10P=",domthres[i], "SNP number =", nrow(meetthres),"\n")} else {
                               meetthres <- dom[!is.na(dom$LogP) & dom$LogP > domthres[i] & dom$LogP <= domthres[i+1],];
                               print(meetthres);cat("----------------------------------------------------------------------------------- -Log10P=",domthres[i], "SNP number =", nrow(meetthres),"\n")}}
      cat("Total SNPs number=",nrow(dom[!is.na(dom$LogP) & dom$LogP > domthres[1],]),"\n", file = "SNP_above_threshold(add).txt", append=T)
      sink()}
    cat("Total number of SNPs = ", numSNP, file="SNP_above_threshold(add).txt", append=T)}
  }