#set workdesk
setwd("~/Analysis/Cortex_thickness/GWAS/")
rm(list=ls())
filename <- "plink.lh_bankssts.assoc.linear"
#.assoc.linear file correspondent to each phenotype without -o option (plink.***.assoc.linear)
files <- list.files("Linear_dominant/")
for(filename in files){
  if(regexpr("^plink.*assoc.linear$",filename) <0){next}
  title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title);
  dom <- read.table(paste("Linear_dominant/",filename,sep=""), header=T);dom <- dom[dom$TEST=="DOM",];
  rec <- read.table(paste("Linear_recessive/",filename,sep=""), header=T);rec <- rec[rec$TEST=="REC",];
  trend <- read.table(paste("Linear_trend/",filename,sep=""), header=T);trend <- trend[trend$TEST=="ADD",];
  trend <- trend[trend$CHR %in% c(1:22,25),]
  minimump <- pmin(dom$P, rec$P, trend$P, na.rm=T)
  dom <- dom[!is.na(minimump),]; rec <- rec[!is.na(minimump),];trend <- trend[!is.na(minimump),];
  minimump <- minimump[!is.na(minimump)];
  data <- rbind(dom[dom$P==minimump,],rec[rec$P==minimump,],trend[trend$P==minimump,])
  data2 <- data[data$P<1e-5 & !is.na(data$P),];
  data2 <- data2[sort.list(data2$P, decreasing=F),]
  cat(paste("\n",title,"\n", sep=""),file = "SNP_list(-LogP>5).txt", append=T)
  sink("SNP_list(-LogP>5).txt",append=T)
  print(data2)
  sink()}