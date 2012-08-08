#set workdesk
setwd("~/Analysis/Cortex_thickness/GWAS/tmp")
rm(list=ls())

library("ggplot2")
library("RColorBrewer")

freq <- read.table("plink.frq", header = T)

#.assoc.linear file correspondent to each phenotype without -o option (plink.***.assoc.linear)
files <- list.files("adjusted_by_age/Linear_recessive/tmp/")
for(filename in files){
  if(regexpr("^plink.*assoc.linear$",filename) <0){next}
  title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title);
  data <- read.table(paste("adjusted_by_age/Linear_recessive/tmp/",filename,sep=""), header=T);data <- data[data$TEST=="REC",];
  data$LogP <- -log10(data$P);
  data$BP <- data$BP/10000;
  MAFthres <- c(0, 0.01, 0.02, 0.03, 0.04, 0.05)
  for(MAF in MAFthres){
    freqthres <- freq[freq$MAF>=MAF,c(2,5)];
    d <- merge(freqthres, data, by="SNP", x.all=T);
    d <- d[order(d$CHR,d$BP),];yline <- -log10(0.05/nrow(d));
    for(i in 1:21){d[d$CHR==i+1,"BP"]<-d[d$CHR==i+1,"BP"]+max(d[d$CHR==i,"BP"])};
    d[d$CHR==25,"BP"]<-d[d$CHR==25,"BP"]+max(d[d$CHR==22,"BP"]);
    d$tick <- c(rep(0,length(nrow(d))));
    for(i in c(1:22,25)){
      d[d$CHR==i,"tick"] <- min(d[d$CHR==i,"BP"])+(max(d[d$CHR==i,"BP"])-min(d[d$CHR==i,"BP"]))/2};
    d[d$CHR==25,"CHR"] <- "XY";
    d$CHR <- as.factor(d$CHR);
    colours <- rep(c(brewer.pal(n = 7, name = "Set1")),5);
    manhattan <- ggplot(data=d, aes(BP, LogP, color=CHR)) +
      geom_point(size=1) +
      scale_colour_manual(values = colours) +
      opts(legend.position = "none") +
      scale_x_continuous(name = "Chromosome", breaks = unique(d$tick), labels = unique(d$CHR)) +
      ylab("-LogP") +
      opts(title = paste("Distribution of marker effects (",title,")", " MAF >=", MAF, sep="")) +
      geom_hline(yintercept=yline)
    ggsave(paste(title,"(MAF",MAF,").png",sep=""), width = 10, height = 5, dpi = 300)}}
