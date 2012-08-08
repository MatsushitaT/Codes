#set workdesk
setwd("~/Analysis/Cortex_thickness/GWAS/adjusted_by_age/Linear_dominant/")
rm(list=ls())

library("ggplot2")
library("RColorBrewer")

#.assoc.linear file correspondent to each phenotype without -o option (plink.***.assoc.linear)
files <- list.files()
for(filename in files){
  if(regexpr("^plink.*assoc.linear$",filename) <0){next}
  title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title);
  d <- read.table(filename, header=T);d <- d[d$TEST=="DOM",];
  d$LogP <- -log10(d$P);
  d$BP <- d$BP/10000; yline <- c(-log10(0.05/nrow(dom)),-log10(0.05/nrow(dom)/13),-log10(0.05/nrow(dom)/34));
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
    opts(title = paste("Distribution of marker effects (",title,")", sep="")) +
    geom_hline(yintercept=yline)
  ggsave(paste(title,".pdf",sep=""), width = 15, height = 7.5)}
