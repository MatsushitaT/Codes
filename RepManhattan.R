#set workdesk
setwd("~/Analysis/Cortex_thickness/")
rm(list=ls())

library("ggplot2")
library("RColorBrewer")

#.assco file correspondent to each phenotype without -o option (plink.***.qassoc)
files <- list.files()
for(filename in files){
if(regexpr("^plink.*qassoc$",filename) <0){next}
  title <- gsub("plink.","",filename);title <- gsub(".qassoc","",title);
  d <- read.table(filename, header=T);
  d$LogP <- -log10(d$P);
  d$BP <- d$BP/10000;
  d2 <- d[d$LogP>5 & !is.na(d$LogP),]
  d2 <- d2[sort.list(d2$LogP, decreasing=T),]
  cat(paste("\n",title,"\n", sep=""),file = "results.txt", append=T)
  sink("results.txt",append=T)
  print(d2)
  sink()
  for(i in 1:25){d[d$CHR==i+1,"BP"]<-d[d$CHR==i+1,"BP"]+max(d[d$CHR==i,"BP"])};
  d$tick <- c(rep(0,length(nrow(d))));
  for(i in 1:26){
    d[d$CHR==i,"tick"] <- min(d[d$CHR==i,"BP"])+(max(d[d$CHR==i,"BP"])-min(d[d$CHR==i,"BP"]))/2};
  d[d$CHR==23,"CHR"] <- "X"; d[d$CHR==24,"CHR"] <- "Y"; d[d$CHR==25,"CHR"] <- "XY"; d[d$CHR==26,"CHR"] <- "MT";
  d$CHR <- as.factor(d$CHR);
  colours <- rep(c(brewer.pal(n = 7, name = "Set1")),5);
  manhattan <- ggplot(data=d, aes(BP, LogP, color=CHR)) +
    geom_point() +
    scale_colour_manual(values = colours) +
    opts(legend.position = "none") +
    scale_x_continuous(name = "Chromosome", breaks = unique(d$tick), labels = unique(d$CHR)) +
    ylab("-LogP") +
    opts(title = paste("Distribution of marker effects (",title,")", sep=""))
ggsave(paste(title,".png",sep=""), width = 10, height = 5, dpi = 300)}
