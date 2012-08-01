#set workdesk
setwd("~/Analysis/Cortex_thickness/GWAS/Linear_trend/")
rm(list=ls())

library("ggplot2")
library("RColorBrewer")

#.assoc.linear file correspondent to each phenotype without -o option (plink.***.assoc.linear)
files <- list.files()
for(filename in files){
if(regexpr("^plink.*assoc.linear$",filename) <0){next}
  title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title);
  d <- read.table(filename, header=T);d <- d[d$TEST=="ADD",];
  d$LogP <- -log10(d$P);
  d2 <- d[d$LogP>5 & !is.na(d$LogP),];
  d2 <- d2[sort.list(d2$LogP, decreasing=T),]
  d$BP <- d$BP/10000;
  cat(paste("\n",title,"\n", sep=""),file = "results.txt", append=T)
  sink("results.txt",append=T)
  print(d2)
  sink()
  d3 <- d[!is.na(d$P),]
  QQ<-rank(d3$P)/(nrow(d3)+1)
  png(file=paste(title,"_QQ.png",sep=""),width=800, height=800)
  plot(-log(QQ,10), d3$LogP, xla="expected -logP values", ylab= "observed -logP values", main=paste("QQ-plot: ", title, sep=""))
  abline(a=0, b=1)
  dev.off()
  for(i in 1:25){d[d$CHR==i+1,"BP"]<-d[d$CHR==i+1,"BP"]+max(d[d$CHR==i,"BP"])};
  d$tick <- c(rep(0,length(nrow(d))));
  for(i in 1:26){
    d[d$CHR==i,"tick"] <- min(d[d$CHR==i,"BP"])+(max(d[d$CHR==i,"BP"])-min(d[d$CHR==i,"BP"]))/2};
  d[d$CHR==23,"CHR"] <- "X"; d[d$CHR==24,"CHR"] <- "Y"; d[d$CHR==25,"CHR"] <- "XY"; d[d$CHR==26,"CHR"] <- "MT";
  d$CHR <- as.factor(d$CHR);
  colours <- rep(c(brewer.pal(n = 7, name = "Set1")),5);
  manhattan <- ggplot(data=d, aes(BP, LogP, color=CHR)) +
    geom_point(size=1) +
    scale_colour_manual(values = colours) +
    opts(legend.position = "none") +
    scale_x_continuous(name = "Chromosome", breaks = unique(d$tick), labels = unique(d$CHR)) +
    ylab("-LogP") +
    opts(title = paste("Distribution of marker effects (",title,")", sep=""))
ggsave(paste(title,".png",sep=""), width = 10, height = 5, dpi = 300)}
