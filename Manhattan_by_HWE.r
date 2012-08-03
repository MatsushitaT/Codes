#set workdesk
setwd("~/Analysis/Cortex_thickness/GWAS/Linear_trend/")
rm(list=ls())

library("ggplot2")
library("RColorBrewer")

hwe <- read.table("plink.hwe", header = T)
freq <- read.table("plink.frq", header = T)
hwe <- hwe[hwe$TEST=="ALL",]

#.assoc.linear file correspondent to each phenotype without -o option (plink.***.assoc.linear)
files <- list.files()
for(filename in files){
  if(regexpr("^plink.*assoc.linear$",filename) <0){next}
  title <- gsub("plink.","",filename);title <- gsub(".assoc.linear","",title);
  data <- read.table(filename, header=T);data <- data[data$TEST=="ADD",];
  data$LogP <- -log10(data$P);
  data$BP <- data$BP/10000;
  hardythres <- c(1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10)
  for(hardy in hardythres){
    hwethres <- hwe[hwe$P>=hardy,c(2,9)];
    d <- merge(hwethres, data, by="SNP", x.all=T);
    d <- d[order(d$CHR,d$BP),]
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
      opts(title = paste("Distribution of marker effects (",title,")", " HWE test P >=", hardy, sep=""))
    ggsave(paste(title,"(HWE",hardy,").png",sep=""), width = 10, height = 5, dpi = 300)}}