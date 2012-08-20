rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")
library("ggplot2")


load("cortical.thickness.rdata")
clinical.data <- read.table("20120807_CoriticalThickness.csv", sep=";", header=T)
clinical.data <- clinical.data[,c("GSKID","gender","ageatexam","diseasecourse","diseaseduration")]
clinical.data[,"GSKID"] <- paste("EPIC",clinical.data[,"GSKID"], sep="")
data <- merge(cortical.thickness,clinical.data,by.x="FID",by.y="GSKID", all.x=T)
data <- cbind(data,median_lt=apply(data[,3:36], 1, median), median_rt=apply(data[,37:70], 1, median))
data <- cbind(data, median=apply(data[3:70], 1, median))
cor.test(data$ageatexam,data$diseaseduration,use="pair",method="spearman")
cat("Statistical analysis: median and disease duration by side. \n\n", file="med_disdur.txt", append=T)
summary(lm(lh_rostralanteriorcingulate ~ gender*ageatexam, data ))
sink("med_disdur.txt", append=T)

cat("\n Corralation between individual median and disease duration. \n")
cat("\n Linear regression \n\n Left median \n")
ans_lt <- lm(median_lt~diseaseduration,data)
print(summary(ans_lt))
cat("\n Right median \n")
ans_rt <- lm(median_rt~diseaseduration,data)
print(summary(ans_rt))

cat("\n Spearman's correlation coefficient \n\n Left median \n")
corr_lt <- cor.test(data$median_lt, data$diseaseduration, use="pair", method="spearman")
print(corr_lt)
cat("\n Right median \n")
corr_rt <- cor.test(data$median_rt, data$diseaseduration, use="pair", method="spearman")
print(corr_rt)
sink()

cat("Statistical analysis: Men vs. Women in each region. \n\n", file="region_dur.txt", append=T)

sink("region_dur.txt", append=T)
cat("\n Corralation between thickness in each region and disease duration. \n")
cat("\n Linear regression \n\n")
region_dur_lm <- apply(data[,c(3:70,75:76)],2,function(x)lm(x~diseaseduration,data))
for(i in 1:70){cat(colnames(data[,c(3:70,75:76)])[i],"\n");print(summary(region_dur_lm[[i]]))}

cat("\n Spearman's correlation coefficient \n\n")
region_dur_cor <- apply(data[,c(3:70,75:76)],2,function(x)cor.test(x, data$ageatexam, use="pair", method="spearman"))
for(i in 1:70){cat(colnames(data[,c(3:70,75:76)])[i],"\n");print(region_dur_cor[[i]])}
sink()

table <- cbind(
  t(sapply(region_dur_lm[c(1:34,69)], function(lm)
    unlist(list(lh.dur=summary(lm)$coefficients["diseaseduration",c(1,4)])))),
  t(sapply(region_dur_lm[c(1:34+34,70)], function(lm)
    unlist(list(rh.dur=summary(lm)$coefficients["diseaseduration",c(1,4)])))),
  t(sapply(region_dur_cor[c(1:34,69)],function(cor)
    c(lh.dur.cor.p=unname(cor$p.value),lh.dur.cor.rho=unname(cor$estimate)))),
  t(sapply(region_dur_cor[c(1:34+34,70)],function(cor)
    c(rh.dur.cor.p=unname(cor$p.value),rh.dur.cor.rho=unname(cor$estimate))))
)
rownames(table)<-gsub("_lt","",rownames(table))
rownames(table)<-gsub("lh_","",rownames(table))
write.csv(table, "diseaseduration_table.csv")


plotdata <- data.frame(median=c(data$median_lt, data$median_rt),
                       side=c(rep("left",nrow(data)),rep("right",nrow(data))),
                       duration=data$diseaseduration,
                       gender=data$gender)

plotdata <- data.frame(plotdata, gender_side=paste(plotdata$gender,plotdata$side, sep="+"))
plotdata2 <- plotdata[!is.na(plotdata$gender),]

median_dur_plot <- ggplot(plotdata2, aes(x=duration,y=median,colour=gender_side))+
  geom_point()+
  stat_smooth(method="lm",se=FALSE)+
  ylab("Median of cortical thickness")+
  xlab("Disease duration")
ggsave("median_duration.pdf",width=10, height=10)


for(i in 3:70){if(i==3){region.data <- data.frame(thickness=data[,i],data[,c(71,74)],region=colnames(data)[i])} else {
  region.data <- rbind(region.data,data.frame(thickness=data[,i],data[,c(71,74)],region=colnames(data)[i]))}}

region.data2 <- region.data[!is.na(region.data$gender),]

#specify a number of division (n <-)
n <- 9
for(t in 1:((68 %/% n)+1)){
  regionplot <- region.data[region.data$region %in% colnames(data[,3:70])[(n*(t-1)+1):(t*n)],]
  region_age <- ggplot(regionplot,aes(diseaseduration,thickness))+
    geom_point()+
    stat_smooth(method="lm",se=FALSE)+
    xlab("Disease duration")+
    facet_wrap(~region)
  ggsave(paste("region_duration",t,".pdf",sep=""), width=12, height=10)}