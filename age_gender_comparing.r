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

cat("Statistical analysis: Men vs. Women in individual median by side. \n\n", file="gender_age.txt", append=T)

sink("gender_age.txt", append=T)
cat("Median's mean in each hemisphere by gender.\n")
cat("\n Left hemisphere \n")
print(gender_mean_lt <- tapply(data$median_lt,data$gender, mean))
cat("\n Right hemisphere \n")
print(gender_mean_rt <- tapply(data$median_rt,data$gender, mean))


cat("\n Lt. medaian: Wilcoxon test by gender \n\n")
print(wilcox.test(median_lt~gender,data))
cat("\n Lt. medaian: t-test by gender \n\n")
print(t.test(median_lt~gender,data))


cat("\n Rt. medaian: Wilcoxon test \n\n")
print(wilcox.test(median_rt~gender,data))
cat("\n Rt. medaian: t-test by gender \n\n")
print(t.test(median_rt~gender,data))

cat("\n Corralation between individual median and age at examination. \n")
cat("\n Linear regression \n\n Left median \n")
ans_lt <- lm(median_lt~ageatexam,data)
print(summary(ans_lt))
cat("\n Right median \n")
ans_rt <- lm(median_rt~ageatexam,data)
print(summary(ans_rt))

cat("\n Spearman's correlation coefficient \n\n Left median \n")
corr_lt <- cor.test(data$median_lt, data$ageatexam, use="pair", method="spearman")
print(corr_lt)
cat("\n Right median \n")
corr_rt <- cor.test(data$median_rt, data$ageatexam, use="pair", method="spearman")
print(corr_rt)
sink()

cat("Statistical analysis: Men vs. Women in each region. \n\n", file="region_gender_age.txt", append=T)

sink("region_gender_age.txt", append=T)
cat("Mean in each region by gender.\n")
print(region_gender_mean <- apply(data[,3:70],2,function(x)tapply(x,data$gender, mean)))

cat("\n Wilcoxon test in each region by gender \n\n")
print(region_gender_wilcoxon <- apply(data[,c(3:70,75:76)],2,function(x)wilcox.test(x~gender,data)))
cat("\n t-test in each region by gender \n\n")
print(region_gender_ttest <- apply(data[,c(3:70,75:76)],2,function(x)t.test(x~gender,data)))

cat("\n Corralation between thickness in each region and age at examination. \n")
cat("\n Linear regression \n\n")
region_age_lm <- apply(data[,c(3:70,75:76)],2,function(x)lm(x~ageatexam,data))
for(i in 1:68){cat(colnames(data)[i+2],"\n");print(summary(region_age_lm[[i]]))}

cat("\n Spearman's correlation coefficient \n\n")
region_age_cor <- apply(data[,c(3:70,75:76)],2,function(x)cor.test(x, data$ageatexam, use="pair", method="spearman"))
for(i in 1:68){cat(colnames(data)[i+2],"\n");print(region_age_cor[[i]])}
sink()

table <- cbind(
  t(sapply(region_gender_ttest[c(1:34,69)], function(tt)
  c(lh.t.test.p=tt$p.value, lh.diff=unname(diff(tt$estimate))))),
               lh.wilcox.p=sapply(region_gender_wilcoxon[c(1:34,69)],"[[","p.value"),
               t(sapply(region_gender_ttest[c(1:34+34,70)], function(tt)
                 c(rh.t.test.p=tt$p.value, rh.diff=unname(diff(tt$estimate))))),
               rh.wilcox.p=sapply(region_gender_wilcoxon[c(1:34+34,70)],"[[","p.value"),
               t(sapply(region_age_lm[c(1:34,69)], function(lm)
                 unlist(list(lh.age=summary(lm)$coefficients["ageatexam",c(1,4)])))),
               t(sapply(region_age_lm[c(1:34+34,70)], function(lm)
                 unlist(list(rh.age=summary(lm)$coefficients["ageatexam",c(1,4)])))),
               t(sapply(region_age_cor[c(1:34,69)],function(cor)
                 c(lh.age.cor.p=unname(cor$p.value),lh.age.cor.rho=unname(cor$estimate)))),
               t(sapply(region_age_cor[c(1:34+34,70)],function(cor)
                 c(rh.age.cor.p=unname(cor$p.value),rh.age.cor.rho=unname(cor$estimate))))
)
rownames(table)<-gsub("_lt","",rownames(table))
rownames(table)<-gsub("lh_","",rownames(table))
write.csv(table, "age_gender_region_table.csv")


plotdata <- data.frame(median=c(data$median_lt, data$median_rt),
                       side=c(rep("left",nrow(data)),rep("right",nrow(data))),
                       age_at_exam=data$ageatexam,
                       gender=data$gender)

plotdata <- data.frame(plotdata, gender_side=paste(plotdata$gender,plotdata$side, sep="+"))
plotdata2 <- plotdata[!is.na(plotdata$gender),]
median_age_plot <- ggplot(plotdata2, aes(x=age_at_exam,y=median,colour=gender_side))+
               geom_point()+
               stat_smooth(method="lm",se=FALSE)+
               ylab("Median of cortical thickness")+
               xlab("Age at Exam.")
ggsave("median_age2.pdf",width=10, height=10)


median_gender <- ggplot(plotdata2, aes(gender,median))+
                geom_boxplot(aes(fill=gender))+geom_jitter(alpha=0.2)+
                facet_wrap(~side)
ggsave("median_gender.pdf", width=8, height=8)


for(i in 3:70){if(i==3){region.data <- data.frame(thickness=data[,i],data[,71:72],region=colnames(data)[i])} else {
region.data <- rbind(region.data,data.frame(thickness=data[,i],data[,71:72],region=colnames(data)[i]))}}

region.data2 <- region.data[!is.na(region.data$gender),]

#specify a number of division (n <-)
n <- 10
for(t in 1:((68 %/% n)+1)){
  regionplot <- region.data2[region.data2$region %in% colnames(data[,3:70])[(n*(t-1)+1):(t*n)],]
  region_gender <- ggplot(regionplot, aes(gender,thickness))+
    geom_boxplot(aes(fill=gender))+
    geom_jitter(alpha=0.2)+
    opts(legend.position = "none") +
    facet_wrap(~region, ncol=n)
  ggsave(paste("region_gender",t,".pdf",sep=""), width=15, height=8)}


#specify a number of division (n <-)
n <- 9
for(t in 1:((68 %/% n)+1)){
  regionplot <- region.data[region.data$region %in% colnames(data[,3:70])[(n*(t-1)+1):(t*n)],]
  region_age <- ggplot(regionplot,aes(ageatexam,thickness))+
  geom_point()+
  stat_smooth(method="lm",se=FALSE)+
  xlab("Age at exam.")+
  facet_wrap(~region)
ggsave(paste("region_age",t,".pdf",sep=""), width=12, height=10)}