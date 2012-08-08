rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")
library("ggplot2")


load("cortical.thickness.rdata")
clinical.data <- read.table("20120807_CoriticalThickness.csv", sep=";", header=T)
clinical.data <- clinical.data[,c("GSKID","gender","ageatexam","diseasecourse","diseaseduration")]
clinical.data[,"GSKID"] <- paste("EPIC",clinical.data[,"GSKID"], sep="")
data <- merge(cortical.thickness,clinical.data,by.x="FID",by.y="GSKID", all.x=T)
data <- cbind(data,median_lt=apply(data[,3:36], 1, median), median_rt=apply(data[,37:70], 1, median))
data2 <- data[!is.na(data$gender),]

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
sink()

cat("Statistical analysis: Men vs. Women in each region. \n\n", file="region_gender_age.txt", append=T)

sink("region_gender_age.txt", append=T)
cat("Mean in each region by gender.\n")
print(region_gender_mean <- apply(data[,3:70],2,function(x)tapply(x,data$gender, mean)))

cat("\n Wilcoxon test in each region by gender \n\n")
print(region_gender_wilcoxon <- apply(data[,3:70],2,function(x)wilcox.test(x~gender,data)))
cat("\n t-test in each region by gender \n\n")
print(region_gender_ttest <- apply(data[,3:70],2,function(x)t.test(x~gender,data)))

cat("\n Corralation between thickness in each region and age at examination. \n")
region_age_lm <- apply(data[,3:70],2,function(x)lm(x~ageatexam,data))
for(i in 1:68){cat(colnames(data)[i+2],"\n");print(summary(region_age_lm[[i]]))}
sink()

plotdata <- data.frame(median=c(data$median_lt, data$median_rt),
                       side=c(rep("left",nrow(data)),rep("right",nrow(data))),
                       age_at_exam=data$ageatexam,
                       gender=data$gender)


median_age_plot <- ggplot(plotdata)+
               geom_point(aes(x=age_at_exam,y=median,colour=side))+
               stat_abline(intercept=ans_lt[[1]][1],slope=ans_lt[[1]][2], colour="red")+
               stat_abline(intercept=ans_rt[[1]][1],slope=ans_rt[[1]][2], colour="green")+
               ylab("Median of cortical thickness")+
               xlab("Age at Exam.")
ggsave("median_age.pdf",width=10, height=10)


plotdata2 <- plotdata[!is.na(plotdata$gender),] 
median_gender <- ggplot(plotdata2, aes(gender,median))+
                geom_boxplot(aes(fill=gender))+geom_jitter(alpha=0.2)+
                facet_wrap(~side)
ggsave("median_gender.pdf", width=7, height=14)

for(i in 3:36){if(i==3){region.data_lt <- data.frame(thickness=data[,i],data[,71:72],region=colnames(data)[i])} else {
region.data_lt <- rbind(region.data_lt,data.frame(thickness=data[,i],data[,71:72],region=colnames(data)[i]))}}

for(i in 37:70){if(i==37){region.data_rt <- data.frame(thickness=data[,i],data[,71:72],region=colnames(data)[i])} else {
  region.data_rt <- rbind(region.data_rt,data.frame(thickness=data[,i],data[,71:72],region=colnames(data)[i]))}}

region.data_lt2 <- region.data_lt[!is.na(region.data_lt$gender),]
region.data_rt2 <- region.data_rt[!is.na(region.data_rt$gender),]

region_gender_lt <- ggplot(region.data_lt2, aes(gender,thickness))+
                 geom_boxplot(aes(fill=gender))+
                 facet_wrap(~region)
ggsave("region_gender_lt.pdf", width=12, height=10)


region_gender_rt <- ggplot(region.data_rt2, aes(gender,thickness))+
  geom_boxplot(aes(fill=gender))+
  facet_wrap(~region)
ggsave("region_gender_rt.pdf", width=12, height=10)

region_age_lt <- ggplot(region.data_lt,aes(ageatexam,thickness))+
                 geom_point()+
                 stat_smooth(method="lm",se=FALSE)+
                 xlab("Age at exam.")+
                 facet_wrap(~region)
ggsave("region_age_lt.pdf", width=12, height=10)

region_age_rt <- ggplot(region.data_rt,aes(ageatexam,thickness))+
  geom_point()+
  stat_smooth(method="lm",se=FALSE)+
  xlab("Age at exam.")+
  facet_wrap(~region)
ggsave("region_age_rt.pdf", width=12, height=10)


