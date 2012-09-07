rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")
library("ggplot2")


load("Final_dataset/cortical.thickness.rdata")
clinical.data <- read.table("20120807_CoriticalThickness.csv", sep=";", header=T)
clinical.data <- clinical.data[,c("GSKID","gender","ageatexam","diseasecourse","diseaseduration")]
clinical.data[,"GSKID"] <- paste("EPIC",clinical.data[,"GSKID"], sep="")
cortical.thickness <- as.data.frame(cortical.thickness)
cortical.thickness$GSKID <- rownames(cortical.thickness)
data <- merge(cortical.thickness,clinical.data, all.x=T, all.y=F)
data <- cbind(data,lh_median=apply(data[,2:35], 1, median), rh_median=apply(data[,36:69], 1, median))
data <- cbind(data, disease_onset_age=data[,71]-data[,73])

lh <- c(2:35,74); rh <- c(2:35+34,75)

##Statistical analysis: Men vs. Women in each region.
##Wilcoxon test in each region by gender
region_gender_wilcoxon <- apply(data[,c(lh,rh)],2,function(x)wilcox.test(x~gender,data))

##t-test in each region by gender \n\n")
region_gender_ttest <- apply(data[,c(lh,rh)],2,function(x)t.test(x~gender,data))

##Corralation between thickness in each region and age at examination. \n")
##Linear regression
region_age_lm <- apply(data[,c(lh,rh)],2,function(x)lm(x~ageatexam,data))

##Spearman's correlation coefficient
region_age_cor <- apply(data[,c(lh,rh)],2,function(x)cor.test(x, data$ageatexam, use="pair", method="spearman"))

##Corralation between thickness in each region and onset age.
##Linear regression
region_onset_lm <- apply(data[,c(lh,rh)],2,function(x)lm(x~disease_onset_age,data))

##Thickness versus age of onset, corrected for age
region_onset_adjust_age_lm <- apply(data[,c(lh,rh)],2,function(x)lm(x~disease_onset_age+ageatexam,data))

table <- cbind(
  t(sapply(region_gender_ttest[c(1:35)], function(tt)
  c(lh.t.test.p=tt$p.value, lh.diff=unname(diff(tt$estimate))))),
               lh.wilcox.p=sapply(region_gender_wilcoxon[c(1:35)],"[[","p.value"),
               t(sapply(region_gender_ttest[c(1:35+35)], function(tt)
                 c(rh.t.test.p=tt$p.value, rh.diff=unname(diff(tt$estimate))))),
               rh.wilcox.p=sapply(region_gender_wilcoxon[c(1:35+35)],"[[","p.value"),
               t(sapply(region_age_lm[c(1:35)], function(lm)
                 unlist(list(lh.age=summary(lm)$coefficients["ageatexam",c(1,4)])))),
               t(sapply(region_age_lm[c(1:35+35)], function(lm)
                 unlist(list(rh.age=summary(lm)$coefficients["ageatexam",c(1,4)])))),
               t(sapply(region_age_cor[c(1:35)],function(cor)
                 c(lh.age.cor.rho=unname(cor$estimate), lh.age.cor.p=unname(cor$p.value)))),
               t(sapply(region_age_cor[c(1:35+35)],function(cor)
                 c(rh.age.cor.rho=unname(cor$estimate),rh.age.cor.p=unname(cor$p.value)))),
               t(sapply(region_onset_lm[c(1:35)], function(lm)
                 unlist(list(lh.onset=summary(lm)$coefficients["disease_onset_age",c(1,4)])))),
               t(sapply(region_onset_lm[c(1:35+35)], function(lm)
                 unlist(list(rh.onset=summary(lm)$coefficients["disease_onset_age",c(1,4)])))),
               t(sapply(region_onset_adjust_age_lm[c(1:35)], function(lm)
                 unlist(list(lh.onset_adjusted_by_age=summary(lm)$coefficients["disease_onset_age",c(1,4)])))),
               t(sapply(region_onset_adjust_age_lm[c(1:35+35)], function(lm)
                 unlist(list(rh.onset_adjusted_by_age=summary(lm)$coefficients["disease_onset_age",c(1,4)]))))
)

rownames(table)<-gsub("lh_","",rownames(table))
write.csv(table, "Final_dataset/age_gender_onset_region_table_case.csv")


plotdata <- data.frame(median=c(data$median_lt, data$median_rt),
                       r_a_cingulate=c(data$lh_rostralanteriorcingulate, data$rh_rostralanteriorcingulate),
                       m_orb_f=c(data$lh_medialorbitofrontal, data$rh_medialorbitofrontal),
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

r_a_cingulate_age_plot <- ggplot(plotdata2, aes(x=age_at_exam,y=r_a_cingulate,colour=gender_side))+
  geom_point()+
  stat_smooth(method="lm",se=FALSE)+
  ylab("Rostal anterior cingulate")+
  xlab("Age at Exam.")
ggsave("r_a_cingurate_age.pdf",width=10, height=10)

m_orb_f_age_plot <- ggplot(plotdata2, aes(x=age_at_exam,y=m_orb_f,colour=gender_side))+
  geom_point()+
  stat_smooth(method="lm",se=FALSE)+
  ylab("Medial orbitofrontal")+
  xlab("Age at Exam.")
ggsave("m_orb_frontal_age.pdf",width=10, height=10)



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