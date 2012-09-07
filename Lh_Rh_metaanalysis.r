rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")
library("metafor")

load("cortical.thickness.rdata")
clinical.data <- read.table("20120807_CoriticalThickness.csv", sep=";", header=T)
clinical.data <- clinical.data[,c("GSKID","gender","ageatexam","diseasecourse","diseaseduration")]
clinical.data[,"GSKID"] <- paste("EPIC",clinical.data[,"GSKID"], sep="")
data <- merge(cortical.thickness,clinical.data,by.x="FID",by.y="GSKID", all.x=T)
data <- cbind(data,median_lt=apply(data[,3:36], 1, median), median_rt=apply(data[,37:70], 1, median))
data <- cbind(data, median=apply(data[3:70], 1, median))
data <- cbind(data, disease_onset_age=data[,72]-data[,74])


lh_F_M <- summary(lm(median_lt~ageatexam*gender,data))
rh_F_M <- summary(lm(median_rt~ageatexam*gender,data))

rma(yi = c(lh_F_M$coefficients["ageatexam:genderM",1],rh_F_M$coefficients["ageatexam:genderM",1]),
    sei = c(lh_F_M$coefficients["ageatexam:genderM",2],rh_F_M$coefficients["ageatexam:genderM",2]),
    method="FE")


