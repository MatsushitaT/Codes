rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")
library("ggplot2")

##Make dataset of cases
load("Final_dataset/cortical.thickness.rdata")
clinical.data <- read.table("20120807_CoriticalThickness.csv", sep=";", header=T)
clinical.data <- clinical.data[,c("GSKID","gender","ageatexam","diseasecourse","diseaseduration")]
clinical.data[,"GSKID"] <- paste("EPIC",clinical.data[,"GSKID"], sep="")
cortical.thickness <- as.data.frame(cortical.thickness)
cortical.thickness$GSKID <- rownames(cortical.thickness)
data <- merge(cortical.thickness,clinical.data, all.x=T, all.y=F)
data <- cbind(data,lh_median=apply(data[,2:35], 1, median), rh_median=apply(data[,36:69], 1, median), phenotype="CASE")
data <- data[,-72:-73]
data <- data[,c(1,74,70,71,2:69,72:73)]
colnames(data)[1] <- "ID"


##Make dataset of controls
thick.ctl <- read.csv("ctrls_4.5_thk.csv")
thick.ctl$fsid. <- gsub("'","",thick.ctl$fsid.)
thick.ctl$TYPE. <- gsub("'","",thick.ctl$TYPE.)
thick.ctl$Gender. <- gsub("'","",thick.ctl$Gender.)

ctl.clinical <- thick.ctl[,1:4]
row.names(thick.ctl) <- thick.ctl[,1]
colnames(ctl.clinical) <- gsub("\\.","",colnames(ctl.clinical))
thick.ctl <- thick.ctl[,c(-1:-5, -40,-75:-76)]
colnames(thick.ctl)[1:34] <- gsub("\\.","",colnames(thick.ctl)[1:34])
colnames(thick.ctl)[1:34] <- paste("lh_",colnames(thick.ctl)[1:34],sep="")
colnames(thick.ctl)[35:68] <- gsub("\\.\\.1","",colnames(thick.ctl)[35:68])
colnames(thick.ctl)[35:68] <- paste("rh_",colnames(thick.ctl)[35:68],sep="")


### Truncated to 3sd range
thick.ctl <- apply(thick.ctl,2,function(x){
  m <- mean(x,na.rm=T)
  sd <- sd(x,na.rm=T)
  pmin(pmax(x,round(m-3*sd,2)),round(m+3*sd,2))})

thick.ctl <- cbind(ctl.clinical, thick.ctl)

thick.ctl$lh_median <- apply(thick.ctl[,5:38], 1, median)
thick.ctl$rh_median <- apply(thick.ctl[,5:38+34],1,median)
colnames(thick.ctl)[1:4] <- colnames(data)[1:4]
thick.ctl[thick.ctl$gender=="Male","gender"] <- "M"
thick.ctl[thick.ctl$gender=="Female", "gender"] <- "F"

##Round Age
thick.ctl$ageatexam <- floor(thick.ctl$ageatexam)

##Confirm the order of regions
sum(colnames(thick.ctl) == colnames(data))

##Combind case and control
total.table <- rbind(data, thick.ctl)

rm(clinical.data, cortical.thickness, ctl.clinical, data,thick.ctl)

##Exclude individual age < 1
total.table <- total.table[-which.min(total.table$ageatexam),]


lh <- c(5:38,73); rh <- c(5:38+34, 74)

## Case vs. Control in each region (t-test)
thickness_phenotype_ttest <- apply(total.table[,c(lh,rh)],2,function(x)t.test(x~phenotype,total.table))

## Case vs. Control in each region (Wilcoxon)
thickness_phenotype_wilcoxon <- apply(total.table[,c(lh,rh)],2,function(x)wilcox.test(x~phenotype,total.table))

##Correlation between  
thickness_age_adjust_by_pheno <- apply(total.table[,c(lh,rh)],2,function(x)lm(x~ageatexam+phenotype,total.table))
thickness_age_adjust_by_pheno_gender <- apply(total.table[,c(lh,rh)],2,function(x)lm(x~ageatexam+phenotype+gender,total.table))

##x~ageatexam*gender+phenotype==x~ageatexam+gender+phenotype+age*gender
thickness_age_agegender <- apply(total.table[,c(lh,rh)],2,function(x)lm(x~ageatexam*gender,total.table))
thickness_age_agepheno <- apply(total.table[,c(lh,rh)],2,function(x)lm(x~ageatexam*phenotype,total.table))
thickness_age_adjust_by_pheno_gender_agegender <- apply(total.table[,c(lh,rh)],2,function(x)lm(x~ageatexam*gender+phenotype,total.table))
thickness_age_adjust_by_pheno_gender_agepheno <- apply(total.table[,c(lh,rh)],2,function(x)lm(x~ageatexam*phenotype+gender,total.table))

table.casecontrol <- cbind(
  t(sapply(thickness_phenotype_ttest[c(1:35)], function(tt)
    c("lh.diff(Control-Case)"=unname(diff(tt$estimate)),"lh.t.test.p"=tt$p.value))),
  lh.wilcox.p=sapply(thickness_phenotype_wilcoxon[c(1:35)],"[[","p.value"),
  t(sapply(thickness_phenotype_ttest[c(1:35+35)], function(tt)
    c("rh.diff(Control-Case)"=unname(diff(tt$estimate)),rh.t.test.p=tt$p.value))),
  rh.wilcox.p=sapply(thickness_phenotype_wilcoxon[c(1:35+35)],"[[","p.value"),
  t(sapply(thickness_age_adjust_by_pheno[1:35], function(lm)
    unlist(list("lh.age~thickness_adjusted_by_phenotype"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno[c(1:35+35)], function(lm)
    unlist(list("rh.age~thickness_adjusted_by_phenotype"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender[c(1:35)], function(lm)
    unlist(list("lh.age~thickness_adjusted_by_phenotype_gender"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender[c(1:35+35)], function(lm)
    unlist(list("rh.age~thickness_adjusted_by_phenotype_gender"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender[c(1:35)], function(lm)
    unlist(list("lh.control_case_adjusted_by_age_gender"=summary(lm)$coefficients["phenotypeCTRL",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender[c(1:35+35)], function(lm)
    unlist(list("rh.control_case_adjusted_by_age_gender"=summary(lm)$coefficients["phenotypeCTRL",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender[c(1:35)], function(lm)
    unlist(list("lh.M_F_adjusted_by_age_phenotype"=summary(lm)$coefficients["genderM",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender[c(1:35+35)], function(lm)
    unlist(list("rh.M_F_adjusted_by_age_phenotype"=summary(lm)$coefficients["genderM",c(1,4)])))),
  t(sapply(thickness_age_agegender[c(1:35)], function(lm)
    unlist(list("lh.slope_female"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_agegender[c(1:35)], function(lm)
    unlist(list("lh.slope_difference(M-F)"=summary(lm)$coefficients["ageatexam:genderM",c(1,4)])))),
  t(sapply(thickness_age_agegender[c(1:35+35)], function(lm)
    unlist(list("rh.slope_female(M-F)"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_agegender[c(1:35+35)], function(lm)
    unlist(list("rh.slope_difference(M-F)"=summary(lm)$coefficients["ageatexam:genderM",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender_agegender[c(1:35)], function(lm)
    unlist(list("lh.slope_female_adjusted_by_phenotype(M-F)"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender_agegender[c(1:35)], function(lm)
    unlist(list("lh.slope_difference_adjusted_by_phenotype(M-F)"=summary(lm)$coefficients["ageatexam:genderM",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender_agegender[c(1:35+35)], function(lm)
    unlist(list("rh.slope_female_adjusted_by_phenotype(M-F)"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender_agegender[c(1:35+35)], function(lm)
    unlist(list("rh.slope_difference_adjusted_by_phenotype(M-F)"=summary(lm)$coefficients["ageatexam:genderM",c(1,4)])))),
  t(sapply(thickness_age_agepheno[c(1:35)], function(lm)
    unlist(list("lh.slope_case"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_agepheno[c(1:35)], function(lm)
    unlist(list("lh.slope_difference(Ctl-Case)"=summary(lm)$coefficients["ageatexam:phenotypeCTRL",c(1,4)])))),
  t(sapply(thickness_age_agepheno[c(1:35+35)], function(lm)
    unlist(list("rh.slope_case"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_agepheno[c(1:35+35)], function(lm)
    unlist(list("rh.slope_difference(Ctl-Case)"=summary(lm)$coefficients["ageatexam:phenotypeCTRL",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender_agepheno[c(1:35)], function(lm)
    unlist(list("lh.slope_case_adjusted_by_gender"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender_agepheno[c(1:35)], function(lm)
    unlist(list("lh.slope_difference_adjusted_by_gender(Ctl-Case)"=summary(lm)$coefficients["ageatexam:phenotypeCTRL",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender_agepheno[c(1:35+35)], function(lm)
    unlist(list("rh.slope_case_adjusted_by_gender"=summary(lm)$coefficients["ageatexam",c(1,4)])))),
  t(sapply(thickness_age_adjust_by_pheno_gender_agepheno[c(1:35+35)], function(lm)
    unlist(list("rh.slope_difference_adjusted_by_gender(Ctl-Case)"=summary(lm)$coefficients["ageatexam:phenotypeCTRL",c(1,4)]))))
)

rownames(table.casecontrol)<-gsub("lh_","",rownames(table.casecontrol))
write.csv(table.casecontrol, "Final_dataset/age_gender_region_table(case_control).csv")

total.table <- total.table[,c(1:4,c(rbind(1:34,1:34+34))+4,73:74)]

plotdata <- NULL
for(i in 5:ncol(total.table)){
  plotdata<-rbind(plotdata, data.frame(Region = rep(colnames(total.table)[i],nrow(total.table)),
                                       Thickness = total.table[,i], total.table[,c(1:4)]))}

plotdata$phenotype_gender <- paste(plotdata$phenotype, plotdata$gender, sep="+")

#specify a number of division (n <-)
n <- 4
for(t in 1:((70 %/% n)+1)){
  regionplot <- plotdata[plotdata$Region %in% colnames(total.table[,c(5:74)])[(n*(t-1)+1):(t*n)],]
  region_age <- ggplot(regionplot,aes(ageatexam,Thickness, colour=phenotype_gender))+
    geom_point()+
    stat_smooth(method="lm",se=F)+
    xlab("Age at exam.")+
    facet_wrap(~Region)
  ggsave(paste("Final_dataset/Age_thickness_gender_case_ctl/region_age_case_ctl",t,".pdf",sep=""), width=15, height=15)}