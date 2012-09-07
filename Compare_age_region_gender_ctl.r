rm(list=ls())
setwd("~/Analysis/Cortex_thickness/")

thick.ctl <- read.csv("ctrls_4.5_thk.csv")
thick.ctl$fsid. <- gsub("'","",thick.ctl$fsid.)
thick.ctl$TYPE. <- gsub("'","",thick.ctl$TYPE.)
thick.ctl$Gender. <- gsub("'","",thick.ctl$Gender.)

#thick.ctl$"lh_mean(unknown)" <- apply(thick.ctl[,5:39],1,mean)
#thick.ctl$"rh_mean(unknown)" <- apply(thick.ctl[,40:74],1,mean)
#thick.ctl$lh_mean <- apply(thick.ctl[,6:39],1,mean)
#thick.ctl$rh_mean <- apply(thick.ctl[,41:74],1,mean)
#thick.ctl$"lh_median(unknown)" <- apply(thick.ctl[,5:39],1,median)
#thick.ctl$"rh_median(unknown)" <- apply(thick.ctl[,40:74],1,median)
#thick.ctl$lh_median <- apply(thick.ctl[,6:39],1,median)
#thick.ctl$rh_median <- apply(thick.ctl[,41:74],1,median)

ctl.clinical <- thick.ctl[,1:4]
row.names(thick.ctl) <- thick.ctl[,1]
colnames(ctl.clinical) <- gsub("\\.","",colnames(ctl.clinical))
thick.ctl <- thick.ctl[,c(-1:-5, -40,-75:-76)]
colnames(thick.ctl)[1:34] <- gsub("\\.","",colnames(thick.ctl)[1:34])
colnames(thick.ctl)[1:34] <- paste("lh_",colnames(thick.ctl)[1:34],sep="")
colnames(thick.ctl)[35:68] <- gsub("\\.\\.1","",colnames(thick.ctl)[35:68])
colnames(thick.ctl)[35:68] <- paste("rh_",colnames(thick.ctl)[35:68],sep="")


## Truncated to 3sd range
thick.ctl <- apply(thick.ctl,2,function(x){
  m <- mean(x,na.rm=T)
  sd <- sd(x,na.rm=T)
  pmin(pmax(x,round(m-3*sd,2)),round(m+3*sd,2))})

thick.ctl <- cbind(ctl.clinical, thick.ctl)

thick.ctl$lh_median <- apply(thick.ctl[,5:38], 1, median)
thick.ctl$rh_median <- apply(thick.ctl[,5:38+34],1,median)

lh <- c(5:38,73); rh <- c(5:38+34,74);

##Statistical analysis: Men vs. Women in each region.
region_gender_wilcoxon <- apply(thick.ctl[,c(lh,rh)],2,function(x)wilcox.test(x~Gender,thick.ctl))

##t-test in each region by gender
region_gender_ttest <- apply(thick.ctl[,c(lh,rh)],2,function(x)t.test(x~Gender,thick.ctl))

##Corralation between thickness in each region and age at examination.

##Exclude individual with extreme age data
thick.ctl2 <- thick.ctl[-which.min(thick.ctl$Age),]
## Linear
region_age_lm <- apply(thick.ctl2[,c(lh,rh)],2,function(x)lm(x~Age,thick.ctl2))

##Spearman's correlation coefficient
region_age_cor <- apply(thick.ctl2[,c(lh,rh)],2,function(x)cor.test(x, thick.ctl2$Age, use="pair", method="spearman"))

table <- cbind(
  t(sapply(region_gender_ttest[c(1:35)], function(tt)
    c(lh.t.test.p=tt$p.value, lh.diff=unname(diff(tt$estimate))))),
  lh.wilcox.p=sapply(region_gender_wilcoxon[c(1:35)],"[[","p.value"),
  t(sapply(region_gender_ttest[c(1:35+35)], function(tt)
    c(rh.t.test.p=tt$p.value, rh.diff=unname(diff(tt$estimate))))),
  rh.wilcox.p=sapply(region_gender_wilcoxon[c(1:35+35)],"[[","p.value"),
  t(sapply(region_age_lm[c(1:35)], function(lm)
    unlist(list(lh.age=summary(lm)$coefficients["Age",c(1,4)])))),
  t(sapply(region_age_lm[c(1:35+35)], function(lm)
    unlist(list(rh.age=summary(lm)$coefficients["Age",c(1,4)])))),
  t(sapply(region_age_cor[c(1:35)],function(cor)
    c(lh.age.cor.rho=unname(cor$estimate),lh.age.cor.p=unname(cor$p.value)))),
  t(sapply(region_age_cor[c(1:35+35)],function(cor)
    c(rh.age.cor.rho=unname(cor$estimate), rh.age.cor.p=unname(cor$p.value))))
)
rownames(table)<-gsub("lh_","",rownames(table))
write.csv(table, "Final_dataset/age_gender_region_table_control.csv")

