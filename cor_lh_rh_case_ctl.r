rm(list=ls())

#set workdesk
setwd("~/Analysis/Cortex_thickness/Final_dataset/")

##Load thickness data. Including individual with age < 1 in controls
## All data is truncated to mean+- 3sd range
load("cortical.thickness.rdata")
load("control_thickness.rdata")
cortical.thickness <- as.data.frame(cortical.thickness)
thick.ctl <- as.data.frame(thick.ctl)

##add median and phenotype collumn
cortical.thickness <- cbind(cortical.thickness,
                            lh_median=apply(cortical.thickness[,1:34],1,median),
                            rh_median=apply(cortical.thickness[,1:34+34],1,median),
                            phenotype="CASE")

thick.ctl <- cbind(thick.ctl,
                  lh_median=apply(thick.ctl[,1:34],1,median),
                  rh_median=apply(thick.ctl[,1:34+34],1,median),
                  phenotype="CONTROL")


##confirm collumn name
sum(colnames(cortical.thickness)==colnames(thick.ctl))


##combine table
total.table <- rbind(cortical.thickness,thick.ctl)
total.table <- total.table[,c(71,1:34,69,1:34+34,70)]

##compare correlations between same region in both hemispheres between case and control
for(i in 2:36){
  tmp <- summary(lm(total.table[,i]~total.table[,i+35]*phenotype, total.table))
  if(i == 2){
    cor_lh_rh_case_ctl <- data.frame(case=t(tmp$coefficients[2,c(1,4)]), 
                                     control_case=t(tmp$coefficients[4,c(1,4)]))} else {
                                       cor_lh_rh_case_ctl <- rbind(cor_lh_rh_case_ctl,
                                                                   data.frame(case=t(tmp$coefficients[2,c(1,4)]), 
                                                                              control_case=t(tmp$coefficients[4,c(1,4)])))}
                                     }
row.names(cor_lh_rh_case_ctl) <- gsub("lh_", "", colnames(total.table[,2:36]))

write.csv(cor_lh_rh_case_ctl,file="cor_lh_rh_case_ctl.csv")

for(i in 2:36){
  hemiplot <- total.table[,c(1,i,i+35)]
  region <- gsub("lh_","",colnames(total.table)[i])
  colnames(hemiplot) <- c("phenotype", "Lh","Rh")
  plot <- ggplot(hemiplot,aes(Rh,Lh,colour=phenotype))+
    geom_point()+
    opts(title = region)+
    stat_smooth(method="lm",se=F)
  ggsave(paste("cor_lh_rh_case_ctl",i,".pdf",sep=""), width=10, height=10)}