rm(list=ls())
setwd("~/Analysis/Cortex_thickness/Final_dataset/")

load("control_thickness(no3sd).rdata")
load("cortical.thickness(no3sd,inclued_outliner).rdata")

cortical.thickness$Phenotype <- "CASE"
thick.ctl$Phenotype <- "CONTROL"
total.table <- rbind(cortical.thickness,thick.ctl)
total.table <- cbind(total.table,
                     lh_median=apply(total.table[,1:34],1,median),
                     rh_median=apply(total.table[,1:34+34],1,median))

tmp <- total.table[, c(rbind(1:34, 1:34+34),70,71,69)]

plotdata <- NULL
for(i in 1:(ncol(tmp)-1)){
  plotdata<-rbind(plotdata, data.frame(Region = rep(names(tmp)[i],nrow(tmp)),
                                       Thickness = tmp[,i],
                                       ID = rownames(tmp),
                                       Phenotype = tmp$Phenotype,
                                       sort = i))}

plotdata$Region_Case <- paste(plotdata$Region,plotdata$Phenotype,sep="_")

ploting <- ggplot(plotdata, aes(reorder(Region_Case,sort), Thickness))
ploting + geom_boxplot(aes(fill =Region_Case)) +
  opts(legend.position = "none") +
  xlab("Region (case-control, lh-rh)")+
  opts(title="Thickness (no 3sd truncated, including outliners)")+
  opts(axis.text.x=theme_text(angle=90, hjust=1))
ggsave("Thickness_no3sd_include_outliner.png", width = 15, height = 8, dpi = 400)

rm(cortical.thickness,plotdata,thick.ctl, tmp, total.table, ploting, i)

load("cortical.thickness.rdata")
load("control_thickness.rdata")
cortical.thickness <- as.data.frame(cortical.thickness)
thick.ctl <- as.data.frame(thick.ctl)
cortical.thickness$Phenotype <- "CASE"
thick.ctl$Phenotype <- "CONTROL"
total.table <- rbind(cortical.thickness,thick.ctl)
total.table <- cbind(total.table,
                     lh_median=apply(total.table[,1:34],1,median),
                     rh_median=apply(total.table[,1:34+34],1,median))

tmp <- total.table[, c(rbind(1:34, 1:34+34),70,71,69)]

plotdata <- NULL
for(i in 1:(ncol(tmp)-1)){
  plotdata<-rbind(plotdata, data.frame(Region = rep(names(tmp)[i],nrow(tmp)),
                                       Thickness = tmp[,i],
                                       ID = rownames(tmp),
                                       Phenotype = tmp$Phenotype,
                                       sort = i))}

plotdata$Region_Case <- paste(plotdata$Region,plotdata$Phenotype,sep="_")

ploting <- ggplot(plotdata, aes(reorder(Region_Case,sort), Thickness))
ploting + geom_boxplot(aes(fill =Region_Case)) +
  opts(legend.position = "none") +
  xlab("Region (case-control, lh-rh)")+
  opts(title="Thickness (3sd truncated, excluding outliners)")+
  opts(axis.text.x=theme_text(angle=90, hjust=1))
ggsave("Thickness_3sd_exclude_outliner.png", width = 15, height = 8, dpi = 400)

rm(cortical.thickness,plotdata,thick.ctl, tmp, total.table, ploting, i)