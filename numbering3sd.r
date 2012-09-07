setwd("~/Analysis/Cortex_thickness/")
load("cortical.thickness(not3sd).rdata")
cortical.thickness <- cortical.thickness[rownames(cortical.thickness)!="EPIC416" & rownames(cortical.thickness)!="EPIC293",]
test <- apply(cortical.thickness[,1:68], 2, function(x){
  m <- mean(x,na.rm=T)
  sd <- sd(x,na.rm=T)
  x[x > m+3*sd|x < m-3*sd]})

test

EPlist <- names(test[[1]])
for(i in 2:68){EPlist <- c(EPlist,names(test[[i]]))}
table(EPlist)