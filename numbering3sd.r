test <- apply(cortical.thickness[,3:70], 2, function(x){
  m <- mean(x,na.rm=T)
  sd <- sd(x,na.rm=T)
  x[x > m+3*sd|x < m-3*sd]})
EPlist <- names(test[[1]])
for(i in 2:68){EPlist <- c(EPlist,names(test[[i]]))}
table(EPlist)