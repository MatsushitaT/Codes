plot.linear.model<-function(xx,yy,...){
plot(xx,yy,...)
model<-lm(yy~xx)
lines(xx[!is.na(xx) & !is.na(yy)],model$fit)
title(main=paste("Slope (beta) =",model[[1]][2]))
title(main=paste("\n\nAnova p-value =",anova(model)[[5]][1]))
}








#xx<-c(NA,1:5,6:2)
#yy<-c(1:10,NA)
#plot.linear.model(xx,yy)