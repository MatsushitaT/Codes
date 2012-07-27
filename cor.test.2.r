## all pairwise correlations done by cor.test()
## written by tm and pk
cor.test.2 <- function(x,meth){
  output <- list()
  for(i in 1:(length(x)-1))
    for(j in (i+1):length(x)){
      output[[names(x)[i]]][[names(x)[j]]] <- cor.test(x[[i]] , x[[j]],method=meth)
    }
  return(unlist(output,recursive=F))
}

## example
#cts <- cor.test.2(list(cortex=1:10, insula=1:10 + rnorm(10),hippocampus=1:10),"pearson" )
#cts
#sapply(cts,"[[","p.value")
