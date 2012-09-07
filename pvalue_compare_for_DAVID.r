rm(list=ls())
#set workdesk
setwd("~/Analysis/Cortex_thickness/GWAS/VEGAS_result/")

load("VEGAScomb.rdata")

##### add blocks to VEGAS output
## set the pvalue threshold
pvalue <- c(0.05, 0.03, 0.01, 0.005, 0.001)

for(p.BlockDefine.threshold in pvalue){
  for(models in names(VEGASES)){
  VEGAS <- VEGASES[[models]]
  VEGAS$block <- NA
  current.block <- 0
  for(current.row in 1:nrow(VEGAS)){
    ## move on if pval is no good
    if(VEGAS[current.row,"Pvalue"] >= p.BlockDefine.threshold) next
    ## if first row, or new block, then bump to next block
    if(current.row == 1 ||
      VEGAS[current.row - 1,"Pvalue"] >= p.BlockDefine.threshold)
      current.block <- current.block + 1
    ## put the block in the matrix
    VEGAS[current.row,"block"] <- current.block
  } ## 20 seconds
  rm(current.row,current.block)
  VEGASES[[models]] <- VEGAS
  }
  if(p.BlockDefine.threshold==0.05){test <- cbind(
    sapply(VEGASES, function(x)max(x$block, na.rm=T)))} else {
  test <- cbind(test,sapply(VEGASES, function(x)max(x$block, na.rm=T)))}
  }

colnames(test) <- pvalue
rownames(test) <- gsub(".out","",rownames(test))

write.csv(test, file="pvalue_comp.csv")
