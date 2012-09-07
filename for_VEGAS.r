rm(list=ls())
#set workdesk
setwd("~/Analysis/Cortex_thickness/GWAS/VEGAS_result/")

##### Read all VEGAS outputs
VEGASES <- list()

for(file in list.files()){
  if(regexpr(".out$",file) <0){next}
  VEGASES[[file]] <- read.table(file, header=T, as.is=T)
}
rm(file)

save(VEGASES, file="VEGAScomb.rdata")

##### Sort VEGAS output by chromosome and position
for(models in rev(names(VEGASES))){
  cat("Now sorting by chromosome:    ",experiment,date(),'\n')
  VEGAS <- VEGASES[[models]]
  #par(mfrow = 1:2)
  #plot(order(VEGAS$Chr,VEGAS$Start,VEGAS$Stop), col = VEGAS$Chr)
  VEGAS <- VEGAS[order(VEGAS$Chr,VEGAS$Start,VEGAS$Stop),]
  #plot(order(VEGAS$Chr,VEGAS$Start,VEGAS$Stop), col = VEGAS$Chr)
  VEGASES[[models]] <- VEGAS
}
rm(VEGAS,models)


##### add blocks to VEGAS output
## set the pvalue threshold
p.BlockDefine.threshold <- 0.05
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
rm(VEGAS,models)

##### select genes with the least p value in blocks including some genes
GENELIST <- list()
for(models in names(VEGASES)){
  VEGAS <- VEGASES[[models]]
  VEGAS <- VEGAS[!is.na(VEGAS$block),]
  write.csv(VEGAS, file=paste("Block_list/",models, "(0.05).csv",sep=""),row.names=F, col.names=T, quote=F)
  dup <- unique(VEGAS[duplicated(VEGAS$block),"block"])
  VEGAS1 <- VEGAS[!(VEGAS$block %in% dup),]
  VEGAS2 <- VEGAS[VEGAS$block %in% dup,]
  VEGAS2 <- do.call("rbind",
                    tapply(rownames(VEGAS2),factor(VEGAS2$block),function(x){
                      temp <- VEGAS2[x,]
                      temp[temp$Pvalue==min(temp[,"Pvalue"]),]
                    }))
  VEGAS <- rbind(VEGAS1,VEGAS2)
  VEGAS <- VEGAS[order(VEGAS$block),]
  GENELIST[[models]] <- VEGAS
}

rm(VEGAS,VEGAS1, VEGAS2, models, dup)


###write gene list for DAVID
for(models in names(GENELIST)){
  VEGAS <- GENELIST[[models]]
  write.table(VEGAS$Gene, file=paste("Gene_list(for_DAVID)/", models, ".txt",sep=""), row.names=F, col.names=F, quote=F)  
}

