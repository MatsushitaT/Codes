rm(list=ls())
load("~/Analysis/Cortex_thickness/GWAS/08212012Big_table.rdata")
cbind(cumsum(rev(table(cut(apply(data[,-1:-7],1,max),9:53/2)))))

data.small <- data[which(apply(data[,-1:-7],1,max) > 12),]
data.small <- data.small[order(data.small$CHR,data.small$BP),]
data.small <- cbind(data.small[1:7],
                    data.small[7 + which(apply(data.small[,-1:-7],2,max) > 12)])

dim(data.small)

data.small2 <- data
data.small2[data.small2$MAF < 0.15 , grep("REC",names(data.small2))] <- 0
cbind(cumsum(rev(table(cut(apply(data.small2[,-1:-7],1,max),9:33/2)))))
data.small2 <- data.small2[which(apply(data.small2[,-1:-7],1,max) > 7),]
data.small2 <- data.small2[order(data.small2$CHR,data.small2$BP),]

data.small2 <- cbind(data.small2[1:7],
                     data.small2[7 + which(apply(data.small2[,-1:-7],2,max) > 7)])
print(data.small2[,1:7])

data.small2 <- data
data.small2[data.small2$MAF < 0.15 , grep("REC",names(data.small2))] <- 0
data.small2[data.small2$MAF < 0.005 , -1:-7] <- 0
cbind(cumsum(rev(table(cut(apply(data.small2[,-1:-7],1,max),9:33/2)))))
data.small2 <- data.small2[which(apply(data.small2[,-1:-7],1,max) > 7),]
dim(data.small2)
data.small2 <- data.small2[order(data.small2$CHR,data.small2$BP),]
data.small2 <- cbind(data.small2[1:7],
                     data.small2[7 + which(apply(data.small2[,-1:-7],2,max) > 7)])

write.csv(unique(data.small2$gene_symbol),file="genelist2.csv")

data.small3 <- data.small2[,grep("sex_and_age", names(data.small2))]

data.small2 <- cbind(data.small2[1:7],data.small3)

length(grep("sex_and_age", names(data.small2)))

plot.data.small2 <- -as.matrix(data.small2[-1:-7])
rownames(plot.data.small2) <- data.small2$gene_symbol
heatmap(t(plot.data.small2))

