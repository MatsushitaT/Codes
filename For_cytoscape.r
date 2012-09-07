rm(list=ls())
setwd("~/Analysis/Cortex_thickness/Final_dataset/")

Rhosif <- read.csv("Correlation_region_rho.csv")

Rhosif <- Rhosif[,c(-1,-3)]
Rhosif$connection <- "DirectedEdge"

##Make sif
## Select Correlation with p < 0.05/561
con.rho <- Rhosif[Rhosif$Cases.p.value<0.05/561,]
con.rho <- con.rho[,c(1,7,2,4)]
write.table(con.rho, file="Cytoscape_data/rho.sif", row.names=F, col.names=F, quote=F)

##Make na file
nodeattri.na <- as.data.frame(table(unlist(con.rho[,c(1,3)])))
nodeattri.na$V3 <- paste(nodeattri.na[,1],nodeattri.na[,2],sep="=")
write.table(nodeattri.na[,3], file="Cytoscape_data/nodeattri.na", row.names=F, col.names=F, quote=F)