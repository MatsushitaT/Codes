set.seed(1980)

######################################################
##### creating randomly generated data ###############
######################################################
## no need to understand this chunk, just run it
tmp <- rnorm(40) +  rep(1:5,each=8) * 5
left.data <- cbind(
  left.A = tmp,
  left.B = tmp + rnorm(40)*4,
  left.C = tmp + 20,
  left.D = tmp + 20 + rnorm(40)*4)
right.data <- left.data + rnorm(length(left.data)) / 10
colnames(right.data) <- gsub("left","right", colnames(right.data))
rm(tmp)

######################################################
##### plot and look at the data ######################
######################################################
## plot # 1
pairs(left.data)
## plot # 2
matplot(left.data,type = 'l',lty=1,lwd=2,col=1:4)
legend('topleft',c("A","B","C","D"),lwd=2,col=1:4)

######################################################
##### Step 1, choose a distance metric ###############
######################################################

##### First distance metric: correlation distance: 
## how to calculate it
dist.cor <- as.dist (1 - cor(left.data))
## look at the distances
round(dist.cor , 2)
## note: small numbers mean "closer", large numbers means "farther"
## note: in this case, A and C are very close!
## that is because they follow the same pattern of peaks and valleys
## this is best seen on plot # 1 above (the pairs() plot)

##### Second distance metric: euclidean distance
## how to calculate it
dist.euc <- dist(t(left.data))
## lok at the distances
round(dist.euc)
## note: small numbers mean "closer", large numbers means "farther"
## note: in this case, A is close to B and C is close to D
## that is because they are close to each other.
## this is best seen on plot # 2 above (the matplot() plot)


##### Third distance metric: manhattan distance
## how to calculate it
dist.man <- dist(t(left.data),method = "manhattan")
## lok at the distances
round(dist.man)
## Note: small numbers mean "closer", large numbers means "farther"
## Note: in this case, A is close to B and C is close to D
## that is because they are close to each other.
## Manhattan is actually quite similar to euclidean for this data.
## I didn't have time to create a dataset highlighting the difference
## between manhattan and euclidean.

######################################################
##### Step 2, choose a clustering algorithm ##########
######################################################
##### as an example, let's use the "ward" algorithm with the "euclidean" distance
hc <- hclust(dist.euc, method = 'ward')
plot(hc) ## note that A and B are close, C and D are close
##### as an example, let's use the "ward" algorithm with the "correlation" distance
hc <- hclust(dist.cor, method = 'ward')
plot(hc) ## note that now A and C are close, B and D are further away.
##### here is a loop for plotting all 7 algorithms
##### using only the "correlation" distance
methods <- c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
pdf("my plots.pdf")
for(m in methods)
  plot(hclust(dist.cor,method = m),main = m)
dev.off()
rm(m)
##### Note: open the file "my plots.pdf" to see the results.
##### Note: in this case, you will see that the seven algorithms seem almost identical.
##### That is because this dataset is very small. In your larger dataset, you will notice that there
##### are big differences between the different algorithms.

######################################################
##### how to cluster left and right together #########
######################################################
all.data <- cbind(left.data,right.data)
plot(hclust(dist(1 - cor(all.data))))