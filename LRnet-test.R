library(LRnetST)
data(example)
Y.n=example$Y # data matrix
p<- dim(Y.n)[2] # no. of nodes
true.dir=example$true.dir  #adjacency matrix of the data generating DAG
temp<- LRnetST::hcSC(Y=Y.n,nodeType=rep("c",p), whiteList=NULL, blackList=NULL, tol = 1e-6, scale=TRUE, maxStep = 1000, restart=10, seed = 1,  verbose = FALSE)
adj.temp=temp$adjacency

#(ii) DAG learning by hill climbing: for bootstrap resamples

library(foreach)
library(doParallel)


boot.adj<- LRnetST::hcSC_boot_parallel(Y=Y.n, n.boot=10, nodeType=rep("c",p), whiteList=NULL, blackList=NULL, scale=TRUE, tol = 1e-6, maxStep = 1000, restart=10, seed = 1,  nodeShuffle=TRUE, numThread = 2,verbose = FALSE)


#(iii) Bootstrap aggregation of DAGs learnt from bootstrap resamples

adj.bag=LRnetST::score_shd(boot.adj, alpha = 1, threshold=0)

sum(adj.bag==1&true.dir==0)/sum(adj.bag==1) ## FDR
sum(adj.bag==1&true.dir==1)/sum(true.dir==1) ## Power