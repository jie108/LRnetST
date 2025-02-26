## LRnetST: Learning Directed Acyclic Graphs for ligands and receptors based on Spatial Transcriptomics data

<img src="Fig1A_new.png" width="700" align="center">

<img src="Fig1B_new.png" width="700" align="center">

- [Reference](#Reference)
- [Overview](#Overview)
- [Installation](#Installation)
- [Usage](#Usage)
- [Arguments](#Arguments)
- [Value](#Value)
- [Examples](#Examples)
- [Contributions](#contributions)

## Reference 
--
https://www.biorxiv.org/content/10.1101/2021.08.03.454931v3

## Overview
```
LRnetST: 
contains the R package "LRnetST" for learning directed acycic graphs based on spatial transcriptomics data.

```


## Installation


### Install LRnetST
```
library(devtools)
install_github("jie108/LRnetST", subdir="LRnetST")
```

## Usage

```
hcSC: A function to learn a DAG model for the given ST data with no bootstrap resamples by the hill climbing algorithm for mixture of continuous and binary variables.

LRnetST::hcSC(Y,nodeType, whiteList, blackList, tol, scale, maxStep, restart, seed,  verbose)

hcSC_boot_parallel: A function to learn a DAG model for every bootstrap resmples of the given ST data by the hill climbing algorithm for mixture of continuous and binary variables.

LRnetST::hcSC_boot_parallel(Y, node.type, n.boot, whiteList, blackList,  scale, tol, maxStep, restart, seed, nodeShuffle, bootDensityThre, numThread, verbose)

score_shd: A function to use structural hamming distance to aggregate DAGs. It aggregates an ensemble of DAGs to obtain a DAG that minimizes the overall distance to the ensemble.

LRnetST::score_shd(boot.adj, alpha, threshold, max.step, blacklist, whitelist, verbose)
```


## Arguments

### Arguments for LRnetST::hcSC and LRnetSt::hcSC_boot_parallel
  
| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| Y	       |	           | an n by p data matrix: n – sample size, p – number of variables
| n.boot (only for hc_boot_parallel) |      1       | an integer: the number of bootstrap resamples of the data matrix Y
| node.type  		       |   NULL      | a vector of length equal to the number of variables specifying the type of variable/node type: "c" for continuous and "b" for binary
| whitelist          | NULL   |  a p by p 0-1 matrix: if the (i,j)th-entry is "1", then the edge i–>j will always be included in the DAG during the search
| blacklist	         | NULL    | a p by p 0-1 matrix: if the (i,j)th-entry is "1", then the edge i–>j will be excluded from the DAG during the search
| scale |  TRUE | logical: whether to scale the continuous nodes such l2_norm^2/n=1 (won't change zero pattern)
| tol     |     1e-06     | a scalar: a number to indicate a threshold below which values will be treated as zero
| maxStep		           | 2000    |an integer: the maximum number of search steps of the hill climbing algorithm
| restart | 10 | an integer: number of times to restart the search algorithm after a local optimal is achieved. The purpose is to search for global optimal
|seed| 1 | an integer: seed used for bootstrap restart and bootstrap resampling
| nodeShuffle (ony for hc_boot_parallel) | FALSE | logical: whether to shuffle the order of the variables before DAG learning
|bootDensityThre (only for hcSC_boot_parallel)| 0.1| numeric between (0,1): lower cutoff of columnwise nonzero-entry rate in bootstrap resamples
| numThread (only for hcSC_boot_parallel) | 2 |  an integer for running parallel computation of DAG learning from bootstrap resamples
| verbose		     | FALSE   | logical: whether print the step information



### Arguments for LRnetSt::score_shd
  
| Parameter                 | Default       | Description   |	
| :------------------------ |:-------------:| :-------------|
| boot.adj	       |	           | A p by p by B array, where B is the number of DAGs to be aggregated. It records the adjacency matrices. It may be the output of the "hcSC_boot_parallel" function.
| alpha         | 1          |a positive scalar: alpha defines which member of the gSHD family should be used to aggregate the DAGs. In general, the larger the alpha, the more aggressive of the aggregation, in that less edges are retained leading to smaller FDR and less power
| threshold 	       |	0	     |a scalar: it defines the frequency cut-off value(=(1-threshold)/2), "0" corresponds to cut-off 0.5
| max.step		           | NULL             |This is a legacy parameter and it does not have any effect 
| blacklist	         | NULL             | a p by p 0-1 matrix: if the (i,j)th-entry is "1", then the edge i–>j will be excluded from the DAG during the search
| whitelist          | NULL           |  a p by p 0-1 matrix: if the (i,j)th-entry is "1", then the edge i–>j will always be included in the DAG during the search
| verbose		     |     FALSE     | logical: whether print the step information



## Value

### Value for LRnetST::hcSC

a list of three components

| Object       | Description   |
| :------------------------ | :-------------|
| adjacency	  | adjacency matrix of the learned DAG
| score       | BIC score at each search step
| operations  | a matrix recording the selected operation, addition, deletion or reversal of an edge, at each search step
| deltaMin    | Minimum value of the score change at every step

### Value for LRnetST::hcSC_boot_parallel
an array
| Object       | Description   |
| :------------------------ | :-------------|
| adjacency	  | adjacency matrix of the learned DAG


### Value for LRnetST::score_shd
a matrix 
| Object       | Description   |
| :------------------------ | :-------------|
| adj.matrix	  | adjacency matrix of the learned DAG

  
## Examples
```
(i) DAG learning by hill climbing for mixture of continuous and binary nodes: no bootstrap resample

library(LRnetST)
data(example)
Y.n=example$Y # data matrix
p<- dim(Y.n)[2] # no. of nodes
true.dir=example$true.dir  #adjacency matrix of the data generating DAG
temp<- LRnetST::hcSC(Y=Y.n,nodeType=rep("c",p), whiteList=NULL, blackList=NULL, tol = 1e-6, scale=TRUE, maxStep = 1000, restart=10, seed = 1,  verbose = FALSE)
adj.temp=temp$adjacency

(ii) DAG learning by hill climbing: for bootstrap resamples

library(foreach)
library(doParallel)
boot.adj<- LRnetST::hcSC_boot_parallel(Y=Y.n, n.boot=10, nodeType=rep("c",p), whiteList=NULL, blackList=NULL, scale=TRUE, tol = 1e-6, maxStep = 1000, restart=10, seed = 1,  nodeShuffle=TRUE, numThread = 2,verbose = FALSE)

(iii) Bootstrap aggregation of DAGs learnt from bootstrap resamples

adj.bag=LRnetST::score_shd(boot.adj, alpha = 1, threshold=0)
sum(adj.bag==1&true.dir==0)/sum(adj.bag==1) ## FDR
sum(adj.bag==1&true.dir==1)/sum(true.dir==1) ## Power
```


## Contributions

If you find small bugs, larger issues, or have suggestions, please email the maintainer at <jiepeng108@gmail.com>. Contributions (via pull requests or otherwise) are welcome.
