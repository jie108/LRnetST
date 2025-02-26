## new implementation by Jie Peng on Feb. 2025 (to replace the score_shd function in dagbag package)
## A family of (generalized) SHD (structure hamming distance) to aggregate the DAGs learned on bootstrap resamples 
## parameters: boot.adj: ensemble of DAGs, p by p by nb array; alpha: parameter for structural hamming distances -- how many operations a reversal is counted
## threshold: used to determine frequency cut off = (1-threshold)/2; 
## max.step: maximum number of steps for hill climbing algorithm -- this is a legacy parameter and has no effect now
## blacklist, whitelist: p by p boolean matrices
## return: a p by p aggregated DAG: edge i->j in the DAG iff (i,j) element =1;
score_shd <- function(boot.adj, alpha = 1, threshold = 0, max.step = NULL, blacklist = NULL, whitelist = NULL, verbose = FALSE){

 ##
 p = dim(boot.adj)[1]     ## number of variables
 nb = dim(boot.adj)[3]    ## number of DAGs in the ensemble 

 ## whitelist and blacklist 
  if(is.null(whitelist)){
         whitelist=matrix(0,p,p)
  }else{## check the whitelist is acyclic 
   library(igraph)
   g <- graph_from_adjacency_matrix(whitelist, mode = "directed")
   if(!is_dag(g)){
    stop("the whitelist must correspond to a DAG!")
   }
}
 
  if(is.null(blacklist)){
       blacklist=matrix(0,p,p)
     }

 ## start with the whitelist 
 res = whitelist 

 ## get cur parents set 
 par.set<-lapply(1:p, function(i) which(res[,i]==1)-1) ## off set index by 1 for the c++ acyclic check function call
  
 ## calculate selection frequency and generalized selection frequency 

 sele.freq=apply(boot.adj, c(1,2), mean) ## selection frequency for each edge 
 gen.freq=sele.freq+(1-alpha/2)*t(sele.freq) ## generalized selection frequency for each edge

 ## make a matrix of four columns: (i,j, sele.freq[i,j], gen.freq[i,j]) with rows sorted according to descending value of gen.freq
  # Get all (i, j) indices
  indices <- which(!is.na(sele.freq), arr.ind = TRUE)

  # Create a data frame with i, j, and corresponding sele.freq and gen.freq values
  df <- data.frame(i = indices[,1], j = indices[,2], SF = sele.freq[indices], GSF=gen.freq[indices])

  # Sort by value in descending order of GSF
  df_sorted <- df[order(-df$GSF), ]

  # Only retain edges with GSF>freq.cut
  freq.cut=(1-threshold)/2 
  #df_filtered <- df_sorted[(df_sorted$GSF > 0.5)&(df_sorted$SF>freq.cut), ]
  df_filtered <- df_sorted[df_sorted$GSF > freq.cut, ]
  
  if(nrow(df_filtered)>0){

    for (k in 1:nrow(df_filtered)){##add edges sequentially according to GSF, each time first check for blacklist/whitelist, then check for acyclicity of the resulting graph 
      from.node=df_filtered[k,"i"]
      to.node=df_filtered[k,"j"]
      
      if(verbose){
       print(paste(k, "th operation:", from.node, "->", to.node))
      }
      
     # browser() ##for debugging 
      if((!blacklist[from.node, to.node])&&(!whitelist[from.node, to.node])){## if not in the blacklist (not allowed) and not in the whitelist (already there)
        if(!edgeOnLoop(from.node-1, to.node-1, par.set)){ ## call Rcpp function: if acyclic, then add this edge to current graph and update par.set; Note that c++ indexing starts with 0, so need to -1 for from.node and to.node indices
          res[from.node, to.node]=1
          par.set[[to.node]]= c(par.set[[to.node]], from.node-1) ## off set index by 1 for c++ function call
      }else if(verbose){
          print(paste("not pass acyclic check!"))
      }
    }##end of if
  }## end of for-loop 
 }##end of if 

 return(res)
}