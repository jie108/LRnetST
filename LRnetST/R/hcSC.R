hcSC <- function(Y, nodeType=NULL, whiteList=NULL, blackList=NULL, scale=TRUE, tol = 1e-6, maxStep = 2000L, restart=10L, seed = 1L, verbose = FALSE) {
##Y: n by p data matrix;  nodeType: n-string: either "c" or "b";
##whiteList: p by p logical matrix of whitelisted edges; blackList: p by p logical matrix of blacklisted edges;
##scale: whether to scale the continous nodes such l2_norm^2/n=1 (won't change zero pattern)
##tol: double, lower cutoff of score improvement for HC algorithm; maxStep:integer, maximum number of search steps of HC algorithm; restart: integer, # of times to perform HC search (to break score ties)
##seed: integer, random number generator seed for breaking score tie;
##return: HC learned adjacency matrix on Y

	p=ncol(Y)
	if(is.null(whiteList)){
		whiteList=matrix(FALSE, p,p)


	}else{##check whiteList format 
      if(!is.matrix(whiteList)||!is.logical(whiteList)||length(dim(whiteList))!=2||dim(whiteList)[1]!=p || dim(whiteList)[2]!=p){
      	stop(paste("whiteList must be a ", p,"by", p, " logical matrix!"))
      }   

	}

	if(is.null(blackList)){
		blackList<-matrix(FALSE, p,p)
        diag(blackList)<-TRUE
	}else{##check black list format 
		if(!is.matrix(blackList)||!is.logical(blackList)||length(dim(blackList))!=2||dim(blackList)[1]!=p || dim(blackList)[2]!=p){
      	stop(paste("blackList must be a ", p,"by", p, " logical matrix!"))
      }   
	}	

	if(is.null(nodeType)){
		nodeType=rep("c",p)
	}else{##check noteType format 
       if(!is.vector(nodeType)||!is.character(nodeType)||!all(nodeType=="c"|nodeType=="b")||length(nodeType)!=p){
       		stop(paste("nodeType must be a length ", p, " character vector with elements either c or b!"))
       }

	}
  

  if(scale){
    for(i in 1:p){
     if(nodeType[i]=="c"){##standardize continuous nodes such that l2_norm^2/n=1
      Y[,i]=Y[,i]/sqrt(sum(Y[,i]^2)/n)
     }
    }
  }


    hcSC_(Y, nodeType, whiteList, blackList, tol, maxStep, restart, seed, verbose)
}


hcSC_boot<-function(Y, n.boot=1, nodeType=NULL,  whiteList=NULL, blackList=NULL, scale=TRUE, tol = 1e-6, maxStep = 2000L, restart=10L, seed = 1L,  nodeShuffle=FALSE, bootDensityThre=0.1, verbose = FALSE){
##Y: n by p data matrix;  n.boot: integer, number of boostrap resamples; nodeType: n-string: either "c" or "b";
##whiteList: p by p logical matrix of whitelisted edges; blackList: p by p logical matrix of blacklisted edges;
##scale: whether to scale the continous nodes such l2_norm^2 /n=1 (won't change zero pattern)
##tol: double, lower cutoff of score improvement for HC algorithm; maxStep:integer, maximum number of search steps of HC algorithm; restart: integer, # of times to perform HC search (to break score ties)
##seed: integer, random number generator seed for boostrap resampling, nodeShuffling, breaking score tie;
##nodeShuffle: whether (T) or not (F) to perform node shuffle (to avoid bias due to search order and score ties);
##bootDensityThre: numeric between (0,1), lower cutoff of columnwise nonzero-entry rate in bootstrap resamples;  
##return: HC learned adjacency matrices on bootstrap resamples 

   p=ncol(Y)
   n=nrow(Y)
   density.min=min(apply(Y!=0, 2, mean)) ##lowest columnwise nonzero-entry rate in the data

   ##check argument type 
   if(is.null(whiteList)){
		whiteList=matrix(FALSE, p,p)
	}else{##check whiteList format 
      if(!is.matrix(whiteList)||!is.logical(whiteList)||length(dim(whiteList))!=2||dim(whiteList)[1]!=p || dim(whiteList)[2]!=p){
      	stop(paste("whiteList must be a ", p,"by", p, " logical matrix!"))
      }   

	}

	if(is.null(blackList)){
		blackList<-matrix(FALSE, p,p)
        diag(blackList)<-TRUE
	}else{##check black list format 
		if(!is.matrix(blackList)||!is.logical(blackList)||length(dim(blackList))!=2||dim(blackList)[1]!=p || dim(blackList)[2]!=p){
      	stop(paste("blackList must be a ", p,"by", p, " logical matrix!"))
      }   
	}	

	if(is.null(nodeType)){
		nodeType=rep("c",p)
	}else{##check noteType format 
       if(!is.vector(nodeType)||!is.character(nodeType)||!all(nodeType=="c"|nodeType=="b")||length(nodeType)!=p){
       		stop(paste("nodeType must be a length ", p, " character vector with elements either c or b!"))
       }

	}

  if(bootDensityThre<=0||bootDensityThre>=density.min){
    stop(paste("bootDensityThre must be in between: ", 0,  " and ", density.min, "(lowest columnwise nonzero-entry rate in data)"))       
  }

  if(scale){
    for(i in 1:p){
     if(nodeType[i]=="c"){##standardize continuous nodes such that l2_norm^2/n=1
      Y[,i]=Y[,i]/sqrt(sum(Y[,i]^2)/n)
     }
    }
  }

   ## hc on bootstrap resamples  
   adj=array(NA,c(p,p,n.boot))

   for(i in 1:n.boot){
   print(paste("fit bootstrap sample ", i, "..."))
   set.seed(i*1001+seed)
   Y.pick = boot_dense(Y, bootDensityThre) ##bootstrap resample with each column having at least bootDensityThre percentage of nozero entries

   if(nodeShuffle){##shuffle node order 
   	node.rand=sample(1:p, p, replace=FALSE)  
   	node.index=sort(node.rand,decreasing=F,ind=T)$ix 
   }else{##not shuffle node order
   	node.rand=1:p
   	node.index=1:p
   }

   Y.B=Y.pick[, node.rand]
   node.B=nodeType[node.rand]
   whiteList.B=whiteList[node.rand, node.rand]
   blackList.B=blackList[node.rand, node.rand]

   curRes=hcSC_(Y.B, node.B, whiteList.B, blackList.B, tol, maxStep, restart, i*11+seed, verbose)
   adjRes=curRes$adjacency
   adj[,,i]=adjRes[node.index, node.index]
   }

   return(adj)
}

hcSC_boot_parallel<-function(Y, n.boot=1, nodeType=NULL,  whiteList=NULL, blackList=NULL, scale=TRUE, tol = 1e-6, maxStep = 2000L, restart=10L, seed = 1L,  nodeShuffle=FALSE, bootDensityThre=0.1, numThread=2, verbose = FALSE){
##Y: n by p data matrix;  n.boot: integer, number of boostrap resamples; nodeType: n-string: either "c" or "b";
##whiteList: p by p logical matrix of whitelisted edges; blackList: p by p logical matrix of blacklisted edges;
##scale: whether to scale the continous nodes such l2_norm^2 /n=1 (won't change zero pattern)
##tol: double, lower cutoff of score improvement for HC algorithm; maxStep:integer, maximum number of search steps of HC algorithm; restart: integer, # of times to perform HC search (to break score ties)
##seed: integer, random number generator seed for boostrap resampling, nodeShuffling, breaking score tie;
##nodeShuffle: whether (T) or not (F) to perform node shuffle (to avoid bias due to search order and score ties);
##bootDensityThre: numeric between (0,1), lower cutoff of columnwise nonzero-entry rate in bootstrap resamples;  
##numThread: number of threads to use in parallel computing
##return: HC learned adjacency matrices on bootstrap resamples 

   p=ncol(Y)
   n=nrow(Y)
   density.min=min(apply(Y!=0, 2, mean)) ##lowest columnwise nonzero-entry rate in the data

   ##check argument type 
   if(is.null(whiteList)){
    whiteList=matrix(FALSE, p,p)
  }else{##check whiteList format 
      if(!is.matrix(whiteList)||!is.logical(whiteList)||length(dim(whiteList))!=2||dim(whiteList)[1]!=p || dim(whiteList)[2]!=p){
        stop(paste("whiteList must be a ", p,"by", p, " logical matrix!"))
      }   

  }

  if(is.null(blackList)){
    blackList<-matrix(FALSE, p,p)
        diag(blackList)<-TRUE
  }else{##check black list format 
    if(!is.matrix(blackList)||!is.logical(blackList)||length(dim(blackList))!=2||dim(blackList)[1]!=p || dim(blackList)[2]!=p){
        stop(paste("blackList must be a ", p,"by", p, " logical matrix!"))
      }   
  } 

  if(is.null(nodeType)){
    nodeType=rep("c",p)
  }else{##check noteType format 
       if(!is.vector(nodeType)||!is.character(nodeType)||!all(nodeType=="c"|nodeType=="b")||length(nodeType)!=p){
          stop(paste("nodeType must be a length ", p, " character vector with elements either c or b!"))
       }

  }

  if(bootDensityThre<=0||bootDensityThre>=density.min){
    stop(paste("bootDensityThre must be in between: ", 0,  " and ", density.min, "(lowest columnwise nonzero-entry rate in data)"))       
  }

  if(scale){
    for(i in 1:p){
     if(nodeType[i]=="c"){##standardize continuous nodes such that l2_norm^2/n=1
      Y[,i]=Y[,i]/sqrt(sum(Y[,i]^2)/n)
     }
    }
  }

   ## hc on bootstrap resamples  
   cl <- makeCluster(numThread)
   registerDoParallel(cl)  
   result<-foreach(i=1:n.boot, .combine="list",.multicombine = TRUE, .maxcombine = n.boot, .packages=c("dagbagSC")) %dopar%{
  
     set.seed(i*1001+seed)
     Y.pick = boot_dense(Y, bootDensityThre) ##bootstrap resample with each column having at least bootDensityThre percentage of nozero entries

     if(nodeShuffle){##shuffle node order 
      node.rand=sample(1:p, p, replace=FALSE)  
      node.index=sort(node.rand,decreasing=F,ind=T)$ix 
     }else{##not shuffle node order
      node.rand=1:p
      node.index=1:p
     }

     Y.B=Y.pick[, node.rand]
     node.B=nodeType[node.rand]
     whiteList.B=whiteList[node.rand, node.rand]
     blackList.B=blackList[node.rand, node.rand]

     curRes=hcSC_(Y.B, node.B, whiteList.B, blackList.B, tol, maxStep, restart, i*11+seed, verbose)
     adjRes=curRes$adjacency
     adjRes=adjRes[node.index, node.index]
     adjRes
   }
  stopCluster(cl)
   
   ##return results as a p by p by n.boot array  
   return(array(as.numeric(unlist(result)), dim=c(p, p, n.boot)))
}




###
boot_dense<-function(Y, bootDensityThre){
##input: data matrix Y; bootDensityThre: numeric between [0,1], lower cutoff of columnwise nonzero-entry rate in bootstrap resamples;  
##return: boostrap resample of Y:  assure that for every column at least bootDensityThre percent entries are nonzero  
 n=nrow(Y)
 sparseCol=TRUE
 while(sparseCol){
   s.pick=sample(1:n, n, replace=TRUE) ##resample data
   Y.pick=Y[s.pick,]
   sparseCol=any(apply(Y.pick!=0, 2, mean)<bootDensityThre)  
 }
 return(Y.pick)
}

###
SC_prepare<-function(logCount){
##input: log(count+1) transformed single cell counts 
##output: inputs (Y, nodeType, whiteList, blackList) for hcSC and hcSC_boot functions
n=nrow(logCount)
p=ncol(logCount)
S=matrix(as.numeric(logCount!=0), n,p)  ##n by p gene on/off indicator matrix 
Y=cbind(logCount, S) ## data matrix for hcSC and hcSC_boot: (logCount, logCount>0): n by 2*p
node.type=c(rep("c",p), rep("b", p))

whiteList=matrix(FALSE, 2*p,2*p)
blackList=matrix(FALSE, 2*p, 2*p)
diag(blackList)=TRUE

for (i in 1:p){
whiteList[p+i, i]=TRUE ## logcount_i>0 -> logcount_i whitelisted
blackList[i, p+i]=TRUE ## logcount_i -> logcount_i>0 blacklisted
}

return(list("Y"=Y, "nodeType"=node.type, "whiteList"=whiteList, "blackList"=blackList))

}