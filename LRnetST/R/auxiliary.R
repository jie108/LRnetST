######## auxiliary functions
cal_order<- function(adj_matrix){
### What: get topological order of nodes in directed acyclic graph
### inputs: 
## adj_matrix: a pxp adjacency matrix: adj_matrix[i,j] indicates the presence (1) or absence (0) of edge i--->j
### outputs:
## node_order: p-vector: the topological order of nodes

  p <- nrow(adj_matrix)
  if(nrow(adj_matrix) != ncol(adj_matrix))  
    stop("adjacency matrix is not square matrix!")

  node_order <- NULL
  ind_p <- numeric(p)

  csums <- colSums(adj_matrix)
  npar_index <- which(csums == 0)
  hpar_index <- which(csums != 0)
  ind_p[npar_index] <- 1

  node_order <- c(node_order, npar_index)
  current <- hpar_index
  while(length(current) > 0){
    j <- current[1]
    par_index <- which(adj_matrix[,j] == 1)
    if(all(ind_p[par_index] == 1)){
      node_order <- c(node_order,j)
      current <- current[-1]
      ind_p[j] <- 1
    }
    else{
      current <- current[-1]
      current <- c(current,j)
    }
  }

  return(node_order)
}



moral_graph<-function(adj.matrix){
### What: given a DAG adjacency matrix, get its moral graph
### inputs:
## adj.matrix: a pxp adjacency matrix : adj.matrix[i,j] indicates the presence (1) or absence (0) of edge i->j
### outputs:
## moral.adj.matrix: a pxp symmetric matrix: adjacency matrix of the moral graph

  p=nrow(adj.matrix)
  moral.adj.matrix=adj.matrix

  n.parents=colSums(adj.matrix)
  for(i in 1:p){
    if(n.parents[i]>=2){
      pars=which(adj.matrix[,i]>0)

      for(j in 1:(n.parents[i]-1)){
        for(k in (j+1):n.parents[i]){
          moral.adj.matrix[pars[j],pars[k]]=1
        }#end for
      }#end for
    }#end if
  }#end for

  moral.adj.matrix=((moral.adj.matrix+t(moral.adj.matrix))>0)
  return(moral.adj.matrix)
}


vstructures<-function(adj.matrix){
### What: find vstructures in a DAG adjacency matrix 
### inputs:
## adj.matrix: p by p DAG adjacency matrix
### outputs: 
## 3-column matrix: (par1, child , par2) where par1 < par2
  
  p=ncol(adj.matrix)
  res=NULL
  for(i in 1:p){
    parent.node=which(adj.matrix[,i]!=0)
    n.par=length(parent.node)
    
    if(n.par>1){
      for(j in 1:(n.par-1)){
        for(k in (j+1):n.par){
          if(adj.matrix[parent.node[j], parent.node[k]]==0&&adj.matrix[parent.node[k], parent.node[j]]==0){
           res=rbind(res,c(parent.node[j],i,parent.node[k]))
          }##end if
        }
      }##for loop
    }##end if
  }## end i loop
  
  if(!is.null(res)){
    colnames(res)=c("par1","child","par2")
  }
  return(res)
}

skeleton<-function(adj.matrix){
### What: given a DAG adjacency matrix, return the skeleton graph
###inputs:
## adj.matrix: p by p DAG adjacency matrix
### outputs: 
## skeleton.matrix: p by p skeleton matrix

skeleton.matrix=(adj.matrix+t(adj.matrix))>0
return(skeleton.matrix)	

}


compare.vstructures<-function(target.vstructures,true.vstructures){
### What: compare two v-structures 
### inputs: two 3-column matrices of v-structures from the "vstructures" function
### outputs: a 3-column matrix, the subset of true-vstructure in the target vstructure
  
  corr.v=NULL
  if(!is.null(target.vstructures)){
    target.vstructures=matrix(target.vstructures,ncol=3)
    target.l=nrow(target.vstructures)
    for(i in 1:target.l){
      res=apply(true.vstructures,1,function(x) all(x==target.vstructures[i,]))
      if(any(res==TRUE)){
        corr.v=rbind(corr.v,true.vstructures[which(res==TRUE),])
      }##end if 
    }##end for
  }##end if 
  
  return(corr.v)
}
