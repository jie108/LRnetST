#########################
install.packages("remotes")
remotes::install_github("jie108/LRnetST", subdir="LRnetST")
library(LRnetST)
library(Seurat)

#tumor.mat: pxn1 ST data for tumor spots
#immune.mat<- pxn2 ST data for immune spots
#tumor.loc: position (n1x2) for tumor spot
#immune.loc: position (n2x2) for immune spots
#num.D: number of downsampled data
#th: threshold for aggregating bootstrap DAGs
#missing.rate.cut: missing rate cut-off for gene filtering
#path.ligrec: path for the ligand receptor database (data frame of 2 columns: ligand and receptor)
#path: path where the downsampled data will be saved

#NIM: function to construct newighbor integrated matrix
#final.DAG: function to construct  

NIM<- function(tumor.mat, immune.mat, tumor.loc, immune.loc, missing.rate.cut, path.ligrec, path, num.D){
  
  mat2b<- tumor.mat; mat2g<- immune.mat; Posb<- tumor.loc; Posg<- immune.loc

  mind.list<-list()
  dist.list<- list()
  for(i in 1:dim(Posb)[1]){
    x<- Posb[i,1]
    y<- Posb[i,2]
    
    d<- NULL
    for(j in 1:dim(Posg)[1]){
      xj<- Posg[j,1]
      yj<- Posg[j,2] 
      
      d<- c(d,sqrt((x-xj)^2 + (y-yj)^2))
    }
    mind<- which.min(d)
    mind.list[[i]]<- mind
    dist.list[[i]]<- min(d)
  }
  md<- unlist(dist.list)
  mnd<- unlist(mind.list)
  
  int<- read.csv(path.ligrec, header=T)
  uni<- union(int$Ligand,int$Receptor)
  
  adj<- matrix(NA, length(uni), length(uni))
  rownames(adj)<- uni
  colnames(adj)<- uni
  for(i in 1:length(uni)){
    for(j in 1:length(uni)){
      x<- uni[i]
      y<- uni[j]

      ix<- which(int$Ligand %in% x)
      iy<- which(int$Receptor %in% y)
      if(length(intersect(ix,iy))>0){
        adj[i,j]<- T
      }else{
        adj[i,j]=F
      }
    }
  }
  diag(adj)<-F

  adj.s<-adj
  for(i in 1:dim(adj)[1]){
    for(j in 1:dim(adj)[2]){
      if(adj[i,j]==T){
        adj.s[i,j]<-1
        adj.s[j,i]<-1
      }
    }
  }

  mat2<- cbind(mat2b,mat2g)
  med<- median(apply(mat2,2,sum))

for(d in 1:num.D){
  print(d)
  
  ######## downsampling ############  
  ds.mat2<- SampleUMI(mat2, max.umi = med, upsample = FALSE, verbose = FALSE)
  colnames(ds.mat2)<- colnames(mat2)
  
  ds.mat2b<- ds.mat2[,1:(dim(mat2b)[2])]
  ds.mat2g<- ds.mat2[,(dim(mat2b)[2]+1):(dim(mat2)[2])]
  
  ds.mat2b.mat2g<- rbind(ds.mat2b, ds.mat2g[,as.numeric(mnd)]) ## for visium
  # t<- NULL
  # for(i in 1:(dim(ds.mat2b)[2])){
  #   t<- cbind(t,ds.mat2g[,which(colnames(ds.mat2g) == mind.list[[i]])])
  # }
  # ds.mat2b.mat2g<- rbind(ds.mat2b, t) ## for merfish
  # 
  rownames(ds.mat2b.mat2g)<- c(rownames(mat2), paste0(rownames(mat2), "_n"))
  
  ds.mat2b.mat2g.bound<- ds.mat2b.mat2g[,which(md <= 2)]
  
  rowna<-(apply(ds.mat2b.mat2g.bound, 1, function(x) mean(x==0)))
  ir<- which(rowna < missing.rate.cut)
  foo<- ds.mat2b.mat2g.bound[ir,]
  
  mat4<- log((foo+1),2) ## for visium

  maty<- matrix(0,dim(mat4)[1],dim(mat4)[2])
  for(i in 1:(dim(mat4)[1])){
    maty[i,]<- ifelse(mat4[i,]!=0,1,0)
  }
  rownames(maty)<- paste0("Y_",rownames(mat4))
  
  ################ overlap with lig-rec ###############
  
  uni.unin<- union(uni, paste0(uni, "_n"))
  
  mat4.sub<- mat4[which(rownames(mat4) %in% uni.unin),]  
  
  maty.sub<- maty[which(rownames(mat4) %in% uni.unin),]  
  
  ######## create lig-rec adj ##########
  p1<- length(which(grepl("_n", rownames(mat4.sub))==F))
  p2<- length(which(grepl("_n", rownames(mat4.sub))==T))
  
  gn<- rownames(mat4.sub)[(p1+1):(p1+p2)]; gn<- gsub("_.*$","",gn)
  
  adj.gtogn1<- adj.s[which(rownames(adj.s) %in% (rownames(mat4.sub)[1:p1])),which(colnames(adj.s) %in% gn)]
  
  adj.gtogn<- adj.gtogn1[order(match(rownames(adj.gtogn1), rownames(mat4.sub)[1:p1])), order(match(colnames(adj.gtogn1), gn))]
  
  adj.gntog1<- adj.s[which(rownames(adj.s) %in% gn),which(colnames(adj.s) %in% rownames(mat4.sub)[1:p1])] 
  
  adj.gntog<- adj.gntog1[order(match(rownames(adj.gntog1), gn)), order(match(colnames(adj.gntog1), rownames(mat4.sub)[1:p1]))]
  
  ########### create prior ############
  p<- 2*dim(mat4.sub)[1]
  
  node.type<- c(rep("c",p/2),rep("b",p/2))
  
  w1<- w2<- w3<- w4<- matrix(F,p/2,p/2)
  
  diag(w3)<- T
  
  whiteList<- rbind(cbind(w1,w2), cbind(w3,w4))
  
  blackList=matrix(T, p,p)
  
  ###########
  b1<-blackList[1:(p1), ((p1)+1):(p/2)] ## g -> gn
  b1[adj.gtogn==T]<- F
  blackList[1:(p1), ((p1)+1):(p/2)]<-b1
  
  b11<-blackList[((p1)+1):(p/2),1:(p1)] ## gn -> g
  b11[adj.gntog==T]<- F
  blackList[((p1)+1):(p/2),1:(p1)]<-b11
  
  ###########
  b2<-blackList[1:(p1), ((p/2)+p1+1):(p)] ## g -> ygn
  b2[adj.gtogn==T]<- F
  blackList[1:(p1), ((p/2)+p1+1):(p)]<- b2
  
  b22<-blackList[((p/2)+p1+1):(p),1:(p1)] ## ygn -> g
  b22[adj.gntog==T]<- F
  blackList[((p/2)+p1+1):(p),1:(p1)]<- b22
  
  ###########
  b3<-blackList[(p/2+1):(p/2 + p1), (p1+1):(p/2)] ## yg -> gn
  b3[adj.gtogn==T]<- F
  blackList[(p/2+1):(p/2 + p1), (p1+1):(p/2)]<-b3
  
  b33<-blackList[(p1+1):(p/2),(p/2+1):(p/2 + p1)] ## gn -> yg
  b33[adj.gntog==T]<- F
  blackList[(p1+1):(p/2),(p/2+1):(p/2 + p1)]<-b33
  
  ############
  b4<-blackList[(p/2+1):(p/2 + p1), (p/2 + p1 + 1):(p)] ## yg -> ygn
  b4[adj.gtogn==T]<- F
  blackList[(p/2+1):(p/2 + p1), (p/2 + p1 + 1):(p)]<- b4
  
  b44<-blackList[(p/2+p1+1):(p),(p/2+1):(p/2+p1)] ## ygn -> yg
  b44[adj.gntog==T]<- F
  blackList[(p/2+p1+1):(p),(p/2+1):(p/2+p1)]<- b44
  
  ########### create temp ##############
  temp<- list()
  temp$Y<- t(as.matrix(rbind(mat4.sub,maty.sub)))
  temp$nodeType=node.type
  temp$whiteList=whiteList
  temp$blackList=blackList 
  rownames(temp$whiteList)<-colnames(temp$whiteList)<-colnames(temp$Y)
  rownames(temp$blackList)<-colnames(temp$blackList)<-colnames(temp$Y)
  
  save(temp, file=paste0(path,"/temp_ds",d,".RData"))
}
}

final.DAG<- function(path, num.D, th){
  resSC.boot.array<- list()
  for(i in 1:num.D){
    print(i)
    load(paste0(path,"/temp_ds",i,".RData"))
    
    resSC.boot=hcSC_boot(Y=temp$Y, n.boot=1, nodeType=temp$nodeType, whiteList=temp$whiteList, blackList=temp$blackList, scale=TRUE, tol = 1e-6, maxStep = 3000, restart=10, seed = 1, nodeShuffle = T, verbose = F)
    resSC.boot.array[[i]]<- resSC.boot[,,1]
    
    colnames(resSC.boot.array[[i]])<- colnames(temp$Y)
    rownames(resSC.boot.array[[i]])<- colnames(resSC.boot.array[[i]])
  }
  
  uni.all<- NULL
  for(i in 1:num.D){
    uni.all<- union(uni.all, rownames(resSC.boot.array[[i]]))
  }
  
  resSC.boot.array.ext<- array(F,c(length(uni.all),length(uni.all),num.D))
  for(i in 1:num.D){
    print(i)
    for(k1 in 1:length(uni.all)){
      for(k2 in 1:length(uni.all)){
        x<- uni.all[k1]
        y<- uni.all[k2]
        
        ix<- which(rownames(resSC.boot.array[[i]]) %in% x)
        iy<- which(colnames(resSC.boot.array[[i]]) %in% y)
        
        if(length(ix)>0 & length(iy)>0){
          if(resSC.boot.array[[i]][ix,iy]==T){
            resSC.boot.array.ext[k1,k2,i]<- T
          }
        }
      }
    }
  }
  
  resSC.boot.array.ext.ord<- resSC.boot.array.ext
  for(i in 1:num.D){
    
    inter<- resSC.boot.array.ext[,,i]
    rownames(inter)<- uni.all
    colnames(inter)<- uni.all
    
    ind.y<- which(grepl("Y_",rownames(inter))==T)
    ind.noy<- which(grepl("Y_",rownames(inter))==F)
    
    inter.ord<- inter[c(ind.noy,ind.y),c(ind.noy,ind.y)]
    
    resSC.boot.array.ext.ord[,,i]<- inter.ord
  }
  
  ####### create new whitelist #########
  new.p<- length(uni.all)
  
  w1<- w2<- w3<- w4<- matrix(F,new.p/2,new.p/2)
  
  diag(w3)<- T
  
  new.white<- rbind(cbind(w1,w2), cbind(w3,w4))
  
  SCboot.agg2=score_shd(resSC.boot.array.ext.ord, alpha = 1, threshold=th, max.step = 3000, whitelist = new.white)
  final.net<- SCboot.agg2
  colnames(final.net)<- uni.all[c(ind.noy,ind.y)]
  rownames(final.net)<- colnames(final.net)
  sum(final.net) - (dim(final.net)[1])/2
  
  return(final.net) 
}

###### example to create NIM and run DAG #####
out = NIM(tumor.mat=mat2b, immune.mat=mat2g, tumor.loc=Posb, immune.loc=Posg, missing.rate.cut=0.5, path.ligrec="ligrec.csv",path = "test", num.D=2)

lrnet<- final.DAG(path = "test", num.D=2, th=0.4)
sum(lrnet)

adj.sub<- lrnet
adj.sub[new.white==1]<- 0
sum(adj.sub)

####### For plot #######
rowind<- NULL
colind<- NULL
for(i in 1:dim(adj.sub)[1]){
  rowind<- c(rowind, which(adj.sub[i,]==1))
  colind<- c(colind, which(adj.sub[,i]==1))
}

adj.sub2<- adj.sub[union(rowind,colind),union(rowind,colind)]

indy<- which(grepl("Y_",rownames(adj.sub2))==T)

indn<- which(grepl("_n",rownames(adj.sub2))==T)

obj<- graph_from_adjacency_matrix(adj.sub2, mode = "directed", weighted = TRUE)

V(obj)$color<- rep(NA, dim(adj.sub2)[1])
V(obj)$color[indn]<- "#99CCFF"
V(obj)$color[-indn]<- "#FF9999"

pdf("test/lrnetst.pdf",width = 20, height = 20, useDingbats=T)
plot.igraph(obj,vertex.size=5, edge.length=10, edge.color="black", vertex.label.cex=1.5,vertex.color=V(obj)$color,layout=layout.fruchterman.reingold(obj, niter=10000))
dev.off()
