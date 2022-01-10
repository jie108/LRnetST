
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"


//List RcppEigen and RcppNumerical headers in  include statements along with the depends() attributes to tell R where to find these header files
#include <RcppEigen.h>
#include <RcppNumerical.h>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::plugins(cpp11)]]          


using namespace Numer;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Upper;


typedef Eigen::Map<Eigen::MatrixXd> MapMat; //Eigen templated class map: to have the contents of an R matrix mapped to the contents of the object from the Eigen class
typedef Eigen::Map<Eigen::VectorXd> MapVec;


//check whether an edge is on a loop  in a graph; note: this function can also check whether adding a new edge to the graph will result in a loop/cycle. (loop version)
bool edgeOnLoop(int fromNode, int toNode, const Rcpp::List& parSet){
  //inputs: fromNode: index of the parent node in the edge under testing; toNode: index of the child node in the edge under testing
  //parSet: list of parent nodes for each node in the graph 
  //output: true: if the edge is on a loop; false: if the edge is not on a loop 

  unsigned int p = parSet.size(); //total number of nodes of the graph
  Rcpp::LogicalVector onPath(p, false); //ith element indicates whether node i is added to path 
  Rcpp::IntegerVector path(p,-1); // record indices of nodes on the path by order of addition 
  
  path[0]=fromNode; // start path at the fromNode, trace back to its ancesters until find toNode (return true) or exhaust the path without finding toNode (return false)
  onPath[fromNode] = true;

  int curPos=0; // current position under examination on path
  int lastPos=0; //last position on path 

  while(true){
   if(curPos>lastPos){//exhaust the path without finding the toNode: no cycle   
    return false; 
    }else{
     int curNode=path[curPos]; 

     if(curNode==toNode){//find cycle:  
      return true;
     }else{
     std::vector<int> curPar=parSet[curNode]; 

     int q=curPar.size();
     if(q>0){//examine parents of curNode
       for(int i=0;i<q; i++){
         int j = curPar[i];  //index of ith parent of curNode 
         if(!onPath[j]){//if ith parent (node j) is not on path yet 
          onPath[j] = true; 
          path[++lastPos]=j; //increase lastPos by 1 and then add node j to the path at lastPos  
         } //end if onPath
       } //end for curPar
      }//end if curPar.size()>0
      curPos++;
     }//end if
    }//end if

  }//end while

  return false; 
}




//get ancester indicators for a given node (loop  version)
Rcpp::LogicalVector ancester(int i, const Rcpp::List& parSet){
//input: i -- index of the node under interest; parSet: list of parents set of the graph; 
//output: a logical vector indicating whether each node is node i's ancester or not
//note: (i) node i itself is included in its ancester set; (i) the graph must be acyclic, otherwise lead to infinite loop

  int p=parSet.size();
  Rcpp::LogicalVector output(p, false); //jth element indicates whether it is node i's ancester 
  Rcpp::IntegerVector path(p, -1);  //search path 
  output[i]=true; //include node i in its ancestor set 
  path[0]=i; //start search at node i, trace back to node i's ancesters 
 
  int curPos=0; //current position on path
  int lastPos=0; //the last(end) position on path 

  while(true){
   if(curPos>lastPos){//exhaust the search path,  return result
     return output;
   }else{
     int curNode=path[curPos]; //current node on path
     std::vector<int> curPar=parSet[curNode]; //current parent set
     int q = curPar.size();
     if(q>0){
      for(int k=0;k<q;k++){//kth parent (node j) of current node
        int j = curPar[k]; 
        if(!output[j]){//if node j is not included in the ancester set yet
         output[j] = true; //include node j in the ancester set
         path[++lastPos]=j; //include node j on path  
        }
      }//end for 
     }//end if q>0
     curPos++; 
   }
  }//end while 

   return output;
}


//get descedant indicators for a given node (loop version)
Rcpp::LogicalVector descend(int i, const Rcpp::List& chdSet){
//input: i -- index of the node under interest; chdSet: list of children set of the graph; note: the graph must be acyclic, otherwise lead to infinite loop
//output: a logical vector indicating whether each node is node i's descendant or not
//note: (i) node i itself is included in its descendant set; (i) the graph must be acyclic, otherwise lead to infinite loop
  
  int p=chdSet.size();
  Rcpp::LogicalVector output(p, false); //jth element indicates whether it is node i's descendant 
  Rcpp::IntegerVector path(p, -1);  //search path 
  output[i]=true; //include node i in its descendant set 
  path[0]=i; //start search at node i, trace down to node i's descendants 
 
  int curPos=0; //current position on path
  int lastPos=0; //the last(end) position on path 

  while(true){
   if(curPos>lastPos){//exhaust the search path,  return result
     return output;
   }else{
     int curNode=path[curPos]; //current node on path
     std::vector<int> curChd=chdSet[curNode]; //current chidren set
     int q = curChd.size();
     if(q>0){
      for(int k=0;k<q;k++){//kth child (node j) of current node
        int j = curChd[k]; 
        if(!output[j]){//if node j is not included in the descendant set yet
         output[j] = true; //include node j in the descandant set
         path[++lastPos]=j; //include node j on path  
        }
      }//end for 
     }//end if q>0
     curPos++; 
   }
  }//end while 
  
   return output;
}


bool acyclicCheck(int fromNode, int toNode, Rcpp::List& parSet, Rcpp::String operType){
  //check where an operation maintains acyclicity 
  //output: true -- acyclic; false -- cycle
  if(operType=="add"){
    return !edgeOnLoop(fromNode, toNode, parSet);
  }else if(operType=="reverse"){
    //remove edge fromNote -> toNode
    std::vector<int> parTo=parSet[toNode]; 
    std::vector<int>::iterator it=std::find(parTo.begin(), parTo.end(), fromNode); //find position of fromNode
    if(it!=parTo.end()){// if edge fromNote -> toNode exists, then remove it 
      parTo.erase(it); //temporarily erase fromNode from toNode parent set 
      parSet[toNode]=parTo; //temporarily update parent set 
      bool output=!edgeOnLoop(toNode, fromNode, parSet); // check adding toNode->fromNode

      parTo.push_back(fromNode); //restore parSet
      parSet[toNode]=parTo;  
      
      return output;
    }else{
      Rcpp::Rcout<<"edge does not exist! can not be reversed!";
      return false;
    }
  }else{
    
      Rcpp::Rcout<<"operation must be either add or reverse";
    return false;
  } 
   
}


void acyclicUpdate(const Rcpp::IntegerVector lastOper, const Rcpp::LogicalVector& fromAn, const Rcpp::LogicalVector& fromDe, const Rcpp::LogicalVector& toAn, const Rcpp::LogicalVector& toDe, 
  const Rcpp::IntegerVector oper, Rcpp::List& parSet, Rcpp::LogicalMatrix& acyStatus){
//input: lastoper: 3 x 1 integer vector of last selected operation; fromAn (toAn), fromDe(toDe): px1 logical vectors of ancester and descendant information of the from(to) node of the last selected operation 
//oper: 3x1 integer vector of operation under consideration; parSet: parent sets of current graph; acyStatus: acyclic status to be updated  
//output: none
 
 int i=oper[0], j=oper[1];

 if(lastOper[2]==1){//last selected operation is to add from->to
  
   if(oper[2]==1&&acyStatus(i,j)&&toDe[i]&&fromAn[j]){//if operation under consideration is to add i->j which is also acyclic in last step, and if i is to Node descendant and j is from node ancester, then leads to circle 
      acyStatus(i,j)=false; //need update
   }

   if(oper[2]==3&&acyStatus(j,i)&&toDe[j]&&fromAn[i]){//if operation under consideration is to reverse i->j which is also acyclic in last step, and if j is to Node descendant and i is from  node ancester, then leads to circle 
       if(i!=lastOper[0]||j!=lastOper[1]){//if operation under consideration is not to reverse from->to
        acyStatus(j,i)=false; //need update  
       }
   }

  //no need to update aycStatus for other situations 
 }else if(lastOper[2]==2){//last selected operation is to delete from->to
   
   if(oper[2]==1&&(!acyStatus(i,j))&&toDe[i]&&fromAn[j]){//if operation under consideration is to add i->j which is also cyclic in last step, and if i is to Node descendant and j is from node ancester, then need to re-check 
    acyStatus(i,j)=acyclicCheck(i,j, parSet, "add"); //need to re-check acyclic status
   }

   if(oper[2]==3&&(!acyStatus(j,i))&&toDe[j]&&fromAn[i]){//if operation under consideration is to reverse i->j which is also cyclic in last step, and if j is to Node descendant and i is from  node ancester, then need to re-check
    acyStatus(j,i)=acyclicCheck(i,j, parSet, "reverse");//need to re-check acyclic status
   }
    
  //no need to update aycStatus for other situations 
 }else if(lastOper[2]==3){//lase selected operation is to reverse from->to
   if(oper[2]==1&&acyStatus(i,j)&&fromDe[i]&&toAn[j]){//if operation under consideration is to add i->j which is also acyclic in last step, and if i is from Node descendant and j is to node ancester, then leads to circle 
      acyStatus(i,j)=false; //need update
   }

   if(oper[2]==1&&(!acyStatus(i,j))&&toDe[i]&&fromAn[j]){//if operation under consideration is to add i->j which is also cyclic in last step, and if i is to Node descendant and j is from node ancester, then need to re-check 
    acyStatus(i,j)=acyclicCheck(i,j, parSet, "add"); //need to re-check acyclic status
   }


   if(oper[2]==3&&acyStatus(j,i)&&fromDe[j]&&toAn[i]){//if operation under consideration is to reverse i->j which is also acyclic in last step, and if j is from Node descendant and i is to  node ancester, then leads to circle 
       if(i!=lastOper[1]||j!=lastOper[0]){//if operation under consideration is not to reverse to->from
       acyStatus(j,i)=false; //need update  
       }
   }   

   if(oper[2]==3&&(!acyStatus(j,i))&&toDe[j]&&fromAn[i]){//if operation under consideration is to reverse i->j which is also cyclic in last step, and if j is to Node descendant and i is from  node ancester, then need to re-check
       acyStatus(j,i)=acyclicCheck(i,j, parSet, "reverse");//need to re-check acyclic status
   }

  //no need to update aycStatus for other situations 

 }else{//no last operation info.: this is the first step: need to check ayclic status
   if(oper[2]==1){//if operation under consideration is to add i->j
      acyStatus(i,j)=acyclicCheck(i,j, parSet, "add");
   } 
 
   if(oper[2]==3){//if operation under consideration is to reverse i->j
      acyStatus(j,i)=acyclicCheck(i,j, parSet, "reverse"); 
   }
 }

}


//create a list from a logical matrix: the ith element contains a vector of indices j such that (j,i)==true 
Rcpp::List listAdjPar(const Rcpp::LogicalMatrix& x){
  //input: x: p by p logical matrix of adjacent matrix; x(j,i)== true iff j->i in the graph 
  //output: list of parent sets of nodes in the graph 
  
  int p=x.nrow();
  
  Rcpp::List output;
  for (int i=0;i<p;i++){
    std::vector<int> temp;
    
     for (int j=0; j<p; j++){
       if(x(j,i)==true){
         temp.push_back(j);
       } 
     }
    output.push_back(temp); 
  }
  
  return output;
}


//create a list from a logical matrix: the jth element contains a vector of indices i such that (j,i)==true 
Rcpp::List listAdjChd(const Rcpp::LogicalMatrix& x){
  //input: x: p by p logical matrix of adjacent matrix; x(j,i)== true iff j->i in the graph 
  //output: list of children sets of nodes in the graph 
  
  int p=x.nrow();
  
  Rcpp::List output;
  for (int j=0;j<p;j++){
    std::vector<int> temp;
    
     for (int i=0; i<p; i++){
       if(x(j,i)==true){
         temp.push_back(i);
       } 
     }
    output.push_back(temp); 
  }
  
  return output;
}

//subset a vector by index indicator
Rcpp::NumericVector vecIndSel(const Rcpp::NumericVector& x, const Rcpp::LogicalVector& Ind){
  int len=x.size();
  int lenSel=Rcpp::sum(Ind);
  Rcpp::NumericVector output(lenSel);
  int counter=0; 
  for (int i=0; i<len; i++){
     if(Ind[i]){
       output[counter++]=x[i];
     }
  }

  return output; 
}



//subset a matrix by column indices and add the column of one 
Rcpp::NumericMatrix matColSelAddOne(const Rcpp::NumericMatrix& x, std::vector<int> idx){
//input: x- a matrix; idx: an integer vector - indices of columns to be selected 
//output: a matrix by selecting columns of x corresponding to idx, plus the column of ones
  
  int nrow=x.nrow();
  int ncol=idx.size();
  
  Rcpp::NumericMatrix output(nrow, ncol+1); //initilize a matrix of size nrow x ncol, filled with zero  
  int counter =0; 
  
  if(ncol>0){
    for (int i=0; i<ncol; i++){
      for (int j=0; j<nrow; j++){
        output[counter++]=x[j+idx[i]*nrow];
      }
    }
  }
  
  for (int j=0; j<nrow; j++){//add column of 1 as the last column
    output[counter++]=1;
  }
  
  return output; 
}

//subset a matrix by row indicators   
Rcpp::NumericMatrix matRowSelInd(const Rcpp::NumericMatrix& x,  const Rcpp::LogicalVector& indRow){
//input: x- a matrix; 
//indRow: a logical vector of length x.nrow() corresponding to rows to be selected 
//output: a matrix by selecting  rows of x corresponding to indRow.
  int nrow=x.nrow();
  int ncol=x.ncol();
  int nrowOut=Rcpp::sum(indRow);
  
  Rcpp::NumericMatrix output(nrowOut, ncol); //initilize a matrix of size nrowOut x ncolOut, filled with zero  
  int counter =0; 
  
    for (int i=0; i<ncol; i++){
      for (int j=0; j<nrow; j++){
        if(indRow[j]){
        output[counter++]=x[j+i*nrow];
        }
      }
    }
  
  return output; 

}



//subset a matrix by column index (and add the column of one ) and row indicators   
Rcpp::NumericMatrix matColSelAddOneRowSelInd(const Rcpp::NumericMatrix& x, const std::vector<int>& idCol, const Rcpp::LogicalVector& indRow){
//input: x- a matrix; 
// idCol: an integer vector - indices of columns to be selected; indRow: a logical vector of length x.nrow() corresponding to rows to be selected 
//output: a matrix by selecting  columns of x corresponding to idCol, plus the column of ones, and rows of x corresponding to indRow.
  int nrow=x.nrow();
  int ncolOut=idCol.size();
  int nrowOut=Rcpp::sum(indRow);
  
  Rcpp::NumericMatrix output(nrowOut, ncolOut+1); //initilize a matrix of size nrowOut x ncolOut, filled with zero  
  int counter =0; 
  
  if(ncolOut>0){
    for (int i=0; i<ncolOut; i++){
      for (int j=0; j<nrow; j++){
        if(indRow[j]){
        output[counter++]=x[j+idCol[i]*nrow];
        }
      }
    }
  }
  
  for (int j=0; j<nrow; j++){//add column of 1 as the last column
    if(indRow[j]){
    output[counter++]=1;
    } 
  }
  
  return output; 

}


//delete one column from a matrix 
Rcpp::NumericMatrix delCol(const Rcpp::NumericMatrix& x, int e){//use loop: this is slightly faster than using iterator 
  //inputs: x: a matrix; e: index of the column to be deleted from x 
  //output: a matrix by deleting column e from x
  
  int nrow=x.nrow(), ncol=x.ncol();
  
  Rcpp::NumericMatrix output=Rcpp::no_init_matrix(nrow, ncol-1); 
  int counter=0;
  for (int i=0;i<ncol; i++){
    for (int j=0;j<nrow; j++){
      if(i!=e){
        output[counter++]=x[j+i*nrow];
      }
    }
  }
  
  return output; 
}


//add a column before the first column of a matrix  
Rcpp::NumericMatrix addCol(const Rcpp::NumericMatrix& x, const Rcpp::NumericVector& v){
  //input: x: a matrix; v: a vector 
  //output: (v, x)
  
  int nrow=x.nrow(), ncol=x.ncol();
   
  Rcpp::NumericMatrix output=Rcpp::no_init_matrix(nrow, ncol+1);
  int counter=0;
  for(int i=0; i<(ncol+1); i++){
    for (int j=0; j<nrow; j++){
      if(i==0){
        output[counter++]=v[j];
      }else{
        output[counter++]=x[j+(i-1)*nrow];
      }
    }
  }
  
  return output;
}

//////////////////////////////////////////////////
// //Eigen implementation of  OLS 
// cross product of a matrix
MatrixXd AtA(const MatrixXd& A){
  const int n(A.cols());
  return MatrixXd(n,n).setZero().selfadjointView<Lower>().rankUpdate(A.adjoint()); 
}


//OLS  by LLT decomposition on XtX: = LLt; where L is a lower triangular matrix 
Rcpp::List fastLM2_(const Rcpp::NumericMatrix& x, const Rcpp::NumericVector& y){
  const MapMat X=Rcpp::as<MapMat>(x);
  const MapVec Y=Rcpp::as<MapVec>(y);
  const unsigned int n(X.rows());
  
  const LLT<MatrixXd> llt(AtA(X)); //LLT class object of XtX;
  const VectorXd betahat(llt.solve(X.adjoint()*Y)); //(XtX)^{-1}Xty
  const VectorXd fitted(X*betahat);
  const VectorXd resid(Y-fitted); 
  //const int df(n-p);
  const double s2=resid.squaredNorm()/double(n); 
  const double loglike= -double(n)/2*(log(2*M_PI)+log(s2)+1);
  
  
  return Rcpp::List::create(
    Rcpp::Named("coefficients") = betahat,
    Rcpp::Named("fitted.values") = fitted,
    Rcpp::Named("residuals")=resid,
    Rcpp::Named("loglikelihood")  = loglike); 
}


/////////////////////////////////////////////////
//logistic regression implementation: copied from RcppNumerical 
class LogisticReg: public MFuncGrad
{
private:
  const MapMat X;
  const MapVec Y;
  const int n;
  Eigen::VectorXd xbeta;  // contains X*beta
  Eigen::VectorXd prob;   // contains log(1+exp(X*beta)) and 1/(1+exp(-X*beta))
public:
  LogisticReg(const MapMat x_, const MapVec y_) :
  X(x_),
  Y(y_),
  n(X.rows()),
  xbeta(n),
  prob(n)
  {}
  
  double f_grad(Constvec& beta, Refvec grad)
  {
    // Negative log likelihood
    //   sum(log(1 + exp(X * beta))) - y' * X * beta
    xbeta.noalias() = X * beta;
    const double yxbeta = Y.dot(xbeta);
    // Calculate log(1 + exp(X * beta)), avoiding overflow
    for(int i = 0; i < n; i++)
      prob[i] = R::log1pexp(xbeta[i]);
    const double f = prob.sum() - yxbeta;
    
    // Gradient
    //   X' * (p - y), p = exp(X * beta) / (1 + exp(X * beta))
    //                   = exp(X * beta - log(1 + exp(X * beta)))
    prob = (xbeta - prob).array().exp();
    grad.noalias() = X.transpose() * (prob - Y);
    
    return f;
  }
  
  Eigen::VectorXd current_xb() const { return xbeta; }
  Eigen::VectorXd current_p()  const { return prob; }
};


Rcpp::List fastLR2_(const Rcpp::NumericMatrix& x, const Rcpp::NumericVector& y,
                   Rcpp::NumericVector start,
                   double eps_f, double eps_g, int maxit)
{
  const MapMat xx = Rcpp::as<MapMat>(x);
  const MapVec yy = Rcpp::as<MapVec>(y);
  // Negative log likelihood
  LogisticReg nll(xx, yy);
  // Initial guess
  Rcpp::NumericVector b = Rcpp::clone(start);
  MapVec beta(b.begin(), b.length());
  
  double fopt;
  int status = optim_lbfgs(nll, beta, fopt, maxit, eps_f, eps_g);
  if(status < 0)
    Rcpp::warning("algorithm did not converge");
  
  return Rcpp::List::create(
    Rcpp::Named("coefficients")      = b,
    Rcpp::Named("fitted.values")     = nll.current_p(),
    Rcpp::Named("linear.predictors") = nll.current_xb(),
    Rcpp::Named("loglikelihood")     = -fopt,
    Rcpp::Named("converged")         = (status >= 0)
  );
}


//calculate bic score given design matrix x, response vector y and response type ("c" or "b"); 
//for "c" nodes, fit OLS only one samples with indRow=1
double bicScoreSC(const Rcpp::NumericMatrix& x, const Rcpp::NumericVector& y, Rcpp::String nodeType, const Rcpp::LogicalVector& indRow){
  int p=x.ncol(), n=y.size();
  Rcpp::List res;
  
  if(nodeType=="c"){//for continuous nodes: fit OLS on samples with indRow=1
    Rcpp::NumericMatrix xS=matRowSelInd(x, indRow);
    Rcpp::NumericVector yS=vecIndSel(y, indRow);  
    res=fastLM2_(xS, yS);
    n= yS.size(); 

  }else{//for binary nodes: fit logistic regreassion on all samples 
    Rcpp::NumericVector start(p);
    res=fastLR2_(x, y, start, 1e-8, 1e-5,300);
   // int n=y.size(); 
  }
  
  double loglike=res["loglikelihood"];
  //Rcpp::Rcout<<"sample size: "<<n<<"\n";
  return -2.0*loglike+log(n)*(p-1); //for "c" nodes: sample size n = # of nonzeros
 }

//calculate BIC score given data, node type and parent set for all nodes 
Rcpp::NumericVector bicScoreAllSC(const Rcpp::NumericMatrix& Y, const Rcpp::CharacterVector& nodeType, const Rcpp::List& parSet){
  //inputs: Y: n by p data matrix; nodeTyep: character vector for each node type ("c": continuous; "b": binary)
  //parSet: list of parent set of each node 
  //output: bic score for each node neighborhood 
  int p=Y.ncol();
  Rcpp::NumericVector score(p);  
  Rcpp::LogicalVector curIndRow(p); 

  for(int i=0;i<p;i++){
    Rcpp::NumericVector curY=Y(Rcpp::_,i); //current node 
    Rcpp::NumericMatrix parX=matColSelAddOne(Y, parSet[i]);//current design matrix 

    if(nodeType[i]=="c"){//for "c" node, fit OLS on samples that current node is nonzero  
      curIndRow = (curY!=0);  
    }


    score[i]=bicScoreSC(parX, curY, nodeType[i], curIndRow);
  }

  return score;
}


//update curGraph, curPar and curScore based on the selected operation: curOper
void updateGraph(const Rcpp::IntegerVector curOper, const Rcpp::NumericVector scoreMin, Rcpp::LogicalMatrix& curGraph, Rcpp::List& curPar, Rcpp::List& curChd, Rcpp::NumericVector& curScore){
//inputs: cuOper: 3 x 1 integer vector of selected operation; scoreMin: 2 x 1 numeric vector of scores involved with the selected operation; 
 //curGraph: p x p logical matrix of current graph to be updated; curPar: list of current parent sets to be updated; curChd: list of current children sets to be updated; curScore: p x 1 of current node scores to be updated
 //output: none 

    int j=curOper[0], i=curOper[1];

    if(curOper[2]==1){//selected operation is addition of j->i
      curGraph(j,i)=true; //add j->i to current graph
      std::vector<int> par_i=curPar[i];
      par_i.push_back(j);
      curPar[i]=par_i; //add node j to node i parent set 
       
      std::vector<int> chd_j=curChd[j];
      chd_j.push_back(i);
      curChd[j]=chd_j;//add node i to node j children set

      curScore[i]=scoreMin[0]; //update node i score 
    }

    if(curOper[2]==2){//selected operation is deletion of j->i
      curGraph(j,i)=false; //delete j->i from current graph
      std::vector<int> par_i=curPar[i];
      std::vector<int>::iterator it=std::find(par_i.begin(), par_i.end(), j); //find position of node j in node i parent set
      par_i.erase(it); 
      curPar[i]=par_i; //delete node j from node i parent set 

      std::vector<int> chd_j=curChd[j];
      std::vector<int>::iterator  it2=std::find(chd_j.begin(), chd_j.end(),i);
      chd_j.erase(it2);
      curChd[j]=chd_j; //delete node i from node j children set  

      curScore[i]=scoreMin[0]; //update node i score 
    }

    if(curOper[2]==3){//selected operation is reversal of j->i
      //first delete j->i:
      curGraph(j,i)=false; //delete j->i from current graph
      std::vector<int> par_i=curPar[i];
      std::vector<int>::iterator it=std::find(par_i.begin(), par_i.end(), j); //find position of node j in node i parent set
      par_i.erase(it); 
      curPar[i]=par_i; //delete node j to node i parent set 
      
      std::vector<int> chd_j=curChd[j];
      std::vector<int>::iterator  it2=std::find(chd_j.begin(), chd_j.end(),i);
      chd_j.erase(it2);
      curChd[j]=chd_j; //delete node i from node j children set  

      curScore[i]=scoreMin[0]; //update node i score 

      //then add i->j:
      curGraph(i,j)=true; //add i->j to current graph
      std::vector<int> par_j=curPar[j];
      par_j.push_back(i);
      curPar[j]=par_j; //add node i to node j parent set 

      std::vector<int> chd_i=curChd[i];
      chd_i.push_back(j);
      curChd[i]=chd_i;//add node j to node i children set

      curScore[j]=scoreMin[1]; //update node j score 
    }

}

//comapre current score change (delta_c) with current deltaMin:
//if delta_c<deltaMin-tol, return true: (and subsequently update current best operation and its relates score, as well as deltaMin)
//if flip&&|delta_c-deltaMin|<=tol, then return true with 50% probability;
//if delta_c>deltaMin+tol: return false (no subsequent update of current best operation);  
bool deltaComp(double delta_c, double deltaMin, double tol, bool flip = true){
 if(delta_c<deltaMin-tol){
  return true;
 }else if(delta_c>deltaMin+tol){
  return false;
 }else{
   if(flip&&(rand()+0.0)/(RAND_MAX+1.0)>0.5)
    return true;
 }

  return false; 
}

//main function: hcSC_: consier SC features -- zero inflated: for "c" nodes fit OLS on nonzero samples  
//run hc once: 
Rcpp::List hcSC1(const Rcpp::NumericMatrix& Y, const Rcpp::CharacterVector& nodeType,const Rcpp::LogicalMatrix&  whiteList,
                       const Rcpp::LogicalMatrix&  blackList, double tol=1e-6, unsigned int maxStep=500, unsigned int seed=1, bool flip=true,  bool verbose=false){
  //inputs: Y: n by p data matrix;
  //whiteList: p by p whitelist: edges always includes ; blackList: p x p blacklist: edges excluded;
  //tol: tolerance for score imporvement: stop search when score improvement is less than tol; maxStep: maximum search steps; 
  //seed: seed for random number generator srand(); verbose: whether print out intermediate results 
  //flip: whether (TRUE) or not (FALSE) flip a coin to break ties in score improvement; 
  //output: a list of: adjacency matrix of the final graph, scores of the final graph, selected operation in each step, best score improvement in each step

  srand(seed); // set random number generator seed: used to break score ties in detalComp(); 
  int p=Y.ncol();
  Rcpp::List stepOper; // a list of selected operation at each step
  Rcpp::List stepDelta; // a list of best score improvement at each step
  
  //// initialization step
  //initialize graph by whiteList: update at end of each step 
  Rcpp::LogicalMatrix  curGraph=Rcpp::clone(whiteList); 

  //initialize parent set (neighborhoods) and children set: update at the end of each step
  Rcpp::List curPar=listAdjPar(curGraph);
  Rcpp::List curChd=listAdjChd(curGraph); 
  
  //initialize neighborhood score: upadate at the end of each step: for "c" nodes fit OLS on nonzero samples  
  Rcpp:: NumericVector curScore=bicScoreAllSC(Y, nodeType,curPar); //current score of each neighborhood 
  
  //initialize score change by operations: start at  NA_REAL; update for each potential operation within each step
  Rcpp::NumericMatrix delta(p,p); //record change of score (delta) by an operation  (add, delete or reverse an edge) to the current graph:
  // if neither i->j nor j->i in the current graph, then delta(i,j) and delta(j,i) record  score change by adding i->j, j->i, respectively; 
  // if j->i in the current graph, so that i->j not in, then delta(j,i) records score change by deleting j->i, and delta(i,j) change by reversing j->i; 
  std::fill(delta.begin(), delta.end(), Rcpp::NumericVector::get_na()); //fill with NA 
  
  //initialize acyclic status: start at NA_LOGICAL; update for each potential operation within each step
  Rcpp::LogicalMatrix  acyStatus(p,p); //record whether an operation (add or reverse an edge) to the current graph maintains acyclicity (true) or not (false): 
  // if neither i->j nor j->i in the current graph, then acyclic(i,j) and acyclic(j,i) record status by adding i->j, j->i, respectively;
  // if j->i in the current graph, so that i->j not in, then acyclic(i,j) record status by reversing j->i and acyclic(j,i) = true (since current graph must be acyclic);
  //note: both i->j and j->i in the current graph is not possible. 
  std::fill(acyStatus.begin(), acyStatus.end(), Rcpp::LogicalVector::get_na()); //fill with NA
  
  ////updating steps : update curGraph/curPar, curScore by selected operations (add, delete or reverse an edge); update aycStatus and delta for each operation within each step;
  
  //record selected operation, curOper, and the corresponding score change, deltaMin; add curOper to the list of selected operations, stepOper. 
  Rcpp::IntegerVector curOper(3,-1); //current operation: curOper[2], 1-add, 2-delete, 3-reverse, edge curOper[0]->curOper[1]; oper[2], -1, no movement (starting and stopping search)
  
  Rcpp::IntegerVector lastOper(3,-1); //last selected operation
  Rcpp::LogicalVector fromAn(p, false), fromDe(p, false), toAn(p,false), toDe(p, false); //ancester and descendant indicators for the from and to nodes of the last selected operation
 

  double deltaMin=0.0; //best score improvement (decrease) at each step: reset to -tol at each step
  Rcpp::NumericVector scoreMin(2); //if selected operation is add/delete j->i, then score[0] is the updated score of node i; if the selected operation is reverse j->i, then score[0] is the updates score of node i and score[1] is the updated score of node j 
  
  double score_c;//score of current operation if deletion or addition
  double score_cr; //score of current operation if reversal 
  double delta_c;//score change of current operation

  Rcpp::LogicalVector indRow_i(p); //sample indicator for continuous nodes : in addition and deletion operations
  Rcpp::LogicalVector indRow_j(p); // for reversal operations 

  ////updating steps:
  unsigned int stepCounter=1; //step counter 
  while(stepCounter<=maxStep){
    for (int i=0;i<p;i++){//consider operations invovling node i
       
      std::vector<int> par_i=curPar[i]; //node i parent set 
      int size_i=par_i.size();
      Rcpp::NumericMatrix X_i=matColSelAddOne(Y, par_i); //node i design matrix 
      Rcpp::NumericVector Y_i=Y(Rcpp::_,i);//node i vector 
      
      if(nodeType[i]=="c"){//for "c" node, fit OLS on samples that current node is nonzero  
       indRow_i = (Y_i!=0);  
      }

      bool last_from_i=(lastOper[0]==i); //whether last selected operation involves node i
      bool last_to_i=(lastOper[1]==i); 

      if(size_i>0){//if node i parent set non-empty: consider deletion and reversal of node i parents
         
          for(int k=0; k<size_i; k++){
             int j=par_i[k]; //kth parent of node i
             acyStatus(j,i)=true;//current edge must be acyclic 

             if(!whiteList(j,i)){//if edge j->i  not in the whitelist:  consider deletion and reversal of j->i 
              
                //deletion:  no need to check/update ayclicity; calculate score change; compare with current deltaMin and make updates on deltaMin and curOper if smaller    
                if(lastOper[2]==-1||last_to_i||(lastOper[2]==3&&last_from_i)){//if node i is involved in last selected operation or no last operation: need to recalcuate its score  
                  Rcpp::NumericMatrix X_c= delCol(X_i,k); //node i design matrix after deletion of kth parent (i.e., node j)
                  score_c=bicScoreSC(X_c, Y_i, nodeType[i], indRow_i); //node i score after deletion of kth parent 
                  delta_c=score_c-curScore[i]; 
                  delta(j,i)=delta_c;// record score change due to deletion of j->i;
                }else{//no need to recalcuate node i score: use calculation from last step
                   delta_c=delta(j,i); 
                   score_c=delta_c+curScore[i];
                }


               if(verbose){
               Rcpp::Rcout <<"\n step count:"<<stepCounter <<": delete edge "<< j<<"->"<<i<< ": delta="<<delta_c<<"\n";
               }

                if(deltaComp(delta_c, deltaMin, tol,flip)){//if score change improved over the current best improvement, upadte curOper and deltaMin, scoreMin 
                  curOper[0]=j;//from node
                  curOper[1]=i;//to node
                  curOper[2]=2; //delete j->i improves score beyond the current best improvement  
                  deltaMin=delta_c;
                  scoreMin[0]=score_c; //score of node i after deleting j->i
                }
                
                //reversal:  check ayclicity; calculate score change; compare with current deltaMin and make updates on deltaMin and curOper if smaller
                if(!blackList(i,j)){//if i->j not in blacklist (since j->i in current graph, so i->j must not in current graph): consider reversal of j->i
                
                acyclicUpdate(lastOper, fromAn, fromDe, toAn, toDe, Rcpp::IntegerVector::create(j,i,3), curPar, acyStatus);//check/update acyclic status

                 if(acyStatus(i,j)){

                    bool last_from_j=(lastOper[0]==j);   //whether last selected operation (lastOper) involves node j 
                    bool last_to_j=(lastOper[1]==j); 

                    if(lastOper[2]==-1||last_to_j||(lastOper[2]==3&&last_from_j)){//if node j is involved in the last selected operation or no last operation: need to recalculate its score

                     Rcpp::NumericMatrix X_j=matColSelAddOne(Y, curPar[j]); //node j design matrix 
                     Rcpp::NumericMatrix X_c=addCol(X_j,Y_i);  //add node i vector to the design matrix
                     Rcpp::NumericVector Y_j=Y(Rcpp::_,j);//node j vector 

                     if(nodeType[j]=="c"){//for "c" node, fit OLS on samples that current node is nonzero  
                      indRow_j = (Y_j!=0);  
                     }

                     score_cr=bicScoreSC(X_c, Y_j, nodeType[j], indRow_j); //node j score after adding node i into its parent 
                     delta_c=score_cr-curScore[j];//change due to adding i->j
                     delta_c=delta_c+delta(j,i); //change of reversing j->i = change due to dropping j->i+change due to adding i->j
                     delta(i,j)=delta_c; //record score change due to reversal of j->i
                    
                     }else{//no need to recalculate node j score: use calculation from last step
                     delta_c=delta(i,j)+delta(j,i);
                     score_cr=delta(i,j)+curScore[j];
                     }
                       
                   if(verbose){
                   Rcpp::Rcout <<"\n step count:"<<stepCounter <<": reverse edge "<< j<<"->"<<i<< ": delta="<<delta_c<<"\n";
                   }

                    if(deltaComp(delta_c, deltaMin, tol,flip)){//if score change improved over the current best improvement, upadte curOper and deltaMin, scoreMin 
                      curOper[0]=j;//from node
                      curOper[1]=i;//to node
                      curOper[2]=3; //reverse j->i improves score beyond the current best improvement  
                      deltaMin=delta_c;
                      scoreMin[0]=score_c; //score of node i after deleting j->i
                      scoreMin[1]=score_cr;//score of node j after adding i->j
                    }  

                 } //end of if acyclic for reversal 
              
               } //end of if black for reversal 

            }//end of if white for deletion/reversal 

        }//end of k-loop  for i node parents 
     } //end of if size>0  for i node neighborhood size
      
    //addition operations: to node i parent set
     for (int j=0; j<p; j++){
       if(!curGraph(j,i)&&!curGraph(i,j)&&!blackList(j,i)){// if j->i  and i->j not in current graph and not in blacklist: consider addition of j->i     
        acyclicUpdate(lastOper, fromAn, fromDe, toAn, toDe, Rcpp::IntegerVector::create(j,i,1), curPar, acyStatus);//check/update acyclic status
    
         if(acyStatus(j,i)){ 
         
           if(lastOper[2]==-1||last_to_i||(lastOper[2]==3&&last_from_i)){//if node i is involved in the last selected operation or no last operation: need to recalculate its score 
            Rcpp::NumericVector Y_j=Y(Rcpp::_,j);//node j vector  
            Rcpp::NumericMatrix X_c=addCol(X_i,Y_j);//node i parent matrix after adding j->i
            score_c=bicScoreSC(X_c, Y_i, nodeType[i], indRow_i);//node i score after adding j->i
        
            delta_c=score_c-curScore[i];//score change due to adding j->i
            delta(j,i)=delta_c;//record score change due to adding j->i
            }else{//no need to recalcuate node i score: use calculation from last step
             delta_c=delta(j,i);
             score_c=delta_c+curScore[i];
            }

            if(verbose){
            Rcpp::Rcout <<"\n step count:"<<stepCounter <<": add edge "<< j<<"->"<<i<< ": delta="<<delta_c<<"\n";
            }
             
           
            if(deltaComp(delta_c, deltaMin, tol,flip)){//if score change improved over the current best improvement, upadte curOper and deltaMin, scoreMin 
              curOper[0]=j;//from node
              curOper[1]=i;//to node
              curOper[2]=1; //add j->i improves score beyond the current best improvement  
              deltaMin=delta_c;
              scoreMin[0]=score_c; //score of node i after adding j->i
            }
            
            
          }//end of if ayclic 
        }//end of if  curGraph&&blacklist for addition
      }// end of j-loop       
    }// end of i-loop  
      
//update curGraph/curPar/curScore and continue searching if there is sufficient improvement; otherwise, breakout the while loop

    if(deltaMin<(-tol)){
      //update curGraph, curPar, curChd, curScore based on the selected operation (curOper) 
      updateGraph(curOper, scoreMin, curGraph, curPar, curChd, curScore);
     
      //update lastOper, fromAn, fromDe, toAn, toDe based on the updated graph 
      lastOper=Rcpp::clone(curOper); 
     
      if(verbose){
      Rcpp::Rcout <<"\n step count:"<<stepCounter<<": lastOper"<<lastOper[0]<<lastOper[1]<<lastOper[2]<<"\n";   
      }

      fromAn=ancester(lastOper[0], curPar);
      toAn=ancester(lastOper[1], curPar);
      fromDe=descend(lastOper[0], curChd);
      toDe=descend(lastOper[1], curChd);
      
      //update  stepOper and stepCounter
      stepOper.push_back(lastOper);
      stepDelta.push_back(deltaMin);
      stepCounter++;   

      //reset deltaMin to the tolerance 
      deltaMin=0.0;
     
      }else{//stop searching when there is no sufficient improvement 

      break;
     }//end if deltaMin<-tol

  }//end while stepCounter  

////
  return Rcpp::List::create(Rcpp::Named("adjacency")=curGraph, Rcpp::Named("score")=curScore, Rcpp::Named("operations")=stepOper, Rcpp::Named("deltaMin")=stepDelta);
}

//run hc restart times with different seeds and return the model with the best (smallest) score
//[[Rcpp::export]]
Rcpp::List hcSC_(const Rcpp::NumericMatrix& Y, const Rcpp::CharacterVector& nodeType,const Rcpp::LogicalMatrix&  whiteList,
                       const Rcpp::LogicalMatrix&  blackList, double tol=1e-6, unsigned int maxStep=500,  unsigned int restart=10, unsigned int seed=1, bool verbose=false){

  Rcpp::List res;
  double scoreMin=R_PosInf;   
  bool flip = restart>1; //when there is only one run, no need to randomize ties; this saves a  bit time when there're many ties 


  for (unsigned int i=0; i<restart; i++){
   Rcpp::List curRes=hcSC1(Y, nodeType, whiteList, blackList, tol, maxStep, i*101+seed, flip, verbose);
   Rcpp::NumericVector curScoreV=curRes["score"];
   double curScore=Rcpp::sum(curScoreV);
   //Rcpp::Rcout<<"restart iteration "<<i+1<< ": score= "<<curScore<<"\n";
   if(curScore<scoreMin){
    res=curRes;
    scoreMin=curScore;
   }
  }
  
  return res;
} 

#pragma GCC diagnostic pop
