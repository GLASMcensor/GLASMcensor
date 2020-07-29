# library(foreach)
# library(doParallel)

# library(Rcpp)
# library(igraph)
# library(mvtnorm)
# library(MASS)
# library(Matrix)

# library(huge)
# library(RSpectra)

#' The Personalized PageRank Augmented Screening method for the screening measure
#' @name Net_FS
#' @param A is the adjacency matrix
#' @param Orig_RM the original ranking measures
#' @param d the the parameter balances the original ranking measure and the network information
#' @param r the size of the variables set including all active predictors.
#' @param max_it  maximum iteration loops.
#' @param tol the stop iteration condition
#' @param threshold set the measures smaller than threshold as 0. default 0.
#' @return Net_RM: the network-based ranking measure.
#' @export

Net_FS<-function (A, Orig_RM,d,r,tol=1e-6,max_it=100,threshold=0) 
{
  #####  INPUT #####
  # A: adjacency matrix
  # Orig_RM: original ranking measures
  # d: the parameter balances the original ranking measure and the network information. d should be in [0,1] empirically 0.85.
  # r: the size of the variable set including all active predictors
  
  ### OUTPUT ####
  # Net_RM: the network-based ranking measure
  Orig_RM[Orig_RM<threshold]=0
  p=length(Orig_RM)
  degree=colSums(A)
  D2=as(matrix(1/colSums(A),1,p),"dgCMatrix")
  A_norm=as(A*(matrix(1,p,1)%*%D2),"dgCMatrix")
  rank_final=matrix(0,p,1)
  temp_matrix=A_norm
  Orig_RM[sort(Orig_RM,index.return=TRUE)$ix[1:(p-r)]]=1e-20
  Orig_RM_norm= matrix(Orig_RM/sum(Orig_RM),p,1)
  diff=1
  max_iterate=1
  rank=matrix(1/p,p,1)
  id=degree>1
  G2=Orig_RM_norm
  temp_matrix2=temp_matrix[id,id]
  while (diff>tol&& max_iterate<max_it){
    G2[id]=matrix(temp_matrix2%*%rank[id])
    rank2=(1-d)*Orig_RM_norm+d*G2
    diff=sum(abs(rank-rank2))
    rank=rank2
    max_iterate=max_iterate+1
  }
  rank_final=rank
  return(rank_final)
}


#' The Graphical Laplacian Augmented Screening method for the screening measure
#' @param amatrix is the adjacency matrix
#' @param utility the original ranking measures
#' @param alpha the the parameter balances the original ranking measure and the network information
#' @param max_it  maximum iteration loops.
#' @param threshold set the measures smaller than threshold as 0. default 0.
#' @return r: the Graphical Laplacian Augmented measure.
#' @export
abrw=function(amatrix,utility,alpha,max_it=100,threshold=0)
{
  #uti=utility
  utility[utility<threshold]=0
  utility=utility/sum(utility)
  pdim=dim(amatrix)[1]
  diag(amatrix)=rep(0,pdim)
  dgc_A=as(amatrix, "dgCMatrix")
  dgc_D=as(diag(rowSums(dgc_A),pdim,pdim), "dgCMatrix")
  #A_norm=as(diag(1/degree,pdim,pdim), "dgCMatrix")%*%dgc_a
  iden=as(diag(1,pdim), "dgCMatrix")
  L=dgc_D-dgc_A
  r=matrix(utility,1,pdim)%*%SparseM::solve(alpha*L+iden)
  return(as.vector(r))
 }


#' The KM for censoring time 
#' @param Ti is the vector of the ordered censored Time (ascent)
#' @param censor is the censor indicator for the Ti. 1 means no censor and 0 means censoring.
#' @return weight: The inverse weight of the C time. used in the cMBKRall function.
#' @export
 KMestimate_CTime = function(Ti,censor){
  N = length(Ti)
  risk = N:1
  
  #Ti no same
  risk=(risk-1)/risk
  #risk*censor
  #risk[N]=1
  risk[censor==1]=1
  
  weight=1/cumprod(risk)
  if(weight[N]==Inf){
    weight[N]=0
  }
  return(weight)
}


KMestimate_KM = function(Ti,censor){
  N = length(Ti)
  risk = N:1
  
  #Ti no same
  risk=(risk-1)/risk
  #risk*censor
  #risk[N]=1
  risk[censor==0]=1

  return(cumprod(risk))
  #return(risk)
}
