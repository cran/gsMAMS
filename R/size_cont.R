
#' @title  Calculates the Sample Size for a Clinical Trial
#' @description This function calculates the sample size per arm of a clinical trial for continuous outcome.
#' @param delta0 Standardized effect size in ineffective arm.
#' @param delta1 Standardized effect size in effective arm.
#' @param alpha Type I error.
#' @param beta Type II error.
#' @param K Number of treatment arms.
#' @return A numeric value indicating the sample size per arm.
#' @examples
#' Size_cont(delta0=0.178,delta1=0.545,alpha = 0.05, beta = 0.1, K = 4)
#' @import stats
#' @import mvtnorm
#' @keywords internal
#' @noRd

 
Size_cont=function(delta0,delta1, alpha, beta, K){
  r=1
  delta=c(delta1,rep(delta0,K-1))
  if(K==1){
    z0=qnorm(1-alpha)
    z1=qnorm(1-beta)
    n=ceiling((z0+z1)^2*(1+r)/(r*delta^2)) }
  if (K>=2) {
    Sigma=matrix((1/(1+r)),K,K)
    diag(Sigma)=1
    root=function(c){
      alpha-(1-mvtnorm::pmvnorm(lower=rep(-Inf,K), upper=rep(c, K),
                       mean=rep(0,K), sigma=Sigma)[1])
    }
    c=uniroot(root, lower=0, upper=999)$root
    2
    Sigma11=Sigma[1:(K-1),1:(K-1)]
    Sigma12=Sigma[1:(K-1),K]
    Sigma21=Sigma[K,1:(K-1)]
    A=(-1)*diag(K)
    A[,1]=1; A=rbind(A[-1,], A[1,])
    B=A%*%Sigma%*%t(A)
    root1=function(n){
      mu=sqrt(r*n/(1+r))*delta
      b=as.numeric(A%*%mu)
      int=mvtnorm::pmvnorm(lower = c(rep(0,K-1),c), upper=rep(Inf,K),
                  mean = b, sigma = B)[1]
      1-beta-as.double(int)
    }
    n=ceiling(uniroot(root1, lower=1, upper=999)$root)
  }
   return(n)
}
