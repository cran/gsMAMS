
#' @title  Calculates the Sample Size for a Clinical Trial
#' @description This function calculates the sample size per arm of a clinical trial of ordinal outcome.
#' @param alpha Type I error.
#' @param beta Type II error.
#' @param K Number of treatment arms.
#' @param prob Probability of ordinal outcomes in control group.
#' @param or0 Odds ratio of ineffective treatment group vs control. 
#' @param or Odds ratio of effective treatment group vs control.
#' @return A numeric value indicating the sample size per arm.
#' @examples 
#' Size_ord(prob=c(0.075, 0.182, 0.319, 0.243, 0.015, 0.166),or=3.06, or0=1.32,alpha=0.05,beta=0.1,K=4)
#' @import stats
#' @import mvtnorm
#' @keywords internal
#' @noRd

Size_ord=function(alpha,beta,K,prob,or0, or)
{
  if (length(prob)>2)
  {
    p0=prob
    Q0=cumsum(p0)
    Qk=1/(1+((1-Q0)/Q0)/or)
    pk=c(Qk[1], diff(Qk, lag=1))
    pbar=(p0+K*pk)/(K+1)
    q <- (1 - sum(pbar^3))/3
    sigma <- 1/sqrt(q)
  }
  if (length(prob)==2)
  {
    p0=prob[1]
    p1=(1+(1-p0)/(p0*or))^(-1)
    pbar=c((p0+p1)/2, 1-(p0+p1)/2)
    2
    q <- (1 - sum(pbar^3))/3
    sigma <- 1/sqrt(q)
  }
  delta=c(log(or),rep(log(or0),K-1))
  if (K==1) {
    V=(qnorm(1-alpha)+qnorm(1-beta))^2/log(or)^2
    n=ceiling(2*V/q) ## sample size per arm ###
    return(n)
  }
  if (K>=2) {
    Sigma=matrix(0.5,K,K)
    diag(Sigma)=1
    root=function(c){
      alpha-(1-mvtnorm::pmvnorm(lower=rep(-Inf,K), upper=rep(c, K), mean=rep(0,K), sigma=Sigma)[1])
    }
    c=uniroot(root, lower=0, upper=999)$root
    Sigma11=Sigma[1:(K-1),1:(K-1)]
    Sigma12=Sigma[1:(K-1),K]
    Sigma21=Sigma[K,1:(K-1)]
    A=(-1)*diag(K)
    A[,1]=1; A=rbind(A[-1,], A[1,])
    B=A%*%Sigma%*%t(A)
    root1=function(n){
      mu=sqrt(n/(2*sigma^2))*delta
      b=as.numeric(A%*%mu)
      int=mvtnorm::pmvnorm(lower = c(rep(0,K-1),c), upper=rep(Inf,K), mean = b, sigma = B)[1]
      3
      1-beta-as.double(int)
    }
    n=ceiling(uniroot(root1, lower=1, upper=999)$root)
  }
  return(n)
}
