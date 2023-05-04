
#' @title  Efficacy and Futility boundary values
#' @description This function calculates the upper and lower bound values for MAMS trial.
#' @param alpha Type I error.
#' @param K Number of treatment arms.
#' @param frac Vector of information time at each look.
#' @return A list of three elements: critical.value, lshape(Futility), and ushape(Efficacy).
#' @examples
#' SCPRT(alpha = 0.05, K = 3, frac = c(1/3, 2/3, 1))
#' @import stats
#' @import mvtnorm
#' @keywords internal
#' @noRd
 

SCPRT=function(alpha, K, frac){
  r=1
  J=length(frac)
  Sigma=matrix((1/(1+r)) ,K,K)
  diag(Sigma)=1
  if (J>=10) stop("Limit the maximum number of stages to ten")
  root=function(c){
    int=mvtnorm::pmvnorm(lower=rep(-Inf,K), upper=rep(c, K), mean=rep(0,K), sigma=Sigma)[1]
    alpha-(1-int)
  }
  c=round(uniroot(root, lower=0, upper=999)$root,3)
  a=c(1.645, 2.109, 2.645, 2.953, 3.166, 3.327, 3.456, 3.562, 3.652, 3.729)
  a=a[J]
  l=round((c*frac-sqrt(2*a*frac*(1-frac)))/sqrt(frac),3)
  u=round((c*frac+sqrt(2*a*frac*(1-frac)))/sqrt(frac),3)
  l[length(frac)]<-u[length(frac)]
  ans=list(critical.value=c, lshape=l, ushape=u)
  return(ans)
}
