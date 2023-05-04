
#' @title  Design the clinical trial for continuous outcome
#' @description This function generates the design parameters of a clinical trial for continuous outcome.
#' @param delta0 Standardized effect size in ineffective arm.
#' @param delta1 Standardized effect size in effective arm.
#' @param alpha Type I error.
#' @param beta Type II error.
#' @param K Number of treatment arms.
#' @param frac Vector of fractions for information time at each look.
#' @return List of cumulative sample size for each stage of treatment and control groups along with maximum total sample size of the trial. It also provides efficacy and futility boundaries of the trial.
#' @examples
#' design_cont(delta0=0.178,delta1=0.545,alpha = 0.05, beta = 0.1, K = 4,frac=c(1/2,1))
#' @import stats
#' @export


design_cont<-function(delta0,delta1, alpha, beta, K,frac){
  n<-Size_cont(delta0=delta0,delta1=delta1,alpha = alpha, beta = beta, K = K) 
  mat <- matrix(NA,nrow=2,ncol=length(frac))
  rownames(mat) <- c("Cumulative sample size for treatment group", "Cumulative sample size for control group")
  colnames(mat)<-paste("Stage",1:length(frac))
  mat[1,] <- ceiling(n*frac)
  mat[2,] <- ceiling(n*frac)
  
  mat1<- matrix(NA,nrow=2,ncol=length(frac))
  rownames(mat1) <- c("Lower bound", "Upper bound")
  colnames(mat1)<-paste("Stage",1:length(frac))
  bv<-SCPRT(alpha = alpha,K=K,frac = frac)
  mat1[1,] <- bv$lshape
  mat1[2,] <- bv$ushape
  p=list("Sample size"=mat,"Maximum total sample size for the trial"=(K+1)*n,"Boundary values"=mat1)
  return(p)
  }
