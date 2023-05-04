
#' @title  Design the clinical trial for ordinal outcome
#' @description This function generates the design parameters of a clinical trial for ordinal outcome.
#' @param alpha Type I error.
#' @param beta Type II error.
#' @param K Number of treatment arms.
#' @param prob Probability of ordinal outcomes in control group.
#' @param or0 Odds ratio of ineffective treatment group vs control. 
#' @param or Odds ratio of effective treatment group vs control.
#' @param frac Vector of fractions for information time at each look.
#' @return List of cumulative sample size for each stage of treatment and control groups along with maximum total sample size of the trial. It also provides efficacy and futility boundaries of the trial.
#' @examples 
#' design_ord(0.05,0.1,K=4,c(0.075, 0.182, 0.319, 0.243, 0.015, 0.166),or=3.06,or0=1.32,c(1/2,1))
#' @import stats
#' @export

design_ord=function(alpha,beta,K,prob,or0,or,frac){
  n=Size_ord(alpha,beta,K,prob,or0,or)
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
  
  
  
  
  
  
  
  
  
  