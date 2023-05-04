
#' @title  Design the clinical trial for survival outcome
#' @description This function generates the design parameters of a clinical trial for survival outcome.
#' @param m0 Median survival time of control group.
#' @param alpha Type I error.
#' @param beta Type II error.
#' @param K Number of treatment arms.
#' @param HR0 Hazard ratio of ineffective treatment group vs control.
#' @param HR1 Hazard ratio of effective treatment group vs control.
#' @param ta Accrual time.
#' @param tf Follow-up time.
#' @param kappa Shape parameter (kappa=1 for exponential distribution).
#' @param eta Rate of loss to follow-up.
#' @param frac Vector of fractions for information time at each look.
#' @return List of cumulative number of events for each stage of combined treatment and control groups along with total number of subjects and maximum total number of events for the trial. It also provides efficacy and futility boundaries of the trial.
#' @examples 
#' design_surv(m0=20,HR0=1, HR1=0.65, ta=20,tf=40,alpha=0.05,beta=0.1,K=3,kappa=1,eta=0,frac=c(1/2,1))
#' @import stats
#' @export

design_surv=function(m0,alpha,beta, K,HR0, HR1, ta, tf,kappa,eta,frac){
  n=Size_surv(m0,alpha,beta, K,HR0, HR1, ta, tf,kappa,eta,frac)
  mat <- matrix(NA,nrow=1,ncol=length(frac))
  rownames(mat)<- c("Cumulative number of events for combined treatment & control")
  colnames(mat)<-paste("Stage",1:length(frac))
  mat[1,] <- ceiling(2*n[1]*frac)
     
  mat1<- matrix(NA,nrow=2,ncol=length(frac))
  rownames(mat1) <- c("Lower bound", "Upper bound")
  colnames(mat1)<-paste("Stage",1:length(frac))
  bv<-SCPRT(alpha = alpha,K=K,frac = frac)
  mat1[1,] <- bv$lshape
  mat1[2,] <- bv$ushape
  p=list("Sample size"=mat,"Maximum total number of events for the trial"=(K+1)*n[1],"Total number of subjects required for the trial"=n[2]*(K+1),"Boundary values"=mat1)
  return(p)
  
}
