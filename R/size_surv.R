

#' @title  Calculates the Sample Size for a Clinical Trial
#' @description This function calculates the sample size per arm of a clinical trial for survival outcome.
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
#' @return A numeric value indicating the sample size per arm.
#' @examples 
#' Size_surv(m0=20,HR0=1, HR1=0.65, ta=20,tf=40,alpha=0.05,beta=0.1,K=3,kappa=1,eta=0,frac=c(1/2,1))
#' @import stats
#' @import mvtnorm
#' @keywords internal
#' @noRd

 

Size_surv=function(m0,alpha,beta, K,HR0, HR1, ta, tf,kappa,eta,frac)
{
  if(K<2 | K>5){stop("K should be between 1 and 6.")}
  c=SCPRT(alpha=alpha,K=K,frac=frac)$critical.value
  hr0=HR0              #exp(-delta0)
  hr1=HR1              #exp(-delta1)
  lambda0=log(2)/m0^kappa
  
  if (K==2){
    lambda1=lambda0*hr1
    lambda2=lambda0*hr0
  }
  
  else if(K==3){
    lambda1=lambda0*hr1
    lambda2=lambda0*hr0
    lambda3=lambda0*hr0
  }
  
  else if(K==4) {
    lambda1=lambda0*hr1
    lambda2=lambda0*hr0
    lambda3=lambda0*hr0
    lambda4=lambda0*hr0
  }
  
  else if(K==5) {
    lambda1=lambda0*hr1
    lambda2=lambda0*hr0
    lambda3=lambda0*hr0
    lambda4=lambda0*hr0
    lambda5=lambda0*hr0
  }
  
  
  
  tau=ta+tf
  
  S0=function(t){exp(-lambda0*t^kappa)}
  h0=function(t){kappa*lambda0*t^(kappa-1)}
  
  
  if (K==1){
    S1=function(t){exp(-lambda1*t^kappa)}
    h1=function(t){kappa*lambda1*t(kappa-1)}
  }
  
  else if (K==2){
    S1=function(t){exp(-lambda1*t^kappa)}
    h1=function(t){kappa*lambda1*t^(kappa-1)}
    S2=function(t){exp(-lambda2*t^kappa)}
    h2=function(t){kappa*lambda2*t^(kappa-1)}
  }
  
  else if(K==3){
    S1=function(t){exp(-lambda1*t^kappa)}
    h1=function(t){kappa*lambda1*t^(kappa-1)}
    S2=function(t){exp(-lambda2*t^kappa)}
    h2=function(t){kappa*lambda2*t^(kappa-1)}
    S3=function(t){exp(-lambda3*t^kappa)}
    h3=function(t){kappa*lambda3*t^(kappa-1)}  
  }  
  
  else if(K==4){ 
    S1=function(t){exp(-lambda1*t^kappa)}
    h1=function(t){kappa*lambda1*t^(kappa-1)}
    S2=function(t){exp(-lambda2*t^kappa)}
    h2=function(t){kappa*lambda2*t^(kappa-1)}
    S3=function(t){exp(-lambda3*t^kappa)}
    h3=function(t){kappa*lambda3*t^(kappa-1)}
    S4=function(t){exp(-lambda4*t^kappa)}
    h4=function(t){kappa*lambda4*t^(kappa-1)}
  }
  
  
  else if(K==5){ 
    S1=function(t){exp(-lambda1*t^kappa)}
    h1=function(t){kappa*lambda1*t^(kappa-1)}
    S2=function(t){exp(-lambda2*t^kappa)}
    h2=function(t){kappa*lambda2*t^(kappa-1)}
    S3=function(t){exp(-lambda3*t^kappa)}
    h3=function(t){kappa*lambda3*t^(kappa-1)}
    S4=function(t){exp(-lambda4*t^kappa)}
    h4=function(t){kappa*lambda4*t^(kappa-1)}
    S5=function(t){exp(-lambda5*t^kappa)}
    h5=function(t){kappa*lambda5*t^(kappa-1)}
  }
  
  
  
  G1=function(t){1-punif(t, tf, tau)}
  G2=function(t){exp(-eta*t)}
  
  omega0=1/(K+1)
  
  if (K==2) {omega1=omega2=1/(K+1)}
  
  else if(K==3){omega1=omega2=omega3=1/(K+1)}
  
  else if(K==4){omega1=omega2=omega3=omega4=1/(K+1)}
  
  else if(K==5){omega1=omega2=omega3=omega4=omega5=1/(K+1)}
  
  
  if (K==2){
    f1=function(t){S1(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h1(t))/(omega1*S1(t)+omega0*S0(t))}
    f2=function(t){S2(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h2(t))/(omega2*S2(t)+omega0*S0(t))}
  }
  
  else if(K==3){
    f1=function(t){S1(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h1(t))/(omega1*S1(t)+omega0*S0(t))}
    f2=function(t){S2(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h2(t))/(omega2*S2(t)+omega0*S0(t))}
    f3=function(t){S3(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h3(t))/(omega3*S3(t)+omega0*S0(t))}
  }
  
  else if (K==4){
    
    f1=function(t){S1(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h1(t))/(omega1*S1(t)+omega0*S0(t))}
    f2=function(t){S2(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h2(t))/(omega2*S2(t)+omega0*S0(t))}
    f3=function(t){S3(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h3(t))/(omega3*S3(t)+omega0*S0(t))}
    f4=function(t){S4(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h4(t))/(omega4*S4(t)+omega0*S0(t))}
  }
  
  else if (K==5){
    
    f1=function(t){S1(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h1(t))/(omega1*S1(t)+omega0*S0(t))}
    f2=function(t){S2(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h2(t))/(omega2*S2(t)+omega0*S0(t))}
    f3=function(t){S3(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h3(t))/(omega3*S3(t)+omega0*S0(t))}
    f4=function(t){S4(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h4(t))/(omega4*S4(t)+omega0*S0(t))}
    f5=function(t){S5(t)*S0(t)*G1(t)*G2(t)*(h0(t)-h5(t))/(omega5*S5(t)+omega0*S0(t))}
  }
  
  
  
  
  
  if (K==2) {
    V1=function(t){(omega1*S1(t)*h0(t)+omega0*S0(t)*h1(t))*S1(t)*S0(t)*G1(t)*G2(t)/(omega1*S1(t)+omega0*S0(t))^2}
    V2=function(t){(omega2*S2(t)*h0(t)+omega0*S0(t)*h2(t))*S2(t)*S0(t)*G1(t)*G2(t)/(omega2*S2(t)+omega0*S0(t))^2}
  }
  
  else if(K==3) {
    V1=function(t){(omega1*S1(t)*h0(t)+omega0*S0(t)*h1(t))*S1(t)*S0(t)*G1(t)*G2(t)/(omega1*S1(t)+omega0*S0(t))^2}
    V2=function(t){(omega2*S2(t)*h0(t)+omega0*S0(t)*h2(t))*S2(t)*S0(t)*G1(t)*G2(t)/(omega2*S2(t)+omega0*S0(t))^2}
    V3=function(t){(omega3*S3(t)*h0(t)+omega0*S0(t)*h3(t))*S3(t)*S0(t)*G1(t)*G2(t)/(omega3*S3(t)+omega0*S0(t))^2}  
  }
  
  else if(K==4) {
    
    V1=function(t){(omega1*S1(t)*h0(t)+omega0*S0(t)*h1(t))*S1(t)*S0(t)*G1(t)*G2(t)/(omega1*S1(t)+omega0*S0(t))^2}
    V2=function(t){(omega2*S2(t)*h0(t)+omega0*S0(t)*h2(t))*S2(t)*S0(t)*G1(t)*G2(t)/(omega2*S2(t)+omega0*S0(t))^2}
    V3=function(t){(omega3*S3(t)*h0(t)+omega0*S0(t)*h3(t))*S3(t)*S0(t)*G1(t)*G2(t)/(omega3*S3(t)+omega0*S0(t))^2}
    V4=function(t){(omega4*S4(t)*h0(t)+omega0*S0(t)*h4(t))*S4(t)*S0(t)*G1(t)*G2(t)/(omega4*S4(t)+omega0*S0(t))^2}  
  }
  
  else if(K==5) {
    
    V1=function(t){(omega1*S1(t)*h0(t)+omega0*S0(t)*h1(t))*S1(t)*S0(t)*G1(t)*G2(t)/(omega1*S1(t)+omega0*S0(t))^2}
    V2=function(t){(omega2*S2(t)*h0(t)+omega0*S0(t)*h2(t))*S2(t)*S0(t)*G1(t)*G2(t)/(omega2*S2(t)+omega0*S0(t))^2}
    V3=function(t){(omega3*S3(t)*h0(t)+omega0*S0(t)*h3(t))*S3(t)*S0(t)*G1(t)*G2(t)/(omega3*S3(t)+omega0*S0(t))^2}
    V4=function(t){(omega4*S4(t)*h0(t)+omega0*S0(t)*h4(t))*S4(t)*S0(t)*G1(t)*G2(t)/(omega4*S4(t)+omega0*S0(t))^2}  
    V5=function(t){(omega5*S5(t)*h0(t)+omega0*S0(t)*h5(t))*S5(t)*S0(t)*G1(t)*G2(t)/(omega5*S5(t)+omega0*S0(t))^2}  
  }
  
  
  
  
  
  if (K==2) {
    V12=function(t){S1(t)*S2(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega1*S1(t)+omega0*S0(t))*(omega2*S2(t)+omega0*S0(t)))}
  }
  
  else if (K==3) {
    V12=function(t){S1(t)*S2(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega1*S1(t)+omega0*S0(t))*(omega2*S2(t)+omega0*S0(t)))}
    V13=function(t){S1(t)*S3(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega1*S1(t)+omega0*S0(t))*(omega3*S3(t)+omega0*S0(t)))}
    V23=function(t){S2(t)*S3(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega2*S2(t)+omega0*S0(t))*(omega2*S3(t)+omega0*S0(t)))}
  }
  
  
  else if(K==4) {
    V12=function(t){S1(t)*S2(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega1*S1(t)+omega0*S0(t))*(omega2*S2(t)+omega0*S0(t)))}
    
    V13=function(t){S1(t)*S3(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega1*S1(t)+omega0*S0(t))*(omega3*S3(t)+omega0*S0(t)))}
    
    V14=function(t){S1(t)*S4(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega1*S1(t)+omega0*S0(t))*(omega4*S4(t)+omega0*S0(t)))}
    
    V23=function(t){S2(t)*S3(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega2*S2(t)+omega0*S0(t))*(omega3*S3(t)+omega0*S0(t)))}
    
    V24=function(t){S2(t)*S4(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega2*S2(t)+omega0*S0(t))*(omega4*S4(t)+omega0*S0(t)))}
    
    V34=function(t){S3(t)*S4(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega3*S3(t)+omega0*S0(t))*(omega4*S4(t)+omega0*S0(t)))}
  }
  
  
  else if(K==5) {
    V12=function(t){S1(t)*S2(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega1*S1(t)+omega0*S0(t))*(omega2*S2(t)+omega0*S0(t)))}
    
    V13=function(t){S1(t)*S3(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega1*S1(t)+omega0*S0(t))*(omega3*S3(t)+omega0*S0(t)))}
    
    V14=function(t){S1(t)*S4(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega1*S1(t)+omega0*S0(t))*(omega4*S4(t)+omega0*S0(t)))}
    
    V15=function(t){S1(t)*S5(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega1*S1(t)+omega0*S0(t))*(omega5*S5(t)+omega0*S0(t)))}
    
    V23=function(t){S2(t)*S3(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega2*S2(t)+omega0*S0(t))*(omega3*S3(t)+omega0*S0(t)))}
    
    V24=function(t){S2(t)*S4(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega2*S2(t)+omega0*S0(t))*(omega4*S4(t)+omega0*S0(t)))}
    
    V25=function(t){S2(t)*S5(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega2*S2(t)+omega0*S0(t))*(omega5*S5(t)+omega0*S0(t)))}
    
    V34=function(t){S3(t)*S4(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega3*S3(t)+omega0*S0(t))*(omega4*S4(t)+omega0*S0(t)))}
    
    V35=function(t){S3(t)*S5(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega3*S3(t)+omega0*S0(t))*(omega5*S5(t)+omega0*S0(t)))}
    
    V45=function(t){S4(t)*S5(t)*h0(t)*S0(t)*G1(t)*G2(t)/((omega4*S4(t)+omega0*S0(t))*(omega5*S5(t)+omega0*S0(t)))}  
    
  }
  
  
  
  
  
  
  
  
  g0=function(t){S0(t)*h0(t)*G1(t)*G2(t)}
  
  
  if (K==2){
    g1=function(t){S1(t)*h1(t)*G1(t)*G2(t)}
    g2=function(t){S2(t)*h2(t)*G1(t)*G2(t)}  
  }
  
  else if (K==3){
    g1=function(t){S1(t)*h1(t)*G1(t)*G2(t)}
    g2=function(t){S2(t)*h2(t)*G1(t)*G2(t)}
    g3=function(t){S3(t)*h3(t)*G1(t)*G2(t)}  
  }
  
  else if (K==4){
    g1=function(t){S1(t)*h1(t)*G1(t)*G2(t)}
    g2=function(t){S2(t)*h2(t)*G1(t)*G2(t)}
    g3=function(t){S3(t)*h3(t)*G1(t)*G2(t)}
    g4=function(t){S4(t)*h4(t)*G1(t)*G2(t)}  
  }
  
  else if (K==5){
    g1=function(t){S1(t)*h1(t)*G1(t)*G2(t)}
    g2=function(t){S2(t)*h2(t)*G1(t)*G2(t)}
    g3=function(t){S3(t)*h3(t)*G1(t)*G2(t)}
    g4=function(t){S4(t)*h4(t)*G1(t)*G2(t)}
    g5=function(t){S5(t)*h5(t)*G1(t)*G2(t)}
  }
  
  
  
  
  
  if (K==2) {P=(integrate(g0, 0,tau)$value+integrate(g1, 0,tau)$value+integrate(g2, 0,tau)$value)/(K+1) }
  
  else if (K==3) {P=(integrate(g0, 0,tau)$value+integrate(g1, 0,tau)$value+integrate(g2, 0,tau)$value+
                       integrate(g3, 0,tau)$value)/(K+1) }
  
  else if (K==4){P=(integrate(g0, 0,tau)$value+integrate(g1, 0,tau)$value+integrate(g2, 0,tau)$value+
                      integrate(g3, 0,tau)$value+integrate(g4, 0,tau)$value)/(K+1)}
  
  else if (K==5){P=(integrate(g0, 0,tau)$value+integrate(g1, 0,tau)$value+integrate(g2, 0,tau)$value+
                      integrate(g3, 0,tau)$value+integrate(g4, 0,tau)$value+integrate(g5, 0,tau)$value)/(K+1)}
  
  
  if (K==2){
    mu1=omega1*omega0*integrate(f1, 0, tau)$value
    mu2=omega2*omega0*integrate(f2, 0, tau)$value
  }
  
  else if(K==3){
    
    mu1=omega1*omega0*integrate(f1, 0, tau)$value
    mu2=omega2*omega0*integrate(f2, 0, tau)$value
    mu3=omega3*omega0*integrate(f3, 0, tau)$value
  }
  
  else if(K==4) {
    
    mu1=omega1*omega0*integrate(f1, 0, tau)$value
    mu2=omega2*omega0*integrate(f2, 0, tau)$value
    mu3=omega3*omega0*integrate(f3, 0, tau)$value
    mu4=omega4*omega0*integrate(f4, 0, tau)$value
  }
  
  
  else if(K==5) {
    
    mu1=omega1*omega0*integrate(f1, 0, tau)$value
    mu2=omega2*omega0*integrate(f2, 0, tau)$value
    mu3=omega3*omega0*integrate(f3, 0, tau)$value
    mu4=omega4*omega0*integrate(f4, 0, tau)$value
    mu5=omega5*omega0*integrate(f5, 0, tau)$value
  }
  
  
  if (K==2) {
    sig11=omega1*omega0*integrate(V1, 0,tau)$value
    sig22=omega2*omega0*integrate(V2, 0,tau)$value
  }
  
  else if(K==3){
    sig11=omega1*omega0*integrate(V1, 0,tau)$value
    sig22=omega2*omega0*integrate(V2, 0,tau)$value
    sig33=omega3*omega0*integrate(V3, 0,tau)$value  
  }
  
  else if(K==4) {
    sig11=omega1*omega0*integrate(V1, 0,tau)$value
    sig22=omega2*omega0*integrate(V2, 0,tau)$value
    sig33=omega3*omega0*integrate(V3, 0,tau)$value
    sig44=omega4*omega0*integrate(V4, 0,tau)$value  
  }
  
  else if(K==5) {
    sig11=omega1*omega0*integrate(V1, 0,tau)$value
    sig22=omega2*omega0*integrate(V2, 0,tau)$value
    sig33=omega3*omega0*integrate(V3, 0,tau)$value
    sig44=omega4*omega0*integrate(V4, 0,tau)$value
    sig55=omega5*omega0*integrate(V5, 0,tau)$value
  }
  
  
  
  if (K==2) {
    sig12=omega1*omega2*omega0*integrate(V12, 0,tau)$value
  }  
  
  else if (K==3) {
    sig12=omega1*omega2*omega0*integrate(V12, 0,tau)$value
    sig13=omega1*omega3*omega0*integrate(V13, 0,tau)$value
    sig23=omega2*omega3*omega0*integrate(V23, 0,tau)$value
  }
  
  else if (K==4) {
    sig12=omega1*omega2*omega0*integrate(V12, 0,tau)$value
    sig13=omega1*omega3*omega0*integrate(V13, 0,tau)$value
    sig14=omega1*omega4*omega0*integrate(V14, 0,tau)$value
    sig23=omega2*omega3*omega0*integrate(V23, 0,tau)$value
    sig24=omega2*omega4*omega0*integrate(V24, 0,tau)$value
    sig34=omega3*omega4*omega0*integrate(V34, 0,tau)$value
  }
  
  else if (K==5) {
    sig12=omega1*omega2*omega0*integrate(V12, 0,tau)$value
    sig13=omega1*omega3*omega0*integrate(V13, 0,tau)$value
    sig14=omega1*omega4*omega0*integrate(V14, 0,tau)$value
    sig15=omega1*omega5*omega0*integrate(V15, 0,tau)$value
    sig23=omega2*omega3*omega0*integrate(V23, 0,tau)$value
    sig24=omega2*omega4*omega0*integrate(V24, 0,tau)$value
    sig25=omega2*omega5*omega0*integrate(V25, 0,tau)$value
    sig34=omega3*omega4*omega0*integrate(V34, 0,tau)$value
    sig35=omega3*omega5*omega0*integrate(V35, 0,tau)$value
    sig45=omega4*omega5*omega0*integrate(V45, 0,tau)$value
  }
  
  
  
  
  
  
  if (K==2){
    rho12=sig12/(sqrt(sig11)*sqrt(sig22))
    rho<-rho12
  }
  
  else if(K==3) {
    rho12=sig12/(sqrt(sig11)*sqrt(sig22))
    rho13=sig13/(sqrt(sig11)*sqrt(sig33))
    rho23=sig23/(sqrt(sig22)*sqrt(sig33))
    rho=c(rho12,rho13,rho23)
  }
  
  else if(K==4) {
    rho12=sig12/(sqrt(sig11)*sqrt(sig22))
    rho13=sig13/(sqrt(sig11)*sqrt(sig33))
    rho14=sig14/(sqrt(sig11)*sqrt(sig44))
    rho23=sig23/(sqrt(sig22)*sqrt(sig33))
    rho24=sig24/(sqrt(sig22)*sqrt(sig44))
    rho34=sig34/(sqrt(sig33)*sqrt(sig44))
    rho=c(rho12,rho13,rho14,rho23,rho24,rho34)
  }
  
  
  else{
    rho12=sig12/(sqrt(sig11)*sqrt(sig22))
    rho13=sig13/(sqrt(sig11)*sqrt(sig33))
    rho14=sig14/(sqrt(sig11)*sqrt(sig44))
    rho15=sig15/(sqrt(sig11)*sqrt(sig55))
    rho23=sig23/(sqrt(sig22)*sqrt(sig33))
    rho24=sig24/(sqrt(sig22)*sqrt(sig44))
    rho25=sig25/(sqrt(sig22)*sqrt(sig55))
    rho34=sig34/(sqrt(sig33)*sqrt(sig44))
    rho35=sig35/(sqrt(sig33)*sqrt(sig55))
    rho45=sig45/(sqrt(sig44)*sqrt(sig55))
    rho=c(rho12,rho13,rho14,rho15,rho23,rho24,rho25,rho34,rho35,rho45)
  }
  
  
  
  
  
  
  
  # Sigma=matrix(c(1,rho12,rho13,rho14, rho12, 1, rho23, rho24, rho13, rho23, 1,
  #                rho34, rho14, rho24, rho34, 1), K,K)
  
  Sigma<-diag(1,K,K)
  
  Sigma[lower.tri(Sigma)]<-rho
  
  Sigma<-Sigma+ t(Sigma)+diag(-1,K,K)
  
  
  # A=matrix(c(1, -1, 0, 0, 1, 0, -1, 0, 1, 0, 0, -1, 1, 0, 0, 0), K, K, byrow=T)
  
  A<-diag(-1,K-1,K-1)
  A<-rbind(A,rep(0,K-1))
  A<-cbind(rep(1,K),A)
  
  
  B=A%*%Sigma%*%t(A)
  
  if (K==2) {
    mu1.bar=mu1/sqrt(sig11)
    mu2.bar=mu2/sqrt(sig22)
    mu=c(mu1.bar, mu2.bar)
  }
  
  else if(K==3) {
    mu1.bar=mu1/sqrt(sig11)
    mu2.bar=mu2/sqrt(sig22)
    mu3.bar=mu3/sqrt(sig33)  
    mu=c(mu1.bar, mu2.bar, mu3.bar)
  }
  
  else if(K==4) {
    mu1.bar=mu1/sqrt(sig11)
    mu2.bar=mu2/sqrt(sig22)
    mu3.bar=mu3/sqrt(sig33)
    mu4.bar=mu4/sqrt(sig44)
    mu=c(mu1.bar, mu2.bar, mu3.bar, mu4.bar)  
  }
  
  else if(K==5) {
    mu1.bar=mu1/sqrt(sig11)
    mu2.bar=mu2/sqrt(sig22)
    mu3.bar=mu3/sqrt(sig33)
    mu4.bar=mu4/sqrt(sig44)
    mu5.bar=mu5/sqrt(sig55)
    mu=c(mu1.bar, mu2.bar, mu3.bar, mu4.bar,mu5.bar)  
  }
  
  
  
  
  root=function(n){
    b=as.numeric(sqrt(n)*(A%*%mu))
    int=pmvnorm(lower = c(rep(0,K-1),c), upper=rep(Inf,K), mean = b, sigma = B)[1]
    1-beta-as.double(int)
  }
  n=uniroot(root, lower=1, upper=99999)$root
  e=ceiling(n*P)
  return(c(ceiling(e/(K+1)), ceiling(n/(K+1))))
}


