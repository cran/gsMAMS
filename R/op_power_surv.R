 


#' @title Provides operating characteristics of group sequential MAMS trial for survival outcome
#' @description Computes power and other characteristics for group-sequential MAMS trial for survival outcome.
#' @param m0 Median survival time of control group.
#' @param alpha Type I error.
#' @param beta Type II error.
#' @param K Number of treatment arms.
#' @param frac Vector of fractions for information time at each look.
#' @param HR0 Hazard ratio of ineffective treatment group vs control.
#' @param HR1 Hazard ratio of effective treatment group vs control.
#' @param nsim Number of simulations.
#' @param ta Accrual time.
#' @param tf Follow-up time.
#' @param kappa Shape parameter (kappa=1 for exponential distribution).
#' @param eta  Rate of loss to follow-up.
#' @param seed Random seed number.
#' @return A list of power, stage-wise probability of success, stopping probability, probability of futility, average number of events happened per arm, average duration of trial.
#' @import stats
#' @importFrom survival survdiff Surv
#' @examples
#' op_power_surv(20,0.05,0.1,4,c(1/2,1),1,0.74,12,40,20,1,0,12)
#' @export




op_power_surv<-function(m0,alpha, beta,K,frac,HR0, HR1,nsim,ta,tf,kappa,eta,seed){

  if (K<=1) {stop("K should be greater than 1.")}
  if(length(frac)==1){stop("The length of frac should be greater than 1.")}
  j<-length(frac)

  bound<-SCPRT(alpha = alpha,K=K,frac = frac)


  hr0=HR0           #exp(-delta0)
  hr1=HR1           #exp(-delta1)

  lambda0=log(2)/m0^kappa

  lambda<-numeric()

  lambda[1]=lambda0*hr1

  for (i in 2:K){
    lambda[i]<-lambda0*hr0
  }

  scale<-vector()
  scale[1]=1/lambda0^(1/kappa)


  for(i in 2:(K+1) ){
    scale[i]<-1/lambda[i-1]^(1/kappa)
  }


  s2<-numeric(length=j)
  

   

  n=Size_surv(m0=m0, HR0=HR0, HR1=HR1, ta=ta, tf=tf, K=K,beta = beta,alpha=alpha ,kappa=kappa,eta = eta,frac=frac)[2]
  d=Size_surv(m0=m0, HR0=HR0, HR1=HR1, ta=ta, tf=tf, K=K,beta = beta,alpha=alpha,kappa=kappa,eta = eta,frac=frac)[1]

   frac<-frac
  d1<-numeric(length = j)
  for(i in 1:j) {
    d1[i]<-ceiling(d*frac[i])
  }



  a<-data.frame()
  for(i in 1:K){
    a<-rbind.data.frame(a,bound$lshape)
  }

  b<-data.frame()
  for(i in 1:K){
    b<-rbind.data.frame(b,bound$ushape)
  }





  tstar=1/j



  tau=ta+tf
  sp<-numeric()
  pf<-numeric()
  smn<-numeric()
  smd<-numeric()
  dur<-numeric()

  #print(a)
  #print(b)
  set.seed(seed)
  for (e in 1:nsim){

    Q=rep(0,j)
    #ASN=0
    datagen=NULL
    for (i in 1:(K+1)){
      w=rweibull(n, kappa, scale[i])
      u=runif(n, 0, ta) ## generate accrual time
      if (eta !=0)
        g=rexp(n,rate=eta)
      if(eta==0)
        g=tau-u
      x=pmax(0, pmin(w,tau-u,g)) ## observed survival time
      cens = as.numeric(w<pmin((tau-u),g) ) ## censoring indicator
      group=rep((i-1),n)
      datagen[[i]]=cbind(x,cens,group,u)
    }


    data=NULL
    for (i in 1:K){
      dat=rbind.data.frame(datagen[[1]],datagen[[i+1]])
      colnames(dat)=c("time","event","group","accrualtime")
      #dat=data.frame(dat)
      dat$calendar_T=dat$accrualtime+dat$time
      dat=dat[order(dat$calendar_T),]
      dat$cum_event=cumsum(dat$event)
      data[[i]]=dat
      #print(head(data[[i]]))
    }


    loc=matrix(NA,K,j)
    z=matrix(NA,K,j)
    var_result=matrix(NA,K,j)
    N_result=matrix(NA,K+1,j)
    d_result=matrix(NA,K+1,j)
    calendarT_look=matrix(NA,K,j)
    stagensizeT=matrix(NA,K,j)
    stagensizeC=matrix(NA,K,j)
    stagedsizeT=matrix(NA,K,j)
    stagedsizeC=matrix(NA,K,j)


    for (h in (1:j)){
      for (k in (1:K)){
        if (h<j){
           
          dlook=ceiling(tstar*h*2*d)
          loc=min(which(data[[k]]$cum_event==dlook))
          calendarT_look[k,h]=data[[k]]$calendar_T[loc];
          data_compare_part1<-data[[k]][1:loc,]
          data_compare_part2<-data[[k]][(loc+1):(2*n), ]
          data_compare_part2=data_compare_part2[
            data_compare_part2$accrualtime<data_compare_part1$calendar_T[loc],]
          data_compare_part2$event=0
          data_compare_part2$time=data_compare_part1$calendar_T[loc]-data_compare_part2$accrualtime
          #data_compare_part2$accrualtime
          data_compare=rbind(data_compare_part1,data_compare_part2)
          stagensizeT[k,h]=nrow(data_compare[data_compare$group!=0,])
          stagensizeC[k,h]=nrow(data_compare[data_compare$group==0,])

          stagedsizeT[k,h]=nrow(data_compare[data_compare$group!=0 &data_compare$event==1,])
          stagedsizeC[k,h]=nrow(data_compare[data_compare$group==0 &data_compare$event==1,])
        }

        #if (e==5 & h==1){browser()}
        if (h==j){
          calendarT_look[k,h]=tau
          data_compare<-data[[k]]
          stagensizeT[k,h]=nrow(data_compare[data_compare$group!=0,])
          stagensizeC[k,h]=nrow(data_compare[data_compare$group==0,])
          stagedsizeT[k,h]=nrow(data_compare[data_compare$group!=0 & data_compare$event==1,])
          stagedsizeC[k,h]=nrow(data_compare[data_compare$group==0 & data_compare$event==1,])
        }
        
        temp<-survdiff(Surv(time, event)~group,data = data_compare)
        z[k,h]=sign(temp$obs[1]-temp$exp[1])*sqrt(temp$chisq)
        #logranktest(data_compare$time, data_compare$event,data_compare$group)
      }
    }


  m1<-as.numeric()
  for (i in 1:K){
    if (i==1) {m1[i]<-z[1,1]>b[i,1]}

    else {m1[i]<-z[1,1]> z[i,1]}
  }
  pk<-prod(m1)
  if (pk>=1){s2[1]=s2[1]+1}


  w=K
  g<-data.frame(matrix(ncol =w-1, nrow = 0))

  for(q in 2:(length(g)+1)) {
    #j<-3
    p<-numeric(length =(j-1)*3)
    k<-seq(1,25,by=3)[1:(j-1)]
    sk<-list(a=1,b=c(2,3))
    #q<-as.numeric()

    for(i in 2:(j)) {

      if (i==2){
        #browser()
        p[k[i-1]]<-z[q,(i-1)]<a[q,(i-1)]
        p[k[i-1]+1]<-z[q,(i-1)]>a[q,(i-1)] & z[q,(i-1)]<b[q,(i-1)]
        p[k[i-1]+2]<-z[1,(i)]>z[q,(i)]

        g[i,q-1]<-0
        for (o in 1:length(sk)){
          g[i,q-1]<-g[i,q-1]+prod(p[sk[[o]]])
        }
        next
      }

      tk<-sk[[length(sk)]]

      sk[[length(sk)]][length(sk[[length(sk)]])]<-sk[[length(sk)]][length(sk[[length(sk)]])]+1

      sk[[length(sk)+1]]<-tk
      sk[[length(sk)]][length(sk[[length(sk)]])]<-sk[[length(sk)]][length(sk[[length(sk)]])]+2
      sk[[length(sk)]][length(sk[[length(sk)]])+1]<-sk[[length(sk)]][length(sk[[length(sk)]])]+1

      p[k[i-1]]<-z[q,(i-1)]<a[q,(i-1)]
      p[k[i-1]+1]<-z[q,(i-1)]>a[q,(i-1)] & z[q,(i-1)]<b[q,(i-1)]
      p[k[i-1]+2]<-z[1,i]>z[q,i]

      g[i,q-1]<-0
      for (o in 1:length(sk)){
        g[i,q-1]<-g[i,q-1]+prod(p[sk[[o]]])
      }

    }
  }


  sj<-data.frame(matrix(ncol=1, nrow = 0))

  for(q in 1:1) {
    #j<-3
    p<-numeric(length =(j-1)*2)
    k<-seq(1,20,by=2)[1:(j-1)]
    sk<-list(a=c(1,2))
    #q<-as.numeric()

    for(i in 2:(j)) {

      if (i==2){
        #browser()
        # p[k[i-1]]<-z[q,(i-1)]<a[q,(i-1)]
        p[k[i-1]]<-z[q,(i-1)]>a[q,(i-1)] & z[q,(i-1)]<b[q,(i-1)]
        p[k[i-1]+1]<-z[1,(i)]>b[q,(i)]

        sj[i,1]<-prod(p[sk[[length(sk)]]])
        next
      }

      #tk<-sk[[length(sk)]]

      sk[[length(sk)+1]]<-sk[[length(sk)]]
      sk[[length(sk)]][length(sk[[length(sk)]])]<-sk[[length(sk)]][length(sk[[length(sk)]])]+1
      sk[[length(sk)]][length(sk[[length(sk)]])+1]<-sk[[length(sk)]][length(sk[[length(sk)]])]+1

      p[k[i-1]]<-z[q,(i-1)]>a[q,(i-1)] & z[q,(i-1)]<b[q,(i-1)]
      p[k[i-1]+1]<-z[1,(i)]>b[q,(i)]


      sj[i,1]<-prod(p[sk[[length(sk)]]])

    }

  }






  ff<-cbind.data.frame(g,sj)

  pl<-numeric(length =j)

  for (y in 2:j ){
    if(prod(ff[y,])>=1) {
      s2[y]<-s2[y]+1}

  }

  #################################
 

  ############################
  stopprob=rep(0,j)
  probfut=rep(0,j)
  stop=rep(0,j)
  stopF=matrix(NA,K,j)
  stopE=matrix(NA,K,j)
  Fflag=rep(0,j)
  Eflag=rep(0,j)

  for (k in 1:K){
    for (h in 1:j){
      stopF[k,h]=z[k,h]<=a[k,h]
      stopE[k,h]=z[k,h]>=b[k,h]
    }
  }


  for (k in 1:K){
    if (any(stopF[k,]==TRUE)) {
      s<-min(which(stopF[k,]==TRUE))
      stopF[k,s:j]<-TRUE
    }
  }


  for (h in 1:j) {
    if (sum(stopF[,h]==TRUE)==K) {Fflag[h]=1}
    if (sum(stopE[,h]==TRUE)>=1) {Eflag[h]=1}
    stop[h]=Eflag[h]+Fflag[h] # stop is rep (0,J) originally, if we stop at stage 3, we will have (0,0,1,0), thus stop records at which stage we stop for current simulation
  }



  stopstage=min(which(stop>=1)) # the stage we should stop, no matter it is due to efficacy or futility
  stopFstage=ifelse(all(Fflag==0),0,min(which(Fflag==1))) # the stage we  stop due to futility

  if(stopFstage==stopstage) {probfut[stopFstage]=1}

  stopprob[stopstage]=1

  samplesize=rep(0,K+1)
  samplesized=rep(0,K+1)

  for (k in 1:K){

    stopFstage1=ifelse(any(stopF[k,]==TRUE),min(which(stopF[k,]==TRUE)),Inf)

    stopstagefork=min(stopstage,stopFstage1)
    samplesize[k]=stagensizeT[k,stopstagefork]
    samplesized[k]=stagedsizeT[k,stopstagefork]
  }
  samplesize[K+1]=max(stagensizeC[,stopstage])
  samplesized[K+1]=max(stagedsizeC[,stopstage])



  duration=max(calendarT_look[,stopstage])



   
  sp<-rbind.data.frame(stopprob,sp)
  pf<-rbind.data.frame(probfut,pf)
  smn[e]<-sum(samplesize)/(K+1)
  smd[e]<-sum(samplesized)/(K+1)
  dur[e]<-duration



  }
   
  power=round(sum(s2)/nsim,3)
   
   
  s2<-rbind.data.frame(s2/nsim)
  names(s2)<-paste0("look", 1:j)
  names(sp)<-paste0("look", 1:j)
  names(pf)<-paste0("look", 1:j)


  p<-list("Power"=power,"Stagewise Power"=colMeans(s2),"Stopping probability under alternative"=colMeans(sp), "Probability of futility under alternative"=colMeans(pf), "Average number of events happened per arm under alternative"=mean(smd), "Average duration of trial(months)"=mean(dur))
  return(p)

}



