
##Helper function
score<-function(r0,rk,n,data0,datak){
  
  a<-length(data0) #number of categories
  temp=0
  temp2=0
  data0=as.vector(data0)
  if (a>2){
    for (u in 1:a){
      if (u==1)
        temp=temp+datak[u]*(sum(data0[u+1:a],na.rm=T))else if (u==a)
          temp=temp+datak[u]*(0-sum(data0[1:u-1],na.rm=T))else
            temp=temp+datak[u]*(sum(data0[u+1:a],na.rm=T)-sum(data0[1:u-1],na.rm=T))
          
          temp2=temp2+((datak[u]+data0[u])/(sum(datak)+sum(data0)))^3
    }
  }
  if (a==2){
    temp=datak[1]*data0[2]-datak[2]*data0[1]
    temp2=((datak[1]+data0[1])/(sum(datak)+sum(data0)))^3+((datak[2]+data0[2])/(sum(datak)+sum(data0)))^3
  }
  
  
  
  
  sk=1/((rk+r0)*n)*temp
  vk=((rk*r0*n)/(3*(rk+r0)))*(1-temp2)
  
  Z=sk/sqrt(vk)
  return(Z)
  
}
  
 
 
#' @title Provides operating characteristics of group sequential MAMS trial for ordinal outcome under null hypothesis
#' @description Computes FWER and other characteristics for group-sequential MAMS trial for ordinal outcome.
#' @param alpha Type I error.
#' @param beta Type II error.
#' @param K Number of treatment arms.
#' @param frac Vector of fractions for information time at each look.
#' @param or0 Odds ratio of ineffective treatment group vs control.
#' @param or Odds ratio of effective treatment group vs control.
#' @param nsim Number of simulations.
#' @param prob Probability of ordinal outcomes in control group.
#' @param seed Random seed number.
#' @return A list of FWER, stage-wise type I error, average sample size used per arm, stopping probability, probability of futility.
#' @import stats
#' @examples
#' op_fwer_ord(0.05,0.1,4,c(0.5,1),1.32,3.06,15,c(0.075, 0.182, 0.319, 0.243, 0.015, 0.166),13)  
#' @export

op_fwer_ord<-function(alpha,beta,K,frac,or0,or,nsim,prob,seed){
  
  if (K<=1) {stop("K should be greater than 1.")}
  if(length(frac)==1){stop("The length of frac should be greater than 1.")}
  j<-length(frac)
  bound<-SCPRT(alpha = alpha,K=K,frac = frac)
  
  l<-Size_ord(prob=prob, or=or, or0=or0, alpha=alpha, beta=beta, K=K)
  #l<-84
  
  prob1<-prob 
  prob2<-prob
  
  
  sp<-numeric()
  pf<-numeric()
  s2<-numeric(length=j)
  s1<-0
  asn<-0
  frac<-frac
  n<-numeric(length = j)
  for(i in 1:j) {
    n[i]<-ceiling(l*frac[i])
  }
  

  #a<-cbind.data.frame(c(-0.630 ,-0.630 ,-0.630 ,-0.630 ),c(0.437,0.437,0.437,0.437),c(2.161,2.161,2.161,2.161))  
  #b<-cbind.data.frame(c(3.126,3.126,3.126,3.126),c(3.092,3.092,3.092,3.092),c(2.161,2.161,2.161,2.161))
    
  a<-data.frame()
  for(i in 1:K){
    a<-rbind.data.frame(a,bound$lshape)
  }

  b<-data.frame()
  for(i in 1:K){
    b<-rbind.data.frame(b,bound$ushape)
  }


  #print(a)
  #print(b)
  set.seed(seed)
  for (e in 1:nsim){
    
    z<-matrix(NA,K,j)
    Q=rep(0,j)
    #ASN=0
    datagen=NULL
    group=matrix(NA, 2,j)
    for (i in 1:(K+1)){
      if (i==1){
        group<-rmultinom(j, size = n[1], prob = prob) #control group
        mySum = t(apply(group, 1, cumsum))
      }else if (i==2){
        group<-rmultinom(j ,size = n[1], prob = prob1)# group  1
        mySum = t(apply(group, 1, cumsum))
      }else{
        group<-rmultinom(j, size = n[1], prob = prob2)# group  2 to K
        mySum = t(apply(group, 1, cumsum))
      }
      datagen[[i]]=mySum
    }
    
    for (v in (1:K)){
      for (h in (1:j))
        #________________For group k vs. control group
        if (j==1){
          z[v,h]=score(h,h,n[1],datagen[[1]],datagen[[v+1]]) 
        }else{
          z[v,h]=score(h,h,n[1],datagen[[1]][,h],datagen[[v+1]][,h])
        }
      # print(c(j,k))
    }
    
    w=K
    g<-data.frame(matrix(ncol=w, nrow = 0))
    mp<-data.frame(matrix(ncol=w, nrow = 0))
    for(q in 1:(length(g))) {
      #j<-3
      p<-numeric(length =(j-1)*3)
      v<-numeric(length =(j-1)*3)
      k<-seq(1,25,by=3)[1:(j-1)]
      sk<-list(a=1,b=c(2,3))
      #q<-as.numeric()
      
      for(i in 2:(j)) {
        
        if (i==2){
          #browser()
          p[k[i-1]]<-z[q,(i-1)]<a[q,(i-1)]
          p[k[i-1]+1]<-z[q,(i-1)]>a[q,(i-1)] & z[q,(i-1)]<b[q,(i-1)]
          p[k[i-1]+2]<-z[q,(i)]<b[q,(i)]
          
          v[k[i-1]]<-z[q,(i-1)]>b[q,(i-1)]
          v[k[i-1]+1]<-z[q,(i-1)]>a[q,(i-1)] & z[q,(i-1)]<b[q,(i-1)]
          v[k[i-1]+2]<-z[q,(i)]>b[q,(i)]
          
          if(i==j){
            g[1,q]<-0
            for (o in 1:length(sk)){
              g[1,q]<-g[1,q]+prod(p[sk[[o]]])
              mp[o,q]<-prod(v[sk[[o]]])
            }
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
        p[k[i-1]+2]<-z[q,(i)]<b[q,(i)]
        
        v[k[i-1]]<-z[q,(i-1)]>b[q,(i-1)]
        v[k[i-1]+1]<-z[q,(i-1)]>a[q,(i-1)] & z[q,(i-1)]<b[q,(i-1)]
        v[k[i-1]+2]<-z[q,(i)]>b[q,(i)]
        
        if(i==j){
          g[1,q]<-0
          for (o in 1:length(sk)){
            g[1,q]<-g[1,q]+prod(p[sk[[o]]])
            mp[o,q]<-prod(v[sk[[o]]])
          }
        }
      }
    }
    
    if(prod(g[1,])>=1) {s1<-s1+1} 
    
    for (y in 1:j ){
      s2[y]<-s2[y]+ ifelse(sum(mp[y,])>=1,1,0)
      
    }
    
     
    
    #################################
    lp<-as.data.frame(z)
    z<-list()
    for (i in 1:length(lp)){
      z[[i]]<-lp[,i]
    }
    
    for (i in 1: (j-1) ){
      
      
      if (i==1) {
        k=numeric(length = j)
        k[i]<-K
        effn<-sum(z[[i]]>=b[1,i])
        fuln<-sum(z[[i]]<=a[1,i])
      }
      
      else{
        arms2=which((z[[i-1]]<=a[1,i-1])==FALSE)
        k[i]<-length(arms2)
        z[[i]]=z[[i]][arms2]
        effn<-sum(z[[i]]>=b[1,i])
        fuln<-sum(z[[i]]<=a[1,i])
      }
      
      if(effn>=1|fuln==k[i]){
        pq=numeric()
        for (s in c(1:i)){
          if(s==i){pq[s]<-(n[s])*(k[s]+1)
          break}
          if(k[s]==k[s+1]){
            pq[s]<-0
          }
          else if(k[s]>k[s+1]) {
            pq[s]<-(k[s]-k[s+1])*n[s]
          }
          
        }
        asn=asn+sum(pq)
        break
      }
      
      else {
        if((i+1)==(j)){
          arms3=which((z[[i]]<=a[1,i])==FALSE)
          k[i+1]<-length(arms3)
          pq=numeric()
          
          for (s in c(1:j)){
            if(s==j){pq[s]<-(n[s])*(k[s]+1)
            break}
            if(k[s]==k[s+1]){
              pq[s]<-0
            }
            else if(k[s]>k[s+1]) {
              pq[s]<-(k[s]-k[s+1])*n[s]
            }
            
          }
          asn<-asn+sum(pq)
        }
        
      }
      
    }
    ############################
    
    z<-lp
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
    
    
    sp<-rbind.data.frame(stopprob,sp)
    pf<-rbind.data.frame(probfut,pf)
    
    
    
    
    
    
  }
  
  
  names(sp)<-paste0("look", 1:j)
  names(pf)<-paste0("look", 1:j)
  
  fwer=round(1-(s1/nsim),3)
  fwg<-rbind.data.frame(s2/nsim)
  names(fwg)<-paste0("look",1:j)
  asn=round(asn/(nsim*(K+1)),3)
  
  p<-list(FWER=fwer,"Stagewise FWER"=colMeans(fwg),"Stopping probability under null"=colMeans(sp),"Probability of futility under null"=colMeans(pf),"Average sample size used per arm under null"=asn)
  return(p)
   
}

