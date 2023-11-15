require(mvtnorm)
require(expm)
library(glasso)
library(magic)
library(rrcov)
setwd("C:/Users/jimoh.ajadi/OneDrive - KFUPM/Desktop/Research/Dr. Ishaq/New Version-April")
source("Phase I estimates.R")

fGCL=function(n,p,h){
 for(l in 1:length(n)){
  b01=1;
  b02=1;
  bi2=1;
  for(i in 1:p){
   b01=(n-i)*b01;
   bi2=(n-i+2)*bi2;
   b02=(n-i)*b02
  }
  b1=b01/(n-1)^p;
  b2=b02*(bi2-b01)/(n-1)^(2*p)
  UCL=(b1+h*sqrt(b2));
  LCL=max((b1-h*sqrt(b2)),0)
 }
 return(data.frame(LCL,UCL))
}

#GPC
fLgCov=function(n,p,mu,MU,sig,SIG,LCL,UCL){
 rl=1
 x <- rmvnorm(n,mu,sig)
 y=t(solve(sqrtm(SIG))%*%t(x-MU))
 S_shrunk=cov(y)
 SVD=svd(S_shrunk)
 logShat=t(SVD$u%*%diag(log(SVD$d)))%*%SVD$v
 diagS=det(logShat)
 while(diagS<UCL & diagS>LCL){
  rl=rl+1
  x <- rmvnorm(n,mu,sig)
  y=t(solve(sqrtm(SIG))%*%t(x-MU))
  S_shrunk=cov(y)
  SVD=svd(S_shrunk)
  logShat=t(SVD$u%*%diag(log(SVD$d)))%*%SVD$v
  diagS=det(logShat)
 }
 return(rl)
}

#PC1
fLCov=function(n,p,mu,MU,sig,SIG,LCL,UCL){
 rl=1
 x <- rmvnorm(n,mu,sig)
 y=t(solve(sqrtm(SIG))%*%t(x-MU))
 S_shrunk=cov(y)
 SVD=svd(S_shrunk)
 logShat=t(SVD$u%*%log(SVD$d))%*%SVD$v
 diagS=sum(logShat)
 while(diagS<UCL & diagS>LCL){
  rl=rl+1
  x <- rmvnorm(n,mu,sig)
  y=t(solve(sqrtm(SIG))%*%t(x-MU))
  S_shrunk=cov(y)
  SVD=svd(S_shrunk)
  logShat=t(SVD$u%*%log(SVD$d))%*%SVD$v
  diagS=sum(logShat)
 }
 return(rl)
}

#GVC

GVCrl = function(n,p,mu,MU,sig,SIG,L)
{
 CL=fGCL(n,p,L)
 UCL=CL$UCL
 LCL=CL$LCL
 x <- rmvnorm(n,mu,sig)
 y=t(solve(sqrtm(SIG))%*%t(x-MU))
 S = cov(y)
 DetS=det(S)
 rl = 1;
 while (DetS > LCL & DetS < UCL)
 {
  rl = rl + 1;
  x <- rmvnorm(n,mu,sig)
  y=t(solve(sqrtm(SIG))%*%t(x-MU))
  S = cov(y)
  DetS=det(S)
 }
 return(rl)
}


#NTCC
NTCCrl = function(n,p,mu,MU,sig,SIG,UCL,LCL)
{
 rl = 1;
 x <- rmvnorm(n,mu,sig)
 y=t(solve(sqrtm(SIG))%*%t(x-MU))
 S = cov(y)
 Trace=sum(diag(S))
 while (Trace > LCL & Trace < UCL)
 {
  rl = rl + 1;
  x <- rmvnorm(n,mu,sig)
  y=t(solve(sqrtm(SIG))%*%t(x-MU))
  S = cov(y)
  Trace=sum(diag(S))
 }
 return(rl)
}

#PLR
fPLR=function(n,p,mu,MU,sig,SIG,UCL,Rho){
 Cn=UCL-1
 rl=0
 while(Cn<UCL){
  rl=rl+1
  x <- rmvnorm(n,mu,sig)
  y=t(solve(sqrtm(SIG))%*%t(x-MU))
  S=t(y)%*%(y)/nrow(y)
  V=glasso(S,Rho)$wi
  logdet=determinant(V)$modulus 
  Cn=sum(diag(S))+logdet-sum(diag(V%*%S))
 }
 return(rl)
}


# Function to compute the ARL for each chart
fARLMultCharts=function(p,n,m,mu,sig,Rho,simno,ChartType,L){
 RL=rep(0,simno)
 if(ChartType=="NTCC"){
  UCL = qchisq(1-L/2,(n-1)*p)/(n-1);
  LCL = qchisq(L/2,(n-1)*p)/(n-1);
 }else if(ChartType=="PC1"){
  if(p==3){
   qtl=readRDS("Control_Limit_LogCovM_p_3_n_8.Rds")
  }else if(p==5){
   qtl=readRDS("Control_Limit_LogCovM_p_5_n_10.Rds")
  }else if(p==10){
   qtl=readRDS("Control_Limit_LogCovM_p_10_n_12.Rds")
  }
  UCL=quantile(qtl,1-L/2)
  LCL=quantile(qtl,L/2)
 }else if(ChartType=="GPC"){
  if(p==3){
   qtl=readRDS("det_Control_Limit_LogCovM_p_3_n_8.Rds")
  }else if(p==5){
   qtl=readRDS("det_Control_Limit_LogCovM_p_5_n_10.Rds")
  }else if(p==10){
   qtl=readRDS("det_Control_Limit_LogCovM_p_10_n_12.Rds")
  }
  UCL=quantile(qtl,1-L/2)
  LCL=quantile(qtl,L/2)
 }else if(ChartType=="PLR"){
  UCL=L
 }else if(ChartType=="GVC"){
  UCL=L
 }
 for(i in 1:simno){
   estimate=fPhaseI_estimates(m,n,p,theta,magnitude,estimator)
   muvhat=estimate$MU
   sigmamhat=estimate$COV
  if(ChartType=="NTCC"){
   RL[i]=NTCCrl(n,p,mu,muvhat,sig,sigmamhat,UCL,LCL)
  }else if(ChartType=="PC1"){
   RL[i]=fLCov(n,p,mu,muvhat,sig,sigmamhat,LCL,UCL)
  }else if(ChartType=="PLR"){
   RL[i]=fPLR(n,p,mu,muvhat,sig,sigmamhat,UCL,Rho)
  }else if(ChartType=="GVC"){
   RL[i]=GVCrl(n,p,mu,muvhat,sig,sigmamhat,L)
  }else if(ChartType=="GPC"){
    RL[i]=fLgCov(n,p,mu,muvhat,sig,sigmamhat,LCL,UCL)
  }
 }
 ARL=mean(RL)
 return(ARL)
}


# Getting Control limits
MULTCHARTSCL<-function(p,n,m,mu,sig,Rho,simno,ChartType,
                       hleft,hright,precision,arl0){
 hright<-hright
 hleft<-hleft
 h<-(hright+hleft)/2
 arl2<-fARLMultCharts(p,n,m,mu,sig,Rho,simno,ChartType,hright)
 
 arl<-fARLMultCharts(p,n,m,mu,sig,Rho,simno,ChartType,h)
 print(h)
 print(arl)
 while(abs(arl-arl0)>abs(precision)){
  arl1=arl2
  arl2=arl
  hleft=hright
  hright<-h
  h<-hright-(hright-hleft)*(arl2-arl0)/(arl2-arl1)
  arl<-fARLMultCharts(p,n,m,mu,sig,Rho,simno,ChartType,h)
  print(h)
  print(arl)
 }
 return(h)
}

###############################################################
## Initialization

fUCL=function(p,ChartType){
 if(ChartType=="PLR"){
  if(p==3){
   L=2.99
  }else if(p==5){
   L=2.955
  }else if(p==10){
    L=5.9
  }
 }else if(ChartType=="GVC"){
  if(p==3){
   L=3
  }else if(p==5){
   L=3
  }else if(p==10){
    L=6.5
  }
 }else{
  L= 1/920
 }
 return(L)
}

#Input parameter to change 
Charts=c("NTCC","GVC","PLR","GPC","PC1")
estimator="MVE"
ChartType=Charts[4]
sim = 10000
Rho=0.1
arl0<-370
precision<-4
## SIMULATION

#set.seed(111)


start.time <- Sys.time()

#print(output)
#OC1
lp=3#c(3,5)
for(k in 1:length(lp)){
 p=lp[k]
 if(p==3){
  n=8
 }else if(p==5){
  n=10
 }else if(p==10){
   n=12
 }
 mu=rep(0,p)
 vm=c(25,50,75,100,150,200)
 theta=0;
 magnitude=1
 sig.ic=diag(p)
 UCL=fUCL(p,ChartType)
 hleft<-UCL-0.01*UCL
 hright<-UCL+0.01*UCL
 H=c()
 for(j in 1:length(vm)){
  m=vm[j]
  H[j]<-MULTCHARTSCL(p,n,m,mu,sig.ic,Rho,sim,ChartType,hleft,hright,precision,arl0)
  print(c(m,H[j]))
 }
 df=data.frame(vm,H)
 filename = paste(ChartType,"_estimator_",estimator, "_14102023_CL_p_" ,p,"n_",n, ".csv", sep="")
 write.csv(df, filename)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

