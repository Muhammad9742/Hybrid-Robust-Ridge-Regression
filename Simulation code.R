###################################################################
# R-Code for Robust M- Ridge Regression #
# In presence of multicollinearity and y-direction outlier #
# Authors:Muhammad Khan  
###################################################################
set.seed(7860)
library(MASS)

sigma1 = c(0.1,1,2,5)
ssize=c(25,50)
pred=c(4,10)
corr=c(0.80,0.90,0.95,0.99)
N=10000                                 #Simulation runs

MSE.r=matrix(0,17,length(corr))
MSE.s=array(0,dim = c(17,length(corr),length(sigma1)))
MSE.n=array(0,dim = c(17,length(corr),length(sigma1),length(ssize)))
MSE.p=array(0,dim = c(17,length(corr),length(sigma1),length(ssize),length(pred)))

#-----------------------
#Predictors loop
for(b in 1:length(pred)){
  p=pred[b]
  I=diag(p)
  
  #------------------------
  #Sample size loop
  for(a in 1:length(ssize)){
    n=ssize[a]
    
z=matrix(0,n,p)
for (i in 1:p) {
  z[,i]=rnorm(n,0,1)                     #Random no generation   
}
x=matrix(0,n,p)
#--------------------------
#Error variance loop
for(S in 1:length(sigma1)){
  sigma=sigma1[S]
for(j in 1:length(corr)){
  r=corr[j]
  
  for (i in 1:p) {
    x[,i]=sqrt(1-r^2)*z[,i]+r*z[,p]    #Generating x-matrix
  }
  
  x=scale(x,center = TRUE,scale = TRUE)    #standardizing x-matrix
  c=t(x)%*%x
  D=eigen(c)$vectors
  ev=eigen(c)$values
  beta=D[,which.max(ev)]
  Z=x%*%D
  lam=t(Z)%*%Z
  lam=round(lam,5)
  lamda=diag(lam)
  alpha=t(D)%*%beta
  
  #Defining Vectors
  K1=rep(0,N)
  K2=rep(0,N)
  K3=rep(0,N)
  K4=rep(0,N)
  K5=rep(0,N)
  K6=rep(0,N)
  K7=rep(0,N)
  #K8=rep(0,N)
  K14=rep(0,N)
  K15=rep(0,N)
  K16=rep(0,N)
  K17=rep(0,N)
  K18=rep(0,N)
  K19=rep(0,N)
  K20=rep(0,N)

  
  
  #Defining matrices
  
  alphahat.ols=matrix(0,nrow = p,ncol = N)
  alphahatk.RR=matrix(0,nrow = p,ncol = N)
  alphahat.m=matrix(0,nrow = p,ncol = N)
  alphahat.MHKB=matrix(0,nrow = p,ncol = N)
  alphahat.MLW=matrix(0,nrow = p,ncol = N)
  alphahat.MAM1=matrix(0,nrow = p,ncol = N)
  alphahat.MAM2=matrix(0,nrow = p,ncol = N)
  alphahat.MAM3=matrix(0,nrow = p,ncol = N)
  alphahat.MAM4=matrix(0,nrow = p,ncol = N)
  alphahat.MAM5=matrix(0,nrow = p,ncol = N)
 
  
  alphahat.AMHKB  = matrix(0,nrow = p,ncol = N)
  alphahat.AMLW   = matrix(0,nrow = p,ncol = N)
  alphahat.AMAM1  = matrix(0,nrow = p,ncol = N)
  alphahat.AMAM2  = matrix(0,nrow = p,ncol = N)
  alphahat.AMAM3  = matrix(0,nrow = p,ncol = N)
  alphahat.AMAM4  = matrix(0,nrow = p,ncol = N)
  alphahat.AMAM5  = matrix(0,nrow = p,ncol = N)
 
  
  
  #Simulation loop starts here
  
  for (i in 1:N) {
    e1=rnorm(ceiling(0.80*n),0,sigma)
    e2=rnorm(floor(0.20*n),0,sigma+20)
    e=c(e1,e2)
    y=Z%*%alpha+e
    betahat=solve(t(x)%*%x)%*%t(x)%*%y
    yhat=x%*%betahat
    sigmahat=sum((y-yhat)^2)/(n-p)
    alphahat=solve(lam)%*%t(Z)%*%y
    
    #OLS estimator
    
    alphahat.ols[,i]=alphahat
    
    #Ridge Regression estimator of HKB method (1975)
    
    HKB=(p*sigmahat)/(sum(alphahat^2))
    K1[i]=HKB
    alphahatk=solve(t(Z)%*%(Z)+I*K1[i])%*%t(Z)%*%y
    alphahatk=c(alphahatk)
    alphahatk.RR[,i] = alphahatk
    
    
    # M-estimator denoted by alphahat.m of Huber 1981
    
    r.fit=rlm(y~x-1,psi=psi.huber,k=1.345,maxit = 1000)
    alphahat.m[,i]=as.vector(r.fit$coefficients)
    # To obtain A2
    e=r.fit$residuals
    s=r.fit$s
    u=e/s
    w.psi=psi.huber(u)
    u.psi=u*w.psi
    deriv.psi=psi.huber(u,deriv = 1)
    num.a2=((s^2)*sum(u.psi^2))/(n-p)
    denum.a2=(sum(deriv.psi)/n)^2
    A2=num.a2/denum.a2
    
    # Silvapulle (1991) Proposed M-estimator of HKB(1975) & LW(1976) 
    
    MHKB=(p*A2)/sum(alphahat.m[,i]^2)
    K2[i]=MHKB
    alphahat.MHKB[,i]=solve(lam+I*K2[i])%*%lam%*%alphahat.m[,i]
    
    MLW=(p*A2)/sum(ev*alphahat.m[,i]^2)
    K3[i]=MLW
    alphahat.MLW[,i]=solve(lam+I*K3[i])%*%lam%*%alphahat.m[,i]
   
    #### Majid etal.[22]
    
    KM= A2/alphahat.m[,i]^2
    
    
    MAM1= sum(KM)/p
    K4[i]=MAM1   
    alphahat.MAM1[,i]=solve(lam+I*K4[i])%*%lam%*%alphahat.m[,i]
    
    MAM2= (prod(KM))^1/p
    K5[i]=MAM2   
    alphahat.MAM2[,i]=solve(lam+I*K5[i])%*%lam%*%alphahat.m[,i]
    
    MAM3= max(KM)
    K6[i]=MAM3   
    alphahat.MAM3[,i]=solve(lam+I*K6[i])%*%lam%*%alphahat.m[,i]
    
    MAM4= p*sum(KM)
    K7[i]=MAM4   
    alphahat.MAM4[,i]=solve(lam+I*K7[i])%*%lam%*%alphahat.m[,i]
    
    MAM5= diag(KM)
    #K8=MAM5   
    alphahat.MAM5[,i]=solve(lam+I*MAM5)%*%lam%*%alphahat.m[,i]
  
    
   
  
  ######Adjusted estimators proposed by Muhammad Khan
  
  tau = max(ev)/p*min(ev)
  Sigma.adj= sigmahat*tau
  
  AMHKB = (p * Sigma.adj) / sum(alphahat.m[,i]^2)
  K14[i] = AMHKB
  alphahat.AMHKB[,i] = solve(lam + I * K14[i]) %*% lam %*% alphahat.m[,i]
  
  AMLW = (p * Sigma.adj) / sum(ev * alphahat.m[,i]^2)
  K15[i] = AMLW
  alphahat.AMLW[,i] = solve(lam + I * K15[i]) %*% lam %*% alphahat.m[,i]
  
 
  AKM = Sigma.adj / (alphahat.m[,i]^2)
  
  AMAM1 = sum(AKM)/p
  K16[i] = AMAM1
  alphahat.AMAM1[,i] = solve(lam + I * K16[i]) %*% lam %*% alphahat.m[,i]
  
  AMAM2 = (prod(AKM))^1/p
  K17[i] = AMAM2
  alphahat.AMAM2[,i] = solve(lam + I * K17[i]) %*% lam %*% alphahat.m[,i]
  
  AMAM3 = max(AKM)
  K18[i] = AMAM3
  alphahat.AMAM3[,i] = solve(lam + I * K18[i]) %*% lam %*% alphahat.m[,i]
  
  AMAM4 = p * sum(AKM)
  K19[i] = AMAM4
  alphahat.AMAM4[,i] = solve(lam + I * K19[i]) %*% lam %*% alphahat.m[,i]
  
  AMAM5 = diag(AKM)
  K20[i] = AMAM5
  alphahat.AMAM5[,i] = solve(lam + I * K20[i]) %*% lam %*% alphahat.m[,i]

  
  
  
  
 
  
   } #End of simulation loop for N=5000 runs.
  
  MSE.ols=sum((alphahat.ols-c(alpha))^2)/N
  MSE.RR=sum((alphahatk.RR-c(alpha))^2)/N
  MSE.m=sum((alphahat.m-c(alpha))^2)/N
  MSE.MHKB=sum((alphahat.MHKB-c(alpha))^2)/N
  MSE.MLW=sum((alphahat.MLW-c(alpha))^2)/N
  MSE.MAM1    = sum((alphahat.MAM1 - c(alpha))^2) / N
  MSE.MAM2    = sum((alphahat.MAM2 - c(alpha))^2) / N
  MSE.MAM3    = sum((alphahat.MAM3 - c(alpha))^2) / N
  MSE.MAM4    = sum((alphahat.MAM4 - c(alpha))^2) / N
  MSE.MAM5    = sum((alphahat.MAM5 - c(alpha))^2) / N
  
  
  MSE.AMHKB  = sum((alphahat.AMHKB - c(alpha))^2) / N
  MSE.AMLW   = sum((alphahat.AMLW - c(alpha))^2) / N
  MSE.AMAM1  = sum((alphahat.AMAM1 - c(alpha))^2) / N
  MSE.AMAM2  = sum((alphahat.AMAM2 - c(alpha))^2) / N
  MSE.AMAM3  = sum((alphahat.AMAM3 - c(alpha))^2) / N
  MSE.AMAM4  = sum((alphahat.AMAM4 - c(alpha))^2) / N
  MSE.AMAM5  = sum((alphahat.AMAM5 - c(alpha))^2) / N
  
 
 
  
  
  
  MSE=c(MSE.ols,MSE.RR,MSE.m,MSE.MHKB,MSE.MLW, MSE.MAM1, MSE.MAM2, MSE.MAM3, MSE.MAM4, MSE.MAM5, MSE.AMHKB, MSE.AMLW, 
        MSE.AMAM1, MSE.AMAM2, MSE.AMAM3, MSE.AMAM4, MSE.AMAM5)
  MSE=round(MSE,5)
  as.matrix(MSE)
  MSE.r[,j]=MSE
  
} #End of correlation loop
  MSE.s[,,S]=MSE.r
} #End of error variance loop
MSE.n[,,,a]=MSE.s
  }  #End of sample size loop
  MSE.p[,,,,b]=MSE.n
}   #End of predictors loop


col.names=c("0.80","0.90","0.95","0.99")
row.names=c("OLS","RR","HUB.M","MHKB","MLW","MAM1","MAM2","MAM3","MAM4","MAM5","AMHKB","AMLW",
            "AMAM1","AMAM2","AMAM3","AMAM4","AMAM5")
matrix.names=c("0.1","1","2","5")
matrix.names2=c("25","50")
matrix.names3=c("4","10")
dimnames(MSE.p)=list(row.names,col.names,matrix.names,matrix.names2,matrix.names3)
write.csv(MSE.p,file = "E:/PHD/csv data1 20%.csv")
