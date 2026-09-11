rm(list=ls())
set.seed(1990)
library(MASS)
p=4                                      #No of explanatory variables
I=diag(p)
setwd("E:/PHD/Real data")      #To change working directory
data=read.csv("Tobacco.csv",header=TRUE)
x1=data$x1
x2=data$x2
x3=data$x3
x4=data$x4
n<- length(x1)
x<-cbind(x1,x2,x3,x4)
x=scale(x,center = TRUE,scale = TRUE)   #Scaling x,,, standardizing x and y
y=data$y
y=(y-mean(y))/sd(y)



#y=scale(x,center = TRUE,scale = TRUE)
eigen(cor(x))
ev=eigen(cor(x))$values
vec=eigen(cor(x))$vectors
which.max(ev)
CN=max(ev)/min(ev)        #Computing condition number
CI=sqrt(max(ev)/min(ev))
beta=vec[,which.max(ev)]
C=t(x)%*%(x)
D=eigen(C)$vectors
Z=x%*%D
lam=t(Z)%*%Z
lam=round(lam,4)
lamda=diag(lam)
alpha=t(D)%*%beta

betahat=solve(t(x)%*%x)%*%t(x)%*%y
yhat=x%*%betahat
sigmahat=(sum((y-yhat)^2))/(n-p)
alphahat=solve(lam)%*%t(Z)%*%y
sigmahat2 = (sum((y - yhat)^2)) / (n - p)
#y[5]=y[5]+10*sigmahat
#y[15]=y[15]-5*sigmahat
#y[20]=y[20]+20*sigmahat
#y[35]=y[35]-750*sigmahat
#y[40]=y[40]+1000*sigmahat

#Estimation of ridge parameter
#OLS estimator
alphahat.ols=alphahat
alphahat.ols=c(alphahat.ols)

# HKB (1976)
HKB=(p*sigmahat)/(sum(alphahat.ols^2))
K1=HKB
alphahatk1=solve(lam+I*K1)%*%lam%*%alphahat.ols
alphahatk1=c(alphahatk1)
alphahatk.RR = alphahatk1

# M-estimator denoted by alphahat.m of Huber 1981

r.fit=rlm(y~x-1,psi=psi.huber,k=1.345,maxit = 1000)
alphahat.m=as.vector(r.fit$coefficients)
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

MHKB=(p*A2)/sum(alphahat.m^2)
K2=MHKB
alphahatk2=solve(lam+I*K2)%*%lam%*%alphahat.m
alphahatk2=c(alphahatk2)
alphahatk.MHKB = alphahatk2

MLW=(p*A2)/sum(ev*alphahat.m^2)
K3=MLW
alphahatk3=solve(lam+I*K3)%*%lam%*%alphahat.m
alphahatk3=c(alphahatk3)
alphahatk.MLW = alphahatk3

#### Majid etal.[22]

KM= (A2)/alphahat.m^2

MAM1= sum(KM)/p
K4=MAM1  
alphahatk4=solve(lam+I*K4)%*%lam%*%alphahat.m
alphahatk4=c(alphahatk4)
alphahatk.MAM1 = alphahatk4

MAM2= (prod(KM))^1/p
K5=MAM2   
alphahatk5=solve(lam+I*K5)%*%lam%*%alphahat.m
alphahatk5=c(alphahatk5)
alphahatk.MAM2 = alphahatk5

MAM3= max(KM)
K6=MAM3   
alphahatk6=solve(lam+I*K6)%*%lam%*%alphahat.m
alphahatk6=c(alphahatk6)
alphahatk.MAM3 = alphahatk6

MAM4= p*sum(KM)
K7=MAM4   
alphahatk7=solve(lam+I*K7)%*%lam%*%alphahat.m
alphahatk7=c(alphahatk7)
alphahatk.MAM4 = alphahatk7

K8 = matrix(0, nrow = p, ncol = p)
MAM5 = diag(KM)
K8 = diag(MAM5)

# Assuming K8 is a diagonal matrix (as per your previous code)
K8 = diag(KM)  # Initialize K8 as a diagonal matrix

# Get the number of rows and columns in the matrix (should be equal in a diagonal matrix)
n_rows = nrow(K8)
n_cols = ncol(K8)

# Ensure that the matrix is indeed square (has equal rows and columns)
if (n_rows != n_cols) {
  stop("K8 is not a square matrix. It should be a diagonal matrix.")
}

# Determine the size of the diagonal
diag_size = n_rows

# Initialize a variable to store the diagonal value
K8_value = NULL

# Extract the diagonal value manually
for (i in 1:diag_size) {
  K8_value=K8[i, i]
}


K8=K8_value

print(K8)


alphahatk8=solve(lam+I*K8)%*%lam%*%alphahat.m
alphahatk8=c(alphahatk8)
alphahatk.MAM5 = alphahatk8

tau = max(ev)/p*min(ev)
Sigma.adj= sigmahat*tau

AMHKB = (p * Sigma.adj) / sum(alphahat.m^2)
K14 = AMHKB
alphahatk14=solve(lam+I*K14)%*%lam%*%alphahat.m
alphahatk14=c(alphahatk14)
alphahatk.AMHKB = alphahatk14

AMLW = (p * Sigma.adj) / sum(ev * alphahat.m^2)
K15 = AMLW
alphahatk15=solve(lam+I*K15)%*%lam%*%alphahat.m
alphahatk15=c(alphahatk15)
alphahatk.AMLW = alphahatk15

AKM = Sigma.adj / (alphahat.m^2)

AMAM1 = sum(AKM)/p
K16 = AMAM1
alphahatk16=solve(lam+I*K16)%*%lam%*%alphahat.m
alphahatk16=c(alphahatk16)
alphahatk.AMAM1 = alphahatk16

AMAM2 = (prod(AKM))^1/p
K17 = AMAM2
alphahatk17=solve(lam+I*K17)%*%lam%*%alphahat.m
alphahatk17=c(alphahatk17)
alphahatk.AMAM2 = alphahatk17

AMAM3 = max(AKM)
K18 = AMAM3
alphahatk18=solve(lam+I*K18)%*%lam%*%alphahat.m
alphahatk18=c(alphahatk18)
alphahatk.AMAM3 = alphahatk18

AMAM4 = p * sum(AKM)
K19 = AMAM4
alphahatk19=solve(lam+I*K19)%*%lam%*%alphahat.m
alphahatk19=c(alphahatk19)
alphahatk.AMAM4 = alphahatk19

K20 = matrix(0, nrow = p, ncol = p)
AMAM5 = diag(AKM)
K20 = diag(AMAM5)

# Assuming K8 is a diagonal matrix (as per your previous code)
K20 = diag(AKM)  # Initialize K8 as a diagonal matrix

# Get the number of rows and columns in the matrix (should be equal in a diagonal matrix)
n_rows1 = nrow(K20)
n_cols1 = ncol(K20)

# Ensure that the matrix is indeed square (has equal rows and columns)
if (n_rows1 != n_cols1) {
  stop("K20 is not a square matrix. It should be a diagonal matrix.")
}

# Determine the size of the diagonal
diag_size1 = n_rows1

# Initialize a variable to store the diagonal value
K20_value = NULL

# Extract the diagonal value manually
for (i in 1:diag_size1) {
  K20_value=K20[i, i]
}


K20=K20_value

print(K20)

alphahatk20=solve(lam+I*K20)%*%lam%*%alphahat.m
alphahatk20=c(alphahatk20)
alphahatk.AMAM5 = alphahatk20


comb.coeff=rbind(alphahat.ols,alphahatk.RR,alphahat.m,alphahatk.MHKB,alphahatk.MLW,alphahatk.MAM1,
                 alphahatk.MAM2,alphahatk.MAM3,alphahatk.MAM4,alphahatk.MAM5,alphahatk.AMHKB,alphahatk.AMLW,alphahatk.AMAM1,
                 alphahatk.AMAM2,alphahatk.AMAM3,alphahatk.AMAM4,alphahatk.AMAM5)
comb.coeff=round(comb.coeff,8)

###########################################################################################################
#MSE of all estimators
#MSE of OLS and RR

K<-c(0,K1) 
a=rep(0,length(K)); b=rep(0,length(K))                            #to obtain MSE of ridge estimators
for(i in 1:length(K)){
  
  a[i]=sum(ev/(ev+K[i])^2)
  b[i]=sum(((K[i]^2)*(alphahat.ols^2))/(ev+K[i])^2)
  
}
MSE1=sigmahat*a+b
#MSE of M and Ridge M-estimators
Km=c(0,K2,K3,K4,K5,K6,K7,K8)
a.m=rep(0,length(Km)); b.m=rep(0,length(Km))
for (i in 1:length(Km)) {
  a.m[i]=sum(((ev)*A2)/((ev+Km[i])^2))
  b.m[i]=sum(((Km[i]^2)*(alphahat.m^2))/(ev+Km[i])^2)
}
MSE2=a.m+b.m

#MSE of adj  Ridge M-estimators
Km1=c(K14,K15,K16,K17,K18,K19,K20)
a.m1=rep(0,length(Km1)); b.m1=rep(0,length(Km1))
for (i in 1:length(Km1)) {
  a.m1[i]=sum(((ev))/((ev+Km1[i])^2))
  b.m1[i]=sum(((Km1[i]^2)*(alphahat.m^2))/(ev+Km1[i])^2)
}
MSE3=Sigma.adj*a.m1+b.m1

MSE=c(MSE1,MSE2,MSE3)
MSE=as.matrix(MSE)
MSE=round(MSE,6)
K.hat=c(K1,K2,K3,K4,K5,K6,K7,K8,K14,K15,K16,K17,K18,K19,K20)
write.csv(K.hat,file = "D:/PHD/K.hat.csv")
row.names(MSE)=c("OLS","RR","HUB.M","MHKB","MLW","MAM1","MAM2","MAM3","MAM4","MAM5","AMHKB","AMLW",
                 "AMAM1","AMAM2","AMAM3","AMAM4","AMAM5")
colnames(MSE)=c("MSE")
write.csv(MSE,file = "E:/PHD/MSEreallife.csv")
write.csv(comb.coeff,file = "E:/PHD/coefficents.csv")

cor.x=round(cor(x),4)
cor.y=round(cor(x,y),4)
cor.xy=cbind(cor.x,cor.y)
write.csv(cor.xy,file = "E:/PHD/correlation.csv")
CN
CI
ev
xnam=paste0("x", 1:p) 
fmla=as.formula(paste("y ~ ", paste(xnam, collapse= "+"))) 
model=lm(fmla)

#Checking outliers
#Influence measures
inflm.model= influence.measures(model)
which(apply(inflm.model$is.inf, 1, any))      #outliers
summary(inflm.model)                          #Summary of outlying observations

library(corrplot)
# Plot the correlation matrix
corrplot(cor.x, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black", 
         number.cex = 0.7, col = colorRampPalette(c("blue", "white", "red"))(200))
MSE
min(MSE)
# ==============================================================
#  Studentized Residuals Method (Y-direction Outlier Detection)
# ==============================================================

# Compute Studentized Residuals
stud_res <- rstudent(model)

# Plot Studentized Residuals
plot(stud_res, type = "h", col = "darkgreen",
     main = "Studentized Residuals for y-direction Outlier Detection",
     ylab = "Studentized Residuals",
     xlab = "Observation Index")

# Add threshold reference lines
abline(h = c(-2, 2), col = "red", lty = 2)      # Common cutoff ±2
abline(h = c(-3, 3), col = "blue", lty = 3)     # Extreme cutoff ±3

# Highlight outliers beyond |2|
outliers_stud <- which(abs(stud_res) > 2)

# Label outlying points
text(x = outliers_stud,
     y = stud_res[outliers_stud],
     labels = outliers_stud,
     pos = 3, col = "red", cex = 0.8)

# Print the detected Y-direction outliers
cat("Observations identified as Y-direction outliers (|Studentized Residual| > 2):\n")
print(outliers_stud)



