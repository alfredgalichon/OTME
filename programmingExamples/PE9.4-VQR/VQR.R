#################################################
###  Optimal Transport Methods in  Economics  ###
#################################################
##########       Alfred Galichon       ##########
#################################################
##########  Princeton University Press ##########
#################################################
#
#
#################################################
##      PE 9.3: Vector Quantile Regression     ##
#################################################
# Reference:
# G. Carlier, V. Chernozhukov, A. Galichon (2016): 
# "Vector Quantile Regression: an optimal transport approach". 
# Annals of Statistics.  
###############################################################
###############################################################
#
###############################################################
######## I. Library initialization and data loading  ##########
###############################################################

library('lattice')
library('Matrix')
library('slam')
library('matrixcalc')
library('gurobi')

load_VQRData<-function(folderName){
  VQRDataPath  <- paste0(getwd(),"/Data/",folderName)
  Xprov 	<- as.matrix(read.csv(paste0(VQRDataPath,"/X.txt"),sep="\t", header=FALSE))
  Y		<- as.matrix(read.csv(paste0(VQRDataPath,"/Y.txt"),sep="\t", header=FALSE))
  n		<- dim(Xprov)[1]
  if (n != dim(Y)[1]) {stop("Numbers of observations of X and Y do not match")}
  d		<- dim(Y)[2]
  r		<- dim(Xprov)[2]+1
  X		<- cbind(matrix(1,n,1),Xprov)
  
  return(list(X=X,Y=Y,n=n,d=d,r=r,step=step))
}

###############################################################
################# II. Core computation ########################
###############################################################

VQRTp<-function(X,Y,U,mu,nu){
# Monge-Kantorovich (transportation version)
# computes VQR via primal program
# (c) by Guillaume Carlier, Victor Chernozhukov and Alfred Galichon

n	<-dim(Y)[1]
d	<-dim(Y)[2]
r	<-dim(X)[2]
m	<-dim(U)[1]
if ((n != dim(X)[1]) |( d != dim(U)[2] )) {stop("wrong dimensions")}
xbar	<- t(nu)%*% X

c	<- -t(matrix(U %*% t(Y),nrow=n*m))
# c <- t(-kronecker(Y,U) %*% matrix(diag(1,d),nrow=d*d)) ### TO BE REMOVED
A1 <- kronecker(sparseMatrix(1:n,1:n),matrix(1,1,m))
A2 <- kronecker(t(X),sparseMatrix(1:m,1:m))
f1 <- matrix(t(nu),nrow=n)
f2 <- matrix(mu %*% xbar,nrow=m*r)
e <- matrix(1,m*n,1)
A <- rbind2(A1,A2)
f <- rbind2(f1,f2)
pi_init <- matrix(mu %*% t(nu),nrow=m*n)

############### LP SOLVING PHASE ###############
result		<- gurobi (list(A=A,obj=c,modelsense="min",rhs=f,ub=e,sense="=",start=pi_init), params=NULL ) 
if (result$status=="OPTIMAL") {pivec <- result$x; Lvec <- t(result$pi) } else {stop("optimization problem with Gurobi")}

#############################################

pi 		<-matrix(pivec,nrow=m)
L1vec 	<-Lvec[1:n]
L2vec 	<-Lvec[(n+1):(n+m*r)]
L1		<-matrix(L1vec,1)
L2		<-matrix(L2vec,m)

psi		<- -t(L1)
b	    <- -L2
val		<- matrix.trace(t(U) %*% pi %*% Y)

#############################################


 return(list(pi=pi,psi=psi,b=b, val=val))
}

###############################################################
################### II. 1-dimensional VQR #####################
###############################################################
ComputeBeta1D <- function( mu,b){
  m <-dim(mu)[1]
  D <-diag(1,m);
  for (i in 1:(m-1)) {D[i+1,i] <- (-1)}
  beta<-diag(c(1/mu))%*% D %*% b
  return(beta)
}
###############################################################
################### III. 2-dimensional VQR #####################
###############################################################

ComputeBetaEtAl2D<-function(b_prov,T,U_prov,pi_prov,step){
  
  m_prov  <- dim(b_prov)[1]
  r		<- dim(b_prov)[2]
  n		<- dim(pi_prov)[2]
  nonzind	<- which((U_prov[,1] != 0) & (U_prov[,2] != 0))
  U 		<- cbind(U_prov[nonzind,1],U_prov[nonzind,2])
  m		<- dim(U)[1]
  mu		<- matrix(1/m,m,1)
  beta1		<- matrix(0,m,r)
  beta2		<- matrix(0,m,r)
  for (k in 1:r) { theGrad <- grad2D(b_prov[,k],T,U_prov,step); beta1[,1] <- theGrad$D1bk[nonzind,1]; beta2[,k]	<-beta1[,1] <- theGrad$D2bk[nonzind,1]}
  pi		<- matrix(0,m,n)
  for (i in 1:n) {pi[,i]	<- pi_prov[nonzind,i] }
  b 		<- b_prov[nonzind]
  
  return(list(beta1=beta1,beta2=beta2,U=U,m=m,mu=mu,pi=pi,b=b))
}

###############################################################
###############################################################
grad2D<-function(f,T,U,step){
  EPS  <- 0.0001
  fact  <- 10/step
  l	<- length(T)
  m	<- dim(U)[1]
  D1	<- matrix(0,m,1)
  D2	<- matrix(0,m,1)
  for (i1 in 2:l){
    u1	<- T[i1]
    for (i2 in 2:l) {
      u2		<- T[i2]
      j		<- which((fact*U[,1] + U[,2])==(fact*u1+u2))
      jprecx	<- which(abs((fact*U[,1] + U[,2])-(fact*(u1-step)+u2))<EPS)
      jprecy	<- which(abs((fact*U[,1] + U[,2])-(fact*u1+(u2-step)))<EPS)
      if ((length(jprecx)!=1)|(length(jprecy)!=1)) {stop("Problem in grad2D")}
      D1[j]		<- (f[j]-f[jprecx])/step
      D2[j]		<- (f[j]-f[jprecy])/step
    }
  }
  
  
  return(list(D1bk=D1,D2bk=D2))
}
###############################################################
###############################################################
PlotFixedx2D<-function(xeval,beta1,beta2,lT){
  yeval=yhat2D(xeval,beta1,beta2)
  persp(z=matrix(yeval$y1s,nrow=(lT-1)),theta=30, phi=30, expand=0.6,col='lightblue', shade=0.75, ltheta=120,ticktype='detailed')
  dev.new()
  persp(z=matrix(yeval$y2s,nrow=(lT-1)),theta=30, phi=30, expand=0.6,col='lightblue', shade=0.75, ltheta=120,ticktype='detailed')
  
}



###############################################################
###############################################################
prepareU2D<-function(T){
  lT  <- length(T)
  m_prov<-lT*lT
  U_prov<- matrix(0,m_prov,2)
  
  for (i in 1:lT) {for (j in 1:lT) {U_prov[i+(j-1)*lT,1] <- T[i]; U_prov[i+(j-1)*lT,2] <- T[j] }}
  mu_prov=matrix(1/m_prov,m_prov,1)
  
  return(list(m=m_prov,U=U_prov,mu=mu_prov))
}

###############################################################
###############################################################
yhat2D<-function(xs,beta1,beta2){
  m  <- dim(beta1)[1]
  r	<- dim(beta2)[2]
  p	<- dim(xs)[1]
  
  y1s	<- matrix(0,p,m)
  y2s	<- matrix(0,p,m)
  for (i in 1:p) { for (j in 1:m) { y1s[i,j] <- beta1[j,] %*% xs[i,] ; y2s[i,j] <- beta2[j,] %*% xs[i,] }}
  return(list(y1s=y1s,y2s=y2s))
}
###############################################################
###############################################################

