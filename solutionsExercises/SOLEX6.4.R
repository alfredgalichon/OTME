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
##          Solution to Exercise 6.4           ##
#################################################
#
library('adagio')

TOL<-1E-9
set.seed(8237)
nr<-40
na<-50
n<-nr*na
d<-2
nr2displ<-7
na2displ<-6
p<-2 #Lp norm

#############################################################
#############################################################
### Auxiliary functions
#############################################################
#############################################################
clayton <-function(theta,n,d,sigma){
  #  y is an nxd matrix of n draws from an d dimensional clayton distribution with parameter theta
  #  theta=0 is uniform, theta=infinity is perfect correlation;  uses Marshall-Olkin's method
  
  samplegamma<-rgamma(n,theta,1)
  sampleunif<-runif(d*n)
  
  
  gamma <- samplegamma %*% array(1,c(1,d))
  y <- (1-log(matrix(sampleunif,nrow=n)) /gamma) ^(-1/theta)
  y<-qnorm(y)
  y<-kronecker(matrix(sigma,nrow=1),matrix(1,n,1))*y
  return(y)
}

#############################################################
sqdist<-function(a,b){
  return(abs(a-b)^p)
}
#############################################################
rearrange<-function(themat,perm){
  nlines = dim(themat)[1]
  ncol = dim(themat)[2]
  newmat = matrix(0,nlines,ncol)
  for (i in 1:nlines)
  {
    newmat[i,] = themat[perm[i],]
    
  }
  return(newmat)
}

#############################################################
#############################################################
nDraws = 100 # change this into 1000 to answe the question
ys = matrix(0,n,2)

for (k in 1:nDraws)
{
  rads<-(1:nr)/nr
  angles<-t((1:na)*2*pi/na)
  xhor<-t(cos(angles))%*%rads
  xver<-t(sin(angles))%*%rads
  x<-cbind(matrix( xhor,ncol=1,byrow=F),matrix(xver,ncol=1,byrow=F))
  y<-clayton(1/3,n,d,c(1, 2))
  C<-outer(x[,1],y[,1],"sqdist")+outer(x[,2],y[,2],"sqdist")
  print(system.time(res <- assignment(100*C)))
  perm<-  res$perm
  
  ys = ys + rearrange(y,perm) / nDraws  
}


plot(x)
plot(ys)

irad1<-which(  abs((x[,1]^2+x[,2]^2)-1)< TOL )
#jrad1<-perm[irad1]
jrad1<-irad1
plot(ys[jrad1,])


rads2plot<-quantile(rads,probs=c((1:nr2displ)/nr2displ),type=1)
angles2plot<-quantile(angles,probs=c((1:na2displ)/nr2displ),type=1)


for (rad in rads2plot) {
  irad<-which(  abs((x[,1]^2+x[,2]^2)-rad^2)< TOL )
  #jrad<-perm[irad]
  jrad<-irad
  #lines(x[irad,],col=rgb(rad, 0,0))
  todraw<-ys[jrad,]
  todraw<-rbind(todraw,todraw[1,])
  lines(todraw,col=rgb(rad, 0,0))
}

for (angle in angles2plot) {
  iangles<-which(  (abs(sqrt(x[,1]^2+x[,2]^2)*cos(angle)-x[,1])+abs(sqrt(x[,1]^2+x[,2]^2)*sin(angle)-x[,2]))< TOL )
  #jangles<-perm[iangles]
  jangles<-iangles
  #lines(x[irad,],col=rgb(rad, 0,0))
  todraw<-ys[jangles,]
  lines(todraw,col=rgb(0, 0,1))
}  