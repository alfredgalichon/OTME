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
##         PE 6.1: Vector Quantiles            ##
#################################################
#
# Reference: 
# Chernozhukov, V., Galichon, A., Hallin, M., and Henry, M. (2016): 
# "Monge-Kantorovich Depth, Quantiles, Ranks and Signs". Annals of Statistics, forthcoming.
#
#################################################
# Requires package 'adagio' 
#################################################

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
#############################################################
rads<-(1:nr)/nr
angles<-t((1:na)*2*pi/na)
xhor<-t(cos(angles))%*%rads
xver<-t(sin(angles))%*%rads
x<-cbind(matrix( xhor,ncol=1,byrow=F),matrix(xver,ncol=1,byrow=F))
y<-clayton(1/3,n,d,c(1, 2))
C<-outer(x[,1],y[,1],"sqdist")+outer(x[,2],y[,2],"sqdist")
print(system.time(res <- assignment(100*C)))
perm<-  res$perm

plot(x)
plot(y)

irad1<-which(  abs((x[,1]^2+x[,2]^2)-1)< TOL )
jrad1<-perm[irad1]
plot(y[jrad1,])


rads2plot<-quantile(rads,probs=c((1:nr2displ)/nr2displ),type=1)
angles2plot<-quantile(angles,probs=c((1:na2displ)/nr2displ),type=1)


for (rad in rads2plot) {
  irad<-which(  abs((x[,1]^2+x[,2]^2)-rad^2)< TOL )
  jrad<-perm[irad]
  
  #lines(x[irad,],col=rgb(rad, 0,0))
  todraw<-y[jrad,]
  todraw<-rbind(todraw,todraw[1,])
  lines(todraw,col=rgb(rad, 0,0))
}

  for (angle in angles2plot) {
    iangles<-which(  (abs(sqrt(x[,1]^2+x[,2]^2)*cos(angle)-x[,1])+abs(sqrt(x[,1]^2+x[,2]^2)*sin(angle)-x[,2]))< TOL )
    jangles<-perm[iangles]
    
    #lines(x[irad,],col=rgb(rad, 0,0))
    todraw<-y[jangles,]
    lines(todraw,col=rgb(0, 0,1))
  }  
