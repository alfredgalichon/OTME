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
##          Solution to Exercise 6.3           ##
#################################################
#
#
#################################################
# Requires packages 'transport' and 'geometry'
#################################################
library('transport')
library('geometry')
#library('animation')
SEED <- 55
MAX_ITER <- 10000
PREC <- 5E-2
EPS <- 5E-2
set.seed(SEED)
nCells <- 10
x_j <- runif(nCells)
y_j <- runif(nCells)
q <- rep(1/nCells,nCells)


anisotropic_power_diagram <- function(x_j,y_j,w_j,scale)
{
  pd1 <- power_diagram(x_j*scale[1],y_j*scale[2],w_j,rect=c(0,1,0,1))
  pd1$sites$xi <- pd1$sites$xi / scale[1]
  pd1$sites$eta <- pd1$sites$eta / scale[2]
  return(pd1)
}


computeOptPwd <- function(scale)
{
  w <- rep(0,nCells)
  t <- 1
  pwd <- anisotropic_power_diagram(x_j,y_j,w,scale)
  cont <- TRUE
  demand <- rep(0,nCells)
  while ((cont==TRUE) && (t<MAX_ITER))
  {
    for (j in 1:nCells)
    {
      cellj <-pwd$cells[[j]]
      demand[j] <- polyarea(cellj[,1],cellj[,2])
    }
    if (max(abs(demand-q))<PREC/nCells) 
    {cont<-FALSE} 
    else 
    {
      t<-t+1
      w <- w - EPS * (demand - q)
      pwd <- anisotropic_power_diagram(x_j,y_j,w,scale)
      
    }}
  return(pwd)
}



demand <- rep(0,nCells)

par(mfrow=c(3,2))
#for (k in 1:49)
for (k in c(1,10,20,30,40,49)) 
{
  print(k)
  theta <- k* pi / 100
  pwd <- computeOptPwd(c(cos(theta),sin(theta)))
  plot(pwd,weights=FALSE)
}



