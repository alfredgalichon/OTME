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
##          Solution to Exercise 4.3           ##
#################################################
#
#
#################################################
# Requires packages 'transport' and 'geometry'
#################################################

require('Matrix')
require('gurobi')

# optimal assignment algorithm, see L1
optAssign<- function (Phi,p,q)
{
  N<-dim(Phi)[1]
  M<-dim(Phi)[2]
  
  c <- t(matrix(Phi,ncol=1))
  A1 <- kronecker(matrix(1,1,M),sparseMatrix(1:N,1:N))
  A2 <- kronecker(sparseMatrix(1:M,1:M),matrix(1,1,N))
  d1 <- matrix(p,ncol=1)
  d2 <- matrix(t(q),ncol=1)
  A <- rbind2(A1,A2)
  d <- rbind2(d1,d2)
  pi_init <- matrix(p %*% t(q) , ncol=1)
  
  result    <- gurobi ( list(A=A,obj=c,modelsense="max",rhs=d,sense="=",start=pi_init), params=list(OutputFlag=0) ) 
  if (result$status=="OPTIMAL") {
    pi <- matrix(result$x,nrow=N)
    u <- result$pi[1:N] + result$pi[N+1]
    v <- result$pi[(N+1):(N+M)] - result$pi[N+1]
    val <- result$objval
  }
  else {stop("optimization problem with Gurobi")}
  return(list(val=val,pi=pi,u=u,v=v))
}


y <- c(0.1, 0.3, 0.6, 1)
q <- c(0.25, 0.25, 0.25, 0.25)
nbDraws <- 1E4

# determination in closed-form
M <- length(y)
z <- y[2:M]
z0 <- y[1]
maxzz <- pmax(z %*% matrix(1,1,M-1), matrix(1,M-1,1) %*% t(z))
maxzz0 <- pmax(z,z0)
maxzz0Mat <- matrix(maxzz0,ncol=1) %*% matrix(1,1,M-1)
A <- maxzz - maxzz0Mat - t(maxzz0Mat) + z0
b <- maxzz0 - z0
c <- z0 / 2
qtilde <- q[2:M]
Aq <- c(A %*% matrix(qtilde,ncol=1))
W_cf <- sum(qtilde * Aq) / 2 + sum(qtilde*b) + c
v_cf <-  c(0,Aq + b) 

# determination by simulation
Xsim <- runif(nbDraws)
PhiSim <- matrix(Xsim,ncol=1) %*% matrix(y,nrow=1)
pSim <- rep(1/nbDraws,nbDraws)
resSim <- optAssign(PhiSim,pSim,q)
W_sim <- resSim$val
v_sim <- resSim$v
# plot u(x) as found by simulation
plot(Xsim,resSim$u)

# compare values found by closed-form and simulation
print(W_cf)
print(W_sim)
print(v_cf)
print(v_sim)

