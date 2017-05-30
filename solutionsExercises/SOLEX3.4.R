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
###   PE 3.1: The Optimal Assignment Problem  ###
###           via Linear Programming          ###
#################################################
#
#
#################################################
# REQUIRES GUROBI, available from gurobi.com
#################################################

require('Matrix')
require('gurobi')
nbpoints=100
nbthetas = 50
seed=777
corrX = 0
corrY = 0.8
p = rep(1/nbpoints,nbpoints)
q = rep(1/nbpoints,nbpoints)
d = c(p,q) 
A1 = kronecker(matrix(1,1,nbpoints),sparseMatrix(1:nbpoints,1:nbpoints))
A2 = kronecker(sparseMatrix(1:nbpoints,1:nbpoints),matrix(1,1,nbpoints))
A = rbind2(A1,A2)

set.seed(seed)
corrMatX = matrix(c(1,corrX,0,sqrt(1-corrX^2)),nrow=2)
corrMatY = matrix(c(1,corrY,0,sqrt(1-corrY^2)),nrow=2)
Xs = matrix(runif(2*nbpoints),nrow=nbpoints) %*% corrMatX
Ys = matrix(runif(2*nbpoints),nrow=nbpoints) %*% corrMatY
Phi1 = Xs[,1]%*%t(Ys[,1])
Phi2 = Xs[,2]%*%t(Ys[,2])

Cs = matrix(0,nbthetas+1,2)
for (k in 1: (nbthetas+1) )
{
  theta = 2*pi*k/nbthetas

  Phi = cos(theta)*Phi1 + sin(theta)*Phi2
  result   = gurobi ( list(A=A,obj=c(Phi),modelsense="max",rhs=d,sense="="), params=list(OutputFlag=0) ) 
  if (result$status=="OPTIMAL") {
    mu = matrix(result$x,nrow=nbpoints)
  } else {stop("optimization problem with Gurobi.") }
  C1 = sum(mu*Phi1)
  C2 = sum(mu*Phi2)
  Cs[k,]=c(C1,C2)
}
plot(Cs)
lines(Cs,col=rgb(0, 0,1))


