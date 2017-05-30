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
require('clue')
n=1000
seed=777

set.seed(seed)
ptm <- proc.time()
Phi=matrix(runif(n*n),nrow=n)
p = rep(1,n)
q = rep(1,n)

N = dim(Phi)[1]
M = dim(Phi)[2]

c = c(Phi)

A1 = kronecker(matrix(1,1,M),sparseMatrix(1:N,1:N))
A2 = kronecker(sparseMatrix(1:M,1:M),matrix(1,1,N))
A = rbind2(A1,A2)

d = c(p,q) 

result   = gurobi ( list(A=A,obj=c,modelsense="max",rhs=d,sense="="), params=list(OutputFlag=0) ) 
if (result$status=="OPTIMAL") {
  pi = matrix(result$x,nrow=N)
  u = result$pi[1:N]
  v = result$pi[(N+1):(N+M)]
  val = result$objval
} else {stop("optimization problem with Gurobi.") }

time <- proc.time()-ptm
print(paste0('Time elapsed (Gurobi) =', time[1], 's.')) 
print(paste0("Value of the problem (Gurobi) =",val))
print(u[1:10])
print(v[1:10])


set.seed(seed)
ptm = proc.time()
j_is = solve_LSAP(Phi,maximum=TRUE)
pih = sparseMatrix(1:n,j_is)

val = sum(pih*Phi)
time = proc.time()-ptm
print(paste0('Time elapsed (Hungarian algorithm) = ', time[1], 's.')) 
print(paste0("Value of the problem (Hungarian algorithm) = ",val))

