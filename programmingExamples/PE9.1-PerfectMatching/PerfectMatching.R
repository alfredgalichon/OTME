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
###         PE 9.1: Perfect matching          ###
#################################################
#
#
#################################################
# REQUIRES GUROBI, available from gurobi.com
#################################################

require('Matrix')
require('gurobi')
nbX=10
nbY=7
seed=777
set.seed(seed)
Phi= - matrix(sample(0:1,size=nbX*nbY,replace=TRUE),nrow=nbX) # simulate a matrix of 0 and -1
p =  rep(1/nbX,nbX)
q = rep(1/nbY,nbY)


c = c(Phi)

A1 = kronecker(matrix(1,1,nbY),sparseMatrix(1:nbX,1:nbX))
A2 = kronecker(sparseMatrix(1:nbY,1:nbY),matrix(1,1,nbX))
A = rbind2(A1,A2)

d = c(p,q) 

result   = gurobi ( list(A=A,obj=c,modelsense="max",rhs=d,sense="="), params=list(OutputFlag=0) ) 
if (result$status=="OPTIMAL") {
  pi = matrix(result$x,nrow=nbX)
  u = result$pi[1:nbX]
  v = result$pi[(nbX+1):(nbX+nbY)]
  val = result$objval
} else {stop("optimization problem with Gurobi.") }

print("u:")
print(u)
print("v=")
print(v)
print(paste0("Value of the problem = ", val))
if (val==0) {print("The value of the problem is zero, thus there is a perfect matching.")
  } else {print(("The value of the problem is positive, hence there is no perfect matching."))}

