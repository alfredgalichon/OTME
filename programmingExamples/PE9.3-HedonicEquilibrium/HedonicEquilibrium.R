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
##         PE 9.3: Hedonic equilibrium         ##
#################################################

library('gurobi')
library('Matrix')

nbX = 3
nbY = 2
nbZ = 2
p = rep(1/nbX,nbX)
q = rep(1/nbY,nbY)
alpha = matrix(runif(nbX*nbZ),nrow=nbX)
gamma = matrix(runif(nbZ*nbY),nrow=nbZ)


nbNodes = nbX+nbZ+nbY
nbArcs = nbX*nbZ + nbZ*nbY
n = c(-p,rep(0,nbZ),q)
Phi = c(alpha,gamma)



# construct node-incidence matrix:
fullGradientMatrix = function(nrow,ncol)
{
narcs = ncol*nrow
nnodes = ncol+nrow
sources = rep(0,narcs)
targets = rep(0,narcs)
for (y in 1:ncol)
{for (x in 1:nrow)
{ sources[(y-1)*nrow+x] = x 
  targets[(y-1)*nrow+x] = nrow+y
}}
return(sparseMatrix(i=1:narcs,j=sources,dims=c(narcs,nnodes),x=-1) + sparseMatrix(i=1:narcs,j=targets,dims=c(narcs,nnodes),x=1))
}


Nabla1=cbind2(fullGradientMatrix(nbX,nbZ),sparseMatrix(i=c(),j=c(),dims=c(nbX*nbZ,nbY)))
Nabla2= cbind2(sparseMatrix(i=c(),j=c(),dims=c(nbY*nbZ,nbX)),fullGradientMatrix(nbZ,nbY))
Nabla = rbind2(Nabla1,Nabla2)

# solve LP via Gurobi
result = gurobi ( list(A=t(Nabla),obj=Phi,modelsense='max',rhs=n,sense='='), params=NULL)
pi_xz = matrix(result$x[1:(nbX*nbZ)],nrow = nbX)
pi_zy = matrix(result$x[(nbX*nbY+1):(nbX*nbZ+nbZ*nbY)],nrow = nbZ)
u_x =  - result$pi[1:nbX]
w_z = - result$pi[(nbX+1):(nbX+nbZ)]
v_y = result$pi[(nbX+nbZ):(nbX+nbZ+nbY)]
print(pi_xz)
print(pi_zy)
print(u_x)
print(w_z)
print(v_y)

