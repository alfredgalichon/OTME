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
###     PE 6.1: The Iterative Proportional    ###
###         Fitting Procedure (IPFP)          ###
#################################################
#
#
seed = 777
tol=1E-12
maxiter = 500


nbX=10
nbY=8


set.seed(seed)

Phi=matrix(runif(nbX*nbY),nrow=nbX)
sigma=1
p=rep(1/nbX,nbX)
q=rep(1/nbY,nbY)

IX=rep(1,nbX)
IY=rep(1,nbY)
tIY=t(IY)
f = p %*% tIY
g = IX %*% t(q)
ptm=proc.time()
v=rep(0,nbY)
cont = TRUE
iter = 0
while(cont)
{
  iter = iter+1
   u = sigma*log(apply(g * exp( ( Phi - IX %*% t(v) ) / sigma ),1,sum))
   vnext = sigma*log(apply(f * exp( ( Phi - u %*% tIY ) / sigma ),2,sum))
   error = max(abs(apply(g * exp( ( Phi - IX %*% t(vnext) - u %*% tIY ) / sigma ),1,sum)-1))
  #print(error)
   if( (error<tol) | (iter >= maxiter)) {cont=FALSE}
   v=vnext
}
pi = f * g * exp( ( Phi - IX %*% t(v) - u %*% tIY ) / sigma )
val = sum(pi*Phi) - sigma* sum(pi*log(pi))
time <- proc.time()-ptm
print(paste0('Time elapsed (IPFP) =', time[1], 's.')) 
if (iter >= 500 ) 
{print('maximum number of iterations reached')
} else {
  print(paste0("Converged in ",iter, " steps."))
  print(paste0("Value of the problem (IPFP) =",val))
  nrow = min(10,nbX)
  ncol = min(10,nbY)
  print("u=")
  print(u[1:nrow])
  print("v=")
  print(v[1:ncol])
  print("pi=")
  print(pi[1:nrow,1:ncol])
}


