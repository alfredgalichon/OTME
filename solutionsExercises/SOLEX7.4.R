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
##          Solution to Exercise 7.4           ##
#################################################
#
#
seed = 777
tol=1E-12


nbX=10
nbY=8


set.seed(seed)

Phi=matrix(runif(nbX*nbY),nrow=nbX)
Phi[1,1]=10 # there is an outlier value in Phi which causes the IPFP to break down below temperature sigma = 0.1


ipfp = function(sigma = 1, maxiter = 1000) # this is the usual IPFP
{
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
    v=vnext - mean(vnext)
  }
  pi = f * g * exp( ( Phi - IX %*% t(v) - u %*% tIY ) / sigma )
  val = sum(pi*Phi) - sigma* sum(pi*log(pi))
  time <- proc.time()-ptm
  print(paste0('')) 
  print(paste0('sigma = ', sigma, ':')) 
  print(paste0('Time elapsed (IPFP) =', time[1], 's.')) 
  if (iter >= maxiter ) 
  {print('maximum number of iterations reached')
  } else {
    print(paste0("Converged in ",iter, " steps."))
    print(paste0("Value of the problem (IPFP) =",val))
    nrow = min(10,nbX)
    ncol = min(10,nbY)
    return(val)
  }
  
}


expm = function(t) (ifelse(t<50,exp(t),(t-49)*exp(50)))  # this is the modified exponential function, which becomes linear when the argument is over 5-

ipfpm = function(sigma = 1, maxiter = 1000) # this is the modified ipfp, which uses the expm function
{
  p=rep(1/nbX,nbX)
  q=rep(1/nbY,nbY)
  
  IX=rep(1,nbX)
  IY=rep(1,nbY)
  tIY=t(IY)
  f = p %*% tIY
  g = IX %*% t(q)
  ptm=proc.time()
  u=rep(0,nbX)
  vnext=v=rep(0,nbY)
  
  cont = TRUE
  iter = 0
  while(cont)
  {
    iter = iter+1
    
    for(i in 1:nbX){
      root_fn <- function(z) (sum(q * expm( ( Phi[i,] - z - v ) / sigma )) -1 )
      u[i] = uniroot(root_fn,c(-1E30,1E30), tol = 1e-300)$root
    }
    for(j in 1:nbY){
      root_fn <- function(z) (sum(p * expm( ( Phi[,j] - z - u ) / sigma )) -1 )
      vnext[j] = uniroot(root_fn, c(-1E30,1E30),tol = 1e-300)$root
    }
    
    error = max(abs(apply(g * expm( ( Phi - IX %*% t(vnext) - u %*% tIY ) / sigma ),1,sum)-1))
    if( (error<tol) | (iter >= maxiter)) {cont=FALSE}
    v=vnext - mean(vnext)
  }
  pi = f * g * expm( ( Phi - IX %*% t(v) - u %*% tIY ) / sigma )
  logpi = log(pi)
  infs = which(logpi == -Inf)
  logpi[infs]=0
  val = sum(pi*Phi) - sigma* sum(pi*logpi)
  time <- proc.time()-ptm
  print(paste0('')) 
  print(paste0('sigma = ', sigma, ':')) 
  print(paste0('Time elapsed (IPFP) =', time[1], 's.')) 
  if (iter >= maxiter ) 
  {print('maximum number of iterations reached')
  } else {
    print(paste0("Converged in ",iter, " steps."))
    print(paste0("Value of the problem (IPFP) =",val))
    nrow = min(10,nbX)
    ncol = min(10,nbY)
    return(val)
  }
  
}




val = ipfp(sigma = 1)
val = ipfp(sigma = 0.1)
val = ipfp(sigma = 0.01)
val = ipfpm(sigma = 0.01)


