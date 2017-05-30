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
##          Solution to Exercise 5.6           ##
#################################################
#
#
BR <- function(v){
  vnext = v
  thesum = sum(exp(-v))
  for (j in (1:length(v)))
  {
    aggr =  thesum- exp(-v[j])
    thef <- function (z) (exp(-z) / (z-1) - aggr)
    vnext[j]=uniroot(thef,c(1,1E50))$root
  }
  return(vnext)
}

M = 10
PREC = 1E-10
vcur = 1+runif(n = M)
cont = TRUE
iter = 0
while (cont)
{
  print(vcur)
  iter = iter+1
  vnext = BR(vcur)
  if (sum(abs(vnext-vcur))<PREC*sum(abs(vcur)))
  {cont = FALSE}
  vcur = vnext
  
}
print(vcur)
print(rep(M /(M-1),M))
