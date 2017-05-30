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
##          Solution to Exercise 8.5           ##
#################################################
#
# The following line from programming example 8.1 

library('gurobi')
library('Matrix')

originNode <- 84 #saint-germain des pres
destinationNode<- 116 #trocadero

T = 50 # number of time periods

thePath = getwd()
arcs = as.matrix(read.csv(paste0(thePath,"/arcs.csv"),sep=";", header=FALSE)) # loads the data
nbBasisNodes = max(arcs[,1]) 
nbBasisArcs = dim(arcs)[1]
nbNodes = nbBasisNodes * T
nbArcs = nbBasisArcs * (T-1)



NablaSpaceP = sparseMatrix(i=1:nbBasisArcs,j=arcs[,2],dims=c(nbBasisArcs,nbBasisNodes),x=1)
NablaSpaceN = sparseMatrix(i=1:nbBasisArcs,j=arcs[,1],dims=c(nbBasisArcs,nbBasisNodes),x=1) 
NablaSpace = NablaSpaceP -NablaSpaceN
PhiSpace <- -arcs[,3] # construct (minus) distance matrix


originNode0 = originNode
destinationNodeT = nbBasisNodes * (T-1) + destinationNode 

n = rep(0,nbNodes) # construct vector of exiting flow
n[c(originNode0,destinationNodeT)]=c(-1,1)


Nabla= sparseMatrix(i=1,j=1,dims=c(nbArcs,nbNodes),x=0)
Phi = rep(0,nbArcs)
for (t in 1:(T-1))
  for (arc in 1:nbBasisArcs)
{
    origin = which(NablaSpace[arc,]==-1) 
    destination = which(NablaSpace[arc,]==1)
    neworigin = origin + nbBasisNodes*(t-1)
    newdestination = destination + nbBasisNodes*t
    newarc = arc + (t-1)*nbBasisArcs
    Nabla[newarc,neworigin] = -1 
    Nabla[newarc,newdestination] = 1
    Phi[newarc] = PhiSpace[arc]
}

# solve LP via Gurobi
result = gurobi ( list(A=t(Nabla),obj=Phi,modelsense='max',rhs=n,sense='=',start=matrix(0,nbArcs,1)), params=NULL)
