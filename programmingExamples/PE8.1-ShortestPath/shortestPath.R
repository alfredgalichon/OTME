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
###           PE 8.1: Shortest path           ###
###           via linear programming          ###
#################################################

library('gurobi')
library('Matrix')

#originNode <- 85 #odeon
#originNode <- 228  #billancourt
#originNode <- 82 #saint-placide
#originNode <- 336 #fontenay-sous-bois
#destinationNode<-84 #saint-germain des pres
#destinationNode<-42 #belleville
#destinationNode <- 7 # charles-de-gaulle etoile
originNode <- 84 #saint-germain des pres
destinationNode<- 116 #trocadero

thePath = getwd()
arcs = as.matrix(read.csv(paste0(thePath,"/arcs.csv"),sep=";", header=FALSE)) # loads the data
namesNodes = as.matrix(read.csv(paste0(thePath,"/names.csv"),sep=";", header=FALSE)) # loads the data
nbNodes = max(arcs[,1])
nbArcs = dim(arcs)[1]

n = rep(0,nbNodes) # construct vector of exiting flow
n[c(originNode,destinationNode)]=c(-1,1)

# construct node-incidence matrix:
Nabla =  sparseMatrix(i=1:nbArcs,j=arcs[,1],dims=c(nbArcs,nbNodes),x=-1) + sparseMatrix(i=1:nbArcs,j=arcs[,2],dims=c(nbArcs,nbNodes),x=1)


Phi <- -arcs[,3] # construct (minus) distance matrix

# solve LP via Gurobi
result = gurobi ( list(A=t(Nabla),obj=Phi,modelsense='max',rhs=n,sense='=',start=matrix(0,nbArcs,1)), params=NULL)
pi = result$x
distance = -result$objval


# deduce minimal distance path:
cont = TRUE
i = originNode
writeLines(paste0(namesNodes[i]," (#", i,")"))
eqpath = which(pi>0)
rank = 0
while(cont)
{ 
  rank = rank+1
  leavingi = which(Nabla[,i]==-1)
  a = intersect(eqpath,leavingi)[1]
  j = which(Nabla[a,]==1)[1]
  writeLines(paste0(rank,": ", namesNodes[j]," (#", j,")"))
  i = j
  if(j==destinationNode) {cont<-FALSE}  
}