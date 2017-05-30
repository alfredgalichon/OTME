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
##          Solution to Exercise 8.6           ##
#################################################
#
# The following line from programming example 8.1 
result = gurobi ( list(A=t(Nabla),obj=Phi,modelsense='max',rhs=n,sense='=',start=matrix(0,nbArcs,1)), params=NULL)
# should be modified into
result = gurobi ( list(A=t(Nabla),obj=Phi,modelsense='max',rhs=n,sense=ifelse(n<0,'>',ifelse(n>0,'<','=')),start=matrix(0,nbArcs,1)), params=NULL)
