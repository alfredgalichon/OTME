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
##          Solution to Exercise 4.4           ##
#################################################
a = 2
b = 1

cdf_P = punif
quantile_P = qunif

cdf_P = pnorm
quantile_Q = qnorm

Phi = function(x,y) (x^a*y^b)
dPhi_dx = function(x,y) (a*x^(a-1)*y^b)
dPhi_dy = function(x,y) (b*x^a*y^(b-1))


Tx = function (x) (quantile_Q(cdf_P(x)))
Tinvy = function (y) (quantile_P(cdf_Q(y)))

ux = function(x) (integrate(f = function(z) (dPhi_dx(z,Tx(z))),lower = 0,upper = x )$value)
vy = function(y) (integrate(f = function(z) (dPhi_dy(Tinvy(z),z)),lower = Tx(0),upper = y )$value
                  ) - Phi(0,Tx(0))

nbIndiv = 100
ranks = (1:nbIndiv) / nbIndiv
wages = rep(0,nbIndiv)

for (k in 1:nbIndiv)
{wages[k] = ux(quantile_P(ranks[k])) }



num = 0
denom = 0
for (i in 1:nbIndiv)
{for (j in 1:nbIndiv)
{
  num = num+ abs(wages[i] - wages[j])
}
  denom = denom + wages[i]
}

gini = num / (2*nbIndiv*denom)

print(gini)