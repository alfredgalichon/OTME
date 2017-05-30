source(paste0(getwd(),"/VQR.R"))
inputs<- load_VQRData("Engel1D")
X<-inputs$X;Y<-inputs$Y;n<-inputs$n;d<-inputs$d;r<-inputs$r
step<- 0.1

nu	<- matrix(1/n,n,1)
T	<- seq(0,1,by=step)
m	<- length(T)
U	<- matrix(T,m,1)
mu	<- matrix(1/m, m,1)

sols	<- VQRTp( X,Y,U,mu,nu ) 
pi<-sols$pi; psi<-sols$psi; b<-sols$b; val<-sols$val
beta		<- ComputeBeta1D( mu,b )
x_quantiles1<- quantile(X[,2],probs=c(0,0.25,0.5,0.75,1))
Nb_display 	<- length(x_quantiles1)
x_quantiles <- cbind(matrix(1,Nb_display,1),matrix(x_quantiles1,Nb_display,1))

yhat_VQR	<- x_quantiles%*%t(beta)

xaxis		<- T[2:(m-1)]
xaxis 	<- rbind(T[2:(m-1)],T[2:(m-1)],T[2:(m-1)],T[2:(m-1)],T[2:(m-1)])
plot(xaxis,yhat_VQR[,2:(m-1)],type='n')
for(i in 1:Nb_display){lines(T[2:(m-1)],yhat_VQR[i,2:(m-1)],col=rgb(0,0.2*i,0))}


