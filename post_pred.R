Hr=c(1:24)
Day=c(0:15)
New.dat=data.frame(Hr=as.factor(Hr),Day=Day,Day2=Day^2)
process.model=~Hr+Day+Day2
New.X=model.matrix(process.model,New.dat)
cur.hr=1

New.X=matrix(0,16,26)
New.X[,1]=1
New.X[,12]=1
New.X[,25]=c(0:15)
New.X[,26]=(c(0:15))^2

#calculate haulout predictions and var-cov matrix
source('c:/users/paul.conn/r_eclipse_work/HO_spotted/R/post_predict_haulout.r')
Psi=post_predict(Tau=out$tau,BetaZ=out$process,Q=out$Q,DM.eta=out$DM.eta,Pred.X=New.X,Zt.minus.Ez=out$E.Zt.minus.Ez,n=100)

if(p.z==1){
	MCMC.p00 <- pnorm(0,matrix(c(1,0,0,0),1,n.y.par)%*%as.matrix(Y.par),1,lower.tail=FALSE)
	MCMC.p01 <- pnorm(0,matrix(c(1,0,1,0),1,n.y.par)%*%as.matrix(Y.par),1,lower.tail=FALSE)
	MCMC.p10 <- pnorm(0,matrix(c(1,1,0,0),1,n.y.par)%*%as.matrix(Y.par),1,lower.tail=FALSE)
	MCMC.p11 <- pnorm(0,matrix(c(1,1,1,1),1,n.y.par)%*%as.matrix(Y.par),1,lower.tail=FALSE)
}
else{
	MCMC.p00 <- pnorm(0,matrix(c(1,0),1,n.y.par)%*%as.matrix(Y.par),1,lower.tail=FALSE)
	MCMC.p01 <- pnorm(0,matrix(c(1,1),1,n.y.par)%*%as.matrix(Y.par),1,lower.tail=FALSE)
	MCMC.p10 <- pnorm(0,matrix(c(1,0),1,n.y.par)%*%as.matrix(Y.par),1,lower.tail=FALSE)
	MCMC.p11 <- pnorm(0,matrix(c(1,1),1,n.y.par)%*%as.matrix(Y.par),1,lower.tail=FALSE)
}	



par(mfrow=c(3,2))
traceplot(out1$p00)
traceplot(out1$p01)
traceplot(out1$p10)
traceplot(out1$p11)
traceplot(out1$psi)
traceplot(out1$tau.eta)

