# Filename: run_HO.R
#
# 
# Author: paul.conn Jul 7, 2011
###############################################################################

require(coda)
#source('c:/users/paul.conn/r_eclipse_work/HO_spotted/R/sim_data_obs_process.r')
#load("c:/users/paul.conn/r_eclipse_work/HOspotted/HOspotted/data/Spotted_Data_Stacked.Rdata")

n.indiv=100
n.obs=100
q.window=1
tau=20
P=c(.9,.8,.2,.1)
Dat=sim_data_probit(n.indiv,n.obs,q.window=2,tau,P)
	
detection.model=~Z*Lag
process.model=~1
prior=list(a.eta=1,b.eta=0.005)
thin=10
control <- list(burnin=100/thin, iter=1100/thin, thin=thin)
q.window=1

#run analysis
out=analyze_temp_occ_probit(detection.model, process.model, Dat, prior, control, q.window)

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
Psi=post_predict_haulout(Tau=out$tau,BetaZ=out$process,Q=out$Q,DM.eta=out$DM.eta,Pred.X=New.X,Zt.minus.Ez=out$E.Zt.minus.Ez,n=100)



par(mfrow=c(3,2))
traceplot(out1$p00)
traceplot(out1$p01)
traceplot(out1$p10)
traceplot(out1$p11)
traceplot(out1$psi)
traceplot(out1$tau.eta)

