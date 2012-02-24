#' Main function for running MCMC simulations for temporal occupancy analysis
#'
#' @param detection.model an object of class 'formula' (as in 'lm') for the detection process.  In addition
#'  to covariates in "Dat" one has access to the keywords "Lag" and "Z"; "Lag" denotes using the observation
#'  in the previous time period as a predictor, while "Z" denotes using the current occupancy state as a covariate.
#' @param process.model   an object of class 'formula' for factors affecting haul-out probability
#' @param Dat   		  a data frame giving the data to be analyzed; the first column should
#'  be a response variable (0=missing data, 1=not occupied, 2=occupied).  Another column should be
#'  named "Indiv" and include unique individual identifiers (numbers are fine).  Remaining columns consist
#'  of explanatory variables specified by the analyst (columns should be named to correspond to 
#'  given formulas for detection and process models).
#' @param prior   a list giving prior gamma parameter for tau (e.g., list(a.eta=1,b.eta=0.005))
#' @param control a list specifying MCMC options (e.g., list(burnin=100, iter=1000, thin=1)
#' @param q.window an integer specifying the number of consecutive hours that are 'replicates' for eta;
#'  this number must be >1 for parameter identifiability; the default is 2
#' @return returns a list with the following slots: \code{Dat} is the Data data frame object, possibly with amended
#'  columns to represent "Lag" and "Z" effects, \code{Q} is the matrix representing the temporal autocorrelation
#'  structure used in the analysis,\code{DM.Y} returns design matrix for detection model for final iteration of MCMC,
#'  \code{DM.eta} returns the design matrix used to scale up to the number of observations
#'  (this is needed when q.window>1), \code{E.Zt.minus.Ez} returns posterior predicitons of residuals needed for making posthoc posterior
#'  predictions, \code{det} returns posterior samples for detection probability parameters,
#'  \code{process} returns posterior samples for occupancy probability parameters, 
#'  \code{tau.eta} returns posterior samples for the precision parameter associated with correlated random effects,
#'  \code{G.m} Goodness-of-fit component for posterior predictive loss
#'  \code{P.m} Variance component (penalty) for posterior predictive loss
#'  \code{D.m} Overall posterior predictive loss score
#' @export
#' @import Matrix
#' @keywords haul-out, mcmc, temporal occupancy
#' @examples 
#' #warning: the following takes a few minutes
#' n.indiv=100
#' n.obs=100
#' tau=20
#' P=c(.9,.8,.2,.1)
#' Dat=sim_data_probit(n.indiv,n.obs,q.window=2,tau,P)
#' detection.model=~Z*Lag
#' process.model=~1
#' prior=list(a.eta=1,b.eta=0.005)
#' thin=10
#' control <- list(burnin=100/thin, iter=1100/thin, thin=thin)
#' q.window=1
#' out=analyze_temp_occ_probit(detection.model, process.model, Dat, prior, control, q.window)
#' 
#' #compute some posterior predictions
#' Hr=c(1:24)
#' Day=10
#' New.dat=data.frame(Hr=as.factor(Hr),Day=rep(Day,24),Day2=rep(Day^2,24))
#' process.model=~Hr+Day+Day2
#' New.X=model.matrix(process.model,New.dat)
#' #now, matrix multiply New.X and t(out$process) to get posterior predictive distribution at desired design points!
#' @author Paul Conn
analyze_temp_occ_probit <-
	function(detection.model, process.model, Dat, prior, control, q.window=2){
	memory.limit(size=8000)
	require(truncnorm)
	require(coda)
	require(Matrix)
	
	Y<-(as.numeric(Dat[,1])>0)		  #y is now 0 (not observed), 1 (observed)
	Z<-(as.numeric(Dat[,1]))
	I.na=which(Z==0)
	n.na=length(I.na)
	Z[I.na]=NA
	Z[-I.na]=Z[-I.na]-1    #z now is 0 (in water), 1 (hauled out), NA (unknown)
	n.obs <- length(Y)
	Lag=c(0,Y[1:(n.obs-1)])  #First detections aren't evaluated so this should work
	
	Dat=as.data.frame(cbind(Y,Z,Dat,Lag))
	rm(Z)
	if("Z"%in%all.vars(detection.model))p.z=1  #indicator for whether true state used as covariate
	#Matrix construction, etc...
	Xy=Matrix(model.matrix(detection.model,Dat))
	Xz=Matrix(model.matrix(process.model,Dat))
	
	Indiv.ids=unique(Dat[,"Indiv"])
	n.indiv = length(Indiv.ids)
	n.y.par=ncol(Xy)
	n.z.par=ncol(Xz)
	Y.par <- rep(0,n.y.par)   
	Z.par <- rep(0,n.z.par)
	get_n_indiv<-function(id,ID)length(which(Dat[,"Indiv"]==id))
	N.obs.ind=unlist(lapply(Indiv.ids,get_n_indiv)) #holds number of hours per individual
	First.obs=N.obs.ind	#condition on first observation so p=1
	First.obs[1]=1
	for(iind in 2:n.indiv)First.obs[iind]=First.obs[iind-1]+N.obs.ind[iind-1]
	n.eta=sum(ceiling(N.obs.ind/q.window))
	
	## construct inverse variance matrix for random walk process as a 
	#  block diagonal 
	id=Indiv.ids[1]
	N.obs.eta=N.obs.ind
	N.obs.eta[1]=ceiling(length(which(Dat[,"Indiv"]==id))/q.window)
	M1=get_Q_1indiv(N.obs.eta[1])
	Mat.list=list(M1=M1)
	if(n.indiv>1){
		for(iind in 2:n.indiv){
			id=Indiv.ids[iind]
			N.obs.eta[iind]=ceiling(length(which(Dat[,"Indiv"]==id))/q.window)
			Tmp=get_Q_1indiv(N.obs.eta[iind])
			eval(parse(text=paste(paste("M",iind,sep=''),'=Tmp')))
			eval(parse(text=paste(paste("M",iind,sep=''),'=list(',paste("M",iind,sep=''),'=',paste("M",iind,sep=''),')')))
			eval(parse(text=paste("Mat.list=c(Mat.list,",paste("M",iind,sep=''),")")))
			eval(parse(text=paste("rm(",paste("M",iind,sep=''),")"))) #clean up
		}
	}
	Q=bdiag(Mat.list)
	rm(Mat.list) #clean up
	
	
	#construct a DM matrix for eta that maps estimated eta's to responses
	DM.eta<-Matrix(0,n.obs,n.eta)
	cur.row=1
	cur.col=1
	for(iind in 1:n.indiv){
		for(ipar in 1:floor(N.obs.ind[iind]/q.window)){
			DM.eta[cur.row:(cur.row+q.window-1),cur.col]=1
			cur.row=cur.row+q.window
			cur.col=cur.col+1
		}
		remainder=N.obs.ind[iind]%%q.window
		if(remainder>0){
			DM.eta[cur.row:(cur.row+remainder-1),cur.col]=1
			cur.row=cur.row+remainder
			cur.col=cur.col+1
		}
	}
	cat(find("t"))
	t.DM.eta=t(DM.eta)
	cross.DM.eta=crossprod(DM.eta,DM.eta)
	
	## Priors
	a.eta <- prior$a.eta  #a.eta, b.eta priors for scalar multiplier tau.z on Q for process
	b.eta <- prior$b.eta
	
	iter <- control$iter*control$thin
	burnin <- control$burnin*control$thin
	storage.iter=control$iter-control$burnin
	
	MCMC.det <- matrix(nrow=storage.iter,ncol=ncol(Xy))
	colnames(MCMC.det)<-colnames(Xy)
	MCMC.proc <- matrix(nrow=storage.iter,ncol=ncol(Xz))
	colnames(MCMC.proc) <- colnames(Xz)
	MCMC.meanz <- matrix(nrow=storage.iter,ncol=1)
	MCMC.Zt.minus.Ez <- matrix(nrow=storage.iter,ncol=n.obs) #keep for posterior predictions
	MCMC.tau.z <- matrix(nrow=storage.iter,ncol=1)
	Pred.sum0 <- rep(0,n.obs)  #holds running sum of 'NA' predictions for calculations of error loss criterion
	Pred.sum1 <- rep(0,n.obs)  #holds running sum of 'wet' predictions 
	colnames(MCMC.tau.z) <- "tau.eta"
	
	V.z.inv <- crossprod(Xz,Xz) 
	A.t <- Matrix(0,n.eta,n.indiv)
	curpl=1
	for(i in 1:n.indiv){
		A.t[curpl:(curpl+N.obs.eta[i]-1),i]=1
		curpl=curpl+N.obs.eta[i]
	}
	A=t(A.t)
	Zeroes=rep(0,n.obs)
	I.n=Diagonal(n.eta)
	n.obs.tau=n.eta-n.indiv #there is one less degree of freedom due to centering

	Eta.par <- rnorm(n.eta,0,.1)  #initial random walk effects on process		
	Eta=DM.eta%*%Eta.par

	if(p.z==1){
	  Dat.temp=Dat
	  Dat.temp[,"Z"]=rep(0,n.obs)
	  Xy0=Matrix(model.matrix(detection.model,Dat.temp))
	  Dat.temp[,"Z"]=rep(1,n.obs)
	  Xy1=Matrix(model.matrix(detection.model,Dat.temp))
    }
		
	cat("\nBeginning MCMC routine ...\n")
	st <- Sys.time()
	
	for(i in 1:iter){
		#cat(iter)
		Ez.q <- as.numeric(Xz%*%Z.par+Eta) #expectation on probit scale
		#print(mean(Ez.q[crap]))
		
		if(p.z==1){
		
			Ey0.q <- as.numeric(Xy0%*%Y.par) #probit scale predictions when z=0
			Ey1.q <- as.numeric(Xy1%*%Y.par) #probit scale predictions when z=1
			
			#Update missing states using bayes rule
			P.zeq1 <- pnorm(Zeroes, Ez.q, 1, lower.tail=FALSE)  
			P.yeq0.zeq0 <- pnorm(Zeroes, Ey0.q, 1)
			P.yeq0.zeq1 <- pnorm(Zeroes, Ey1.q, 1)
			
			P.z=P.yeq0.zeq1*P.zeq1
			P.z=P.z/(P.z+P.yeq0.zeq0*(1-P.zeq1))
		}
		else P.z=pnorm(Zeroes,Ez.q,1,lower.tail=FALSE)
		
		Dat[I.na,"Z"] <- rbinom(n.na, 1, P.z[I.na])
		
		#Update y.tilde
	    if(p.z==1)Xy=Matrix(model.matrix(detection.model,Dat)) #DM is dynamic if Z included in formula
		Ey.q<-as.numeric(Xy%*%Y.par)
		Y.tilde <- rtruncnorm(n.obs-n.indiv, a=ifelse(Y[-First.obs]==0,-Inf,0), b=ifelse(Y[-First.obs]==0,0,Inf), Ey.q[-First.obs], 1)
		
		#Update z.tilde
		Z.tilde <- rtruncnorm(n.obs, a=ifelse(Dat[,"Z"]==0,-Inf,0), b=ifelse(Dat[,"Z"]==0,0,Inf), Ez.q, 1)
		
		#Update detection parameters
		V.y.inv <- crossprod(Xy[-First.obs,],Xy[-First.obs,]) 
		M.y <- solve(V.y.inv, crossprod(Xy[-First.obs,],Y.tilde))
		Y.par <- M.y + solve(chol(V.y.inv), rnorm(n.y.par,0,1))
		
		#Update process parameters
		M.z <- solve(V.z.inv, crossprod(Xz,Z.tilde-Eta))#+crossprod(Q.g,mu.g))
		Z.par <- M.z + solve(chol(V.z.inv), rnorm(n.z.par,0,1))
		
		#Update tau
		tau.z <- rgamma(1, n.obs.tau/2 + a.eta, as.numeric(crossprod(Eta.par, Q %*% Eta.par)/2) + b.eta)
		
		#Update Eta
		Zt.minus.Ez=Z.tilde-Xz%*%Z.par
		V.eta.inv <- cross.DM.eta + tau.z*Q
		M.eta <- solve(V.eta.inv, t.DM.eta%*%Zt.minus.Ez)
		Eta.par <- M.eta + solve(chol(V.eta.inv), rnorm(n.eta,0,1))
		#center using eq 2.30 of Rue and Held		
		Eta.par=Eta.par-V.eta.inv %*% A.t %*% solve(A %*% V.eta.inv %*% A.t,A%*%Eta.par)
		Eta=DM.eta%*%Eta.par
				
		if(i>burnin & i%%control$thin==0){
			MCMC.det[(i-burnin)/control$thin,] <- as.numeric(Y.par)
			MCMC.proc[(i-burnin)/control$thin,] <- as.numeric(Z.par)
			MCMC.Zt.minus.Ez[(i-burnin)/control$thin,] <- as.numeric(Zt.minus.Ez)
			MCMC.tau.z[(i-burnin)/control$thin] <- as.numeric(tau.z)
			#posterior predictive loss stuff
			z.rep <- (rnorm(n.obs,as.numeric(Xz%*%Z.par+Eta),1)>0)
			y.rep<-(rnorm(n.obs,as.numeric(Xy%*%Y.par),1)>0)
			y.rep=y.rep*(z.rep+1)			
			Pred.sum0=Pred.sum0+(y.rep==0)
			Pred.sum1=Pred.sum1+(y.rep==1)
		}
		if(i==15){
			tpi <- as.numeric(difftime(Sys.time(), st, units="secs"))/15
			ttc <- round((iter-15)*tpi/3600, 2)
			cat("\nApproximate time till completion: ", ttc, " hours\n")
		}
		if(100*(i/iter) >= 10 & (100*(i/iter))%%10==0) cat("\n", 100*(i/iter), "% completed\n")
	}	

	#calculate posterior predictive loss statistic - don't use first observations since these are conditioned on
	y.obs.0=(Dat[,1]==0)
	y.obs.1=(Dat[,1]==1)
	y.obs.2=(Dat[,1]==2)
	p0.tilde=Pred.sum0/storage.iter
	p1.tilde=Pred.sum1/storage.iter
	p2.tilde=1-p0.tilde-p1.tilde
	G.m=sum((p0.tilde[-First.obs]-y.obs.0[-First.obs])^2)+sum((p1.tilde[-First.obs]-y.obs.1[-First.obs])^2)+sum((p2.tilde[-First.obs]-y.obs.2[-First.obs])^2)
	P.m=sum(p0.tilde[-First.obs]*(1-p0.tilde[-First.obs]))+sum(p1.tilde[-First.obs]*(1-p1.tilde[-First.obs]))+sum(p2.tilde[-First.obs]*(1-p2.tilde[-First.obs]))
	D.m=G.m+P.m
	
	MCMC.det <- mcmc(MCMC.det)
	MCMC.proc <- mcmc(MCMC.proc)
	MCMC.Zt.minus.Ez
	tau.eta <- mcmc(MCMC.tau.z)
	out <- list(Dat=Dat,Q=Q,DM.eta=DM.eta,DM.Y=Xy,E.Zt.minus.Ez=MCMC.Zt.minus.Ez,det=MCMC.det,process=MCMC.proc,tau.eta=tau.eta,G.m=G.m,P.m=P.m,D.m=D.m)
	return(out)
}


