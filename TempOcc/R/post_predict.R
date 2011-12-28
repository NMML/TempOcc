#' Generate posterior predictions for occupancy at prespecified sets of covariates
#'  after MCMC is completed
#' @param Tau a real vector of posterior samples for tau (as obtained via MCMC)
#' @param BetaZ a matrix giving posterior samples for occupancy regression coefficients, Beta.  Each column
#'  references a separate beta value.
#' @param Pred.X a design matrix giving the coefficient values for which one wants predictions (one prediction would 
#'  correspond to a row vector
#' @param Q	 an inverse variance component specifying the temporal structure associated with the occupancy model
#' @param DM.eta	the design matrix applied to the Eta random effects to scale up to the level of the data.  This depends on the "window" assumed in the estimation model.
#' @param Zt.minus.Ez	Posterior predictions of residuals (these are outputs from the main MCMC simulation)
#' @param n integer giving desired number of Monte Carlo replicates for posterior predictions.  Defaults to 1000
#' @return Returns a list composed of the following: \code{Psi.mean}: a vector of posterior
#'  mean predictions for occupancy, \code{Psi.median} a vector of posterior median predictions, 
#'  \code{Psi.VC} an empirical variance-covariance matrix for posterior predictions, \code{Psi.lower} a vector of lower
#'  90\% confidence limits for predictions,and \code{Psi.upper} a vector of upper 90\% limits for occupancy predictions
#' @export
#' @keywords CAR, RW1, Q matrix
#' @author Paul Conn
post_predict<-function(Tau,BetaZ,Pred.X,Q,DM.eta,Zt.minus.Ez,n=1000){
	memory.limit(size=8000)
	require(MASS)
	n.eta=nrow(Q)
	n.preds=nrow(Pred.X)
	n.mcmc=length(Tau)
	cross.DM.eta=crossprod(DM.eta,DM.eta)
	t.DM.eta=t(DM.eta)
	
	Pred=matrix(0,n.preds,n) 
	for(i in 1:n){
		print(i)
		i.samp=sample(c(1:n.mcmc),1)
		V.eta.inv <- Matrix(cross.DM.eta + Tau[i.samp]*Q)
		M.eta <- solve(V.eta.inv, t.DM.eta%*%Zt.minus.Ez[i.samp,])
		Eta.par <- M.eta + solve(chol(V.eta.inv), rnorm(n.eta,0,1))
		Eta=sample(Eta.par,n.preds,replace=1)
		for(j in 1:n.preds){
			Pred[j,i]=1-pnorm(0,Pred.X[j,]%*%BetaZ[i.samp,]+Eta[j],1)
		}
	}
	Psi.mean=apply(Pred,1,'mean')
	Psi.median=apply(Pred,1,'median')
	quant95<-function(X)quantile(X,probs=0.95)
	quant05<-function(X)quantile(X,probs=0.05)
	Psi.lower=apply(Pred,1,'quant05')
	Psi.upper=apply(Pred,1,'quant95')
	Sigma=matrix(0,n.preds,n.preds)
	diag(Sigma)=apply(Pred,1,'var')
	if(n.preds>1){
	  for(i in 2:n.preds){
		  for(j in 1:(i-1)){
			  Sigma[i,j]=cov(Pred[i,],Pred[j,])
		  }
	  }
  	}
	out=list(Psi.mean=Psi.mean,Psi.median=Psi.median,Psi.VC=Sigma,Psi.lower=Psi.lower,Psi.upper=Psi.upper)
	return(out)
}

