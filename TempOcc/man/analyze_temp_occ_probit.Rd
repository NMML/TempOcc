\name{analyze_temp_occ_probit}
\alias{analyze_temp_occ_probit}
\title{Main function for running MCMC simulations for temporal occupancy analysis}
\usage{
  analyze_temp_occ_probit(detection.model, process.model,
    Dat, prior, control, q.window = 2)
}
\arguments{
  \item{detection.model}{an object of class 'formula' (as
  in 'lm') for the detection process.  In addition to
  covariates in "Dat" one has access to the keywords "Lag"
  and "Z"; "Lag" denotes using the observation in the
  previous time period as a predictor, while "Z" denotes
  using the current occupancy state as a covariate.}

  \item{process.model}{an object of class 'formula' for
  factors affecting haul-out probability}

  \item{Dat}{a data frame giving the data to be analyzed;
  the first column should be a response variable (0=missing
  data, 1=not occupied, 2=occupied).  Another column should
  be named "Indiv" and include unique individual
  identifiers (numbers are fine).  Remaining columns
  consist of explanatory variables specified by the analyst
  (columns should be named to correspond to given formulas
  for detection and process models).}

  \item{prior}{a list giving prior gamma parameter for tau
  (e.g., list(a.eta=1,b.eta=0.005))}

  \item{control}{a list specifying MCMC options (e.g.,
  list(burnin=100, iter=1000, thin=1)}

  \item{q.window}{an integer specifying the number of
  consecutive hours that are 'replicates' for eta; this
  number must be >1 for parameter identifiability; the
  default is 2}
}
\value{
  returns a list with the following slots: \code{Dat} is
  the Data data frame object, possibly with amended columns
  to represent "Lag" and "Z" effects, \code{Q} is the
  matrix representing the temporal autocorrelation
  structure used in the analysis,\code{DM.Y} returns design
  matrix for detection model for final iteration of MCMC,
  \code{DM.eta} returns the design matrix used to scale up
  to the number of observations (this is needed when
  q.window>1), \code{E.Zt.minus.Ez} returns posterior
  predicitons of residuals needed for making posthoc
  posterior predictions, \code{det} returns posterior
  samples for detection probability parameters,
  \code{process} returns posterior samples for occupancy
  probability parameters, \code{tau.eta} returns posterior
  samples for the precision parameter associated with
  correlated random effects, \code{G.m} Goodness-of-fit
  component for posterior predictive loss \code{P.m}
  Variance component (penalty) for posterior predictive
  loss \code{D.m} Overall posterior predictive loss score
}
\description{
  Main function for running MCMC simulations for temporal
  occupancy analysis
}
\examples{
#warning: the following takes a few minutes
n.indiv=100
n.obs=100
tau=20
P=c(.9,.8,.2,.1)
Dat=sim_data_probit(n.indiv,n.obs,q.window=2,tau,P)
detection.model=~Z*Lag
process.model=~1
prior=list(a.eta=1,b.eta=0.005)
thin=10
control <- list(burnin=100/thin, iter=1100/thin, thin=thin)
q.window=1
out=analyze_temp_occ_probit(detection.model, process.model, Dat, prior, control, q.window)

#compute some posterior predictions
Hr=c(1:24)
Day=10
New.dat=data.frame(Hr=as.factor(Hr),Day=rep(Day,24),Day2=rep(Day^2,24))
process.model=~Hr+Day+Day2
New.X=model.matrix(process.model,New.dat)
#now, matrix multiply New.X and t(out$process) to get posterior predictive distribution at desired design points!
}
\author{
  Paul Conn
}
\keyword{haul-out,}
\keyword{mcmc,}
\keyword{occupancy}
\keyword{temporal}

