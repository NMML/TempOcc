pkgname <- "TempOcc"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('TempOcc')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("analyze_temp_occ_probit")
### * analyze_temp_occ_probit

flush(stderr()); flush(stdout())

### Name: analyze_temp_occ_probit
### Title: Main function for running MCMC simulations for temporal
###   occupancy analysis
### Aliases: analyze_temp_occ_probit
### Keywords: haul-out, mcmc, occupancy temporal

### ** Examples

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



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
