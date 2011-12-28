# Filename: run_HO.R
#
# 
# Author: paul.conn Jul 7, 2011
###############################################################################
load("z:/conn/haulout/HO_temp_occ.Rdata")
Dat=cbind(HO2[,5],HO2[,c(1,2,3,4,7,8,9,10)]) #put 'Dry' first
Dat[which(Dat[,2]=="B"),2]="S"
Dat[,2]=factor(Dat[,2])
Dat=cbind(Dat,Dat[,"Day.since.March1"]^2)
Dat[,"Hr"]=as.factor(Dat[,"Hr"])
levels(Dat[,2])=c("R","S")
Dat[,"Sex"]=as.factor(as.character(Dat[,"Sex"]))
colnames(Dat)[1]="Dry"
colnames(Dat)[8]="Day"
colnames(Dat)[10]="Day2"

rm(HO2)
detection.model=~Instr*Lag*Z  #"Lag" and "Z" are keywords - "Lag" denotes observation previous time step, "Z" denotes occupancy status
process.model=~(Species+Sex+Age+Day)^2+Day2+Day2:Sex+Day2:Age+Day2:Species+Hr
prior=list(a.eta=1,b.eta=0.005)
thin=10
control <- list(burnin=100/thin, iter=1100/thin, thin=thin)
p.z=1
q.window=2

#run analysis
out=analyze_temp_occ_probit(detection.model, process.model, Dat, prior, control, q.window)

save(out,"out_InstrLagZ.Rdata")

