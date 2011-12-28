#' Simulate temporal occupancy data
#' @param n.indiv integer giving number of sites (individuals)
#' @param n.obs integer giving number of observation per site (assumed equal for every individual)
#' @param q.window is an integer giving the number of replicates for each eta (must be >1; default is 2)
#' @param tau gives the precision (1/variance) of the CAR process
#' @param P is a vector with 4 detection probabilities: (1) state_t=not occupied, detected at t-1,
#'  (2) state_t = occupied, detected at t-1, (3) state_t=not occupied, not detected at t-1, and 
#'  (4) state_t = occupied, not detected at t-1
#' @return Stacked data.frame object with "Response" as the first column (0=missing, 1=Not occupied, 2=Occupied), "Hr" (timestep) as second column, and "Ind" (individual identifier) as third column
#' @export
#' @keywords occupancy, probit model, simulation
#' @author Paul Conn
sim_data_probit<-function(n.indiv,n.obs,q.window=2,tau,P){
	
	max.t=n.obs+1  # +1 needed because time starts at t=1 for indexing simplicity
	sd=sqrt(1/tau)
	
	
	p1=P[1]  #detection probability | state=0, obs_t-1>0 (in water)
	p2=P[2]  #detection probability | state=1, obs_t-1>0 (hauled out) 
	p3=P[3]	#detection probability | state=0, obs_t-1=0
	p4=P[4]  #detectoin probability | state=1, obs_t-1=0
	
	Response=matrix(0,n.indiv,max.t-1)
	
	Q=get_Q_1indiv((max.t-1)/q.window)
	
	#construct a DM matrix for eta that maps estimated eta's to responses
	n.eta=n.obs/q.window
	DM.eta<-matrix(0,n.obs,n.eta)
	cur.row=1
	cur.col=1
	for(ipar in 1:floor(n.obs/q.window)){
		DM.eta[cur.row:(cur.row+q.window-1),cur.col]=1
		cur.row=cur.row+q.window
		cur.col=cur.col+1
	}
	remainder=n.obs%%q.window
	if(remainder>0){
		DM.eta[cur.row:(cur.row+remainder-1),cur.col]=1
		cur.row=cur.row+remainder
		cur.col=cur.col+1
	}
	
	## perform for each animal
	#initial state
	for(iind in 1:n.indiv){
		Eta.par=rrw(tau*Q)
		Eta=DM.eta%*%Eta.par
		for(ihr in 1:n.obs){
			Response[iind,ihr]=(Eta[ihr]+rnorm(1,0,1))>0
		}
	}
	
	
	###observation model
	Y=round(Response+1)
	Response=round(Response)
	for(iind in 1:n.indiv){
		for(ihr in 2:(max.t-1)){
			if(Y[iind,ihr-1]>0 & Response[iind,ihr]==0)Y[iind,ihr]=rbinom(1,1,p1)*Y[iind,ihr]
			if(Y[iind,ihr-1]>0 & Response[iind,ihr]==1)Y[iind,ihr]=rbinom(1,1,p2)*Y[iind,ihr]
			if(Y[iind,ihr-1]==0 & Response[iind,ihr]==0)Y[iind,ihr]=rbinom(1,1,p3)*Y[iind,ihr]
			if(Y[iind,ihr-1]==0 & Response[iind,ihr]==1)Y[iind,ihr]=rbinom(1,1,p4)*Y[iind,ihr]		
		}
	}
	
	### post processing to stack data
	#Final data in form:  Response, hour, individual
	Stacked=matrix(0,1,3)
	First.ti=rep(0,n.indiv)
	Last.ti=First.ti
	for(i in 1:n.indiv){
		Cur.dat=Y[i,]
		I.pos=which(Cur.dat>0)
		First.ti[i]=min(I.pos)
		Last.ti[i]=max(I.pos)
		Stacked=rbind(Stacked,cbind(Y[i,First.ti[i]:Last.ti[i]],First.ti[i]:Last.ti[i],rep(i,Last.ti[i]-First.ti[i]+1)))
	}
	Stacked=Stacked[-1,]
	colnames(Stacked)=c("Response","Hr","Indiv")
	Out=as.data.frame(Stacked)

	return(Out)
}
