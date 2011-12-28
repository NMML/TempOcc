#' Plot temporal occupancy data for inidivual "sites" using ggplot2
#' @param Dat matrix object with number of rows equal to the number of sites (individuals) and the number of columns equal to the number of occasions
#'  (analagous to an encounter history).  NAs give missing, 0 gives not hauled out, 1 gives hauled out
#' @param Time	a vector giving times that each occasion is associated with; for evenly spaced occasions this could simply be c(1:ncol(Dat))
#' @param Ind	a vector of individual IDs (can be a character or integer vector)
#' @param plot.ncol   an integer specifying the number of columns desired in the figure
#' @param DayMoYear	  a vector giving the day, month, and year, of each element in Time (formatted using the date.mmddyy function in the 'date' package)
#' @param Label.times  a vector giving the elements of DayMoYear to plot as labels on the graph
#' @param First	an integer vector giving the indices of the first observations of animals in the study; assumed to be in same order as "Ind"
#' @param Last  an integer vector givign the indices of the last observations of animals
#' @param min.bar.length relative number of hours each haulout "block" takes up on the graph
#' @export
#' @keywords plot, temporal occupancy
#' @author Paul Conn
plot_temp_occ <- function(Dat,Time,Ind,plot.ncol=1,DayMoYear,Label.times,First,Last,min.bar.length=1){
	require(ggplot2)
	require(date)
	#Dat[which(Dat[,1]==0),1]=NA
	#Dat[,1]=Dat[,1]+0.1	 	#relative length of rectangles (NAs are blank)
	Dat=(Dat+0.1)/1.1	 	#relative length of rectangles (NAs are blank)
	Dat[which(is.na(Dat)==1)]=0
	
	Active=rep(0,length(Time)*length(Ind))
	for(i in 1:length(Ind)){
		Active[((i-1)*length(Time)+First[i]):((i-1)*length(Time)+Last[i])]=.5
	}
	#turn into data frame for ggplot2
	H.df=data.frame(
		indiv=rep(Ind,each=length(Time)),
		date=rep(Time,length(Ind)),
		Active=Active,
		H=as.vector(t(Dat))
	)
	
	H.df2=data.frame(
		indiv=Ind,
		First=First,
		Last=Last
	)
	crap=H.df[which(H.df[,"indiv"]==15),]
	myplot=ggplot(H.df2, aes(xmin=First,xmax=Last,ymin=-.5,ymax=.5)) + geom_rect(fill="gray")
	myplot=myplot+geom_rect(data=H.df,aes(xmin=date,xmax=date+min.bar.length,ymin=-H,ymax=H),fill="black")  #5 used here because length=1 doesn't show up!
	
	myplot=myplot + facet_wrap(~ indiv,ncol=plot.ncol)
	myplot=myplot + scale_x_continuous(breaks=Label.times,labels=DayMoYear[Label.times]) 
	myplot=myplot + scale_y_continuous(breaks=NA) + labs(x="Time",y=NULL)
	
	myplot
}


 

