# SIMULATE AN ICAR PROCESS (x is the precision matrix = tau*Q)
#' Simulate an ICAR Process
#' @param x a precision matrix = tau*Q
#' @return returns a vector of spatially correlated random effects of length = to the number of rows of x
#' @export
#' @keywords CAR, simulation
#' @author Devin Johnson
rrw <- function(x){
	v <- eigen(x, TRUE)
	val.inv <- sqrt(ifelse(v$values>sqrt(.Machine$double.eps), 1/v$values, 0))
	P <- v$vectors
	sim <- P%*%diag(val.inv)%*%rnorm(dim(x)[1], 0, 1)
	X <- rep(1,length(sim))
	if(sum(val.inv==0)==2) X <- cbind(X, 1:length(sim))  #extra col for centering RW2
	sim <- sim-X%*%solve(crossprod(X), crossprod(X,sim))
	return(sim)
}

#' Formulate an RW1 Q matrix (precision matrix = tau*Q) for one individual
#' @param n number of timesteps an individual is available to be observed
#' @return returns an (n x n) Q matrix
#' @export
#' @keywords CAR, RW1, Q matrix
#' @author Paul Conn
get_Q_1indiv<-function(n){
	require(Matrix)
	Tmp.mat2=diag(n-1)
	Tmp.mat3=cbind(Tmp.mat2,rep(0,n-1))
	Tmp.mat3=rbind(rep(0,n),Tmp.mat3)
	Tmp.mat4=cbind(rep(0,n-1),Tmp.mat2)
	Tmp.mat4=rbind(Tmp.mat4,rep(0,n))
	Tmp.mat=-Tmp.mat3-Tmp.mat4
	diag(Tmp.mat)=-apply(Tmp.mat,1,"sum")
	Matrix(Tmp.mat)		
}

#' Formulate an RW2 Q matrix (precision matrix = tau*Q) for one individual
#' @param n number of timesteps an individual is available to be observed
#' @export
#' @keywords CAR, RW2, Q matrix
#' @author Paul Conn
get_Q_1indiv_RW2<-function(n){
	Tmp.mat1=diag(n-2)
	Tmp.mat1=cbind(Tmp.mat1,rep(0,n-2),rep(0,n-2))
	Tmp.mat1=rbind(rep(0,n),rep(0,n),Tmp.mat1)
	Tmp.mat2=diag(n-2)
	Tmp.mat2=cbind(rep(0,n-2),rep(0,n-2),Tmp.mat2)
	Tmp.mat2=rbind(Tmp.mat2,rep(0,n),rep(0,n))
	Tmp.mat3=diag(n-1)*-4
	Tmp.mat3[1]=-2
	Tmp.mat3[n-1,n-1]=-2
	Tmp.mat4=cbind(Tmp.mat3,rep(0,n-1))
	Tmp.mat4=rbind(rep(0,n),Tmp.mat4)
	Tmp.mat5=cbind(rep(0,n-1),Tmp.mat3)
	Tmp.mat5=rbind(Tmp.mat5,rep(0,n))
	Tmp.mat=Tmp.mat1+Tmp.mat2+Tmp.mat4+Tmp.mat5
	diag(Tmp.mat)=-apply(Tmp.mat,1,"sum")
	Matrix(Tmp.mat)
}

#' Inverse logit function
#' @param x value on the logit scale
#' @return returns logit^(-1) (x)
#' @export
#' @keywords logit
#' @author Paul Conn
expit <-
		function(x){1/(1+exp(-x))}


