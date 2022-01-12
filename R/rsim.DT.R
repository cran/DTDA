


rsim.DT<-function(n, tau, model=NULL){

if(model==1){


X <- runif(n, 0, 1)
U <- runif(n,-tau, 1)
V <- U + tau
for (i in 1:n){

	while (U[i] > X[i] | V[i] < X[i]){
	X[i] <- runif(1, 0, 1)
	U[i] <- runif(1, -tau, 1)
	V[i] <- U[i] + tau
	}

}


T<-table(U<= X&X<=V)
if(sum(T[names(T)=="TRUE"])!=length(X) |sum(T[names(T)=="TRUE"]==0)){
  stop("Condition of double truncation is violated","\n")
}

}
if (model==2){


X <- runif(n, 0, 1)
U <- (runif(n, 0, 1)^2) * (1 + tau) - tau
V <- U + tau
for (i in 1:n){
while (U[i] > X[i] | V[i] < X[i]){
X[i] <- runif(1, 0, 1)
U[i] <- (runif(1, 0, 1)^2) * (1 + tau) - tau
V[i] <- U[i] + tau
}
}
T<-table(U<= X&X<=V)
if(sum(T[names(T)=="TRUE"])!=length(X) |sum(T[names(T)=="TRUE"]==0)){
  stop("Condition of double truncation is violated","\n")
}



}

return(invisible(as.data.frame(cbind(X,U,V))))
}

















