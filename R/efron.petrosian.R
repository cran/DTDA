efron.petrosian <-
	function(X, U=NA, V=NA, wt=NA, error=NA,
	nmaxit=NA , boot=TRUE, B=NA, alpha=NA,
		display.F=FALSE, display.S=FALSE){

 # set truncation both side by default
    trunc <- "both"

   # analize if we have a NA and remove then
    if (all(is.na(U)) == TRUE & !all(is.na(V)) == TRUE) {
      trunc <- "right"
      cat("case U=NA","\n")
	  cat("warning: data on one truncation limit are missing; the result is based on an iterative algorithm, but an explicit-form NPMLE (Lynden-Bell estimator) exists","\n")


      if (any(is.na(V)) == TRUE | any(is.na(X)) == TRUE) {
        navec <- c(which(is.na(X)), which(is.na(V)))
        X <- X[-navec]
        V <- V[-navec]
      }
    }



	 # check for NA in the resulting vector
    if (all(is.na(V)) == TRUE & !all(is.na(U)) == TRUE) {
      trunc <- "left"
      cat("case V=NA","\n")
cat("warning: data on one truncation limit are missing; the result is based on an iterative algorithm, but an explicit-form NPMLE (Lynden-Bell estimator) exists","\n")

      if (any(is.na(U)) == TRUE | any(is.na(X)) == TRUE) {
        navec <- c(which(is.na(X)), which(is.na(U)))
        X <- X[-navec]
        U <- U[-navec]
      }
    }
	
if (all(is.na(V)) == TRUE & all(is.na(U)) == TRUE) {
cat("case U=NA and V=NA","\n")
stop("warning: data on at least one truncation limit (left or right) is required","\n")
}	
	
	# if trunccation is set both sides prepare the arrays
    if (trunc == "both") {
      if (any(is.na(U)) == TRUE |
          any(is.na(V)) == TRUE | any(is.na(X)) == TRUE) {
        navec <- c(which(is.na(X)), which(is.na(U)), which(is.na(V)))
        X <- X[-navec]
        U <- U[-navec]
        V <- V[-navec]
      }
    }

# check for some na in wight matrices
    if (is.na(wt) == TRUE)
      wt <- rep(1 / length(X), times = length(X))

    D <- cbind(X, U, V) # concatenate them into D

    if (all(is.na(V)) == TRUE)
      D[, 3] <- rep(max(D[, 1]) + 1, length(X))

    if (all(is.na(U)) == TRUE) {
      D[, 2] <- rep(min(D[, 1]) - 1, length(X))
      D[, 1] <- D[, 1]
      D[, 3] <- D[, 3]
    }


	C<- matrix(0, nrow = nrow(D), ncol = ncol(D))
    EE <- matrix(0, nrow = nrow(D), ncol = ncol(D))

    ########################################
    # ordenaçao or ordering zone
    ########################################

      ord <- order(D[, 1], method = "auto")
      C[, 1] <- sort(D[, 1], method = "auto")
      C[, 2:ncol(D)] <-
      D[ord, 2:ncol(D)] # use the ordered indeces for ordering

	T<-table(C[,2]<= C[,1]&C[,1]<=C[,3])
    if(sum(T[names(T)=="TRUE"])!=length(X) |sum(T[names(T)=="TRUE"]==0)){
      stop("Condition of double truncation is violated","\n")
    }

    if (is.na(error) == TRUE)
      error <- 1e-6


	au <- outer(C[, 1], C[, 2], ">=")
    av <- outer(C[, 1], C[, 3], "<=")
    auu <- outer(C[, 2], C[, 2], "<=") * 1L

    J <- au * av
    SC <- colSums(J)
    ad <- length(which(SC == 1))
        if (ad != 0) {
      #if (my_debug)
        cat("Warning. Non-uniqueness or no existence of the NPMLE", "\n")
    }

# now we proceed after good sample
    JI <- t(J)
    f0 <- matrix(data = 1 / nrow(C),
                 ncol = 1,
                 nrow = nrow(C))

    f <- f0
    S0 <- 1
    if (is.na(nmaxit) == TRUE)
      nmaxit <- 100000

    iter <- 0
    while (S0 > error | iter > nmaxit) {
      iter <- iter + 1
      if (iter > nmaxit)
        stop("Default number of iterations not enough for convergence")

F0<-JI%*%f
IF0<-1/F0
If1<-J%*%IF0
f<-1/If1
if(sum(f)!=1)f<-f/sum(f)

S0<-max(abs(f-f0))
f0<-f
}

F0<-JI%*%f
mult <- tabulate(match(C[,1],unique(C[,1])))
if(sum(mult)==length(unique(C[,1]))){
Fval <- (f*mult)}
if(sum(mult)>length(unique(C[,1]))){
weigth<-f[!duplicated(C[,1])]
Fval<- (weigth*mult)}

x<-unique(C[,1])
events<-sum(mult)
n.event<-mult
#f<-Fval
FF<-cumsum(Fval)
FFF<-cumsum(f)
Sob<-1-FF+Fval
Sob[Sob < 1e-12] <- 0
    Sob0<-1-FFF+f
	Sob0[Sob0<1e-12]<-0


if (boot==TRUE){
if (is.na(B)==TRUE) B<-500

B <- 500 # this is harcoded


        if (B < 40) {
          cat("Warning. Number of replicates less than 40", "\n")
          cat("Confidence bands cannot be computed", "\n")
        }


if(trunc=="both"|trunc=="left"){
        ########################################


        ########################################
        # BOOSTRAP FORK FOR BITH SIDE TRUNCATION
        ########################################
########################################
          # PARALELL BOOTSTARP BOTH
          ########################################
  if(!requireNamespace("foreach")){
    install.packages("foreach")

  }

  if(!requireNamespace("doParallel")){
    install.packages("doParallel")

  }


			cores=detectCores()
			if (cores[1] > 1) cl <- makeCluster(cores[1]-1, type = "PSOCK")
			else cl <- makeCluster(cores[1], type = "PSOCK")
			registerDoParallel(cl)
          ########################################
          # Start of parallel state both 500 iteractions B = 500
          final_boot <- foreach(i=1:B, .combine='rbind') %dopar% {
            # dont forget to put code repeat

            ## Preallocation Local to the core (no B indexing)
            #preallocate out internal matrices and variables (each core must have its own copy) no B

            n_sampling_tries <- 0
            # container for bootstrap operations
			M1b <- matrix(0, nrow = nrow(C), ncol = ncol(C)) # usado dentro do loop
            M2b <- matrix(0, nrow = nrow(C), ncol = ncol(C)) # usado dentro do loop
			M_IF0<-matrix(0,nrow=nrow(C))
			M_IF01<-matrix(0,,nrow=nrow(C))
			M_IF0Sob<-matrix(0,nrow=nrow(C))
			ind<-seq(1,nrow(C),by=1)
			indbb<-seq(1,nrow(C),by=1)
			indbb1<-seq(1,nrow(C),by=1)



			# repeat sample if condition not met
            repeat{
              # get the sample
              indb <- sample(ind, nrow(C), replace = TRUE)
              M1b <- C[indb,]

              ########################################
              # ordenaçao or ordering zone
              ########################################
              ord <- order(M1b[, 1], method = "auto")
              M2b[, 1] <- sort(M1b[, 1], method = "auto")

              M2b[, 2:ncol(M1b)] <- M1b[ord, 2:ncol(M1b)]

              # preallocate
			  Aun<-unique(M2b)
				Jb <- matrix(data = 0,ncol = nrow(M2b), nrow = nrow(M2b))

              aub <- outer(M2b[, 1], M2b[, 2], ">=")
			  aubun <- outer(Aun[, 1], Aun[, 2], ">=")
              avb <- outer(M2b[, 1], M2b[, 3], "<=")
			  avbun<-outer(Aun[, 1], Aun[, 3], "<=")
              auub <- outer(M2b[, 2], M2b[, 2], "<=") * 1L
              Jb <- aub * avb
			  Jbun<-aubun*avbun
              SCb <-colSums(Jbun)
				SCb2<-rowSums(Jbun)
              adb <- length(which(SCb == 1))
			  adbb <- length(which(SCb2 == 1))

              n_sampling_tries <- n_sampling_tries + 1 # counter for number of tries

            # if the condition is set, we break the repat, otherwise repeat sampling  process
            if (adb == 0 & adbb==0)
              break
            else
              cat("Warning. Non-uniqueness or no existence of the NPMLE - New Round","\n")
            }
            ## End os sample repeat

            # now we proceed after good sample
            JIb <- t(Jb)
			# preallocation
            f0b <- matrix(data = wt, ncol = 1, nrow = nrow(M2b))
            f1b <- f0b

            S0b <- 1
            iterb <- 0

            # some inner loop
            while (S0b > error | iterb > nmaxit) {
              iterb <- iterb + 1

              if (iterb > nmaxit)
                stop("Default number of iterations not enough for convergence")


				F0b<-JIb%*%f1b
				IF0b<-1/F0b
				If1b<-Jb%*%IF0b
				f1b<-1/If1b
				if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
				S0b<-max(abs(f1b-f0b))
				f0b<-f1b
			}


			ff0b<-numeric(nrow(C))
			for(i in 1:nrow(C)){
			indbb1<-(C[,1]==C[i,1])
			pos1<-min(which(indbb1==TRUE))
			if(pos1==1){
			ff0b[indbb1]<-sum(f1b[indbb1])
			}
			if(pos1>1){
			ff0b[indbb1]<-sum(f1b[indbb1])
			}
			}


			FF0b<-numeric(nrow(C))
			for(i in 1:nrow(C)){
			indbb<-(C[,1]==C[i,1])
			pos<-min(which(indbb==TRUE))
			if(pos==1){
			FF0b[indbb]<-sum(f1b[indbb])}
			if(pos>1){
			FF0b[indbb]<-sum(f1b[indbb])+FF0b[pos-1]
			}
			}

			Sobb <- 1 - FF0b + ff0b
            Sobb[Sobb < 1e-12] <- 0


			# SAVE start saving stufs (before indexed to b index)
            M_IF0 <- as.vector(FF0b)
            M_IF01 <- as.vector(f1b)
            M_IF0Sob <- as.vector(Sobb)

            ########################################
            # ordenaçao or ordering zone
            ########################################
            ordA <- order(M2b[, 2], method = "auto")
            ord3 <- sort(M2b[, 2], method = "auto")

            ordB <- order(M2b[, 3], method = "auto")
            ord4 <- sort(M2b[, 3], method = "auto")

            # end ordenaçao or ordering zone
            ########################################

			return(list( M_IF0, M_IF01, M_IF0Sob, n_sampling_tries))

          # End of parallel state both
          ########################################
          }
          #stop cluster
          stopCluster(cl)

          # computation of the B aggegated results r comes as a list concatenated as row (thsi enables any size lists)

          # preallocate the original ones as matrix a,d put the list indexing
          ############################################
          n_sampling_tries <- 0
          M_IF0 <- matrix(0, nrow = nrow(C), ncol = B)
          M_IF01 <- matrix(0, nrow = nrow(C), ncol = B)
          M_IF0Sob <- matrix(0, nrow = nrow(C), ncol = B)
          stderror<-matrix(0, nrow = nrow(C), ncol = 1)
		  BootRepeat<-matrix(0,ncol=B)
          ###########################################

			###########################################

          	for(b in 1:B){
			M_IF0[,b]<-as.numeric(unlist(final_boot[b,1]))
			M_IF01[,b]<-as.numeric(unlist(final_boot[b,2]))
			M_IF0Sob[,b]<-as.numeric(unlist(final_boot[b,3]))
			BootRepeat[,b]<-as.numeric((unlist(final_boot[b,4])))
			}

		  stderror<-apply(M_IF0,1,sd)
		  M_IF0_sort<-matrix(0,nrow=nrow(C),ncol=B)
			for(i in 1:nrow(M_IF0_sort)){
			M_IF0_sort[i,]<-sort(M_IF0[i,])
			}
		if(is.na(alpha)==TRUE) alpha<-0.05
		lowerF<-M_IF0_sort[,floor(alpha*B/2)]
		upperF<-M_IF0_sort[,floor((1-alpha/2)*B)]

		M_IF0_sort1<-matrix(0,nrow=nrow(C),ncol=B)
		for(i in 1:nrow(M_IF0_sort1)){
		M_IF0_sort1[i,]<-sort(M_IF0Sob[i,])
		}
		lowerS<-M_IF0_sort1[,floor(alpha*B/2)]
		upperS<-M_IF0_sort1[,floor((1-alpha/2)*B)]


			} ##end trunc both and left


		if(trunc=="right"){

		########################################
          # PARALELL BOOTSTARP BOTH
          ########################################
		  if(!requireNamespace("foreach")){
		    install.packages("foreach")

		  }

		  if(!requireNamespace("doParallel")){
		    install.packages("doParallel")

		  }


			cores=detectCores()
			if (cores[1] > 1) cl <- makeCluster(cores[1]-1, type = "PSOCK")
			else cl <- makeCluster(cores[1], type = "PSOCK")
			registerDoParallel(cl)
          ########################################
          # Start of parallel state both 500 iteractions B = 500
          final_boot <- foreach(i=1:B, .combine='rbind') %dopar% {
            # dont forget to put code repeat

            ## Preallocation Local to the core (no B indexing)
            #preallocate out internal matrices and variables (each core must have its own copy) no B

            n_sampling_tries <- 0
			 M1b <- matrix(0, nrow = nrow(C), ncol = ncol(C)) # usado dentro do loop
            M2b <- matrix(0, nrow = nrow(C), ncol = ncol(C)) # usado dentro do loop

			M_IF0<-matrix(0,nrow=nrow(C))
			M_IF01<-matrix(0,nrow=nrow(C))
			M_IF0Sob<-matrix(0,nrow=nrow(C))
			ind<-seq(1,nrow(C),by=1)
			indbb<-seq(1,nrow(C),by=1)
			indbbb<-seq(1,nrow(C),by=1)

			# repeat sample if condition not met
            repeat{
              # get the sample
              indb <- sample(ind, nrow(C), replace = TRUE)
              M1b <- C[indb,]

              ########################################
              # ordenaçao or ordering zone
              ########################################
              ord <- order(M1b[, 1], method = "auto")
              M2b[, 1] <- sort(M1b[, 1], method = "auto")

              M2b[, 2:ncol(M1b)] <- M1b[ord, 2:ncol(M1b)]

			              # preallocate
						  Aun<-unique(M2b)
              Jb <- matrix(data = 0,ncol = nrow(M2b), nrow = nrow(M2b))

              aub <- outer(M2b[, 1], M2b[, 2], ">=")
			  aubun <- outer(Aun[, 1], Aun[, 2], ">=")
              avb <- outer(M2b[, 1], M2b[, 3], "<=")
			  avbun<-outer(Aun[, 1], Aun[, 3], "<=")
              auub <- outer(M2b[, 2], M2b[, 2], "<=") * 1L
              Jb <- aub * avb
			  Jbun<-aubun*avbun
              SCb <-colSums(Jbun)
				SCb2<-rowSums(Jbun)
              adb <- length(which(SCb == 1))
			  adbb <- length(which(SCb2 == 1))
              n_sampling_tries <- n_sampling_tries + 1 # counter for number of tries

            # if the condition is set, we break the repat, otherwise repeat sampling  process
            if (adb == 0 & adbb==0)
              break
            else
              cat("Warning. Non-uniqueness or no existence of the NPMLE - New Round","\n")
            }
            ## End os sample repeat

            # now we proceed after good sample
            JIb <- t(Jb)
			# preallocation
            f0b <- matrix(data = wt, ncol = 1, nrow = nrow(M2b))
            f1b <- f0b

            S0b <- 1
            iterb <- 0

            # some inner loop
            while (S0b > error | iterb > nmaxit) {
              iterb <- iterb + 1

              if (iterb > nmaxit)
                stop("Default number of iterations not enough for convergence")
			F0b<-JIb%*%f1b
			IF0b<-1/F0b
			If1b<-Jb%*%IF0b
			f1b<-1/If1b
			if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
			S0b<-max(abs(f1b-f0b))
			f0b<-f1b
			}



			FF0b<-numeric(nrow(C))
			for(i in 1:nrow(C)){
			indbb<-(C[,1]==C[i,1])
			pos<-min(which(indbb==TRUE))
			if(pos==1){
			FF0b[indbb]<-sum(f1b[indbb])}
			if(pos>1){
			FF0b[indbb]<-sum(f1b[indbb])+FF0b[pos-1]
			}
			}


			fb<-numeric(nrow(C))
			for(i in 1:nrow(C)){
			indbbb<-(C[,1]==C[i,1])
			pos1<-min(which(indbbb==TRUE))
			if(pos1==1){
			fb[indbbb]<-sum(f1b[indbbb])}
			if(pos1>1){
			fb[indbbb]<-sum(f1b[indbbb])
			}
			}


			Sobb<-1-FF0b+fb
			Sobb[Sobb<1e-12]<-0
			FF0b[FF0b<1e-12]<-0


			# SAVE start saving stufs (before indexed to b index)
            M_IF0 <- as.vector(FF0b)
            M_IF01 <- as.vector(f1b)
            M_IF0Sob <- as.vector(Sobb)


			#M_IF0[,b]<-as.vector(FF0b)
			#M_IF01[,b]<-as.vector(fb)
			#M_IF0Sob[,b]<-as.vector(Sobb)


			########################################
            # ordenaçao or ordering zone
            ########################################
            ordA <- order(M2b[, 2], method = "auto")
            ord3 <- sort(M2b[, 2], method = "auto")

            ordB <- order(M2b[, 3], method = "auto")
            ord4 <- sort(M2b[, 3], method = "auto")

			# variales to be outputed of the clusters (paralell)

          return(list( M_IF0, M_IF01, M_IF0Sob,n_sampling_tries))

          # End of parallel state both
          ########################################
          }
          #stop cluster
          stopCluster(cl)

          # computation of the B aggegated results r comes as a list concatenated as row (thsi enables any size lists)

          # preallocate the original ones as matrix a,d put the list indexing
          ############################################
          n_sampling_tries <- 0
          M_IF0 <- matrix(0, nrow = nrow(C), ncol = B)
          M_IF01 <- matrix(0, nrow = nrow(C), ncol = B)
          M_IF0Sob <- matrix(0, nrow = nrow(C), ncol = B)
          stderror<-matrix(0, nrow = nrow(C), ncol = 1)
		  BootRepeat<-matrix(0,ncol=B)
          ###########################################

          ###########################################
          # Now convert back to the indices from the array list of
          # TODO
          ###########################################

          	for(b in 1:B){
			M_IF0[,b]<-as.numeric(unlist(final_boot[b,1]))
			M_IF01[,b]<-as.numeric(unlist(final_boot[b,2]))
			M_IF0Sob[,b]<-as.numeric(unlist(final_boot[b,3]))
			BootRepeat[,b]<-as.numeric((unlist(final_boot[b,4])))
			}

		  stderror<-apply(M_IF0,1,sd)

		  M_IF0_sort<-matrix(0,nrow=nrow(C),ncol=B)
			for(i in 1:nrow(M_IF0_sort)){
			M_IF0_sort[i,]<-sort(M_IF0[i,])
			}
		if(is.na(alpha)==TRUE) alpha<-0.05
		lowerF<-M_IF0_sort[,floor(alpha*B/2)]
		upperF<-M_IF0_sort[,floor((1-alpha/2)*B)]

		M_IF0_sort1<-matrix(0,nrow=nrow(C),ncol=B)
		for(i in 1:nrow(M_IF0_sort1)){
		M_IF0_sort1[i,]<-sort(M_IF0Sob[i,])
		}
		lowerS<-M_IF0_sort1[,floor(alpha*B/2)]
		upperS<-M_IF0_sort1[,floor((1-alpha/2)*B)]


} ###end trunc rigth






		} ##############################
        # end of bootstrap true
        ##############################



	if(boot==TRUE){

	if(trunc=="both"|trunc=="left"){
	x<-unique(C[,1])
	if((display.F==TRUE)&(display.S==TRUE)){

	dev.new()
	par(mfrow=c(1,2))


plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="CDF", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperF[i],C[i+1,1],upperF[i], lty=3)
segments(C[i+1,1],upperF[i],C[i+1,1],upperF[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],lowerF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],lowerF[i],C[i+1,1],lowerF[i], lty=3)
segments(C[i+1,1],lowerF[i],C[i+1,1],lowerF[i+1],lty=3)
}


plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob0))



for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}

segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperS[i],C[i+1,1],upperS[i], lty=3)
segments(C[i+1,1],upperS[i],C[i+1,1],upperS[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)



for(i in 1:(length(C[,1])-1)){
segments(C[i,1],lowerS[i],C[i+1,1],lowerS[i], lty=3)
segments(C[i+1,1],lowerS[i],C[i+1,1],lowerS[i+1],lty=3)
}


} ####display true and true



if((display.S==FALSE)&(display.F==TRUE)){
#x<-C[,1]
dev.new()
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="CDF", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperF[i],C[i+1,1],upperF[i], lty=3)
segments(C[i+1,1],upperF[i],C[i+1,1],upperF[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(x),1,lty=3)
segments(C[,1][1],0,C[,1][1],lowerF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[i,1],lowerF[i],C[i+1,1],lowerF[i], lty=3)
segments(C[i+1,1],lowerF[i],C[i+1,1],lowerF[i+1],lty=3)
}
} ###display FALSE AND TRUE


if((display.F==FALSE)&(display.S==TRUE)){
dev.new()
plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob0))



for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}

segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[i,1],upperS[i],C[i+1,1],upperS[i], lty=3)
segments(C[i+1,1],upperS[i],C[i+1,1],upperS[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[i,1],lowerS[i],C[i+1,1],lowerS[i], lty=3)
segments(C[i+1,1],lowerS[i],C[i+1,1],lowerS[i+1],lty=3)
}
}#### DISPLAY FALSE AND TRUE

} ###end trunc both and left

if(trunc=="right"){


if((display.F==TRUE)&(display.S==TRUE)){
x<-unique(C[,1])
dev.new()
par(mfrow=c(1,2))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="CDF", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FFF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(max(C[,1]),1,max(C[,1]),max(upperF), lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],upperF[i],C[,1][i+1],upperF[i], lty=3)
segments(C[,1][i+1],upperF[i],C[,1][i+1],upperF[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(max(C[,1]),1,max(C[,1]),max(lowerF), lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],lowerF[i],C[,1][i+1],lowerF[i], lty=3)
segments(C[,1][i+1],lowerF[i],C[,1][i+1],lowerF[i+1],lty=3)
}


plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob0))



for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}

segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],upperS[i],C[,1][i+1],upperS[i], lty=3)
segments(C[,1][i+1],upperS[i],C[,1][i+1],upperS[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],lowerS[i],C[,1][i+1],lowerS[i], lty=3)
segments(C[,1][i+1],lowerS[i],C[,1][i+1],lowerS[i+1],lty=3)
}

} ##display TRUE AND TRUE

if((display.F==FALSE)&(display.S==TRUE)){
x<-unique(C[,1])
dev.new()
par(mfrow=c(1,1))
plot(x,Sob,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob0))




for(i in 1:(length(x)-1)){
segments(x[i],Sob[i],x[i+1],Sob[i])
segments(x[i+1],Sob[i],x[i+1],Sob[i+1])
}

segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],upperS[i],C[,1][i+1],upperS[i], lty=3)
segments(C[,1][i+1],upperS[i],C[,1][i+1],upperS[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)


for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],lowerS[i],C[,1][i+1],lowerS[i], lty=3)
segments(C[,1][i+1],lowerS[i],C[,1][i+1],lowerS[i+1],lty=3)
}
}###DISPLAY FALSE AND TRUE

if((display.F==TRUE)&(display.S==FALSE)){
x<-unique(C[,1])
dev.new()
par(mfrow=c(1,1))
plot(x,FF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="CDF", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FFF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FF[i],x[i+1],FF[i])
segments(x[i+1],FF[i],x[i+1],FF[i+1])
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],upperF[i],C[,1][i+1],upperF[i], lty=3)
segments(C[,1][i+1],upperF[i],C[,1][i+1],upperF[i+1],lty=3)
}


segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),0,C[,1][1],0,lty=3)
segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
segments(C[,1][1],0,C[,1][1],lowerF[1], lty=3)

for(i in 1:(length(C[,1])-1)){
segments(C[,1][i],lowerF[i],C[,1][i+1],lowerF[i], lty=3)
segments(C[,1][i+1],lowerF[i],C[,1][i+1],lowerF[i+1],lty=3)
}

}###DISPLAY TRUE AND FALSE

} ##end trunc rigth
} ###end BOOT  true


if (boot==FALSE|B<40){


if(trunc=="both"|trunc=="left"){

if((display.F==TRUE)&(display.S==TRUE)){
x<-C[,1]
dev.new()
par(mfrow=c(1,2))
plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="CDF", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FFF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FFF[i],x[i+1],FFF[i])
segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
}




plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="X",ylab="")
segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob0))
segments(x0 = min(x),
         x1 = min(x),
         y0 = Sob0[1],
         y1 = 1) 


for(i in 1:(length(x)-1)){
segments(x[i],Sob0[i],x[i+1],Sob0[i])
segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
}

}



if((display.S==FALSE)&(display.F==TRUE)){
x<-C[,1]
dev.new()
plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="CDF", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FFF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FFF[i],x[i+1],FFF[i])
segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
}

}

if((display.F==FALSE)&(display.S==TRUE)){
x<-C[,1]
dev.new()
plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="X",ylab="")
segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob0))
segments(x0 = min(x),
         x1 = min(x),
         y0 = Sob0[1],
         y1 = 1) 


for(i in 1:(length(x)-1)){
segments(x[i],Sob0[i],x[i+1],Sob0[i])
segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
}



}
}## END trunc both and left

if(trunc=="right"){

Sob0[Sob0<1e-12]<-0
FFF[FFF<1e-12]<-0


if((display.F==TRUE)&(display.S==TRUE)){
x<-C[,1]
dev.new()
par(mfrow=c(1,2))
plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="CDF", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FFF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FFF[i],x[i+1],FFF[i])
segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
}

plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="X",ylab="")
segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob0))
segments(x0 = min(x),
         x1 = min(x),
         y0 = Sob0[1],
         y1 = 1) 


for(i in 1:(length(x)-1)){
segments(x[i],Sob0[i],x[i+1],Sob0[i])
segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
}



 }

if((display.F==FALSE)&(display.S==TRUE)){
x<-C[,1]
dev.new()

par(mfrow=c(1,1))
plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
segments(max(x),0,max(x),min(Sob0))
segments(x0 = min(x),
         x1 = min(x),
         y0 = Sob0[1],
         y1 = 1) 


for(i in 1:(length(x)-1)){
segments(x[i],Sob0[i],x[i+1],Sob0[i])
segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
}
 }

if((display.F==TRUE)&(display.S==FALSE)){
x<-C[,1]
dev.new()
par(mfrow=c(1,1))
plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="CDF", xlab="X",ylab="")

segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
segments(x[1],0,x[1],FFF[1])

for(i in 1:(length(x)-1)){
segments(x[i],FFF[i],x[i+1],FFF[i])
segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
}
 }

}###end trunc rigth

}##end boot false

cat("n.iterations",iter,"\n")
cat("S0",S0,"\n")
cat("events",events,"\n")


if(boot==TRUE){
 cat("B",B,"\n")
 cat("alpha",alpha,"\n")

#summary<-cbind("time"=x,"n.event"=mult,"density"=round(f,5),"cumulative.df"=round(FF,5),"survival"=round(Sob,5))

#colnames(summary)<-c("time","n.event","density", "cumulative.df", "survival")
#rownames(summary)<-rep("",times=length(x))
#print(summary,digits=5, justify="left")

return(invisible(list(n.iterations=iter, events=events, B=B, alpha=alpha,time=C[,1], n.event=mult, density=round(as.vector(f),5), cumulative.df=round(FFF,5), survival=
round(as.vector(Sob0),5), truncation.probs=round(as.vector(F0),5), upper.df=round(upperF,5),lower.df=round(lowerF,5),upper.Sob=round(upperS,5),
lower.Sob=round(lowerS,5), sd.boot=round(stderror,5),boot.repeat=as.vector(BootRepeat))))


}

if(boot==FALSE){

#summary<-cbind("time"=x,"n.event"=mult,"density"=round(f,5),"cumulative.df"=round(FF,5),"survival"=round(Sob,5))

#colnames(summary)<-c("time","n.event","density", "cumulative.df", "survival")
#rownames(summary)<-rep("",times=length(x))
#print(summary,digits=5, justify="left")
return(invisible(list(n.iterations=iter, events=events, time=C[,1], n.event=mult, density=round(as.vector(f),5), cumulative.df=round(FFF,5), survival=
round(as.vector(Sob0),5), truncation.probs=round(as.vector(F0),5))))


}
}

