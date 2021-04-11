shen <-
  function (X,
            U = NA,
            V = NA,
            wt = NA,
            error = NA,
            nmaxit = NA,
			boot = TRUE,
            boot.type = "simple",
            B = NA,
            alpha = NA,
			display.FS=FALSE,
			display.UV=FALSE,
                 plot.joint=FALSE,
				 plot.type=NULL)
			{
    ##############################
    # start of shen
    ##############################



    # set truncation both side by default
    trunc <- "both"

    # analize if we have a NA and remove then
    if (all(is.na(U)) == TRUE) {
      trunc <- "right"
      cat("case U=NA","\n")

      if (any(is.na(V)) == TRUE | any(is.na(X)) == TRUE) {
        navec <- c(which(is.na(X)), which(is.na(V)))
        X <- X[-navec]
        V <- V[-navec]
      }
    }


    # check for NA in the resulting vector
    if (all(is.na(V)) == TRUE) {
      trunc <- "left"
      cat("case V=NA","\n")


      if (any(is.na(U)) == TRUE | any(is.na(X)) == TRUE) {
        navec <- c(which(is.na(X)), which(is.na(U)))
        X <- X[-navec]
        U <- U[-navec]
      }
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


    # preallocation
    C <- matrix(0, nrow = nrow(D), ncol = ncol(D))
    EE <- matrix(0, nrow = nrow(D), ncol = ncol(D))

    ########################################
    # ordenaçao or ordering zone
    ########################################

      ord <- order(D[, 1], method = "auto")
      C[, 1] <- sort(D[, 1], method = "auto")
      C[, 2:ncol(D)] <-
      D[ord, 2:ncol(D)] # use the ordered indeces for ordering

	T<-table(C[,2]<= C[,1]&C[,1]<=C[,3])
	if(length(T)!=1){
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

      F0 <- JI %*% f
      k0 <- ((sum(1 / F0)) ^ (-1)) * (1 / F0)
      if (sum(k0) != 1)
        k0 <- k0 / sum(k0)

      k <- k0
      K0 <- J %*% k
      f <- ((sum(1 / K0)) ^ (-1)) * (1 / K0)
      if (sum(f) != 1)
        f <- f / sum(f)

      S0 <-
        max(abs(f - f0)) # get the max from the absoliute differences in the array
      f0 <- f
      k0 <- k
    }

    F0 <- JI %*% f
    mult <- tabulate(match(C[, 1], unique(C[, 1])))

    if (sum(mult) == length(unique(C[, 1])))
      Fval <- (f * mult)

    if (sum(mult) > length(unique(C[, 1]))) {
      weigth <- f[!duplicated(C[, 1])]
      Fval <- (weigth * mult)
    }

    x <- unique(C[, 1])
    events <- sum(mult)
    n.event <- mult
    FF <- cumsum(Fval) # do cumsum computations
    FFF<-cumsum(f) # do cumsum without ties
    Sob <- 1 - FF + Fval
    Sob[Sob < 1e-12] <- 0
    Sob0<-1-FFF
	Sob[Sob<1e-12]<-0
	Sob0[Sob0<1e-12]<-0




    ########################################
    # FORK FOR BOTH SIDES TRUNCATION
    ########################################
    if (trunc == "both") {


      indbbb <- seq(1, nrow(C), by = 1)
      indbbu <- seq(1, nrow(C), by = 1)
      indbbv <- seq(1, nrow(C), by = 1)
      kMUV <- cbind(C[, 2], C[, 3], k)

      ordUV <- order(kMUV[, 1])

      kMUV[, 1] <- sort(kMUV[, 1])
      kMUV[, 2] <- kMUV[ordUV, 2]
      kMUV[, 3] <- kMUV[ordUV, 3]
      kuv <- numeric(nrow(C))

      kMU <- cbind(C[, 2], k)

      ########################################
      # ordenaçao or ordering zone
      ########################################
      ordU <- order(kMU[, 1], method = "auto")
      kMU[, 1] <- sort(kMU[, 1], method = "auto")
      ########################################


      kMU[, 2] <- kMU[ordU, 2]
      kk0u <- numeric(nrow(C))
      kk0u <- kMU[, 2]

      UU <- unique(kMU[, 1])
		fU<-cumsum(kk0u)
      kMV <- cbind(C[, 3], k)

      ########################################
      # ordenaçao or ordering zone
      ########################################
      ordV <- order(kMV[, 1], method = "auto")
      kMV[, 1] <- sort(kMV[, 1], method = "auto")
      ########################################

      kMV[, 2] <- kMV[ordV, 2]
      kk0v <- numeric(nrow(C))
		kk0v <- kMV[, 2]


      VV <- unique(kMV[, 1])
	  fV<-cumsum(kk0v)

      ########################################
      # ordenaçao or ordering zone
      ########################################

        ordUV <- order(kMUV[, 1], method = "auto")



      EE[, 1] <- C[ordUV, 1]


      ########################################
      # ordenaçao or ordering zone
      ########################################

	  Gf<-matrix(data=k,ncol=ncol(J),nrow=nrow(J), byrow=T)
		Gff<-J*Gf
		biasf<-rowSums(Gff)


      KK <- matrix(data = 0,
                   ncol = nrow(C),
                   nrow = nrow(C))


	  	auk <- outer(kMU[, 1], kMUV[, 1], ">=")
    		avk <- outer(kMV[, 1], kMUV[, 2], ">=")
		Jk<-auk*avk
	 	 Gfk<-matrix(data=kMUV[,3],ncol=ncol(Jk),nrow=nrow(Jk), byrow=T)
		Gffk<-Jk*Gfk
		KK<-apply(Gffk,1,cumsum)



     if(plot.joint==TRUE){
      if(plot.type=="image"){
        sU<-sort(U)
        sV<-sort(V)
        if(any(diff(sU)==0)& all(diff(sV)!=0)){
          stepU<-diff(range(U))/(length(U)*length(unique(U)))
          image(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution")
          filled.contour(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
          filled.contour(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
        }

        if(all(diff(sU)!=0)& any(diff(sV)==0)){
          stepV<-diff(range(V))/(length(V)*length(unique(V)))
          image(sort(U),sort(V+runif(length(U),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution")
          filled.contour(sort(U),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
          filled.contour(sort(U),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
        }

        if(any(diff(sU)==0)& any(diff(sV)==0)){
          stepU<-diff(range(U))/(length(U)*length(unique(U)))
          stepV<-diff(range(V))/(length(V)*length(unique(V)))
          image(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(U),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution")
          filled.contour(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
          filled.contour(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
        }

        if(all(diff(sU)!=0)& all(diff(sV)!=0)){
          image(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution")
          filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
          filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
        }}

      if(plot.type=="persp"){
        fcol<-topo.colors(10)[cut(KK[2:nrow(C),2:nrow(C)],10,include.lowest=TRUE)]
        sU<-sort(U)
        sV<-sort(V)
        if(any(diff(sU)==0)& all(diff(sV)!=0)){
          stepU<-diff(range(U))/(length(U)*length(unique(U)))
          persp(x=sort(U+runif(length(U),0,stepU)),y=sort(V),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
        }
        if(all(diff(sU)!=0)& any(diff(sV)==0)){
          stepV<-diff(range(V))/(length(V)*length(unique(V)))
          persp(x=sort(U),y=sort(V+runif(length(V),0,stepV)),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
        }
        if(any(diff(sU)==0)& any(diff(sV)==0)){
          stepU<-diff(range(U))/(length(U)*length(unique(U)))
          stepV<-diff(range(V))/(length(V)*length(unique(V)))
          persp(x=sort(U+runif(length(U),0,stepU)),y=sort(V+runif(length(V),0,stepV)),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
        }
        if(all(diff(sU)!=0)& all(diff(sV)!=0)){
          persp(x=sort(U),y=sort(V),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
        }
      }


    }
    if(plot.joint==FALSE){
    }

      # END FORK FOR BOTH SIDE TRUNCATION
      ########################################
    }




	########################################
    # FORK FOR LEFT SIDE TRUNCATION
    ########################################
    if (trunc == "left") {
      indbbb <- seq(1, nrow(C), by = 1)
      indbb9 <- seq(1, nrow(C), by = 1)
      kMU <- cbind(C[, 2], k)

      ########################################
      # ordenaçao or ordering zone
      ########################################
      ordU <- order(kMU[, 1], method = "auto")
      kMU[, 1] <- sort(kMU[, 1], method = "auto")

      kMU[, 2] <- kMU[ordU, 2]
      EE[, 1] <- X[ordU]

      # pre allocation
      KKG <- matrix(data = 0,
                    ncol = 1,
                    nrow = nrow(EE))



      multG <- tabulate(match(EE[, 1], unique(EE[, 1])))

	  Juleft <- outer(C[, 1], C[, 2], ">=")

	  Gfleft<-matrix(data=kMU[,2],ncol=ncol(Juleft),nrow=nrow(Juleft), byrow=T)
		Gffleft<-Juleft*Gfleft
		biasf<-rowSums(Gffleft)



      ########################################
      # ordenaçao or ordering zone
      ########################################
      ordc <- order(unique(EE[, 1]), method = "auto")

      #biasf <- fGval[ordc]
      kk0b <- numeric(nrow(kMU))


      #############PODEREI APAGAR ISTO TUDO######################


        kk0b<-kMU[,2]
        UU <- unique(kMU[, 1])
    	  ku<-kk0b
      	  fU<-cumsum(kk0b)

    # END FORK FOR LEFT SIDE TRUNCATION
    ########################################
    }


	########################################
    # FORK FOR RIGHT SIDE TRUNCATION
    ########################################
    if (trunc == "right") {
      indbbb <- seq(1, nrow(C), by = 1)
      indbb9 <- seq(1, nrow(C), by = 1)
      kMV <- cbind(C[, 3], k)

      ########################################
      # ordenaçao or ordering zone
      ########################################
      ordV <- order(kMV[, 1], method = "auto")
      kMV[, 1] <- sort(kMV[, 1], method = "auto")


      kMV[, 2] <- kMV[ordV, 2]
      EE[, 1] <- X[ordV]

      # pre allocation
      KKG <- matrix(data = 0, ncol = 1, nrow = nrow(EE))


      multG <- tabulate(match(EE[, 1], unique(EE[, 1])))

	Jvrigth <- outer(C[, 1], C[, 3], "<=")

	  Gfrigth<-matrix(data=kMV[,2],ncol=ncol(Jvrigth),nrow=nrow(Jvrigth), byrow=T)
		Gffrigth<-Jvrigth*Gfrigth
		biasf<-rowSums(Gffrigth)

      ########################################
      # ordenaçao or ordering zone
      ########################################
      ordc <- order(unique(EE[, 1]), method = "auto")

            kk0b <- numeric(nrow(kMV))

      	kk0b<-kMV[,2]


         VV <- unique(kMV[, 1])
      kv <- kk0b
      fV <- cumsum(kv)

    # END FORK FOR RIGHT SIDE TRUNCATION
    ########################################
    }



	########################################
    # FORK FOR BOOTSTRAP
    ########################################
    if (boot == TRUE) {

      if (boot.type == "simple") {
        if (is.na(B) == TRUE)

		 B <- 500 # this is harcoded


        if (B < 40) {
          cat("Warning. Number of replicates less than 40", "\n")
          cat("Confidence bands cannot be computed", "\n")
        }

        ########################################


        ########################################
        # BOOSTRAP FORK FOR BITH SIDE TRUNCATION
        ########################################
        if (trunc == "both") {


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

            M_IF0 <- matrix(0, nrow = nrow(C)) # as vector
            M_IF01 <- matrix(0, nrow = nrow(C)) # as vector
            M_IF0Sob <- matrix(0, nrow = nrow(C)) # as vector

            kMUb <- matrix(0, nrow = nrow(C))
            kMVb <- matrix(0, nrow = nrow(C))
            M_fU <- matrix(0, nrow = nrow(C))
            M_fV <- matrix(0, nrow = nrow(C))
            M3b <- matrix(0, nrow = nrow(C))
            M4b <- matrix(0, nrow = nrow(C))
            k1bU <- matrix(0, nrow = nrow(C))
            k1bV <- matrix(0, nrow = nrow(C))

            ind <- seq(1, nrow(C), by = 1)
            indbb1 <- seq(1, nrow(C), by = 1)
            indbb7 <- seq(1, nrow(C), by = 1)

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
              Jb <- matrix(data = 0,ncol = nrow(M2b), nrow = nrow(M2b))

              aub <- outer(M2b[, 1], M2b[, 2], ">=")
              avb <- outer(M2b[, 1], M2b[, 3], "<=")
              auub <- outer(M2b[, 2], M2b[, 2], "<=") * 1L
              Jb <- aub * avb
              SCb <- apply(Jb, 2, "sum")
              adb <- length(which(SCb == 1))

              n_sampling_tries <- n_sampling_tries + 1 # counter for number of tries

            # if the condition is set, we break the repat, otherwise repeat sampling  process
            if (adb == 0)
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

              F0b <- JIb %*% f1b
              k0b <- ((sum(1 / F0b)) ^ (-1)) * (1 / F0b)
              if (sum(k0b) != 1)
                k0b <- k0b / sum(k0b)

              k1b <- k0b
              K0b <- Jb %*% k1b
              f1b <- ((sum(1 / K0b)) ^ (-1)) * (1 / K0b)

              if (sum(f1b) != 1)
                f1b <- f1b / sum(f1b)

              S0b <- max(abs(f1b - f0b))
              f0b <- f1b
              k0b <- k1b
            }

            ff0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb1 <- (C[, 1] == C[i, 1])
              pos1 <- min(which(indbb1 == TRUE))

              if (pos1 >= 1)
                ff0b[indbb1] <- sum(f1b[indbb1])

            }

            FF0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb7 <- (C[, 1] == C[i, 1])
              pos0 <- min(which(indbb7 == TRUE))

              if (pos0 == 1)
                FF0b[indbb7] <- sum(f1b[indbb7])

              else if (pos0 > 1)
                FF0b[indbb7] <- sum(f1b[indbb7]) + FF0b[pos0 - 1]
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


            M3b <- ord3
            M4b <- ord4
            k1bU <- k1b[ordA]
            k1bV <- k1b[ordB]


            # check this was two loops
            for (i in 1:nrow(C)) {
              for (j in 1:nrow(C)) {
                if (M3b[i] <= kMU[j, 1])
                  M_fU[j] <- sum(k1bU[1:i])

                if (M4b[i] <= kMV[j, 1])
                  M_fV[j] <- sum(k1bV[1:i])
              }
            }

			#aUb<-outer(M3b, kMU[,1] "<=")

            # variales to be outputed of the clusters (paralell)

          return(list( M_IF0, M_IF01, M_IF0Sob, M_fU, M_fV,n_sampling_tries))

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
          M_fU <- matrix(0, nrow = nrow(C), ncol = B)
          M_fV <- matrix(0, nrow = nrow(C), ncol = B)
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
			M_fU[,b]<-as.numeric(unlist(final_boot[b,4]))
			M_fV[,b]<-as.numeric(unlist(final_boot[b,5]))
			BootRepeat[,b]<-as.numeric((unlist(final_boot[b,6])))
			}

		  stderror<-apply(M_IF0,1,sd)

          # this is same as before
          M_IF0_sort <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_IF0_sort)) {
            M_IF0_sort[i,] <- sort(M_IF0[i,])
          }

          if (is.na(alpha) == TRUE)
            alpha <- 0.05

          lowerF <- M_IF0_sort[, floor(alpha * B / 2)]
          upperF <- M_IF0_sort[, floor((1 - alpha / 2) *
                                         B)]
          M_IF0_sort1 <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_IF0_sort1)) {
            M_IF0_sort1[i,] <- sort(M_IF0Sob[i,])
          }

          lowerS <- M_IF0_sort1[, floor(alpha * B / 2)]
          upperS <- M_IF0_sort1[, floor((1 - alpha / 2) * B)]
          M_fU_sort <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_fU_sort)) {
            M_fU_sort[i,] <- sort(M_fU[i,])
          }

          lowerU <- M_fU_sort[, floor(alpha * B / 2)]
          upperU <- M_fU_sort[, floor((1 - alpha / 2) * B)]
          M_fV_sort <- matrix(0, nrow = nrow(C), ncol = B)

          for (i in 1:nrow(M_fV_sort)) {
            M_fV_sort[i,] <- sort(M_fV[i,])
          }

          lowerV <- M_fV_sort[, floor(alpha * B / 2)]
          upperV <- M_fV_sort[, floor((1 - alpha / 2) * B)]

        }
        ##############################
        # end of bootstrap trucn BOTH
        ##############################




	if (trunc == "left") {




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

            M_IF0 <- matrix(0, nrow = nrow(C)) # as vector
            M_IF01 <- matrix(0, nrow = nrow(C)) # as vector
            M_IF0Sob <- matrix(0, nrow = nrow(C)) # as vector

            kMUb <- matrix(0, nrow = nrow(C))
            kMVb <- matrix(0, nrow = nrow(C))
            M_fU <- matrix(0, nrow = nrow(C))
            M_fV <- matrix(0, nrow = nrow(C))
            M3b <- matrix(0, nrow = nrow(C))
            M4b <- matrix(0, nrow = nrow(C))
            k1bU <- matrix(0, nrow = nrow(C))
            k1bV <- matrix(0, nrow = nrow(C))

            ind <- seq(1, nrow(C), by = 1)
            indbb1 <- seq(1, nrow(C), by = 1)
            indbb7 <- seq(1, nrow(C), by = 1)


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
              Jb <- matrix(data = 0,ncol = nrow(M2b), nrow = nrow(M2b))

              aub <- outer(M2b[, 1], M2b[, 2], ">=")
              avb <- outer(M2b[, 1], M2b[, 3], "<=")
              auub <- outer(M2b[, 2], M2b[, 2], "<=") * 1L
              Jb <- aub * avb
              SCb <- apply(Jb, 2, "sum")
              adb <- length(which(SCb == 1))

              n_sampling_tries <- n_sampling_tries + 1 # counter for number of tries

            # if the condition is set, we break the repat, otherwise repeat sampling  process
            if (adb == 0)
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

              F0b <- JIb %*% f1b
              k0b <- ((sum(1 / F0b)) ^ (-1)) * (1 / F0b)
              if (sum(k0b) != 1)
                k0b <- k0b / sum(k0b)

              k1b <- k0b
              K0b <- Jb %*% k1b
              f1b <- ((sum(1 / K0b)) ^ (-1)) * (1 / K0b)

              if (sum(f1b) != 1)
                f1b <- f1b / sum(f1b)

              S0b <- max(abs(f1b - f0b))
              f0b <- f1b
              k0b <- k1b
            }

            ff0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb1 <- (C[, 1] == C[i, 1])
              pos1 <- min(which(indbb1 == TRUE))

              if (pos1 >= 1)
                ff0b[indbb1] <- sum(f1b[indbb1])

            }

            FF0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb7 <- (C[, 1] == C[i, 1])
              pos0 <- min(which(indbb7 == TRUE))

              if (pos0 == 1)
                FF0b[indbb7] <- sum(f1b[indbb7])

              else if (pos0 > 1)
                FF0b[indbb7] <- sum(f1b[indbb7]) + FF0b[pos0 - 1]
            }

            Sobb <- 1 - FF0b + ff0b
            Sobb[Sobb < 1e-12] <- 0


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


            M3b <- ord3
            M4b <- ord4
            #k1bU <- k1b[ordU]
			k1bU <- k1b[ordA]
            #k1bV <- k1b[ordV]



			# check this was two loops
            for (i in 1:nrow(C)) {
              for (j in 1:nrow(C)) {
                if (M3b[i] <= kMU[j, 1])
                  M_fU[j] <- sum(k1bU[1:i])
              }
            }



            # variales to be outputed of the clusters (paralell)

          return(list( M_IF0, M_IF01, M_IF0Sob, M_fU,n_sampling_tries))

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
          M_fU <- matrix(0, nrow = nrow(C), ncol = B)
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
			M_fU[,b]<-as.numeric(unlist(final_boot[b,4]))
			BootRepeat[,b]<-as.numeric((unlist(final_boot[b,5])))
			}

          stderror<-apply(M_IF0,1,sd)

          # this is same as before
          M_IF0_sort <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_IF0_sort)) {
            M_IF0_sort[i,] <- sort(M_IF0[i,])
          }

          if (is.na(alpha) == TRUE)
            alpha <- 0.05

          lowerF <- M_IF0_sort[, floor(alpha * B / 2)]
          upperF <- M_IF0_sort[, floor((1 - alpha / 2) *
                                         B)]
          M_IF0_sort1 <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_IF0_sort1)) {
            M_IF0_sort1[i,] <- sort(M_IF0Sob[i,])
          }

          lowerS <- M_IF0_sort1[, floor(alpha * B / 2)]
          upperS <- M_IF0_sort1[, floor((1 - alpha / 2) * B)]
          M_fU_sort <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_fU_sort)) {
            M_fU_sort[i,] <- sort(M_fU[i,])
          }

          lowerU <- M_fU_sort[, floor(alpha * B / 2)]
          upperU <- M_fU_sort[, floor((1 - alpha / 2) * B)]

        }
        ##############################
        # end of bootstrap trucn LEFT
        ##############################



	if (trunc == "right") {




	########################################
          # PARALELL BOOTSTARP RIGHT
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
            M2bb<-matrix(0,nrow=nrow(C),ncol=ncol(C))
            M_IF0 <- matrix(0, nrow = nrow(C)) # as vector
            M_IF01 <- matrix(0, nrow = nrow(C)) # as vector
            M_IF0Sob <- matrix(0, nrow = nrow(C)) # as vector


            kMVb <- matrix(0, nrow = nrow(C))
			M_fV <- matrix(0, nrow = nrow(C))
			M4b <- matrix(0, nrow = nrow(C))
			k1bV <- matrix(0, nrow = nrow(C))

            ind <- seq(1, nrow(C), by = 1)
            indbb1 <- seq(1, nrow(C), by = 1)
            indbb7 <- seq(1, nrow(C), by = 1)





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
              Jb <- matrix(data = 0,ncol = nrow(M2b), nrow = nrow(M2b))

              aub <- outer(M2b[, 1], M2b[, 2], ">=")
              avb <- outer(M2b[, 1], M2b[, 3], "<=")
              auub <- outer(M2b[, 2], M2b[, 2], "<=") * 1L
              Jb <- aub * avb
              SCb <- apply(Jb, 2, "sum")
              adb <- length(which(SCb == 1))

              n_sampling_tries <- n_sampling_tries + 1 # counter for number of tries

            # if the condition is set, we break the repat, otherwise repeat sampling  process
            if (adb == 0)
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

              F0b <- JIb %*% f1b
              k0b <- ((sum(1 / F0b)) ^ (-1)) * (1 / F0b)
              if (sum(k0b) != 1)
                k0b <- k0b / sum(k0b)

              k1b <- k0b
              K0b <- Jb %*% k1b
              f1b <- ((sum(1 / K0b)) ^ (-1)) * (1 / K0b)

              if (sum(f1b) != 1)
                f1b <- f1b / sum(f1b)

              S0b <- max(abs(f1b - f0b))
              f0b <- f1b
              k0b <- k1b
            }

            ff0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb1 <- (C[, 1] == C[i, 1])
              pos1 <- min(which(indbb1 == TRUE))

              if (pos1 >= 1)
                ff0b[indbb1] <- sum(f1b[indbb1])

            }

            FF0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb7 <- (C[, 1] == C[i, 1])
              pos0 <- min(which(indbb7 == TRUE))

              if (pos0 == 1)
                FF0b[indbb7] <- sum(f1b[indbb7])

              else if (pos0 > 1)
                FF0b[indbb7] <- sum(f1b[indbb7]) + FF0b[pos0 - 1]
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


            ordB <- order(M2b[, 3], method = "auto")
            ord4 <- sort(M2b[, 3], method = "auto")

            # end ordenaçao or ordering zone
            ########################################



            M4b <- ord4

	k1bV <- k1b[ordB]

            # check this was two loops
########################################################
            for (i in 1:nrow(C)) {
              for (j in 1:nrow(C)) {
                if (M4b[i] <= kMV[j, 1])
                  M_fV[j] <- sum(k1bV[1:i])
              }
            }



            # variales to be outputed of the clusters (paralell)

          return(list( M_IF0, M_IF01, M_IF0Sob, M_fV,n_sampling_tries))

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
          M_fV <- matrix(0, nrow = nrow(C), ncol = B)
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
			M_fV[,b]<-as.numeric(unlist(final_boot[b,4]))
			BootRepeat[,b]<-as.numeric((unlist(final_boot[b,5])))
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


        M_fV_sort<-matrix(0,nrow=nrow(C),ncol=B)
        for(i in 1:nrow(M_fV_sort)){
          M_fV_sort[i,]<-sort(M_fV[i,])
        }
        lowerV<-M_fV_sort[,floor(alpha*B/2)]
        upperV<-M_fV_sort[,floor((1-alpha/2)*B)]

        }

        ##############################
        # end of bootstrap trucn RIGTH
        ##############################



	}##end of simple



	if (boot.type == "obvious") {

	 if (is.na(B)==TRUE) B<-500

	if(B<40){
          cat("Warning. Number of replicates less than 40","\n")
          cat("Confidence bands cannot be computed","\n")
        }

	if (is.na(alpha) == TRUE)
            alpha <- 0.05


	if(trunc=="both"){


		########################################
    # PARALELL BOOTSTARP BOTH
    ########################################
    #library(foreach)
    #library(doParallel)
    #cores=detectCores()
    #cl <- makeCluster(3, type = "PSOCK") # ihnerith all from master


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
    ########################################


    final_boot <- foreach(i=1:B, .combine='rbind') %dopar% {
      ## Preallocation Local to the core (no B indexing)
      #preallocate out internal matrices and variables (each core must have its own copy)

		# dont forget to put code repeat

            ## Preallocation Local to the core (no B indexing)
            #preallocate out internal matrices and variables (each core must have its own copy) no B

            n_sampling_tries <- 0
			wt <- rep(1 / length(C[,1]), times = length(C[,1]))
			error <- 1e-6
			nmaxit<-10000000
            # container for bootstrap operations

            # container for bootstrap operations

        W1<-as.vector(f)
        W2<-as.vector(k)
        ind<-1:nrow(C)
        M_IF0<-matrix(0,nrow=nrow(C))
        M_IF01<-matrix(0,nrow=nrow(C))
        M_IF0Sob<-matrix(0,nrow=nrow(C))
        kMUb<-matrix(0,nrow=nrow(C))
        kMVb<-matrix(0,nrow=nrow(C))
        M_fU<-matrix(0,nrow=nrow(C))
        M_fV<-matrix(0,nrow=nrow(C))
        M3b<-matrix(0,nrow=nrow(C))
        M4b<-matrix(0,nrow=nrow(C))
        k1bU<-matrix(0,nrow=nrow(C))
        k1bV<-matrix(0,nrow=nrow(C))

        indbb<-seq(1,nrow(C),by=1)
        indbb1<-seq(1,nrow(C),by=1)

	repeat {

		  DB<-matrix(0,nrow=nrow(C),ncol=3)
          DB[,1]<-sample(C[,1],size=nrow(C),replace=TRUE,prob=W1)
          indUV<-sample(ind,size=nrow(C),replace=TRUE,prob=W2)
          DB[,2:3]<-C[indUV,2:3]

		  log1<-(DB[,1]>=DB[,2])
          log2<-(DB[,1]<=DB[,3])
          indlog<-(log1*log2==0)

          while(sum(indlog)>0){
            DB[indlog,1]<-sample(C[,1],size=sum(indlog),replace=TRUE,prob=W1)
            indUV<-sample(ind,size=sum(indlog),replace=TRUE,prob=W2)
            DB[indlog,2:3]<-C[indUV,2:3]
            log1<-(DB[,1]>=DB[,2])
            log2<-(DB[,1]<=DB[,3])
            indlog<-(log1*log2==0)
          }

		  ordB<-order(DB[,1],method = "auto")
          DBB<-matrix(0,nrow=nrow(C),ncol=3)
          DBB[,1]<-sort(DB[,1],method = "auto")
          DBB[,2:ncol(DB)]<-DB[ordB,2:ncol(DBB)]


          Jb<-matrix(data=0,ncol=nrow(DBB),nrow=nrow(DBB))
          aub <- outer(DBB[, 1], DBB[, 2], ">=")
              avb <- outer(DBB[, 1], DBB[, 3], "<=")
              auub <- outer(DBB[, 2], DBB[, 2], "<=") * 1L
              Jb <- aub * avb
              SCb <- apply(Jb, 2, "sum")
              adb <- length(which(SCb == 1))

              n_sampling_tries <- n_sampling_tries + 1 # counter for number of tries

            # if the condition is set, we break the repat, otherwise repeat sampling  process
            if (adb == 0)
              break
            else
              cat("Warning. Non-uniqueness or no existence of the NPMLE - New Round","\n")
            }
            ## End os sample repeat

          # now we proceed after good sample
            JIb <- t(Jb)

            # preallocation
            f0b <- matrix(data = wt, ncol = 1, nrow = nrow(DBB))
            f1b <- f0b

            S0b <- 1
            iterb <- 0

            # some inner loop
            while (S0b > error | iterb > nmaxit) {
              iterb <- iterb + 1

              if (iterb > nmaxit)
                stop("Default number of iterations not enough for convergence")

              F0b <- JIb %*% f1b
              k0b <- ((sum(1 / F0b)) ^ (-1)) * (1 / F0b)
              if (sum(k0b) != 1)
                k0b <- k0b / sum(k0b)

              k1b <- k0b
              K0b <- Jb %*% k1b
              f1b <- ((sum(1 / K0b)) ^ (-1)) * (1 / K0b)

              if (sum(f1b) != 1)
                f1b <- f1b / sum(f1b)

              S0b <- max(abs(f1b - f0b))
              f0b <- f1b
              k0b <- k1b
            }

            ff0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb1 <- (C[, 1] == C[i, 1])
              pos1 <- min(which(indbb1 == TRUE))

              if (pos1 >= 1)
                ff0b[indbb1] <- sum(f1b[indbb1])

            }

            FF0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb7 <- (C[, 1] == C[i, 1])
              pos0 <- min(which(indbb7 == TRUE))

              if (pos0 == 1)
                FF0b[indbb7] <- sum(f1b[indbb7])

              else if (pos0 > 1)
                FF0b[indbb7] <- sum(f1b[indbb7]) + FF0b[pos0 - 1]
            }

            Sobb <- 1 - FF0b + ff0b
            Sobb[Sobb < 1e-12] <- 0


			Sobb <- 1 - FF0b + ff0b
            Sobb[Sobb < 1e-12] <- 0

            # SAVE start saving stufs (before indexed to b index)
            M_IF0 <- as.vector(FF0b)
            M_IF01 <- as.vector(f1b)
            M_IF0Sob <- as.vector(Sobb)

			ordA<-order(DBB[,2],method = "auto")
          ordB<-order(DBB[,3],method = "auto")
          ord3<-sort(DBB[,2],method = "auto")
          ord4<-sort(DBB[,3],method = "auto")



			M3b <- ord3
            M4b <- ord4
            #k1bU <- k1b[ordU]
            #k1bV <- k1b[ordV]
			k1bU <- k1b[ordA]
			k1bV <- k1b[ordB]


            # check this was two loops
            for (i in 1:nrow(C)) {
              for (j in 1:nrow(C)) {
                if (M3b[i] <= kMU[j, 1])
                  M_fU[j] <- sum(k1bU[1:i])

                if (M4b[i] <= kMV[j, 1])
                  M_fV[j] <- sum(k1bV[1:i])
              }
            }



    # variales to be outputed of the clusters (paralell)

          return(list( M_IF0, M_IF01, M_IF0Sob, M_fU, M_fV,n_sampling_tries))

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
          M_fU <- matrix(0, nrow = nrow(C), ncol = B)
          M_fV <- matrix(0, nrow = nrow(C), ncol = B)
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
			M_fU[,b]<-as.numeric(unlist(final_boot[b,4]))
			M_fV[,b]<-as.numeric(unlist(final_boot[b,5]))
			BootRepeat[,b]<-as.numeric((unlist(final_boot[b,6])))
			}

          stderror<-apply(M_IF0,1,sd)

	# this is same as before
          M_IF0_sort <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_IF0_sort)) {
            M_IF0_sort[i,] <- sort(M_IF0[i,])
          }



          lowerF <- M_IF0_sort[, floor(alpha * B / 2)]
          upperF <- M_IF0_sort[, floor((1 - alpha / 2) *
                                         B)]
          M_IF0_sort1 <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_IF0_sort1)) {
            M_IF0_sort1[i,] <- sort(M_IF0Sob[i,])
          }

          lowerS <- M_IF0_sort1[, floor(alpha * B / 2)]
          upperS <- M_IF0_sort1[, floor((1 - alpha / 2) * B)]
          M_fU_sort <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_fU_sort)) {
            M_fU_sort[i,] <- sort(M_fU[i,])
          }

          lowerU <- M_fU_sort[, floor(alpha * B / 2)]
          upperU <- M_fU_sort[, floor((1 - alpha / 2) * B)]
          M_fV_sort <- matrix(0, nrow = nrow(C), ncol = B)

          for (i in 1:nrow(M_fV_sort)) {
            M_fV_sort[i,] <- sort(M_fV[i,])
          }

          lowerV <- M_fV_sort[, floor(alpha * B / 2)]
          upperV <- M_fV_sort[, floor((1 - alpha / 2) * B)]

        }

        ##############################
        # end of bootstrap trucn BOTH
        ##############################


	if (trunc == "left") {

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
			wt <- rep(1 / length(C[,1]), times = length(C[,1]))
			error <- 1e-6
			nmaxit<-10000000
            # container for bootstrap operations

            # container for bootstrap operations

		W1<-as.vector(f)
        W2<-as.vector(k)
        ind<-1:nrow(C)
        M_IF0<-matrix(0,nrow=nrow(C))
        M_IF01<-matrix(0,nrow=nrow(C))
        M_IF0Sob<-matrix(0,nrow=nrow(C))
        kMUb<-matrix(0,nrow=nrow(C))
        M_fU<-matrix(0,nrow=nrow(C))
        M3b<-matrix(0,nrow=nrow(C))
        k1bU<-matrix(0,nrow=nrow(C))

        indbb<-seq(1,nrow(C),by=1)
        indbb1<-seq(1,nrow(C),by=1)

	repeat {

		  DB<-matrix(0,nrow=nrow(C),ncol=3)
          DB[,1]<-sample(C[,1],size=nrow(C),replace=TRUE,prob=W1)
          indUV<-sample(ind,size=nrow(C),replace=TRUE,prob=W2)
          DB[,2:3]<-C[indUV,2:3]

		  log1<-(DB[,1]>=DB[,2])
          log2<-(DB[,1]<=DB[,3])
          indlog<-(log1*log2==0)

          while(sum(indlog)>0){
            DB[indlog,1]<-sample(C[,1],size=sum(indlog),replace=TRUE,prob=W1)
            indUV<-sample(ind,size=sum(indlog),replace=TRUE,prob=W2)
            DB[indlog,2:3]<-C[indUV,2:3]
            log1<-(DB[,1]>=DB[,2])
            log2<-(DB[,1]<=DB[,3])
            indlog<-(log1*log2==0)
          }

		  ordB<-order(DB[,1],method = "auto")
          DBB<-matrix(0,nrow=nrow(C),ncol=3)
          DBB[,1]<-sort(DB[,1],method = "auto")
          DBB[,2:ncol(DB)]<-DB[ordB,2:ncol(DBB)]


          Jb<-matrix(data=0,ncol=nrow(DBB),nrow=nrow(DBB))
          aub <- outer(DBB[, 1], DBB[, 2], ">=")
              avb <- outer(DBB[, 1], DBB[, 3], "<=")
              auub <- outer(DBB[, 2], DBB[, 2], "<=") * 1L
              Jb <- aub * avb
              SCb <- apply(Jb, 2, "sum")
              adb <- length(which(SCb == 1))

              n_sampling_tries <- n_sampling_tries + 1 # counter for number of tries

            # if the condition is set, we break the repat, otherwise repeat sampling  process
            if (adb == 0)
              break
            else
              cat("Warning. Non-uniqueness or no existence of the NPMLE - New Round","\n")
            }
            ## End os sample repeat

          # now we proceed after good sample
            JIb <- t(Jb)

            # preallocation
            f0b <- matrix(data = wt, ncol = 1, nrow = nrow(DBB))
            f1b <- f0b

            S0b <- 1
            iterb <- 0

            # some inner loop
            while (S0b > error | iterb > nmaxit) {
              iterb <- iterb + 1

              if (iterb > nmaxit)
                stop("Default number of iterations not enough for convergence")

              F0b <- JIb %*% f1b
              k0b <- ((sum(1 / F0b)) ^ (-1)) * (1 / F0b)
              if (sum(k0b) != 1)
                k0b <- k0b / sum(k0b)

              k1b <- k0b
              K0b <- Jb %*% k1b
              f1b <- ((sum(1 / K0b)) ^ (-1)) * (1 / K0b)

              if (sum(f1b) != 1)
                f1b <- f1b / sum(f1b)

              S0b <- max(abs(f1b - f0b))
              f0b <- f1b
              k0b <- k1b
            }

            ff0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb1 <- (C[, 1] == C[i, 1])
              pos1 <- min(which(indbb1 == TRUE))

              if (pos1 >= 1)
                ff0b[indbb1] <- sum(f1b[indbb1])

            }

            FF0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb7 <- (C[, 1] == C[i, 1])
              pos0 <- min(which(indbb7 == TRUE))

              if (pos0 == 1)
                FF0b[indbb7] <- sum(f1b[indbb7])

              else if (pos0 > 1)
                FF0b[indbb7] <- sum(f1b[indbb7]) + FF0b[pos0 - 1]
            }

            Sobb <- 1 - FF0b + ff0b
            Sobb[Sobb < 1e-12] <- 0


			Sobb <- 1 - FF0b + ff0b
            Sobb[Sobb < 1e-12] <- 0

            # SAVE start saving stufs (before indexed to b index)
            M_IF0 <- as.vector(FF0b)
            M_IF01 <- as.vector(f1b)
            M_IF0Sob <- as.vector(Sobb)

			ordA <- order( DBB[, 2], method = "auto")
            ord3 <- sort( DBB[, 2], method = "auto")




            # end ordenaçao or ordering zone
            ########################################


            M3b <- ord3
            k1bU <- k1b[ordA]


           # check this was two loops
            for (i in 1:nrow(C)) {
              for (j in 1:nrow(C)) {
                if (M3b[i] <= kMU[j, 1])
                  M_fU[j] <- sum(k1bU[1:i])
              }
            }



            # variales to be outputed of the clusters (paralell)

          return(list( M_IF0, M_IF01, M_IF0Sob, M_fU,n_sampling_tries))

          # End of parallel state both
          ########################################
          }
          #stop cluster
          stopCluster(cl)


	  # computation of the B aggegated results r comes as a list concatenated as row (thsi enables any size lists)

          # preallocate the original ones as matrix a,d put the list indexing

	      n_sampling_tries <- 0
          M_IF0 <- matrix(0, nrow = nrow(C), ncol = B)
          M_IF01 <- matrix(0, nrow = nrow(C), ncol = B)
          M_IF0Sob <- matrix(0, nrow = nrow(C), ncol = B)
          M_fU <- matrix(0, nrow = nrow(C), ncol = B)
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
			M_fU[,b]<-as.numeric(unlist(final_boot[b,4]))
			BootRepeat[,b]<-as.numeric((unlist(final_boot[b,5])))
			}

          stderror<-apply(M_IF0,1,sd)

          # this is same as before
          M_IF0_sort <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_IF0_sort)) {
            M_IF0_sort[i,] <- sort(M_IF0[i,])
          }

          if (is.na(alpha) == TRUE)
            alpha <- 0.05

          lowerF <- M_IF0_sort[, floor(alpha * B / 2)]
          upperF <- M_IF0_sort[, floor((1 - alpha / 2) *
                                         B)]
          M_IF0_sort1 <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_IF0_sort1)) {
            M_IF0_sort1[i,] <- sort(M_IF0Sob[i,])
          }

          lowerS <- M_IF0_sort1[, floor(alpha * B / 2)]
          upperS <- M_IF0_sort1[, floor((1 - alpha / 2) * B)]
          M_fU_sort <- matrix(0, nrow = nrow(C), ncol = B)
          for (i in 1:nrow(M_fU_sort)) {
            M_fU_sort[i,] <- sort(M_fU[i,])
          }

          lowerU <- M_fU_sort[, floor(alpha * B / 2)]
          upperU <- M_fU_sort[, floor((1 - alpha / 2) * B)]

        }
        ##############################
        # end of bootstrap trucn LEFT
        ##############################

	if (trunc == "right") {

	########################################
          # PARALELL BOOTSTARP RIGHT
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
			wt <- rep(1 / length(C[,1]), times = length(C[,1]))
			error <- 1e-6
			nmaxit<-10000000
            # container for bootstrap operations

		W1<-as.vector(f)
        W2<-as.vector(k)
        ind<-1:nrow(C)
        M_IF0<-matrix(0,nrow=nrow(C))
        M_IF01<-matrix(0,nrow=nrow(C))
        M_IF0Sob<-matrix(0,nrow=nrow(C))

        kMVb<-matrix(0,nrow=nrow(C))
        M_fV<-matrix(0,nrow=nrow(C))
        M4b<-matrix(0,nrow=nrow(C))
        k1bV<-matrix(0,nrow=nrow(C))

        indbb<-seq(1,nrow(C),by=1)
        indbb1<-seq(1,nrow(C),by=1)

	repeat {

		  DB<-matrix(0,nrow=nrow(C),ncol=3)
          DB[,1]<-sample(C[,1],size=nrow(C),replace=TRUE,prob=W1)
          indUV<-sample(ind,size=nrow(C),replace=TRUE,prob=W2)
          DB[,2:3]<-C[indUV,2:3]

		  log1<-(DB[,1]>=DB[,2])
          log2<-(DB[,1]<=DB[,3])
          indlog<-(log1*log2==0)

          while(sum(indlog)>0){
            DB[indlog,1]<-sample(C[,1],size=sum(indlog),replace=TRUE,prob=W1)
            indUV<-sample(ind,size=sum(indlog),replace=TRUE,prob=W2)
            DB[indlog,2:3]<-C[indUV,2:3]
            log1<-(DB[,1]>=DB[,2])
            log2<-(DB[,1]<=DB[,3])
            indlog<-(log1*log2==0)
          }

		  ordB<-order(DB[,1],method = "auto")
          DBB<-matrix(0,nrow=nrow(C),ncol=3)
          DBB[,1]<-sort(DB[,1],method = "auto")
          DBB[,2:ncol(DB)]<-DB[ordB,2:ncol(DBB)]

          error<-1e-6
          Jb<-matrix(data=0,ncol=nrow(DBB),nrow=nrow(DBB))
          aub <- outer(DBB[, 1], DBB[, 2], ">=")
              avb <- outer(DBB[, 1], DBB[, 3], "<=")
              auub <- outer(DBB[, 2], DBB[, 2], "<=") * 1L
              Jb <- aub * avb
              SCb <- apply(Jb, 2, "sum")
              adb <- length(which(SCb == 1))

              n_sampling_tries <- n_sampling_tries + 1 # counter for number of tries

            # if the condition is set, we break the repat, otherwise repeat sampling  process
            if (adb == 0)
              break
            else
              cat("Warning. Non-uniqueness or no existence of the NPMLE - New Round","\n")
            }
            ## End os sample repeat

          # now we proceed after good sample
            JIb <- t(Jb)

            # preallocation
            f0b <- matrix(data = wt, ncol = 1, nrow = nrow(DBB))
            f1b <- f0b

            S0b <- 1
            iterb <- 0

            # some inner loop
            while (S0b > error | iterb > nmaxit) {
              iterb <- iterb + 1

              if (iterb > nmaxit)
                stop("Default number of iterations not enough for convergence")

              F0b <- JIb %*% f1b
              k0b <- ((sum(1 / F0b)) ^ (-1)) * (1 / F0b)
              if (sum(k0b) != 1)
                k0b <- k0b / sum(k0b)

              k1b <- k0b
              K0b <- Jb %*% k1b
              f1b <- ((sum(1 / K0b)) ^ (-1)) * (1 / K0b)

              if (sum(f1b) != 1)
                f1b <- f1b / sum(f1b)

              S0b <- max(abs(f1b - f0b))
              f0b <- f1b
              k0b <- k1b
            }

            ff0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb1 <- (C[, 1] == C[i, 1])
              pos1 <- min(which(indbb1 == TRUE))

              if (pos1 >= 1)
                ff0b[indbb1] <- sum(f1b[indbb1])

            }

            FF0b <- numeric(nrow(C))

            for (i in 1:nrow(C)) {
              indbb7 <- (C[, 1] == C[i, 1])
              pos0 <- min(which(indbb7 == TRUE))

              if (pos0 == 1)
                FF0b[indbb7] <- sum(f1b[indbb7])

              else if (pos0 > 1)
                FF0b[indbb7] <- sum(f1b[indbb7]) + FF0b[pos0 - 1]
            }

            Sobb <- 1 - FF0b + ff0b
            Sobb[Sobb < 1e-12] <- 0


			Sobb <- 1 - FF0b + ff0b
            Sobb[Sobb < 1e-12] <- 0

            # SAVE start saving stufs (before indexed to b index)
            M_IF0 <- as.vector(FF0b)
            M_IF01 <- as.vector(f1b)
            M_IF0Sob <- as.vector(Sobb)


            ordB <- order( DBB[, 3], method = "auto")
            ord4 <- sort( DBB[, 3], method = "auto")

            # end ordenaçao or ordering zone
            ########################################



            M4b <- ord4
            #k1bV <- k1b[ordV]
			k1bV <- k1b[ordB]

            # check this was two loops
            for (i in 1:nrow(C)) {
              for (j in 1:nrow(C)) {
                if (M4b[i] <= kMV[j, 1])
                  M_fV[j] <- sum(k1bV[1:i])
              }
            }


		# variales to be outputed of the clusters (paralell)

		# variales to be outputed of the clusters (paralell)

          return(list( M_IF0, M_IF01, M_IF0Sob, M_fV,n_sampling_tries))

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
          M_fV <- matrix(0, nrow = nrow(C), ncol = B)
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
			M_fV[,b]<-as.numeric(unlist(final_boot[b,4]))
			BootRepeat[,b]<-as.numeric((unlist(final_boot[b,5])))
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


        M_fV_sort<-matrix(0,nrow=nrow(C),ncol=B)
        for(i in 1:nrow(M_fV_sort)){
          M_fV_sort[i,]<-sort(M_fV[i,])
        }
        lowerV<-M_fV_sort[,floor(alpha*B/2)]
        upperV<-M_fV_sort[,floor((1-alpha/2)*B)]

        }

        ##############################
        # end of bootstrap trucn RIGTH
        ##############################

      }##end of obvious



	}## end of boot

 if (boot==TRUE){

    if(trunc=="both"){

	  x<-C[,1]
	  UU<-sort(C[,2])
	  VV<-sort(C[,3])
      if((display.FS==TRUE)&(display.UV==TRUE)){
        dev.new()
        par(mfrow=c(2,2))

        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
        segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
        segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],upperF[i],C[i+1,1],upperF[i], lty=3)
          segments(C[i+1,1],upperF[i],C[i+1,1],upperF[i+1],lty=3)
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
        segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
        segments(C[,1][1],0,C[,1][1],lowerF[1], lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],lowerF[i],C[i+1,1],lowerF[i], lty=3)
          segments(C[i+1,1],lowerF[i],C[i+1,1],lowerF[i+1],lty=3)
        }


        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
        }

        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),1,C[,1][1],1,lty=3)
        segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
        segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],upperS[i],C[i+1,1],upperS[i], lty=3)
          segments(C[i+1,1],upperS[i],C[i+1,1],upperS[i+1],lty=3)
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),1,C[,1][1],1,lty=3)
        segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
        segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],lowerS[i],C[i+1,1],lowerS[i], lty=3)
          segments(C[i+1,1],lowerS[i],C[i+1,1],lowerS[i+1],lty=3)
        }




        plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

        segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
        segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
        segments(UU[1],0,UU[1],fU[1])

        for(i in 1:(length(UU)-1)){
          segments(UU[i],fU[i],UU[i+1],fU[i])
          segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
        }

        CC<-sort(C[,2])

        segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
        segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
        segments(sort(C[,2])[1],0,sort(C[,2])[1],upperU[1], lty=3)

        for(i in 1:(length(C[,2])-1)){
          segments(CC[i],upperU[i],CC[i+1],upperU[i], lty=3)
          segments(CC[i+1],upperU[i],CC[i+1],upperU[i+1],lty=3)
        }


        segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
        segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
        segments(sort(C[,2])[1],0,sort(C[,2])[1],lowerU[1], lty=3)

        for(i in 1:(length(C[,2])-1)){
          segments(CC[i],lowerU[i],CC[i+1],lowerU[i], lty=3)
          segments(CC[i+1],lowerU[i],CC[i+1],lowerU[i+1],lty=3)
        }


        plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

        segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
        segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
        segments(VV[1],0,VV[1],fV[1])

        for(i in 1:(length(VV)-1)){
          segments(VV[i],fV[i],VV[i+1],fV[i])
          segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
        }


        CCC<-sort(C[,3])

        segments(min(C[,3])-(max(C[,3])-min(C[,3]))/length(C[,3]),0,sort(C[,3])[1],0,lty=3)
        segments(max(C[,3]),1,max(C[,3])+(max(C[,3])-min(C[,3]))/length(C[,3]),1,lty=3)
        segments(sort(C[,3])[1],0,sort(C[,3])[1],upperV[1], lty=3)

        for(i in 1:(length(C[,3])-1)){
          segments(CCC[i],upperV[i],CCC[i+1],upperV[i], lty=3)
          segments(CCC[i+1],upperV[i],CCC[i+1],upperV[i+1],lty=3)
        }


        segments(min(C[,3])-(max(C[,3])-min(C[,3]))/length(C[,3]),0,sort(C[,3])[1],0,lty=3)
        segments(max(C[,3]),1,max(C[,3])+(max(C[,3])-min(C[,3]))/length(C[,3]),1,lty=3)
        segments(sort(C[,3])[1],0,sort(C[,3])[1],lowerV[1], lty=3)

        for(i in 1:(length(C[,3])-1)){
          segments(CCC[i],lowerV[i],CCC[i+1],lowerV[i], lty=3)
          segments(CCC[i+1],lowerV[i],CCC[i+1],lowerV[i+1],lty=3)
        }

      }

      if((display.FS==FALSE)&(display.UV==TRUE)){
        dev.new()
        par(mfrow=c(1,2))
        plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

        segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
        segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
        segments(UU[1],0,UU[1],fU[1])

        for(i in 1:(length(UU)-1)){
          segments(UU[i],fU[i],UU[i+1],fU[i])
          segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
        }

        CC<-sort(C[,2])

        segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
        segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
        segments(sort(C[,2])[1],0,sort(C[,2])[1],upperU[1], lty=3)

        for(i in 1:(length(C[,2])-1)){
          segments(CC[i],upperU[i],CC[i+1],upperU[i], lty=3)
          segments(CC[i+1],upperU[i],CC[i+1],upperU[i+1],lty=3)
        }


        segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
        segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
        segments(sort(C[,2])[1],0,sort(C[,2])[1],lowerU[1], lty=3)

        for(i in 1:(length(C[,2])-1)){
          segments(CC[i],lowerU[i],CC[i+1],lowerU[i], lty=3)
          segments(CC[i+1],lowerU[i],CC[i+1],lowerU[i+1],lty=3)
        }


        plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

        segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
        segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
        segments(VV[1],0,VV[1],fV[1])

        for(i in 1:(length(VV)-1)){
          segments(VV[i],fV[i],VV[i+1],fV[i])
          segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
        }


        CCC<-sort(C[,3])

        segments(min(C[,3])-(max(C[,3])-min(C[,3]))/length(C[,3]),0,sort(C[,3])[1],0,lty=3)
        segments(max(C[,3]),1,max(C[,3])+(max(C[,3])-min(C[,3]))/length(C[,3]),1,lty=3)
        segments(sort(C[,3])[1],0,sort(C[,3])[1],upperV[1], lty=3)

        for(i in 1:(length(C[,3])-1)){
          segments(CCC[i],upperV[i],CCC[i+1],upperV[i], lty=3)
          segments(CCC[i+1],upperV[i],CCC[i+1],upperV[i+1],lty=3)
        }


        segments(min(C[,3])-(max(C[,3])-min(C[,3]))/length(C[,3]),0,sort(C[,3])[1],0,lty=3)
        segments(max(C[,3]),1,max(C[,3])+(max(C[,3])-min(C[,3]))/length(C[,3]),1,lty=3)
        segments(sort(C[,3])[1],0,sort(C[,3])[1],lowerV[1], lty=3)

        for(i in 1:(length(C[,3])-1)){
          segments(CCC[i],lowerV[i],CCC[i+1],lowerV[i], lty=3)
          segments(CCC[i+1],lowerV[i],CCC[i+1],lowerV[i+1],lty=3)
        }
      }
      if((display.FS==TRUE)&(display.UV==FALSE)){
        dev.new()
        par(mfrow=c(1,2))
        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
        segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
        segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],upperF[i],C[i+1,1],upperF[i], lty=3)
          segments(C[i+1,1],upperF[i],C[i+1,1],upperF[i+1],lty=3)
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
        segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
        segments(C[,1][1],0,C[,1][1],lowerF[1], lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],lowerF[i],C[i+1,1],lowerF[i], lty=3)
          segments(C[i+1,1],lowerF[i],C[i+1,1],lowerF[i+1],lty=3)
        }


        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
        }

        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),1,C[,1][1],1,lty=3)
        segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
        segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],upperS[i],C[i+1,1],upperS[i], lty=3)
          segments(C[i+1,1],upperS[i],C[i+1,1],upperS[i+1],lty=3)
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),1,C[,1][1],1,lty=3)
        segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
        segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],lowerS[i],C[i+1,1],lowerS[i], lty=3)
          segments(C[i+1,1],lowerS[i],C[i+1,1],lowerS[i+1],lty=3)
        }

      }


      if((display.FS==FALSE)&(display.UV==FALSE)){


      }
    }





    if(trunc=="left"){

	  X<-C[,1]
	  UU<-sort(C[,2])
      if((display.FS==TRUE)&(display.UV==TRUE)){
        dev.new()
        par(mfrow=c(2,2))
        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="EP estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
        segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
        segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],upperF[i],C[i+1,1],upperF[i], lty=3)
          segments(C[i+1,1],upperF[i],C[i+1,1],upperF[i+1],lty=3)
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
        segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
        segments(C[,1][1],0,C[,1][1],lowerF[1], lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],lowerF[i],C[i+1,1],lowerF[i], lty=3)
          segments(C[i+1,1],lowerF[i],C[i+1,1],lowerF[i+1],lty=3)
        }


        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
        }

        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(C[,1]),1,C[,1][1],1,lty=3)
        segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
        segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],upperS[i],C[i+1,1],upperS[i], lty=3)
          segments(C[i+1,1],upperS[i],C[i+1,1],upperS[i+1],lty=3)
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),1,C[,1][1],1,lty=3)
        segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
        segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],lowerS[i],C[i+1,1],lowerS[i], lty=3)
          segments(C[i+1,1],lowerS[i],C[i+1,1],lowerS[i+1],lty=3)
        }




        plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

        segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
        segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
        segments(UU[1],0,UU[1],fU[1])

        for(i in 1:(length(UU)-1)){
          segments(UU[i],fU[i],UU[i+1],fU[i])
          segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
        }

        CC<-sort(C[,2])

        segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
        segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
        segments(sort(C[,2])[1],0,sort(C[,2])[1],upperU[1], lty=3)

        for(i in 1:(length(C[,2])-1)){
          segments(CC[i],upperU[i],CC[i+1],upperU[i], lty=3)
          segments(CC[i+1],upperU[i],CC[i+1],upperU[i+1],lty=3)
        }


        segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
        segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
        segments(sort(C[,2])[1],0,sort(C[,2])[1],lowerU[1], lty=3)

        for(i in 1:(length(C[,2])-1)){
          segments(CC[i],lowerU[i],CC[i+1],lowerU[i], lty=3)
          segments(CC[i+1],lowerU[i],CC[i+1],lowerU[i+1],lty=3)
        }


      }

      if((display.FS==FALSE)&(display.UV==TRUE)){
        dev.new()
        par(mfrow=c(1,1))
        plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

        segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
        segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
        segments(UU[1],0,UU[1],fU[1])

        for(i in 1:(length(UU)-1)){
          segments(UU[i],fU[i],UU[i+1],fU[i])
          segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
        }

        CC<-sort(C[,2])

        segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
        segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
        segments(sort(C[,2])[1],0,sort(C[,2])[1],upperU[1], lty=3)

        for(i in 1:(length(C[,2])-1)){
          segments(CC[i],upperU[i],CC[i+1],upperU[i], lty=3)
          segments(CC[i+1],upperU[i],CC[i+1],upperU[i+1],lty=3)
        }


        segments(min(C[,2])-(max(C[,2])-min(C[,2]))/length(C[,2]),0,sort(C[,2])[1],0,lty=3)
        segments(max(C[,2]),1,max(C[,2])+(max(C[,2])-min(C[,2]))/length(C[,2]),1,lty=3)
        segments(sort(C[,2])[1],0,sort(C[,2])[1],lowerU[1], lty=3)

        for(i in 1:(length(C[,2])-1)){
          segments(CC[i],lowerU[i],CC[i+1],lowerU[i], lty=3)
          segments(CC[i+1],lowerU[i],CC[i+1],lowerU[i+1],lty=3)
        }


      }

      if((display.FS==TRUE)&(display.UV==FALSE)){
        dev.new()
        par(mfrow=c(1,2))
        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
        segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
        segments(C[,1][1],0,C[,1][1],upperF[1], lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],upperF[i],C[i+1,1],upperF[i], lty=3)
          segments(C[i+1,1],upperF[i],C[i+1,1],upperF[i+1],lty=3)
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),0,C[,1][1],0,lty=3)
        segments(max(C[,1]),1,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),1,lty=3)
        segments(C[,1][1],0,C[,1][1],lowerF[1], lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],lowerF[i],C[i+1,1],lowerF[i], lty=3)
          segments(C[i+1,1],lowerF[i],C[i+1,1],lowerF[i+1],lty=3)
        }


        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
        }

        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),1,C[,1][1],1,lty=3)
        segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
        segments(max(C[,1]),0,max(C[,1]),min(upperS), lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],upperS[i],C[i+1,1],upperS[i], lty=3)
          segments(C[i+1,1],upperS[i],C[i+1,1],upperS[i+1],lty=3)
        }


        segments(min(C[,1])-(max(C[,1])-min(C[,1]))/length(x),1,C[,1][1],1,lty=3)
        segments(max(C[,1]),0,max(C[,1])+(max(C[,1])-min(C[,1]))/length(C[,1]),0,lty=3)
        segments(max(C[,1]),0,max(C[,1]),min(lowerS), lty=3)

        for(i in 1:(length(C[,1])-1)){
          segments(C[i,1],lowerS[i],C[i+1,1],lowerS[i], lty=3)
          segments(C[i+1,1],lowerS[i],C[i+1,1],lowerS[i+1],lty=3)
        }

      }


      if((display.FS==FALSE)&(display.UV==FALSE)){


      }
    }




    if(trunc=="right"){
      x<-C[,1]
	  VV<-sort(C[,3])
      if((display.FS==TRUE)&(display.UV==TRUE)){
        dev.new()
        par(mfrow=c(2,2))

        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
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

        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
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



        plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

        segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
        segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
        segments(VV[1],0,VV[1],fV[1])

        for(i in 1:(length(VV)-1)){
          segments(VV[i],fV[i],VV[i+1],fV[i])
          segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
        }

        segments(min(V)-(max(V)-min(V))/length(V),0,sort(V)[1],0,lty=3)
        segments(max(V),1,max(V)+(max(V)-min(V))/length(V),1,lty=3)
        segments(sort(V)[1],0,sort(V)[1],upperV[1], lty=3)

        CCCC<-sort(V)

        for(i in 1:(length(V)-1)){
          segments(CCCC[i],upperV[i],CCCC[i+1],upperV[i], lty=3)
          segments(CCCC[i+1],upperV[i],CCCC[i+1],upperV[i+1],lty=3)
        }


        segments(min(V)-(max(V)-min(V))/length(V),0,CCCC[1],0,lty=3)
        segments(max(V),1,max(V)+(max(V)-min(V))/length(V),1,lty=3)
        segments(CCCC[1],0,CCCC[1],lowerV[1], lty=3)

        for(i in 1:(length(V)-1)){
          segments(CCCC[i],lowerV[i],CCCC[i+1],lowerV[i], lty=3)
          segments(CCCC[i+1],lowerV[i],CCCC[i+1],lowerV[i+1],lty=3)
        }


      }

      if((display.FS==FALSE)&(display.UV==TRUE)){
        dev.new()
        par(mfrow=c(1,1))
        plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

        segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
        segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
        segments(VV[1],0,VV[1],fV[1])

        for(i in 1:(length(VV)-1)){
          segments(VV[i],fV[i],VV[i+1],fV[i])
          segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
        }

        segments(min(V)-(max(V)-min(V))/length(V),0,sort(V)[1],0,lty=3)
        segments(max(V),1,max(V)+(max(V)-min(V))/length(V),1,lty=3)
        segments(sort(V)[1],0,sort(V)[1],upperV[1], lty=3)

        CCCC<-sort(V)

        for(i in 1:(length(V)-1)){
          segments(CCCC[i],upperV[i],CCCC[i+1],upperV[i], lty=3)
          segments(CCCC[i+1],upperV[i],CCCC[i+1],upperV[i+1],lty=3)
        }


        segments(min(V)-(max(V)-min(V))/length(V),0,CCCC[1],0,lty=3)
        segments(max(V),1,max(V)+(max(V)-min(V))/length(V),1,lty=3)
        segments(CCCC[1],0,CCCC[1],lowerV[1], lty=3)

        for(i in 1:(length(V)-1)){
          segments(CCCC[i],lowerV[i],CCCC[i+1],lowerV[i], lty=3)
          segments(CCCC[i+1],lowerV[i],CCCC[i+1],lowerV[i+1],lty=3)
        }


      }
      if((display.FS==TRUE)&(display.UV==FALSE)){
        dev.new()
        par(mfrow=c(1,2))
        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
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

        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
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

      }

      if((display.FS==FALSE)&(display.UV==FALSE)){

      }
    }

  }











  if (boot==FALSE|B<40){
    if(trunc=="both"){
       x<-C[,1]
	   UU<-sort(C[,2])
	   VV<-sort(C[,3])
      if(plot.joint==TRUE){
        if(plot.type=="image"){
          sU<-sort(U)
          sV<-sort(V)
          if(any(diff(sU)==0)& all(diff(sV)!=0)){
            stepU<-diff(range(U))/(length(U)*length(unique(U)))
            image(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution")
            filled.contour(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
            filled.contour(sort(U+runif(length(U),0,stepU)),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
          }

          if(all(diff(sU)!=0)& any(diff(sV)==0)){
            stepV<-diff(range(V))/(length(V)*length(unique(V)))
            image(sort(U),sort(V+runif(length(U),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution")
            filled.contour(sort(U),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
            filled.contour(sort(U),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
          }

          if(any(diff(sU)==0)& any(diff(sV)==0)){
            stepU<-diff(range(U))/(length(U)*length(unique(U)))
            stepV<-diff(range(V))/(length(V)*length(unique(V)))
            image(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(U),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution")
            filled.contour(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
            filled.contour(sort(U+runif(length(U),0,stepU)),sort(V+runif(length(V),0,stepV)),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
          }

          if(all(diff(sU)!=0)& all(diff(sV)!=0)){
            image(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution")
            filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
            filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
          }}

        if(plot.type=="persp"){
          fcol<-topo.colors(10)[cut(KK[2:nrow(C),2:nrow(C)],10,include.lowest=TRUE)]
          sU<-sort(U)
          sV<-sort(V)
          if(any(diff(sU)==0)& all(diff(sV)!=0)){
            stepU<-diff(range(U))/(length(U)*length(unique(U)))
            persp(x=sort(U+runif(length(U),0,stepU)),y=sort(V),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
          }
          if(all(diff(sU)!=0)& any(diff(sV)==0)){
            stepV<-diff(range(V))/(length(V)*length(unique(V)))
            persp(x=sort(U),y=sort(V+runif(length(V),0,stepV)),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
          }
          if(any(diff(sU)==0)& any(diff(sV)==0)){
            stepU<-diff(range(U))/(length(U)*length(unique(U)))
            stepV<-diff(range(V))/(length(V)*length(unique(V)))
            persp(x=sort(U+runif(length(U),0,stepU)),y=sort(V+runif(length(V),0,stepV)),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
          }
          if(all(diff(sU)!=0)& all(diff(sV)!=0)){
            persp(x=sort(U),y=sort(V),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")
          }
        }


      }
      if(plot.joint==FALSE){
      }




      if((display.FS==TRUE)&(display.UV==TRUE)){
        dev.new()
        par(mfrow=c(2,2))
        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
        }



        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
        }


        plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

        segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
        segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
        segments(UU[1],0,UU[1],fU[1])

        for(i in 1:(length(UU)-1)){
          segments(UU[i],fU[i],UU[i+1],fU[i])
          segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
        }


        plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

        segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
        segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
        segments(VV[1],0,VV[1],fV[1])

        for(i in 1:(length(VV)-1)){
          segments(VV[i],fV[i],VV[i+1],fV[i])
          segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
        }

      }

      if((display.FS==FALSE)&(display.UV==TRUE)){
        dev.new()
        par(mfrow=c(1,2))
        plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

        segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
        segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
        segments(UU[1],0,UU[1],fU[1])

        for(i in 1:(length(UU)-1)){
          segments(UU[i],fU[i],UU[i+1],fU[i])
          segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
        }


        plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

        segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
        segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
        segments(VV[1],0,VV[1],fV[1])

        for(i in 1:(length(VV)-1)){
          segments(VV[i],fV[i],VV[i+1],fV[i])
          segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
        }

      }

      if((display.FS==TRUE)&(display.UV==FALSE)){
        dev.new()
        par(mfrow=c(1,2))
        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
        }



        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
        }
      }


      if((display.FS==FALSE)&(display.UV==FALSE)){


      }
    }


    if(trunc=="left"){
	UU<-sort(C[,2])
      x<-C[,1]
      if((display.FS==TRUE)&(display.UV==TRUE)){
        dev.new()
        par(mfrow=c(2,2))
        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
        }



        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
        }


        plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

        segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
        segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
        segments(UU[1],0,UU[1],fU[1])

        for(i in 1:(length(UU)-1)){
          segments(UU[i],fU[i],UU[i+1],fU[i])
          segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
        }

      }

      if((display.FS==FALSE)&(display.UV==TRUE)){
        dev.new()
        par(mfrow=c(1,1))
        plot(UU,fU,ylim=c(0,1),xlim=c(min(UU)-(max(UU)-min(UU))/length(UU),max(UU)+(max(UU)-min(UU))/length(UU)),type="n",main="Marginal U", xlab="U",ylab="")

        segments(min(UU)-(max(UU)-min(UU))/length(UU),0,UU[1],0)
        segments(max(UU),1,max(UU)+(max(UU)-min(UU))/length(UU),1)
        segments(UU[1],0,UU[1],fU[1])

        for(i in 1:(length(UU)-1)){
          segments(UU[i],fU[i],UU[i+1],fU[i])
          segments(UU[i+1],fU[i],UU[i+1],fU[i+1])
        }
      }

      if((display.FS==TRUE)&(display.UV==FALSE)){
        dev.new()
        par(mfrow=c(1,2))
        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
        }



        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
        }

      }


      if((display.FS==FALSE)&(display.UV==FALSE)){


      }
    }


    if(trunc=="right"){
      x<-C[,1]
	  VV<-sort(C[,3])
      if((display.FS==TRUE)&(display.UV==TRUE)){

        dev.new()
        par(mfrow=c(2,2))
        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
        }


        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
        }


        plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

        segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
        segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
        segments(VV[1],0,VV[1],fV[1])

        for(i in 1:(length(VV)-1)){
          segments(VV[i],fV[i],VV[i+1],fV[i])
          segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
        }

      }

      if((display.FS==FALSE)&(display.UV==TRUE)){
        dev.new()
        par(mfrow=c(1,1))
        plot(VV,fV,ylim=c(0,1),xlim=c(min(VV)-(max(VV)-min(VV))/length(VV),max(VV)+(max(VV)-min(VV))/length(VV)),type="n",main="Marginal V", xlab="V",ylab="")

        segments(min(VV)-(max(VV)-min(VV))/length(VV),0,VV[1],0)
        segments(max(VV),1,max(VV)+(max(VV)-min(VV))/length(VV),1)
        segments(VV[1],0,VV[1],fV[1])

        for(i in 1:(length(VV)-1)){
          segments(VV[i],fV[i],VV[i+1],fV[i])
          segments(VV[i+1],fV[i],VV[i+1],fV[i+1])
        }
      }
      if((display.FS==TRUE)&(display.UV==FALSE)){
        dev.new()
        par(mfrow=c(1,2))
		x<-C[,1]
        plot(x,FFF,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Shen estimator", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),0,x[1],0)
        segments(max(x),1,max(x)+(max(x)-min(x))/length(x),1)
        segments(x[1],0,x[1],FFF[1])

        for(i in 1:(length(x)-1)){
          segments(x[i],FFF[i],x[i+1],FFF[i])
          segments(x[i+1],FFF[i],x[i+1],FFF[i+1])
        }


        plot(x,Sob0,ylim=c(0,1),xlim=c(min(x)-(max(x)-min(x))/length(x),max(x)+(max(x)-min(x))/length(x)),type="n",main="Survival", xlab="Time of interest",ylab="")

        segments(min(x)-(max(x)-min(x))/length(x),1,x[1],1)
        segments(max(x),0,max(x)+(max(x)-min(x))/length(x),0)
        segments(max(x),0,max(x),min(Sob0))


        for(i in 1:(length(x)-1)){
          segments(x[i],Sob0[i],x[i+1],Sob0[i])
          segments(x[i+1],Sob0[i],x[i+1],Sob0[i+1])
        }

      }

      if((display.FS==FALSE)&(display.UV==FALSE)){

      }
    }


  }



	cat("n.iterations",iter,"\n")
  cat("S0",S0,"\n")
  cat("events",events,"\n")


	if(boot==TRUE){

    cat("B",B,"\n")
    cat("alpha",alpha,"\n")
    cat("Boot",boot.type,"\n")


    if(trunc=="both"){
      return(invisible(list(n.iterations=iter,events=events, B=B, alpha=alpha,Boot=boot.type,time=C[,1], n.event=mult, density=round(as.vector(f),5), cumulative.df=round(FFF,5), survival=
                              round(as.vector(Sob0),5), truncation.probs=round(as.vector(F0),5), biasf=round(biasf,5), upper.df=round(upperF,5),lower.df=round(lowerF,5),upper.Sob=round(upperS,5),
                            lower.Sob=round(lowerS,5), density.joint=round(as.vector(k),5),marginal.U=round(fU,5), marginal.V=round(fV,5), upper.fU=round(upperU,5), lower.fU=round(lowerU,5),upper.fV=round(upperV,5), lower.fV=round(lowerV,5),cumulative.joint=round(KK,5), sd.boot=round(stderror,5),Boot.Repeat=as.vector(BootRepeat) )))

    }

    if(trunc=="left"){
      return(invisible(list(n.iterations=iter,events=events, B=B, alpha=alpha,Boot=boot.type,time=C[,1], n.event=mult, density=round(as.vector(f),5), cumulative.df=round(FFF,5), survival=
                              round(as.vector(Sob0),5), truncation.probs=round(as.vector(F0),5), biasf=round(biasf,5), upper.df=round(upperF,5),lower.df=round(lowerF,5),upper.Sob=round(upperS,5),
                            lower.Sob=round(lowerS,5), density.joint=round(as.vector(ku),5),marginal.U=round(fU,5),upper.fU=round(upperU,5), lower.fU=round(lowerU,5),sd.boot=round(stderror,5),Boot.Repeat=as.vector(BootRepeat))))
    }
    if(trunc=="right"){

      return(invisible(list(n.iterations=iter,events=events, B=B, alpha=alpha,Boot=boot.type,time=C[,1], n.event=mult, density=round(as.vector(f),5), cumulative.df=round(FFF,5), survival=
                              round(as.vector(Sob0),5), truncation.probs=round(as.vector(F0),5), biasf=round(biasf,5), upper.df=round(upperF,5),lower.df=round(lowerF,5),upper.Sob=round(upperS,5),
                            lower.Sob=round(lowerS,5), density.joint=round(as.vector(kv),5), marginal.V=round(fV,5), upper.fV=round(upperV,5), lower.fV=round(lowerV,5),sd.boot=round(stderror,5),Boot.Repeat=as.vector(BootRepeat) )))

    }
  }



	 if(boot==FALSE){

    if(trunc=="both"){
      return(invisible(list(n.iterations=iter,events=events, time=C[,1], n.event=mult, density=round(as.vector(f),5), cumulative.df=round(FFF,5), survival=
                              round(as.vector(Sob0),5), truncation.probs=round(as.vector(F0),5), biasf=round(biasf,5), density.joint=round(as.vector(k),5),marginal.U=round(fU,5), marginal.V=round(fV,5), cumulative.joint=round(KK,5) )))

    }


	if(trunc=="left"){
      return(invisible(list(n.iterations=iter,events=events, time=C[,1], n.event=mult, density=round(as.vector(f),5), cumulative.df=round(FFF,5), survival=
                              round(as.vector(Sob0),5), truncation.probs=round(as.vector(F0),5),biasf=round(biasf,5), density.joint=round(as.vector(k),5), marginal.U=round(fU,5))))
    }

	if(trunc=="right"){

      return(invisible(list(n.iterations=iter,events=events, time=x, n.event=mult, density=round(as.vector(f),5), cumulative.df=round(FFF,5), survival=
                              round(as.vector(Sob0),5), truncation.probs=round(as.vector(F0),5), biasf=round(biasf,5), density.joint=round(as.vector(k),5), marginal.V=round(fV,5) )))

    }


	}


	}
