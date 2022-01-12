

densityDT <-function(X,U,V, bw="DPI2", from,to,n, wg=NA){

  CDF <-function(X, U=NA, V=NA, wt=NA, error=NA, nmaxit=NA){

    if(any(is.na(U))==TRUE | any(is.na(V))==TRUE|any(is.na(X))==TRUE){
      navec<-c(which(is.na(X)),which(is.na(U)),which(is.na(V)))
      X<-X[-navec]
      U<-U[-navec]
      V<-V[-navec]

    }
    if(is.na(wt)==TRUE) wt<-rep(1/length(X),times=length(X))

    D<-cbind(X,U,V)
    ord<-order(D[,1])
    C<-matrix(0,nrow=nrow(D),ncol=ncol(D))
    C[,1]<-sort(D[,1])
    C[,2:ncol(D)]<-D[ord,2:ncol(D)]

    T<-table(C[,2]<= C[,1]&C[,1]<=C[,3])
     if(sum(T[names(T)=="TRUE"])!=length(X) |sum(T[names(T)=="TRUE"]==0)){
      stop("Condition of double truncation is violated","\n")
    }



    if(is.na(error)==TRUE) error<-1e-6

    au<-outer(C[,1],C[,2],">=")
    av<-outer(C[,1],C[,3],"<=")
    auu<-outer(C[,2],C[,2],"<=")*1L

    J<-matrix(data=0,ncol=nrow(C),nrow=nrow(C))
    J<-au*av
    SC<-colSums(J)
    ad<-length(which(SC==1))
    if(ad!=0){
      cat("Warning. Non-uniqueness or no existence of the NPMLE","\n")
    }
    JI<-t(J)
    f0<-matrix(data=1/nrow(C),ncol=1,nrow=nrow(C))
    f<-f0
    k<-rep(1,times=length(f))
    S0<-1
    S1<-1
    if(is.na(nmaxit)==TRUE)nmaxit<-100000000000000000
    iter<-0
    while(( S0>error |S1>error )| iter>nmaxit){
      iter<-iter+1
      if (iter>nmaxit) stop("Default number of iterations not enough for convergence")
      F0<-JI%*%f
      k0<-((sum(1/F0))^(-1))*(1/F0)
      if(sum(k0)!=1)k0<-k0/sum(k0)
      K0<-J%*%k0
      f<-((sum(1/K0))^(-1))*(1/K0)
      if(sum(f)!=1)f<-f/sum(f)
      S0<-max(abs(f-f0))
      f0<-f
      S1<-max(abs(k-k0))
      k<-k0

    }

    F0<-JI%*%f
    K0<-J%*%k0

    k<-k0

    mult <- tabulate(match(C[,1],unique(C[,1])))

    x<-unique(C[,1])
    events<-sum(mult)
    n.event<-mult
    FFF<-cumsum(f)
    Sob0<-1-FFF
    Sob0[Sob0<1e-12]<-0


    uu<-unique(C[,2])
    KK<-apply(auu*as.vector(k),2,"sum")
    kMUV<-cbind(C[,2],C[,3],k)

    kuv<-numeric(nrow(C))
    kuv<-sum(k)

    Gf<-matrix(data=k,ncol=ncol(J),nrow=nrow(J), byrow=T)
    Gff<-J*Gf
    biasf<-apply(Gff,1,"sum")
    return(invisible(list(n.iterations=iter,events=events, time=C[,1], n.event=mult, density=as.vector(f), cumulative.df=round(FFF,5), survival=
                            round(as.vector(Sob0),5), truncation.probs=round(as.vector(F0),5), biasf=biasf, density.joint=round(as.vector(kuv),5), k=k , cumulative.dk=round(KK,5), survival0=
                            round(as.vector(Sob0),5))))
  }

  kdens<-function(data,h,xx,wg){
    w<-dnorm(apply(outer(xx,data,"-"),2,function(xx)outer(xx,h,"/")))
    Mw<-matrix(rep(wg,times=(length(xx)*length(h))),ncol=(length(xx)*length(h)),nrow=length(wg))
    sw<-as.vector(w%*%Mw[,1])
    return(1/h*t(matrix(sw,ncol=length(h))))
  }


  CVfunEP=function(h,data, wg,alfanfEP,gfuncEP){
    t1<-sapply(1:length(h), function(i) {integrate(function(xx) kdens(data, h[i],xx,wg)^2,lower=min(xx), upper= max(xx), subdivisions = 2000)$value})
    #t1<-integrate(function(x) kdens(data, h,x,wg)^2,lower=min(x), upper= max(x), abs.tol=1e-4)$value
    r<-sapply(1:length(data),function(i){kdens(data[-i], h,data[i],wg[-i])})
    t2<-sapply(1: nrow(r), function (i) t2<-2*alfanfEP*sum(r[i,]*(gfuncEP)^(-1))/(length(data)-1))
    return(t1-t2)
  }

  trap_rule<-function(c,fg,fh){
    return(c/2 * (ISE_int(fh[1],fg[1]) + ISE_int(fh[length(fh)],fg[length(fg)])) +
             c * sum(ISE_int(fh[2:(length(fh)-1)],fg[2:(length(fg)-1)])) )
  }


  ISE_int<-function(f,fh){
    return((f-fh)^2)
  }


if(any(is.na(U))==TRUE | any(is.na(V))==TRUE|any(is.na(X))==TRUE){
      navec<-c(which(is.na(X)),which(is.na(U)),which(is.na(V)))
      X<-X[-navec]
      U<-U[-navec]
      V<-V[-navec]
	wg<-wg[-navec]

    }



  vxM<-CDF(X,U,V)
  ifelse(is.na(wg)==TRUE, wg<-WgEP<-vxM[[5]], wg)
  ifelse(missing(from)==TRUE, from<-min(X)+0.0001, from)
  ifelse(missing(to)==TRUE, to<-max(X)-0.0001, to)
  ifelse(missing(n)==TRUE, n<-500, n)


    ord<-order(X)
    X<-sort(X)
    U<-U[ord]
    V<-V[ord]


  gfuncEP<-biasfM<-vxM[[9]]
  f2hatM<-vxM[[10]]
  alfanfEP<-((1/length(X))*sum(1/gfuncEP))^(-1)
  gnEP<-sum(gfuncEP^(-2))
  Mdata<-X

  data<-X
  xx<-seq(from+0.0001,to,length=n)
  vec<-xx

  ####################
  #bandwidth selectors
  ####################


  if(is.numeric(bw)){
    densities<-as.vector(kdens(data,bw,xx,wg))
  }


  dx1<-outer(Mdata,Mdata,"-")
  rr1<-outer(1/gfuncEP,1/gfuncEP, "*")

  HPhi<-function(x,h){
    x<-x/h
    res<-(x^6-15*x^4+45*x^2-15)*dnorm(x)
    return(res)
  }

  HPhi2<-function(x,h){
    x<-x/h
    res<-(x^4-6*x^2+3)*dnorm(x)
    return(res)
  }

  wgemp<-1/length(X)
  numfuncEP<-alfanfEP^2*wgemp*gnEP
  meandt<-wgemp*alfanfEP*sum(Mdata*(biasfM)^(-1) )
  standt<-sapply(1:length(X), function(i) standt<-((Mdata[i]-meandt)^2 *(gfuncEP[i])^(-1)))
  standtt<-sum(standt)
  stand<-sqrt(wgemp*alfanfEP*standtt)

  h<-seq(0.001,max(X), by=(max(X)-min(X))/1000)



  if (bw == "DPI2") {

    Psi8NREP<- (105)/(32*pi^(1/2)*stand^9)
    K6zero<--15*(2*pi)^(-1/2)
    g1mEP<-((-2*numfuncEP*K6zero)/(Psi8NREP*length(X)))^(1/9)
    hatqsi6g1EP<-sum(sapply(dx1,HPhi,g1mEP)*rr1)/(alfanfEP^(-2)*length(X)^2*g1mEP^7)
    K4zero<-3*(2*pi)^(-1/2)
    g2mEP<-((-2*numfuncEP*K4zero)/(hatqsi6g1EP*length(X)))^(1/7)
    hatqsi4g2EP<-sum(sapply(dx1,HPhi2,g2mEP)*rr1)/(alfanfEP^(-2)*n^2*g2mEP^5)
    hopt_DPI2<-(((2*pi^(1/2))^(-1)*numfuncEP)/(hatqsi4g2EP*n))^(1/5)

    densitiesDPI2<-as.vector(kdens(data,hopt_DPI2,xx,wg))

  }


  if (bw == "NR") {

    F2<-sapply(1:length(vxM[[6]]), function(i) F2<-(vxM[[6]][i]==vxM[[6]][floor(length(X)*(1-0.25))+1]))
    F1<-sapply(1:length(vxM[[6]]), function(i) F1<-(vxM[[6]][i]==vxM[[6]][floor(length(X)*0.25)+1]))
    IQR1<-vxM[[3]][F2==TRUE]
    IQR2<-vxM[[3]][F1==TRUE]
    IQRm<-IQR1-IQR2
    hoptNR<-(4/3)^(1/5)*min(stand,IQRm/1.349)*(numfuncEP)^(1/5)*length(X)^(-1/5)

    densitiesNR<-as.vector(kdens(data,hoptNR,xx,wg))
  }



  if (bw == "LSCV") {

    Score_EP<-CVfunEP(h, Mdata, wg,alfanfEP,gfuncEP)
    ind_EP<-which.min(Score_EP)
    hopt_LSCV<-h[ind_EP]
    densitiesLSCV<-as.vector(kdens(data,hopt_LSCV,xx,wg))
  }





  if (bw == "DPI1" ) {

    Psi6NREP<- -15/(16*pi^(1/2)*stand^7)
    K4zero<-3*(2*pi)^(-1/2)
    g11mEP<-((-2*numfuncEP*K4zero)/(Psi6NREP*length(X)))^(1/7)
    hatqsi4g21EP<-sum(sapply(dx1,HPhi2,g11mEP)*rr1)/(alfanfEP^(-2)*length(X)^2*g11mEP^5)
    hopt_DPI1<-(((2*pi^(1/2))^(-1)*numfuncEP)/(hatqsi4g21EP*length(X)))^(1/5)

    densitiesDPI1<-as.vector(kdens(data,hopt_DPI1,xx,wg))
  }





  if (bw == "SBoot") {

    B <-200
    h<-seq(0.001,max(X), by=(max(X)-min(X))/300)
    c<-(max(xx)-min(xx))/length(xx)
    vxx<-cbind(X,U,V)
    indOB<-1:nrow(vxx)
    Psi8NREP<- (105)/(32*pi^(1/2)*stand^9)
    K6zero<--15*(2*pi)^(-1/2)
    g1mEP<-((-2*numfuncEP*K6zero)/(Psi8NREP*length(X)))^(1/9)
    hatqsi6g1EP<-sum(sapply(dx1,HPhi,g1mEP)*rr1)/(alfanfEP^(-2)*length(X)^2*g1mEP^7)
    K4zero<-3*(2*pi)^(-1/2)
    g2mEP<-((-2*numfuncEP*K4zero)/(hatqsi6g1EP*length(X)))^(1/7)
    hatqsi4g2EP<-sum(sapply(dx1,HPhi2,g2mEP)*rr1)/(alfanfEP^(-2)*n^2*g2mEP^5)
    hopt_DPI2<-(((2*pi^(1/2))^(-1)*numfuncEP)/(hatqsi4g2EP*n))^(1/5)
    densitiesDPI2<-as.vector(kdens(data,hopt_DPI2,xx,wg))

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

      wt <- rep(1 / length(vxx[,1]), times = length(vxx[,1]))
      error <- 1e-6
      nmaxit<-10000000
       cnt_n_rep <- 0

      vxBOBEP<-matrix(0,nrow=nrow(vxx),ncol=3)
      fEPB<-matrix(data=0,nrow=length(X))
      Mat_ISETB_EP<-matrix(0,nrow=length(h))
      Vec_ISETB_EP<-matrix(0,ncol=1)


      repeat {

        vxBOBEP[,1]<-(max(X)^2-hopt_DPI2^2)^1/2*sample(vxx[,1],size=nrow(vxx),replace=T,prob=as.vector(vxM[[5]])) + (hopt_DPI2)*rnorm(1,0,sqrt(var(vxx[,1])))
        indUVEP<-sample(indOB,size=nrow(vxx),replace=T,prob=vxM[[11]])
        vxBOBEP[,2:3]<-vxx[indUVEP,2:3]

        log1EP<-(vxBOBEP[,1]>=vxBOBEP[,2])
        log2EP<-(vxBOBEP[,1]<=vxBOBEP[,3])
        indlogB<-(log1EP*log2EP==0) #this indicator vector gives the position that do not agree

        while(sum(indlogB)>0){
          vxBOBEP[indlogB,1]<-(max(X)^2-hopt_DPI2^2)^1/2*sample(vxx[,1],replace=T,prob=vxM[[5]],size=sum(indlogB))+ (hopt_DPI2)*rnorm(1,0,sqrt(var(vxx[,1])))
          indUVEP<-sample(indOB,size=sum(indlogB),replace=T,prob=vxM[[11]])
          vxBOBEP[indlogB,2:3]<-vxx[indUVEP,2:3]
          log1EP<-(vxBOBEP[,1]>=vxBOBEP[,2])
          log2EP<-(vxBOBEP[,1]<=vxBOBEP[,3])
          indlogB<-(log1EP*log2EP==0)
        }

        ordB<-order(vxBOBEP[,1])

        vxxBEP<-matrix(0,nrow=nrow(vxBOBEP),ncol=ncol(vxBOBEP))
        vxxBEP[,1]<-sort(vxBOBEP[,1])
        vxxBEP[,2:ncol(vxBOBEP)]<-vxBOBEP[ordB,2:ncol(vxBOBEP)]

        Jb<-matrix(data=0,ncol=nrow(vxxBEP),nrow=nrow(vxxBEP))

        aub <- outer(vxxBEP[,1], vxxBEP[,2], ">=")
        avb <- outer(vxxBEP[,1], vxxBEP[,3], "<=")
        auub <- outer(vxxBEP[,2], vxxBEP[,2], "<=")*1L
        Jb <- aub*avb
        SCb <- apply(Jb,2,"sum")
        adb <- length(which(SCb==1))
        cat("Number of changes in Bootstrap sample",adb,"\n") # test condition to repeat proceess
        cnt_n_rep <- cnt_n_rep +1   # add 1 to counter
        # if the condition is set, we proceeed, otherwise repeat process
        if (adb == 0)
          break
        else
          cat("Warning. Non-uniqueness or no existence of the NPMLE - New Round","\n")

      }
      # now we proceed after good sample
      JIb <- t(Jb)

      # preallocation
      f0b <- matrix(data = wt, ncol = 1, nrow = nrow(vxxBEP))
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

      F0b<-JIb%*%f1b
      K0<-Jb%*%k0b

      wgEPB<-f1b
      fEPB<-kdens(vxxBEP[,1],h,vec,wgEPB)
      fgEP<-densitiesDPI2
      Vec_ISETB_EP<-numeric(length(h))
      Vec_ISETB_EP<-sapply( 1: length(h), function(i) Vec_ISETB_EP[i]<-trap_rule(c,fgEP,fEPB[i,]))
      for(i in 1:dim(fEPB)[1]){
        Mat_ISETB_EP[i,]<-trap_rule(c,fgEP,fEPB[i,])
      }

      # variales to be outputed of the clusters (paralell)

      return(list(Mat_ISETB_EP,  cnt_n_rep ))

      ########################################
      # End of parallel state both
      ########################################
    }
    #stop cluster
    stopCluster(cl)

    # computation of the B aggegated results r comes as a list concatenated as row (thsi enables any size lists)
    # preallocate the original ones as matrix a,d put the list indexing

    Mat_ISETB_EPF<-matrix(0,nrow=length(h),ncol=B) # In each line, we will store Vec_ISE
    BootRepeat<-matrix(0,ncol=B)

    for(b in 1:B){
      Mat_ISETB_EPF[,b]<-as.numeric(unlist(final_boot[b,1]))
      BootRepeat[,b]<-as.numeric(unlist(final_boot[b,2]))

    }

    if(sum(BootRepeat>B)){
      cat("Warning. At least one bootstrap resample was replaced by a new one due to the violation of the existence and/or uniqueness condition for the NPMLE","\n")
    }

    vec_MISETB_EP<-apply(Mat_ISETB_EPF,1,"mean")
    indB_SBoot<-which.min(vec_MISETB_EP)


    hoptB_SBoot<-h[indB_SBoot]

    densitiesSBoot<-as.vector(kdens(data,hoptB_SBoot,xx,wg))

  }




  if(is.numeric(bw)){
    return(invisible(list(x=round(xx,5),y=round(densities,5),bw=bw)))
  }

  if (bw == "NR") {
    return(invisible(list(x=round(xx,5),y=round(densitiesNR,5),bw=hoptNR)))
  }

  if (bw == "LSCV") {
    return(invisible(list(x=round(xx,5),y=round(densitiesLSCV,5),bw=hopt_LSCV)))
  }

  if (bw == "DPI1") {
    return(invisible(list(x=round(xx,5),y=round(densitiesDPI1,5),bw=hopt_DPI1)))
  }


  if (bw == "DPI2") {
	return(invisible(list(x=round(xx,5),y=round(densitiesDPI2,5),bw=hopt_DPI2)))
	}

	if (bw == "SBoot") {
	return(invisible(list(x=round(xx,5),y=round(densitiesSBoot,5),bw=hoptB_SBoot)))
	}


}





