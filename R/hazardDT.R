


hazardDT<-function(X,U,V, bw="LSCV", from,to,n, wg=NA){

 shen0 <-function(X, U=NA, V=NA, wt=NA, error=NA, nmaxit=NA){

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
    if(length(T)!=1){
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
    if(is.na(nmaxit)==TRUE)nmaxit<-10000000
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




CVhazardEP=function(h,data, wg,bias,vxMMEP,alpnEP,biasfEP){
	t1<-sapply(1:length(h), function(i) {integrate(function(xx) (kdens(data, h[i],xx,wg)/sum(1/vxM[[9]]))^2,lower=min(xx), upper= max(xx), subdivisions=2000)$value})
	#t1<-integrate(function(vec) kdens(data, h,vec,wg)^2,lower=min(vec), upper= max(vec), abs.tol=1e-4)$value
	r<-sapply(1:length(data),function(i){kdens(data[-i], h,data[i],wg[-i])/sum(1/vxMMEP[,i][[9]])})
	d1cv<-matrix(data=0, ncol=1, nrow=length(data))
	d2cv<-matrix(data=0, ncol=1, nrow=length(data))

	d1cv[2:length(data)]<-sapply(2:length(data), function(i) d1cv[i]<-sum(1/vxMMEP[[9]][1:i-1]))
	d2cv<-sapply(1:length(data),function(i) d2cv[i]<-sum(1/vxMMEP[,i][[9]]))
	sob0Mcv<-1-d1cv/d2cv
	sob0Mcv[sob0Mcv==0]<-0.0000001
	t2<-sapply(1: nrow(r), function (i) t2<-2*alpnEP*sum(r[i,]*((biasfEP)^(-1)/(sob0Mcv)))/(length(data)-1))

	return(t1-t2)
}




ord<-order(X)
X<-sort(X)
U<-U[ord]
V<-V[ord]
vxx<-cbind(X,U,V)

vxM<-shen0(X,U,V)
bias<-vxM[[9]]
f3hatM<-1/vxM[[9]]
	d1<-matrix(data=0, ncol=1, nrow=length(X))
	d1[2:length(X)]<-sapply(2:length(X), function(i) d1[i]<-sum(1/vxM[[9]][1:i-1]))
	d2<-sum(f3hatM)
	sob0M<-1-d1/d2


if(is.na(wg) == TRUE){


	wg<-wgEP<-f3hatM/sob0M
}

Mdata<-X
data<-X

xx<-seq(from+0.0001,to,length=n)
#xx<-seq(min(x)+0.0001,max(x),length=500)
vec<-xx
h<-seq(0.001,max(X), by=(max(X)-min(X))/500)
####################
  #bandwidth selectors
  ####################


  if(is.numeric(bw)){
    hazard<-as.vector(kdens(data,bw,xx,wg)/sum(1/vxM[[9]]))
  }

else {

vxMMEP<-list()
vxMMEP<-sapply(1:length(X), function(i) vxMMEP<-shen0(vxx[-i,1],vxx[-i,2],vxx[-i,3]))
biasfSEP<-list()

hatpsiEP<-list()
#hatpsiEP<-sapply(1:length(X), function(i) hatpsiEP<-vxMMEP[,i][[11]])
hatpsiEP<-sapply(1:length(X), function(i) hatpsiEP<-vxMMEP[,i][[9]])

aux1EP<-list()
aux1EP<-sapply(1:length(X), function(i) aux1EP<-( vxx[i,1]>=vxx[-i,2])& ( vxx[i,1]<=vxx[-i,3]))*1L

weightsEP<-aux1EP*hatpsiEP
biasfEP<- apply(weightsEP,2,"sum")
alpnEP<-((1/(length(X)-1))*sum(1/biasfEP))^(-1)
####fhaz1<-vxMMEP[,i][[9]]

Score_EP<-CVhazardEP(h, data, wgEP,bias, vxMMEP,alpnEP,biasfEP)
ind_EP<-which.min(Score_EP)

hopt_LSCV<-h[ind_EP]

hazard_LSCV<-as.vector(kdens(data,hopt_LSCV,xx,wgEP)/sum(1/vxM[[9]]))



}



if(is.numeric(bw)){
    return(invisible(list(x=round(xx,5),y=round(hazard,5),bw=bw)))
  }

else{

return(invisible(list(x=round(xx,5),y=round(hazard_LSCV,5),bw=hopt_LSCV)))
  }
}





