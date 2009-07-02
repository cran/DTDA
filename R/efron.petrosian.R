`efron.petrosian` <-
function(X, U=NA, V=NA, wt=NA, error=NA,
 nmaxit=NA , boot=TRUE, B=NA, alpha=NA, 
display.F=FALSE, display.S=FALSE){

trunc<-"double"

if(all(is.na(U))==TRUE){

trunc<-"right"
cat("case U=NA","\n")
if(any(is.na(V))==TRUE|any(is.na(X))==TRUE){
navec<-c(which(is.na(X)),which(is.na(V)))
X<-X[-navec]
V<-V[-navec]
}}
if(all(is.na(V))==TRUE){

trunc<-"left"
cat("case V=NA","\n")
if(any(is.na(U))==TRUE|any(is.na(X))==TRUE){
navec<-c(which(is.na(X)),which(is.na(U)))
X<-X[-navec]
U<-U[-navec]
}}

if(trunc=="double"){
if(any(is.na(U))==TRUE | any(is.na(V))==TRUE|any(is.na(X))==TRUE){

navec<-c(which(is.na(X)),which(is.na(U)),which(is.na(V)))
X<-X[-navec]
U<-U[-navec]
V<-V[-navec]

}}

if(is.na(wt)==TRUE) wt<-rep(1/length(X),times=length(X))
D<-cbind(X,U,V)

if(all(is.na(U))==TRUE){
D[,2]<--D[,3]
D[,1]<--D[,1]
D[,3]<-rep(max(D[,1])+1,length(X))
}
if(all(is.na(V))==TRUE){
D[,3]<-rep(max(D[,1])+1,length(X))
}

ord<-order(D[,1])
C<-matrix(0,nrow=nrow(D),ncol=ncol(D))
C[,1]<-sort(D[,1])
C[,2:ncol(D)]<-D[ord,2:ncol(D)]


if(is.na(error)==TRUE) error<-1/(10*nrow(C))
J<-matrix(data=0,ncol=nrow(C),nrow=nrow(C))
for (k in 1:nrow(C)){
for(j in 1:nrow(C)) {
a1<-min(C[j,2],C[j,3])
     b1<-max(C[j,2],C[j,3])
if(C[k,1]>=a1 & C[k,1]<=b1) J[k,j]<-1
}}

JI<-t(J)
f0<-matrix(data=wt,ncol=1,nrow=nrow(C))
f<-f0
S0<-1
if(is.na(nmaxit)==TRUE)nmaxit<-100
iter<-0
while(S0>error | iter>nmaxit){
iter<-iter+1
cat("iter",iter,"\n")
F0<-JI%*%f
IF0<-1/F0
If1<-J%*%IF0
f<-1/If1
if(sum(f)!=1)f<-f/sum(f)

S0<-max(abs(f-f0))
f0<-f
}

FF<-matrix(data=0, ncol=1, nrow=nrow(C))
for(i in 1:nrow(C)){
FF[i,]<-sum(f[1:i,])
 }
Sob<-matrix(data=0, ncol=1, nrow=nrow(C))
for(j in 1:nrow(C)){
Sob[j,]<-1-FF[j,]+f[j,]
 }

if (boot==TRUE){
if (is.na(B)==TRUE) B<-500

M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
ind<-seq(1,nrow(C),by=1)

for (b in 1:B){
cat("Resample B",b,"\n")
indb<-sample(ind,nrow(C),replace=TRUE)
M1b<-C[indb,]
M2b<-matrix(0,nrow=nrow(M1b),ncol=ncol(M1b))
ord<-order(M1b[,1])
M2b[,1]<-sort(M1b[,1])
M2b[,2:ncol(M1b)]<-M1b[ord,2:ncol(M1b)]

Jb<-matrix(data=0,ncol=nrow(M2b),nrow=nrow(M2b))
for (k in 1:nrow(M2b)){
for(j in 1:nrow(M2b)) {
a2<-min(M2b[j,2],M2b[j,3])
     b2<-max(M2b[j,2],M2b[j,3])
if(M2b[k,1]>=a2 & M2b[k,1]<=b2) Jb[k,j]<-1
}}
JIb<-t(Jb)
f0b<-matrix(data=wt,ncol=1,nrow=nrow(M2b))
f1b<-f0b
 S0b<-1
iterb<-0
while(S0b>error|iterb>nmaxit){
iterb<-iterb+1
F0b<-JIb%*%f1b
IF0b<-1/F0b
If1b<-Jb%*%IF0b
f1b<-1/If1b
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
}
FF0b<-matrix(data=0, ncol=1, nrow=nrow(M2b))
for(i in 1:nrow(M2b)){
FF0b[i,]<-sum(f1b[1:i,])
 }
Sobb<-matrix(data=0, ncol=1, nrow=nrow(M2b))
for(j in 1:nrow(M2b)){
Sobb[j,]<-1-FF0b[j,]+f1b[j,]
}


M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(f1b)
M_IF0Sob[,b]<-as.vector(Sobb)

}

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
}

if(boot==TRUE){

if(trunc=="double"|trunc=="left"){

if((display.F==TRUE)&(display.S==TRUE)){

dev.new()
par(mfrow=c(1,2))
plot(C[,1],FF,main="EP estimator",xlab="Time of interest",ylab="",lty=1,type="l")
 lines(C[,1],upperF,lty=3)
 lines(C[,1],lowerF,lty=3)
plot(C[,1],Sob,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 lines(C[,1],upperS,lty=3)
 lines(C[,1],lowerS,lty=3)
}

if((display.S==FALSE)&(display.F==TRUE)){
dev.new()
plot(C[,1],FF,main="EP estimator",xlab="Time of interest",ylab="",lty=1,type="l")
 lines(C[,1],upperF,lty=3 )
 lines(C[,1],lowerF,lty=3)
}

if((display.F==FALSE)&(display.S==TRUE)){
dev.new()
plot(C[,1],Sob,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 lines(C[,1],upperS,lty=3)
 lines(C[,1],lowerS,lty=3)
}}

if(trunc=="right"){

if((display.F==TRUE)&(display.S==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(-C[,1],FF,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 lines(-C[,1],upperF,lty=3)
 lines(-C[,1],lowerF,lty=3)
}

if((display.F==FALSE)&(display.S==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(-C[,1],FF,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 lines(-C[,1],upperF,lty=3)
 lines(-C[,1],lowerF,lty=3)
}

if((display.F==TRUE)&(display.S==FALSE)){
dev.new()
par(mfrow=c(1,1))
plot(-C[,1],FF,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 lines(-C[,1],upperF,lty=3)
 lines(-C[,1],lowerF,lty=3)
}

}
}
if (boot==FALSE){
upperF<-rep(NA,length(C))
lowerF<-rep(NA,length(C))
upperS<-rep(NA,length(C))
lowerS<-rep(NA,length(C))


if(trunc=="double"|trunc=="left"){

if((display.F==TRUE)&(display.S==TRUE)){

dev.new()
par(mfrow=c(1,2))
plot(C[,1],FF,main="EP estimator",xlab="Time of interest",ylab="",lty=1,type="l")
 plot(C[,1],Sob,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 }

if((display.S==FALSE)&(display.F==TRUE)){
dev.new()
plot(C[,1],FF,main="EP estimator",xlab="Time of interest",ylab="",lty=1,type="l")
 }

if((display.F==FALSE)&(display.S==TRUE)){
dev.new()
plot(C[,1],Sob,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
}}

if(trunc=="right"){

if((display.F==TRUE)&(display.S==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(-C[,1],FF,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 }

if((display.F==FALSE)&(display.S==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(-C[,1],FF,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 }

if((display.F==TRUE)&(display.S==FALSE)){
dev.new()
par(mfrow=c(1,1))
plot(-C[,1],FF,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 }

}

}
return(list(density=f,cumulative.df=FF,truncation.probs=F0,S0=S0,survival=Sob,n.iterations=iter,
Boot="simple",B=B,alpha=alpha,upper.df=upperF,lower.df=lowerF,upper.Sob=upperS,lower.Sob=lowerS))

}

