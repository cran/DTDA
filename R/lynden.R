`lynden` <-
function(X, U=NA, V=NA, error=NA, nmaxit=NA,
 boot=TRUE, B=NA, alpha=NA, display.F=FALSE, 
display.S=FALSE){



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


NJ<-matrix(data=0,ncol=1,nrow=nrow(C))
K<-NJ
for(j in 1:nrow(C)){
for(k in 1:nrow(C)){
if (C[k,2]<=C[j,1]& C[k,1]>=C[j,1]) {NJ[j,1]<-NJ[j,1]+1}
}}
K[,1]<-rep(1,length=nrow(C))


mNJ<-min(NJ[1:nrow(C)-1,1])
if (mNJ!=1) K[,1]<-NJ[,1]
if (mNJ==1){
istar<-max(which(NJ[1:nrow(C)-1,1]==1))
K[1:istar,1]<-2
if(istar<nrow(C)){K[(istar+1):nrow(C),1]<-NJ[(istar+1):nrow(C),1]}
}
NJ<-as.matrix(apply(cbind(NJ,K),1,max),ncol=1)


J<-matrix(data=0,ncol=nrow(C),nrow=nrow(C))
for (k in 1:nrow(C)){
for(j in 1:nrow(C)) {
a1<-min(C[j,2],C[j,3])
b1<-max(C[j,2],C[j,3])
if(C[k,1]>=a1 & C[k,1]<=b1) J[k,j]<-1
}}
h0<-1/NJ


G0<-matrix(data=0,ncol=1,nrow=nrow(C))
G0[1,1]<-1
for(j in 2:nrow(C)){
G0[j,1]<-exp(sum(log(1-h0[1:(j-1),1])))
}


f0<-matrix(data=G0[nrow(C),1],ncol=1,nrow=nrow(C))
for(j in 1:nrow(C)-1){
f0[j,1]<-(G0[j,1]-G0[j+1,1])
}

Gvi0<-matrix(data=0,ncol=1, nrow=nrow(C))
for(j in 1:nrow(C)){
auxi<-0
for (k in 1:nrow(C)){
if (C[j,3]< C[k,1]) {auxi<-sum(f0[k,])+ auxi}
}
Gvi0[j,]<-auxi
}

F0<-t(J)%*%f0
Q0<-Gvi0/F0
h<-h0
for(i in 1:nrow(h0)){
h[i,]<- 1/(NJ[i,]+(J[i,])%*%Q0)
}

S0<-1
if(is.na(error)==TRUE) error<-1/(10*nrow(C))

if(is.na(nmaxit)==TRUE)nmaxit<-100
iter<-0

while(S0>error|iter>nmaxit){
iter<-iter+1
cat("iter",iter,"\n")
G1<-matrix(data=0,ncol=1,nrow=nrow(C))
G1[1,1]<-1
for(j in 2:nrow(C)){
G1[j,1]<-exp(sum(log(1-h[1:(j-1),1])))
}

f<-matrix(data=G1[nrow(C),1],ncol=1,nrow=nrow(C))
for(j in 1:nrow(C)-1){
f[j,1]<-(G1[j,1]-G1[j+1,1])
f[nrow(C),1]<-1-sum(f[1:nrow(C)-1,1])
}


Gvi1<-matrix(data=0,ncol=1, nrow=nrow(C))
for(j in 1:nrow(C)){
auxi<-0
for (k in 1:nrow(C)){
if (C[j,3]< C[k,1]) {auxi<-sum(f[k,])+ auxi}
}
Gvi1[j,]<-auxi
}



F0<-t(J)%*%f
Q0<-Gvi0/F0
for(i in 1:nrow(h)){
h[i,]<- 1/(NJ[i,]+(J[i,])%*%Q0)
}

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

NJb<-matrix(data=0,ncol=1,nrow=nrow(C))
Kb<-NJb
for(j in 1:nrow(C)){
for(k in 1:nrow(C)){
if (M2b[k,2]<=M2b[j,1]& M2b[k,1]>=M2b[j,1]) {NJb[j,1]<-NJb[j,1]+1}
}}
Kb[,1]<-rep(1,length=nrow(C))


mNJb<-min(NJb[1:nrow(C)-1,1])
if (mNJb!=1) {Kb[,1]<-NJb[,1]}
if (mNJb==1){
istarb<-max(which(NJb[1:nrow(C)-1,1]==1))
Kb[1:istarb,1]<-2
if(istarb<nrow(M2b)){Kb[(istarb+1):nrow(M2b),1]<-NJb[(istarb+1):nrow(M2b),1]}
}
NJb<-as.matrix(apply(cbind(NJb,Kb),1,max),ncol=1)


Jb<-matrix(data=0,ncol=nrow(M2b),nrow=nrow(M2b))
for (k in 1:nrow(M2b)){
for(j in 1:nrow(M2b)) {
a2<-min(M2b[j,2],M2b[j,3])
     b2<-max(M2b[j,2],M2b[j,3])
if(M2b[k,1]>=a2 & M2b[k,1]<=b2) Jb[k,j]<-1
}}
h0b<-1/NJb

G0b<-matrix(data=0,ncol=1,nrow=nrow(C))
G0b[1,1]<-1
for(j in 2:nrow(C)){
G0b[j,1]<-exp(sum(log(1-h0b[1:(j-1),1])))
}
f0b<-matrix(data=G0b[nrow(C),1],ncol=1,nrow=nrow(C))
for(j in 1:nrow(C)-1){
f0[j,1]<-(G0b[j,1]-G0b[j+1,1])
}
Gvi0b<-matrix(data=0,ncol=1, nrow=nrow(C))
for(j in 1:nrow(C)){
auxib<-0
for (k in 1:nrow(C)){
if (M2b[j,3]< M2b[k,1]) {auxib<-sum(f0b[k,])+ auxib}
}
Gvi0b[j,]<-auxib
}
F0b<-t(Jb)%*%f0b
Q0b<-Gvi0b/F0b
h1b<-h0b
for(i in 1:nrow(h0b)){
h1b[i,]<- 1/(NJb[i,]+(Jb[i,])%*%Q0b)
}

 S0b<-1
if(is.na(error)==TRUE) error<-1/(10*nrow(C))

if(is.na(nmaxit)==TRUE)nmaxit<-100

iterb<-0
while(S0b>error|iterb>nmaxit){
iterb<-iterb+1
G1b<-matrix(data=0,ncol=1,nrow=nrow(C))
G1b[1,1]<-1
for(j in 2:nrow(C)){
G1b[j,1]<-exp(sum(log(1-h1b[1:(j-1),1])))
}
f1b<-matrix(data=G1b[nrow(C),1],ncol=1,nrow=nrow(C))
for(j in 1:nrow(C)-1){
f1b[j,1]<-(G1b[j,1]-G1b[j+1,1])
f1b[nrow(C),1]<-1-sum(f1b[1:nrow(C)-1,1])
}
Gvi1b<-matrix(data=0,ncol=1, nrow=nrow(C))
for(j in 1:nrow(C)){
auxib<-0
for (k in 1:nrow(C)){
if (M2b[j,3]< M2b[k,1]) {auxib<-sum(f1b[k,])+ auxib}
}
Gvi1b[j,]<-auxib
}
F0b<-t(Jb)%*%f1b
Q0b<-Gvi0b/F0b
for(i in 1:nrow(h1b)){
h1b[i,]<- 1/(NJb[i,]+(Jb[i,])%*%Q0b)
}
S0b<-max(abs(f1b-f0b))
f0b<-f1b


}


FF0b<-cumsum(as.vector(f1b))
Sobb<-numeric(nrow(C))
Sobb<-1-FF0b+as.vector(f1b)


M_IF0[,b]<-FF0b
M_IF01[,b]<-as.vector(f1b)
M_IF0Sob[,b]<-Sobb
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



if (boot==TRUE){

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

if(display.S==FALSE){
dev.new()
plot(C[,1],FF,main="EP estimator",xlab="Time of interest",ylab="",lty=1,type="l")
 lines(C[,1],upperF,lty=3 )
 lines(C[,1],lowerF,lty=3)
}

if(display.F==FALSE){
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

if(display.S==FALSE){
dev.new()
plot(C[,1],FF,main="EP estimator",xlab="Time of interest",ylab="",lty=1,type="l")
 }

if(display.F==FALSE){
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

}

}
return(list(NJ=NJ,density=f,cumulative.df=FF,truncation.probs=F0,hazard=h,S0=S0,
survival=Sob,n.iterations=iter,Boot="simple",B=B,alpha=alpha,upper.df=upperF,
lower.df=lowerF,upper.Sob=upperS,lower.Sob=lowerS))

}

