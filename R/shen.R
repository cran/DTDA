`shen` <-
function(X, U=NA, V=NA, wt=NA, error=NA, nmaxit=NA,
 boot=TRUE, boot.type="simple",B=NA, alpha=NA, 
display.FS=FALSE, display.UV=FALSE,
 plot.joint=FALSE, plot.type=NULL){

trunc<-"both"
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

if(trunc=="both"){
if(any(is.na(U))==TRUE|any(is.na(V))==TRUE|any(is.na(X))==TRUE){
navec<-c(which(is.na(X)),which(is.na(U)),which(is.na(V)))
X<-X[-navec]
U<-U[-navec]
V<-V[-navec]
}}



if(is.na(wt)==TRUE) wt<-rep(1/length(X),times=length(X))

D<-cbind(X,U,V)
if(all(is.na(V))==TRUE)D[,3]<-rep(Inf,times=length(X))
if(all(is.na(U))==TRUE){
D[,2]<--D[,3]
D[,1]<--D[,1]
D[,3]<-rep(Inf,length(X))
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
f0<-matrix(data=1/nrow(C),ncol=1,nrow=nrow(C))
f<-f0
S0<-1
if(is.na(nmaxit)==TRUE)nmaxit<-100
iter<-0
while(S0>error | iter>nmaxit){
iter<-iter+1
cat("iter",iter,"\n")
F0<-JI%*%f
k0<-((sum(1/F0))^(-1))*(1/F0)
if(sum(k0)!=1)k0<-k0/sum(k0)
k<-k0
K0<-J%*%k
f<-((sum(1/K0))^(-1))*(1/K0)
if(sum(f)!=1)f<-f/sum(f)
S0<-max(abs(f-f0))
f0<-f
k0<-k
}
FF<-matrix(data=0, ncol=1, nrow=nrow(C))
for(i in 1:nrow(C)){
FF[i,]<-sum(f[1:i,])
}

Sob<-matrix(data=0, ncol=1, nrow=nrow(D))
for(j in 1:nrow(D)){
Sob[j,]<-1-FF[j,]+f[j,]
}




if(trunc=="both"){

kMUV<-cbind(U,V,k)
kMU<-cbind(U,k)
ordU<-order(kMU[,1])
kMU[,1]<-sort(kMU[,1])
kMU[,2]<-kMU[ordU,2]
kMV<-cbind(V,k)
ordV<-order(kMV[,1])
kMV[,1]<-sort(kMV[,1])
kMV[,2]<-kMV[ordV,2]


fU<-matrix(data=0, nrow=nrow(C),ncol=1)
fV<-matrix(data=0, nrow=nrow(C),ncol=1)

for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(kMU[i,1]<=kMU[j,1])
                        fU[j,]<-sum(kMU[1:i,2])
                        }}

for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(kMV[i,1]<=kMV[j,1])
fV[j,]<-sum(kMV[1:i,2])
}}


KK<-matrix(data=0, ncol=nrow(C), nrow=nrow(C))
cat("Computing joint distribution","\n")
for (i in 1:nrow(C)){
for (j in 1:nrow(C)){
for (l in 1:nrow(C)){
if(kMU[i,1]>=kMUV[l,1]& kMV[j,1]>=kMUV[l,2])
KK[i,j]<-KK[i,j]+kMUV[l,3]
}}}


if(plot.joint==TRUE){
if(plot.type=="image"){
image(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution")
filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
}
if(plot.type=="persp"){
fcol<-topo.colors(10)[cut(KK[2:nrow(C),2:nrow(C)],10,include.lowest=TRUE)]
                   persp(x=sort(U),y=sort(V),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")

}

}
if(plot.joint==FALSE){
}

}

if(trunc=="left"){

kMU<-cbind(U,k)
ordU<-order(kMU[,1])
kMU[,1]<-sort(kMU[,1])
kMU[,2]<-kMU[ordU,2]


fU<-matrix(data=0, nrow=nrow(C),ncol=1)

for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(kMU[i,1]<=kMU[j,1])
                        fU[j,]<-sum(kMU[1:i,2])
                        }}
fV<-rep(NA,length(fU))
KK<-NA
}

if(trunc=="right"){


kMV<-cbind(V,k)
ordV<-order(kMV[,1])
kMV[,1]<-sort(kMV[,1])
kMV[,2]<-kMV[ordV,2]


fV<-matrix(data=0, nrow=nrow(C),ncol=1)


for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(kMV[i,1]<=kMV[j,1])
fV[j,]<-sum(kMV[1:i,2])
}}

fU<-rep(NA,length(fV))
KK<-NA

}


if(boot==TRUE){

if(boot.type=="simple"){
cat("Simple Bootstrap","\n")
if (is.na(B)==TRUE) B<-500

if(trunc=="both"){

M1b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M2b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
kMUb<-matrix(0,nrow=nrow(C),ncol=B)
kMVb<-matrix(0,nrow=nrow(C),ncol=B)
M_fU<-matrix(0,nrow=nrow(C),ncol=B) 
M_fV<-matrix(0,nrow=nrow(C),ncol=B) 
M3b<-matrix(0,nrow=nrow(C),ncol=B)
M4b<-matrix(0,nrow=nrow(C),ncol=B)
k1bU<-matrix(0,nrow=nrow(C),ncol=B) 
k1bV<-matrix(0,nrow=nrow(C),ncol=B)


ind<-seq(1,nrow(C),by=1)

for (b in 1:B){
cat("Resample B",b,"\n")
indb<-sample(ind,nrow(C),replace=TRUE)
M1b<-C[indb,]
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
while(S0b>error){
iterb<-iterb+1
F0b<-JIb%*%f1b
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
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



ordA<-order(M2b[,2])
ordB<-order(M2b[,3])
ord3<-sort(M2b[,2])
ord4<-sort(M2b[,3])
M3b[,b]<-ord3
M4b[,b]<-ord4
k1bU[,b]<-k1b[ordU,]
k1bV[,b]<-k1b[ordV,]




for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M3b[i,b]<=kMU[j,1])
                        M_fU[j,b]<-sum(k1bU[1:i,b])
                        }}

for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M4b[i,b]<=kMV[j,1])
M_fV[j,b]<-sum(k1bV[1:i,b])
}}



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


M_fU_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_fU_sort)){
M_fU_sort[i,]<-sort(M_fU[i,])
}
lowerU<-M_fU_sort[,floor(alpha*B/2)]
upperU<-M_fU_sort[,floor((1-alpha/2)*B)]

M_fV_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_fV_sort)){
M_fV_sort[i,]<-sort(M_fV[i,])
}
lowerV<-M_fV_sort[,floor(alpha*B/2)]
upperV<-M_fV_sort[,floor((1-alpha/2)*B)]

}


if(trunc=="left"){

M1b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M2b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
kMUb<-matrix(0,nrow=nrow(C),ncol=B)

M_fU<-matrix(0,nrow=nrow(C),ncol=B) 
M3b<-matrix(0,nrow=nrow(C),ncol=B)
k1bU<-matrix(0,nrow=nrow(C),ncol=B) 


ind<-seq(1,nrow(C),by=1)

for (b in 1:B){
cat("Resample B",b,"\n")
indb<-sample(ind,nrow(C),replace=TRUE)
M1b<-C[indb,]
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
while(S0b>error){
iterb<-iterb+1
F0b<-JIb%*%f1b
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
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



ordA<-order(M2b[,2])
ord3<-sort(M2b[,2])
M3b[,b]<-ord3
k1bU[,b]<-k1b[ordU,]




for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M3b[i,b]<=kMU[j,1])
                        M_fU[j,b]<-sum(k1bU[1:i,b])
                        }}



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


M_fU_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_fU_sort)){
M_fU_sort[i,]<-sort(M_fU[i,])
}
lowerU<-M_fU_sort[,floor(alpha*B/2)]
upperU<-M_fU_sort[,floor((1-alpha/2)*B)]

upperV<-rep(NA,length(C))
lowerV<-rep(NA,length(C))
}

if(trunc=="right"){

M1b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M2b<-matrix(0,nrow=nrow(C),ncol=ncol(C))
M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
kMVb<-matrix(0,nrow=nrow(C),ncol=B)
M_fV<-matrix(0,nrow=nrow(C),ncol=B) 
M4b<-matrix(0,nrow=nrow(C),ncol=B)
k1bV<-matrix(0,nrow=nrow(C),ncol=B)


ind<-seq(1,nrow(C),by=1)

for (b in 1:B){
cat("Resample B",b,"\n")
indb<-sample(ind,nrow(C),replace=TRUE)
M1b<-C[indb,]
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
while(S0b>error){
iterb<-iterb+1
F0b<-JIb%*%f1b
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
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



ordB<-order(M2b[,3])
ord4<-sort(M2b[,3])

M4b[,b]<-ord4
k1bV[,b]<-k1b[ordV,]





for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M4b[i,b]<=kMV[j,1])
M_fV[j,b]<-sum(k1bV[1:i,b])
}}



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


M_fV_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_fV_sort)){
M_fV_sort[i,]<-sort(M_fV[i,])
}
lowerV<-M_fV_sort[,floor(alpha*B/2)]
upperV<-M_fV_sort[,floor((1-alpha/2)*B)]

upperU<-rep(NA,length(C))
lowerU<-rep(NA,length(C))



}

}




if(boot.type=="obvious"){
cat("Obvious Bootstrap","\n")

if (is.na(B)==TRUE) B<-500

if(trunc=="both"){


W1<-as.vector(f)
W2<-as.vector(k)
ind<-1:nrow(C)
if (is.na(B)==TRUE) B<-500
M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
kMUb<-matrix(0,nrow=nrow(C),ncol=B)
kMVb<-matrix(0,nrow=nrow(C),ncol=B)
M_fU<-matrix(0,nrow=nrow(C),ncol=B) 
M_fV<-matrix(0,nrow=nrow(C),ncol=B) 
M3b<-matrix(0,nrow=nrow(C),ncol=B)
M4b<-matrix(0,nrow=nrow(C),ncol=B)
k1bU<-matrix(0,nrow=nrow(C),ncol=B) 
k1bV<-matrix(0,nrow=nrow(C),ncol=B)


for(b in 1:B){
cat("Resample B",b,"\n")
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

ordB<-order(DB[,1])
DBB<-matrix(0,nrow=nrow(C),ncol=3)
DBB[,1]<-sort(DB[,1])
DBB[,2:ncol(DB)]<-DB[ordB,2:ncol(DBB)]

error<-1/(10*nrow(DBB))
Jb<-matrix(data=0,ncol=nrow(DBB),nrow=nrow(DBB))
for (k in 1:nrow(DBB)){
for(j in 1:nrow(DBB)) {
a2<-min(DBB[j,2],DBB[j,3])
      b2<-max(DBB[j,2],DBB[j,3])
if(DBB[k,1]>=a2 & DBB[k,1]<=b2) Jb[k,j]<-1
}}
JIb<-t(Jb)
f0b<-matrix(data=1/nrow(DBB),ncol=1,nrow=nrow(DBB))
f1b<-f0b
S0b<-1
if(is.na(nmaxit)==TRUE)nmaxit<-100
iterb<-0
while(S0b>error){
iterb<-iterb+1
F0b<-JIb%*%f1b
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
}
FF0b<-matrix(data=0, ncol=1, nrow=nrow(DBB))
for(i in 1:nrow(DBB)){
FF0b[i,]<-sum(f1b[1:i,])
}
Sobb<-matrix(data=0, ncol=1, nrow=nrow(DBB))
for(j in 1:nrow(DBB)){
Sobb[j,]<-1-FF0b[j,]+f1b[j,]
}

M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(f1b)
M_IF0Sob[,b]<-as.vector(Sobb)



ordA<-order(DBB[,2])
ordB<-order(DBB[,3])
ord3<-sort(DBB[,2])
ord4<-sort(DBB[,3])
M3b[,b]<-ord3
M4b[,b]<-ord4
k1bU[,b]<-k1b[ordU,]
k1bV[,b]<-k1b[ordV,]




for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M3b[i,b]<=kMU[j,1])
                        M_fU[j,b]<-sum(k1bU[1:i,b])
                        }}

for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M4b[i,b]<=kMV[j,1])
M_fV[j,b]<-sum(k1bV[1:i,b])
}}



}

M_IF0_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_IF0_sort)){
M_IF0_sort[i,]<-sort(M_IF0[i,])
}
if(is.na(alpha)==TRUE) alpha<-0.05
lowerF<-M_IF0_sort[,floor(alpha*B/2)]
upperF<-M_IF0_sort[,floor((1-alpha/2)*B)]

M_IF0_sort1<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_IF0_sort1)){
M_IF0_sort1[i,]<-sort(M_IF0Sob[i,])
}
lowerS<-M_IF0_sort1[,floor(alpha*B/2)]
upperS<-M_IF0_sort1[,floor((1-alpha/2)*B)]


M_fU_sort<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_fU_sort)){
M_fU_sort[i,]<-sort(M_fU[i,])
}
lowerU<-M_fU_sort[,floor(alpha*B/2)]
upperU<-M_fU_sort[,floor((1-alpha/2)*B)]


M_fV_sort<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_fV_sort)){
M_fV_sort[i,]<-sort(M_fV[i,])
}
lowerV<-M_fV_sort[,floor(alpha*B/2)]
upperV<-M_fV_sort[,floor((1-alpha/2)*B)]



}


if(trunc=="left"){


W1<-as.vector(f)
W2<-as.vector(k)
ind<-1:nrow(C)
if (is.na(B)==TRUE) B<-500
M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)
kMUb<-matrix(0,nrow=nrow(C),ncol=B)
M_fU<-matrix(0,nrow=nrow(C),ncol=B) 

M3b<-matrix(0,nrow=nrow(C),ncol=B)
k1bU<-matrix(0,nrow=nrow(C),ncol=B) 


for(b in 1:B){
cat("Resample B",b,"\n")
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

ordB<-order(DB[,1])
DBB<-matrix(0,nrow=nrow(C),ncol=3)
DBB[,1]<-sort(DB[,1])
DBB[,2:ncol(DB)]<-DB[ordB,2:ncol(DBB)]

error<-1/(10*nrow(DBB))
Jb<-matrix(data=0,ncol=nrow(DBB),nrow=nrow(DBB))
for (k in 1:nrow(DBB)){
for(j in 1:nrow(DBB)) {
a2<-min(DBB[j,2],DBB[j,3])
      b2<-max(DBB[j,2],DBB[j,3])
if(DBB[k,1]>=a2 & DBB[k,1]<=b2) Jb[k,j]<-1
}}
JIb<-t(Jb)
f0b<-matrix(data=1/nrow(DBB),ncol=1,nrow=nrow(DBB))
f1b<-f0b
S0b<-1
if(is.na(nmaxit)==TRUE)nmaxit<-100
iterb<-0
while(S0b>error){
iterb<-iterb+1
F0b<-JIb%*%f1b
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
}
FF0b<-matrix(data=0, ncol=1, nrow=nrow(DBB))
for(i in 1:nrow(DBB)){
FF0b[i,]<-sum(f1b[1:i,])
}
Sobb<-matrix(data=0, ncol=1, nrow=nrow(DBB))
for(j in 1:nrow(DBB)){
Sobb[j,]<-1-FF0b[j,]+f1b[j,]
}

M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(f1b)
M_IF0Sob[,b]<-as.vector(Sobb)



ordA<-order(DBB[,2])
ord3<-sort(DBB[,2])
M3b[,b]<-ord3
k1bU[,b]<-k1b[ordU,]




for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M3b[i,b]<=kMU[j,1])
                        M_fU[j,b]<-sum(k1bU[1:i,b])
                        }}




}

M_IF0_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_IF0_sort)){
M_IF0_sort[i,]<-sort(M_IF0[i,])
}
if(is.na(alpha)==TRUE) alpha<-0.05
lowerF<-M_IF0_sort[,floor(alpha*B/2)]
upperF<-M_IF0_sort[,floor((1-alpha/2)*B)]

M_IF0_sort1<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_IF0_sort1)){
M_IF0_sort1[i,]<-sort(M_IF0Sob[i,])
}
lowerS<-M_IF0_sort1[,floor(alpha*B/2)]
upperS<-M_IF0_sort1[,floor((1-alpha/2)*B)]


M_fU_sort<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_fU_sort)){
M_fU_sort[i,]<-sort(M_fU[i,])
}
lowerU<-M_fU_sort[,floor(alpha*B/2)]
upperU<-M_fU_sort[,floor((1-alpha/2)*B)]

upperV<-rep(NA,length(C))
lowerV<-rep(NA,length(C))





}



if(trunc=="right"){


W1<-as.vector(f)
W2<-as.vector(k)
ind<-1:nrow(C)
if (is.na(B)==TRUE) B<-500
M_IF0<-matrix(0,nrow=nrow(C),ncol=B)
M_IF01<-matrix(0,nrow=nrow(C),ncol=B)
M_IF0Sob<-matrix(0,nrow=nrow(C),ncol=B)

kMVb<-matrix(0,nrow=nrow(C),ncol=B)
 
M_fV<-matrix(0,nrow=nrow(C),ncol=B) 
M4b<-matrix(0,nrow=nrow(C),ncol=B) 
k1bV<-matrix(0,nrow=nrow(C),ncol=B)


for(b in 1:B){
cat("Resample B",b,"\n")
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

ordB<-order(DB[,1])
DBB<-matrix(0,nrow=nrow(C),ncol=3)
DBB[,1]<-sort(DB[,1])
DBB[,2:ncol(DB)]<-DB[ordB,2:ncol(DBB)]

error<-1/(10*nrow(DBB))
Jb<-matrix(data=0,ncol=nrow(DBB),nrow=nrow(DBB))
for (k in 1:nrow(DBB)){
for(j in 1:nrow(DBB)) {
a2<-min(DBB[j,2],DBB[j,3])
      b2<-max(DBB[j,2],DBB[j,3])
if(DBB[k,1]>=a2 & DBB[k,1]<=b2) Jb[k,j]<-1
}}
JIb<-t(Jb)
f0b<-matrix(data=1/nrow(DBB),ncol=1,nrow=nrow(DBB))
f1b<-f0b
S0b<-1
if(is.na(nmaxit)==TRUE)nmaxit<-100
iterb<-0
while(S0b>error | iterb>nmaxit){
iterb<-iterb+1
F0b<-JIb%*%f1b
k0b<-((sum(1/F0b))^(-1))*(1/F0b)
if(sum(k0b)!=1)k0b<-k0b/sum(k0b)
k1b<-k0b
K0b<-Jb%*%k1b
f1b<-((sum(1/K0b))^(-1))*(1/K0b)
if(sum(f1b)!=1)f1b<-f1b/sum(f1b)
S0b<-max(abs(f1b-f0b))
f0b<-f1b
k0b<-k1b
}
FF0b<-matrix(data=0, ncol=1, nrow=nrow(DBB))
for(i in 1:nrow(DBB)){
FF0b[i,]<-sum(f1b[1:i,])
}
Sobb<-matrix(data=0, ncol=1, nrow=nrow(DBB))
for(j in 1:nrow(DBB)){
Sobb[j,]<-1-FF0b[j,]+f1b[j,]
}

M_IF0[,b]<-as.vector(FF0b)
M_IF01[,b]<-as.vector(f1b)
M_IF0Sob[,b]<-as.vector(Sobb)



ordB<-order(DBB[,3])
ord4<-sort(DBB[,3])
M4b[,b]<-ord4
k1bV[,b]<-k1b[ordV,]



for (i in 1:nrow(C)){
for(j in 1:nrow(C)){
if(M4b[i,b]<=kMV[j,1])
M_fV[j,b]<-sum(k1bV[1:i,b])
}}



}
M_IF0_sort<-matrix(0,nrow=nrow(C),ncol=B)
for(i in 1:nrow(M_IF0_sort)){
M_IF0_sort[i,]<-sort(M_IF0[i,])
}
if(is.na(alpha)==TRUE) alpha<-0.05
lowerF<-M_IF0_sort[,floor(alpha*B/2)]
upperF<-M_IF0_sort[,floor((1-alpha/2)*B)]

M_IF0_sort1<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_IF0_sort1)){
M_IF0_sort1[i,]<-sort(M_IF0Sob[i,])
}
lowerS<-M_IF0_sort1[,floor(alpha*B/2)]
upperS<-M_IF0_sort1[,floor((1-alpha/2)*B)]


M_fV_sort<-matrix(0,nrow=nrow(DBB),ncol=B)
for(i in 1:nrow(M_fV_sort)){
M_fV_sort[i,]<-sort(M_fV[i,])
}
lowerV<-M_fV_sort[,floor(alpha*B/2)]
upperV<-M_fV_sort[,floor((1-alpha/2)*B)]

upperU<-rep(NA,length(C))
lowerU<-rep(NA,length(C))


}

}


if(trunc=="both"){

if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(2,2))
plot(C[,1],FF,xlab="Time of interest",ylab="",main="Shen estimator",lty=1,type="l")
 lines(C[,1],upperF,lty=3)
 lines(C[,1],lowerF,lty=3)
plot(C[,1],Sob,type="l",xlab="Time of interest",ylab="",main="Survival",lty=1)
 lines(C[,1],upperS,lty=3)
 lines(C[,1],lowerS,lty=3)
plot(C[,1],fU,type="l",main="Marginal U",xlab="Time of interest",ylab="",lty=1)
 lines(C[,1],upperU,lty=3)
lines(C[,1],lowerU,lty=3)
plot(C[,1],fV,type="l",main="Marginal V",xlab="Time of interest",ylab="",lty=1)
 lines(C[,1],upperV,lty=3)
 lines(C[,1],lowerV,lty=3)
}

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,2))
plot(C[,1],fU,type="l",main="Marginal U",xlab="Time of interest",ylab="",lty=1)
 lines(C[,1],upperU,lty=3)
 lines(C[,1],lowerU,lty=3)
plot(C[,1],fV,type="l",main="Marginal V",xlab="Time of interest",ylab="",lty=1)
 lines(C[,1],upperV,lty=3)
 lines(C[,1],lowerV,lty=3)
}
if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(C[,1],FF,main="Shen estimator",xlab="Time of interest",ylab="",lty=1,type="l")
 lines(C[,1],upperF,lty=3)
 lines(C[,1],lowerF,lty=3)
plot(C[,1],Sob,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 lines(C[,1],upperS,lty=3)
 lines(C[,1],lowerS,lty=3)
}


if((display.FS==FALSE)&(display.UV==FALSE)){


}
}


if(trunc=="left"){

if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(2,2))
plot(C[,1],FF,xlab="Time of interest",ylab="",main="Shen estimator",lty=1,type="l")
 lines(C[,1],upperF,lty=3)
 lines(C[,1],lowerF,lty=3)
plot(C[,1],Sob,type="l",xlab="Time of interest",ylab="",main="Survival",lty=1)
 lines(C[,1],upperS,lty=3)
 lines(C[,1],lowerS,lty=3)
plot(C[,1],fU,type="l",main="Marginal U",xlab="Time of interest",ylab="",lty=1)
 lines(C[,1],upperU,lty=3)
lines(C[,1],lowerU,lty=3)
}

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(C[,1],fU,type="l",main="Marginal U",xlab="Time of interest",ylab="",lty=1)
 lines(C[,1],upperU,lty=3)
 lines(C[,1],lowerU,lty=3)
}
if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(C[,1],FF,main="Shen estimator",xlab="Time of interest",ylab="",lty=1,type="l")
 lines(C[,1],upperF,lty=3)
 lines(C[,1],lowerF,lty=3)
plot(C[,1],Sob,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 lines(C[,1],upperS,lty=3)
 lines(C[,1],lowerS,lty=3)
}


if((display.FS==FALSE)&(display.UV==FALSE)){


}
}


if(trunc=="right"){

if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,2))
plot(-C[,1],FF,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 lines(-C[,1],upperF,lty=3)
 lines(-C[,1],lowerF,lty=3)
plot(-C[,1],fV,type="l",main="Marginal V",xlab="Time of interest",ylab="",lty=1)
 lines(-C[,1],upperV,lty=3)
 lines(-C[,1],lowerV,lty=3)

}

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(-C[,1],fV,type="l",main="Marginal V",xlab="Time of interest",ylab="",lty=1)
 lines(-C[,1],upperV,lty=3)
 lines(-C[,1],lowerV,lty=3)

}
if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(-C[,1],FF,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 lines(-C[,1],upperF,lty=3)
 lines(-C[,1],lowerF,lty=3)

}

if((display.FS==FALSE)&(display.UV==FALSE)){

}
}
}

if (boot==FALSE){
upperF<-rep(NA,length(C))
lowerF<-rep(NA,length(C))
upperS<-rep(NA,length(C))
lowerS<-rep(NA,length(C))
upperU<-rep(NA,length(C))
lowerU<-rep(NA,length(C))
upperV<-rep(NA,length(C))
lowerV<-rep(NA,length(C))

if(trunc=="both"){

if(plot.joint==TRUE){
if(plot.type=="image"){
image(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution")
filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=rainbow)
filled.contour(sort(U),sort(V),KK,xlab="U",ylab="V",main="Joint distribution",color.palette=topo.colors)
}
if(plot.type=="persp"){
fcol<-topo.colors(10)[cut(KK[2:nrow(C),2:nrow(C)],10,include.lowest=TRUE)]
                   persp(x=sort(U),y=sort(V),KK,theta=-135,phi=40,col=fcol,xlab="U",ylab="V",zlab="Joint distribution")

}

}
if(plot.joint==FALSE){
}




if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(2,2))
plot(C[,1],FF,xlab="Time of interest",ylab="",main="Shen estimator",lty=1,type="l")
 plot(C[,1],Sob,type="l",xlab="Time of interest",ylab="",main="Survival",lty=1)
 plot(C[,1],fU,type="l",main="Marginal U",xlab="Time of interest",ylab="",lty=1)
 plot(C[,1],fV,type="l",main="Marginal V",xlab="Time of interest",ylab="",lty=1)
 }

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,2))
plot(C[,1],fU,type="l",main="Marginal U",xlab="Time of interest",ylab="",lty=1)
 plot(C[,1],fV,type="l",main="Marginal V",xlab="Time of interest",ylab="",lty=1)
 }
if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(C[,1],FF,main="Shen estimator",xlab="Time of interest",ylab="",lty=1,type="l")
 plot(C[,1],Sob,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 }


if((display.FS==FALSE)&(display.UV==FALSE)){


}
}


if(trunc=="left"){

if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(2,2))
plot(C[,1],FF,xlab="Time of interest",ylab="",main="Shen estimator",lty=1,type="l")
 plot(C[,1],Sob,type="l",xlab="Time of interest",ylab="",main="Survival",lty=1)
 plot(C[,1],fU,type="l",main="Marginal U",xlab="Time of interest",ylab="",lty=1)
 }

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(C[,1],fU,type="l",main="Marginal U",xlab="Time of interest",ylab="",lty=1)
 }
if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(C[,1],FF,main="Shen estimator",xlab="Time of interest",ylab="",lty=1,type="l")
 plot(C[,1],Sob,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 }


if((display.FS==FALSE)&(display.UV==FALSE)){


}
}


if(trunc=="right"){

if((display.FS==TRUE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,2))
plot(-C[,1],FF,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 plot(-C[,1],fV,type="l",main="Marginal V",xlab="Time of interest",ylab="",lty=1)
 
}

if((display.FS==FALSE)&(display.UV==TRUE)){
dev.new()
par(mfrow=c(1,1))
plot(-C[,1],fV,type="l",main="Marginal V",xlab="Time of interest",ylab="",lty=1)
 
}
if((display.FS==TRUE)&(display.UV==FALSE)){
dev.new()
par(mfrow=c(1,2))
plot(-C[,1],FF,type="l",main="Survival",xlab="Time of interest",ylab="",lty=1)
 
}

if((display.FS==FALSE)&(display.UV==FALSE)){

}
}


}

return(list(density=f,cumulative.df=FF,truncation.probs=F0,S0=S0,survival=Sob,density.joint=kMU[,2],
marginal.U=fU,marginal.V=fV,cummulative.joint=KK,
n.iterations=iter,Boot=boot.type,B=B,alpha=alpha,upper.df=upperF,lower.df=lowerF,upper.Sob=upperS,lower.Sob=lowerS,
upper.fU=upperU,lower.fU=lowerU,upper.fV=upperV,lower.fV=lowerV))
}

