MultiScale.DP.binom = function(n,n0,n1,UniqGene,Gene,Concentration=c(1,1),alpha=1,K=10,NumDraws=2000){
N = rep(0,K)
Probs = matrix(nrow = length(UniqGene), ncol = K)
C = rep(1,length(UniqGene))
V = rep(0,K)
Pi = rep(0,K)
a = rep(1,K)
b = rep(1,K)
pp = rbeta(K,1,1)
pAs = matrix(ncol = NumDraws,nrow = length(UniqGene))
posts = matrix(ncol = NumDraws,nrow = length(Gene))
post = rep(0.5,length(Gene))
Realize = rbinom(length(Gene),1,post)
pA = rep(0,length(UniqGene))
Gamma=matrix(nrow=length(Gene),ncol = 2)
Gamma0=matrix(nrow=length(Gene),ncol = 2)
Gamma1=matrix(nrow=length(Gene),ncol = 2)
pzA.marg = rep(0,length(Gene))
pzNotA.marg = rep(0,length(Gene))
pAvec = rep(0,length(Gene))
nG = rep(0, length(UniqGene))
kG = rep(0, length(UniqGene))
for(i in 1:length(Gene)){
	Gamma[i,] = Concentration+n[i,];
	Gamma0[i,]= Concentration+n0[i,];
	Gamma1[i,] = Concentration+n1[i,];
	pzA.marg[i] = mlbeta(Gamma[i,])-mlbeta(Gamma[i,]-n[i,])
	pzNotA.marg[i] = mlbeta(Gamma0[i,])+mlbeta(Gamma1[i,])-mlbeta(Gamma0[i,]-n0[i,])-mlbeta(Gamma1[i,]-n1[i,]) 
}
logpz = apply(cbind(pzA.marg,pzNotA.marg),1,logSum)
	logpAcond = pzA.marg-logpz;
	post = exp(logpAcond)
for(i in 1:length(UniqGene)) nG[i] = sum(Gene==UniqGene[i])
Indices = list()
for(i in 1:length(UniqGene)) Indices[[i]] = which(Gene==UniqGene[i])

for(j in 1:NumDraws){
Realize = rbinom(length(post),1,post)

#Update Dirichlet parameters
for(k in 1:K) N[k] = sum(C==k)
for(k in 1:(K-1)) V[k] = rbeta(1,1+N[k],alpha+sum(N[(k+1):K]))	
V[K]= 1

##Weights:
Pi[1]=V[1]
for(k in 2:K){ 
	Pi[k] = V[k]*prod(1-V[1:(k-1)])
}

####Find data probs; estimate C
for(k in 1:K){
	Probs[,k] = Pi[k]*dbinom(kG,nG,pp[k])
}
for(i in 1:length(pA)){
	C[i] = c(1:K)[rmultinom(1,1,Probs[i,])==1]
	pA[i] = pp[C[i]]
}
for(i in 1:length(UniqGene)){pAvec[Indices[[i]]] = pA[i]}
pAs[,j] = pA

for(k in 1:K) pp[k] = rbeta(1,1+sum(kG[C==k]), 1+sum(nG[C==k])-sum(kG[C==k]))

for(i in 1:length(UniqGene)){kG[i] = sum(Realize[Indices[[i]]])}
pzA = log(pAvec)+pzA.marg
pzNotA= log(1-pAvec)+pzNotA.marg
logpz = apply(cbind(pzA,pzNotA),1,logSum)
logpAcond = pzA-logpz;
post = exp(logpAcond)
posts[,j] = post

}

if(length(UniqGene)==1)pA.est = mean(pAs[floor(NumDraws/4):NumDraws])
if(length(UniqGene)>1) pA.est = rowMeans(pAs[,floor(NumDraws/4):NumDraws])
post.est = rowMeans(posts[,floor(NumDraws/4):NumDraws])
return(list(pG = pA.est,posts = post.est))
}

####Additional functions:


mbeta <- function(...) { 
  exp(sum(lgamma(c(...)))-lgamma(sum(c(...))))
}
mlbeta <- function(...) { 
  sum(lgamma(c(...)))-lgamma(sum(c(...)))
}
logSum <- function(l){max(l) + log(sum(exp(l-max(l))))}


	


