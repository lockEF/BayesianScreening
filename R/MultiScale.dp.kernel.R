MultiScale.dp.kernel = function(X,Class,UniqGene,Gene,mu,Sigma,Concentration=0.5,alpha=1,KK=10,NumDraws=1000){

N = dim(X)[2]
M = dim(X)[1]
K = length(mu)
fc = array(dim = c(M,N,K))
tau0 = matrix(rep(1/K,M*K),nrow = M, ncol = K)
tau1 = matrix(rep(1/K,M*K),nrow = M, ncol = K)
tau00 = matrix(rep(1/K,M*K),nrow = M, ncol = K)
tau11 = matrix(rep(1/K,M*K),nrow = M, ncol = K)
tau = matrix(rep(1/K,M*K),nrow = M, ncol = K)
Zp = array(dim = c(N,K))
Z = array(dim = c(N,K))
for(m in 1:M)for(k in 1:K) fc[m,,k] = dtnorm(X[m,],mu[k],Sigma[k],0,1)

NN = rep(0,KK)
Probs = matrix(nrow = length(UniqGene), ncol = KK)
C = rep(1,length(UniqGene))
V = rep(0,KK)
Pi = sort(rdirichlet(1,rep(alpha,KK)), decreasing=TRUE)
a = rep(1,KK)
b = rep(1,KK)
pp = rbeta(KK,1,1)
pps = matrix(ncol = NumDraws,nrow=KK)
Pis = matrix(ncol = NumDraws,nrow=KK)
pAs = matrix(ncol = NumDraws,nrow = length(UniqGene))
posts = matrix(ncol = NumDraws,nrow = length(Gene))
post = rep(0.5,length(Gene))
Realize = rbinom(length(Gene),1,post)
pA = rep(0.5,length(UniqGene))
pzA.marg = rep(0,length(Gene))
pzNotA.marg = rep(0,length(Gene))
pAvec = rep(0,length(Gene))
nG = rep(0, length(UniqGene))
kG = rep(0, length(UniqGene))
for(i in 1:length(UniqGene)) nG[i] = sum(Gene==UniqGene[i])
Indices = list()
for(i in 1:length(UniqGene)) Indices[[i]] = which(Gene==UniqGene[i])

for(j in 75:NumDraws){
for(m in 1:M){
  Zp[Class==0,] = t(t(fc[m,Class==0,])*tau0[m,]);
  Zp[Class==1,] = t(t(fc[m,Class==1,])*tau1[m,]);
  Z = t(apply(Zp,1,rmultinom,size=1,n=1))
	n0 = colSums(Z[Class==0,])
	n1 = colSums(Z[Class==1,])
  n=n0+n1
  Gamma = Concentration+n;
	Gamma0 = Concentration+n0;
	Gamma1 = Concentration+n1;
	pzA.marg[m] = mlbeta(Gamma)-mlbeta(Gamma-n)
	pzNotA.marg[m] = mlbeta(Gamma0)+mlbeta(Gamma1)-mlbeta(Gamma0-n0)-mlbeta(Gamma1-n1) 
	tau[m,] = rdirichlet(1,Gamma)
	tau00[m,] =rdirichlet(1,Gamma0)
	tau11[m,] = rdirichlet(1,Gamma1)}
	tau0 = (1-post)*tau00+post*tau
	tau1 = (1-post)*tau11+post*tau
for(i in 1:length(UniqGene)){pAvec[Indices[[i]]] = pA[i]}  ##set gene-level probabilities for all sites
	pAs[,j] = pA
	pzA = log(pAvec)+pzA.marg
	pzNotA= log(1-pAvec)+pzNotA.marg
	logpz = apply(cbind(pzA,pzNotA),1,logSum) ##find sum
	logpAcond = pzA-logpz;
post = exp(logpAcond)
Realize = rbinom(length(post),1,post)
posts[,j] = post
for(i in 1:length(UniqGene)){kG[i] = sum(Realize[Indices[[i]]])}

###Find Data probs, estimate C
for(k in 1:KK){
  Probs[,k] = Pi[k]*dbinom(kG,nG,pp[k])
}
for(i in 1:length(pA)){
  C[i] = c(1:KK)[rmultinom(1,1,Probs[i,])==1]
  pA[i] = pp[C[i]]
}
#Update Dirichlet parameters
 for(k in 1:KK) NN[k] = sum(C==k)
for(k in 1:(KK-1)){
	V[k] = rbeta(1,1+NN[k],alpha+sum(NN[(k+1):KK]))	
}
V[KK]= 1

##Weights:
Pi[1]=V[1]
for(k in 2:KK){ 
	Pi[k] = V[k]*prod(1-V[1:(k-1)])
}

#update pp
for(k in 1:KK) pp[k] = rbeta(1,1+sum(kG[C==k]), 1+sum(nG[C==k])-sum(kG[C==k]))

pps[,j]=pp
Pis[,j]=Pi
}
if(length(UniqGene)==1)pA.est = mean(pAs[floor(NumDraws/4):NumDraws])
if(length(UniqGene)>1) pA.est = rowMeans(pAs[,floor(NumDraws/4):NumDraws])
post.est = rowMeans(posts[,floor(NumDraws/4):NumDraws])
return(list(pG = pA.est,posts = post.est))

}

mbeta <- function(...) { 
  exp(sum(lgamma(c(...)))-lgamma(sum(c(...))))
}
mlbeta <- function(...) { 
  sum(lgamma(c(...)))-lgamma(sum(c(...)))
}
logSum <- function(l){max(l) + log(sum(exp(l-max(l))))}
	


