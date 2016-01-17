DictFitTestFixedP <- function(X,Class,mu,Sigma,pA=0.5,Concentration = 0.5,NumDraws = 1000){
#Concentration = 0.1; NumDraws = 1100; mu= Dict$Mu; Sigma=Dict$Sigma 
#pA=0.5;
N = dim(X)[2]
M = dim(X)[1]
K = length(mu)
fc = array(dim = c(M,N,K))
tau0 = matrix(rep(1/K,M*K),nrow = M, ncol = K)
tau1 = matrix(rep(1/K,M*K),nrow = M, ncol = K)
tau00 = matrix(rep(1/K,M*K),nrow = M, ncol = K)
tau11 = matrix(rep(1/K,M*K),nrow = M, ncol = K)
Zp = array(dim = c(M,N,K))
Z = array(dim = c(M,N,K))
for(m in 1:M)for(k in 1:K) fc[m,,k] = dtnorm(X[m,],mu[k],Sigma[k],0,1)
A = rbinom(M,1,pA)==1;
APvec = matrix(nrow = M, ncol=NumDraws);
tau0vec = array(dim=c(M,NumDraws,K));
tau1vec = array(dim=c(M,NumDraws,K));
for(w in 1:NumDraws){
	for(k in 1:K){ 
			Zp[,Class==0,k] = fc[,Class==0,k]*tau0[,k];
	Zp[,Class==1,k] = fc[,Class==1,k]*tau1[,k];} 
for(m in 1:M){
	for(i in 1:N) Z[m,i,] = rmultinom(1,1,Zp[m,i,]);
	n = colSums(Z[m,,]);
	Gamma = Concentration+n;
	n0 = colSums(Z[m,Class==0,])
	n1 = colSums(Z[m,Class==1,])
	Gamma0= Concentration+n0;
	Gamma1 = Concentration+n1;
	pzA = log(pA)+mlbeta(Gamma)-mlbeta(Gamma-n)
	pzNotA= log(1-pA)+mlbeta(Gamma0)+mlbeta(Gamma1)-mlbeta(Gamma0-n0)-mlbeta(Gamma1-n1)
	logpz = logSum(c(pzA,pzNotA))
	logpAcond = pzA-logpz;
	pAcond = exp(logpAcond)
	A[m] = rbinom(1,1,pAcond)==1;
	tau = rdirichlet(1,Gamma)
	tau00[m,] =rdirichlet(1,Gamma0)
	tau11[m,] = rdirichlet(1,Gamma1)
	tau0[m,] = (1-pAcond)*tau00[m,]+pAcond*tau
	tau1[m,] = (1-pAcond)*tau11[m,]+pAcond*tau
#	if(A[m]){ 
#		tau0[m,] = rdirichlet(1,Gamma);
#		tau1[m,]= tau0[m,];}
#	if(!A[m]){
#		tau0[m,] = rdirichlet(1,Gamma0)
#		tau1[m,] = rdirichlet(1,Gamma1)}
	APvec[m,w] = pAcond;
}
	tau0vec[,w,] = tau0;
	tau1vec[,w,] = tau1;
	#Concentration = ConcEst(rbind(tau00,tau11))
	#pA = rbeta(1,1+sum(A),1+M-sum(A))
}
Results = list(postDraws=APvec,tao0Draws=tau0vec,tao1Draws=tau1vec)
return(Results)
}

logSum <- function(l){max(l) + log(sum(exp(l-max(l))))}

mbeta <- function(...) { 
    exp(sum(lgamma(c(...)))-lgamma(sum(c(...))))
}
mlbeta <- function(...) { 
    sum(lgamma(c(...)))-lgamma(sum(c(...)))
}
