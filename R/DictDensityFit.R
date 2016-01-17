DictDensityFit <- function(X,mu,Sigma,Concentration = 0.1,NumDraws = 1000){
###Estimate  weights based on dictionary densities (single-class model)
N = length(X)
K = length(mu)
fc = matrix(nrow = N, ncol = K)
tau = rep(1/K,K)
Zp = matrix(nrow=N,ncol = K)
Z = matrix(nrow=N,ncol = K)
tauVec = matrix(nrow= NumDraws,ncol=K) 
for(k in 1:K) fc[,k] = dtnorm(X,mu[k],Sigma[k],0,1)
for(w in 1:NumDraws){
	for(k in 1:K) Zp[,k] = fc[,k]*tau[k];
	for(i in 1:N) Z[i,] = rmultinom(1,1,Zp[i,]);
	n = colSums(Z);
	Gamma = Concentration+n;
	tau = rdirichlet(1,Gamma);
	tauVec[w,] = tau;
	#Concentration = ConcEst(tau)
}
return(colMeans(tauVec[floor(NumDraws/5):NumDraws,]))
}


