EstimateTruncatedDictionary <- function(X,K=2, a0 = 0.5,b0 = 0.5,mu0 = 0.5,Concentration = 0.5,NumDraws = 1000){
M = dim(X)[1]
N = dim(X)[2]
Gamma = matrix(rep(1,M*K),nrow = M,ncol = K)
StDev = sd(as.vector(X))
a0 = 1
b0 = StDev^2
mu = seq(from = 1/(2*K), to= 1-1/(2*K), length.out = K)
Sigma = rep(1/(2*K),K)
S =  matrix(nrow = 1,ncol = K)
PostMean =  matrix(nrow = 1,ncol = K)
A =  matrix(nrow = 1,ncol = K)
B = matrix(nrow = 1,ncol = K)
Tau = matrix(nrow = 1,ncol = K) 
Lambda = matrix(nrow = 1,ncol = K) 
Z = array(dim = c(M,N,K))
logF = array(dim = c(M,N,K))
logZ = array(dim = c(M,N,K))
tau = matrix(rep(1/K,M*K),nrow = M, ncol = K)
n = matrix(nrow = M, ncol = K)
Y=X
muVec = matrix(nrow = NumDraws, ncol = K)
SigmaVec = matrix(nrow = NumDraws,ncol = K)
ConcentrationVec = matrix(nrow = NumDraws,ncol = K)
##Run MCMC
for(w in 1:NumDraws){
	for(k in 1:K) logF[,,k] = log(tau[,k])-log(Sigma[k])+(-log(2*pi)-((X-mu[k])/Sigma[k])^2)/2-log(pnorm(1,mu[k],Sigma[k])-pnorm(0,mu[k],Sigma[k]))
	Zp = exp(logF)	
	for(m in 1:M){
 		for(i in 1:N) Z[m,i,] = rmultinom(1,1,Zp[m,i,])
		for(k in 1:K){
 			CDF = ptnorm(X[m,Z[m,,k]==1],mu[k],Sigma[k],0,1)
			Y[m,Z[m,,k]==1] = qnorm(CDF,mu[k],Sigma[k])}
		n[m,] = colSums(Z[m,,])
}
	num = colSums(n);
	for(k in 1:K){
		Ylist = rep(0,num[k])
		for(m in 1:M) if(n[m,k]>0) Ylist[(sum(n[0:(m-1),k])+1):sum(n[1:m,k])] = as.vector(Y[m,Z[m,,k]==1])
		if(num[k]>1){
			S[k] = sd(Ylist)^2
			PostMean[,k] = sum(Ylist)/(num[k]+1)
			B[,k] = b0+0.5*(num[k]*S[k]+num[k]*(mean(Ylist)-mu0)^2/(1+num[k]))}
		if(num[k]==1){
			PostMean[,k] = (mu0+Ylist)/2
			B[,k] = b0+0.5*(Ylist-mu0)^2/2}
		if(num[k]==0){
			PostMean[,k] = mu0
			B[,k] = b0}
		Lambda[k] = 1+num[k]
		A[,k] = a0+num[k]/2
		Tau[,k] = rgamma(1,shape=A[,k],rate=B[,k])
		mu[k] = rnorm(1,PostMean[,k],sqrt(1/(Tau[,k]*Lambda[k])))
		Sigma[k] = sqrt(1/Tau[,k])}
 	     
		for(m in 1:M){
 			Gamma[m,] = Concentration+n[m,];
		    tau[m,] = rdirichlet(1,Gamma[m,])}
		
		muVec[w,]= mu
		SigmaVec[w,] = Sigma
	Concentration = ConcEst(tau)
	ConcentrationVec[w,] = Concentration
} 
#Save dictionary estimates
Mu = colMeans(muVec[NumDraws/5:NumDraws,])
Sigma = colMeans(SigmaVec[NumDraws/5:NumDraws,])
Concentration = colMeans(ConcentrationVec[NumDraws/5:NumDraws,])
return(list(Mu=Mu,Sigma=Sigma,Concentration=Concentration))
}

####Additional functions:

ConcEst <- function(Probs){
##Finds an approximate MLE for Dirichlet parameters using method of moments as described in ROnning (1989); see also
##http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf
Ep = colMeans(Probs)
Varp = apply(Probs,2,'var')
Term = Ep*(1-Ep)/Varp-1
Tot = exp(mean(log(Term[1:(length(Term)-1)])))
Alpha = Ep*Tot
return(Alpha)}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

logSum <- function(l){max(l) + log(sum(exp(l-max(l))))}