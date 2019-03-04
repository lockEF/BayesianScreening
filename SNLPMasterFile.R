# Structured and Non-Local Priors Code from 
#"Bayesian GWAS with Structured and Non-Local Priors by 
#A. Kaplan, M. Fiecas, and E.F. Lock
# February 26th 2019 #

#install "truncnorm package"
#install.packages("truncnorm")
require(truncnorm)

#Two sub-routines are coded in Rcpp 
#install.packages("Rcpp"); install.packages("RcppArmadillo")
require(Rcpp); require(RcppArmadillo); 
#Have to source the .cpp file that should be in the same folder as this file is.
sourceCpp("~/DPFunctionsCPP.cpp")

# Sub Routines for Master Chib#

#Calculates the log-sum to assist in averaging many densities to approximate the marginal likelihood under
#the altnerative #
#l should be a vector of density values 
logSum <- function(l){
  max(l) + log(sum(exp(l-max(l))))}

#The log-Non-Local density up to a proportionality constant (given that the point null is 0) 
#tau: proposed "effect size" for the non-local density
#theta: also alpha_1 in the paper (the non local support input)
#v and k are the other non-local parameters that are needed 
propnonloc <- function(tau, theta, v, k){
  res = log(k) + (v/2)*log(tau) - log(gamma(v/(2*k))) - ((v+1)/2)*log(theta^2) - ((theta^2)/(tau))^(-k)
  return(res)
}

#Metropolis Hastings sampling for the Non Local random variable
#NSIM: number of samples drawn from the Non-Local distribution 
#v, tau, k: Non local parameter values 
#jit: random normal jitter for the acceptance portion of the MH sampling; we set jit to 0.7
MHsamps <- function(NSIM, v, tau, jit, k){
  accept = 0
  theta = 1
  thetastor = numeric(NSIM)
  for (i in 1:NSIM){
    
    thetastar = theta + rnorm(1, 0, jit)
    
    if(propnonloc(tau,thetastar,v, k) - propnonloc(tau,theta,v, k) > log(runif(1))){theta = thetastar; accept = accept+1}else{
      theta = theta}
    thetastor[i] = theta;
  }
  aRATE = accept/NSIM
  res <- list("aRATE" = aRATE, "samps" = thetastor)
}
#aRate: acceptance rate of the MH sampler 
#samps: vector of drawn samples of length NSIM

#Approximating the marginal alternative likelihood
#nu: matrix, with the first column being the number of control individuals without minor allele presence 
# and the second column being the number of control individuals with minor allele presence for a given SNP
# the number of rows should be the number of SNPs one is analyzing 
#na: similar to na in structure, but for the case individuals 
#N = nu + na
#parameters = vector with the input in this order of the Non-Local parameters c( v, tau, k)
#if you have hundreds of thousands of SNPs, large = 1 allows for faster computation by 
#piecing up the evaluations and then reconstructing them for the evaluation of the 
#marginal likelihood 
MarginalAltLike <- function(nu,na,N,parameters,large = 1){
  
  v = parameters[1]; tau = parameters[2]; k = parameters[3]
  nonLocSamp <- MHsamps(10000, v = v, tau = tau, jit= .7, k = k)$samps #alpha 1 nonlocal
  nonlocsamps <- sample(nonLocSamp, 4000) 
  nullsamps <- rnorm(4000)
  probsNu <- pnorm(nullsamps - nonlocsamps)
  probsNa <- pnorm(nullsamps + nonlocsamps)
  #break into groups of about groups of about 10000 each -- helps with memory for ADNI application
  if(large == 1){
    ngroups <- ceiling(dim(nu)[1]/10000)
    groupIndex <- split(c(1:dim(nu)[1]), sample(ngroups,dim(nu)[1],replace=TRUE))
    ApproxMargLike<-rep(0,dim(nu)[1])
    for(i in 1:ngroups){
      #  print(i)
      nu.temp <- nu[groupIndex[[i]],]
      na.temp <- na[groupIndex[[i]],]
      likeNu <- outer(rep(1,4000),log(choose(rowSums(nu.temp),nu.temp[,2])))+outer(log(probsNu),nu.temp[,2])+outer(log(1-probsNu),nu.temp[,1])
      likeNu[is.nan(likeNu)]=0
      likeNa <- outer(rep(1,4000),log(choose(rowSums(na.temp),na.temp[,2])))+outer(log(probsNa),na.temp[,2])+outer(log(1-probsNa),na.temp[,1])
      likeNa[is.nan(likeNa)]=0
      AltB4Integrate <- likeNu+likeNa
      ApproxMargLike[groupIndex[[i]]] = (-log(4000) + apply(AltB4Integrate,2, function(x) (logSum(x))))
      gc()
    }
    P.Alt.Mark <- exp(ApproxMargLike)}else{
      
      likeNu <- sapply(1:nrow(nu), function(x) dbinom(nu[x,2], sum(nu[x,]),prob = probsNu, log = T))
      likeNa <- sapply(1:nrow(na), function(x) dbinom(na[x,2], sum(na[x,]),prob = probsNa, log = T))
      AltB4Integrate <- likeNu+likeNa
      ApproxMargLike = (-log(4000) + apply(AltB4Integrate,2, function(x) (logSum(x))))
      P.Alt.Mark <- exp(ApproxMargLike)
      rm(likeNu,likeNa,AltB4Integrate);gc()}
  return(P.Alt.Mark); 
}


#Initializing parameters and storage before the loop#
#Collect: how many samples one would like to have at the end of the analysis
#thin: thin the chain, collect every thin sample
#Burnin: warm-up period for the sampler to omit from collection 
#Gene: Vector, the length of the number of SNPs, with each entry being a numerical value 
# corresponding to the parent-gene that SNP is annotated for
# for our analyses and simulations we ran as.factor(GENE ANNOTATION VECTOR) 
# to obtain this vector. 
#X.Dat: marker-covariate matrix; number of cols = number of covariates; number of rows= number of SNPs
#colfx: number of columns of X.Dat; if colfx = 0 then the analysis runs without covariate estimation 
#H: number of stick breaks for the Dirichlet Process; we set this to 20
#nu: matrix, with the first column being the number of control individuals without minor allele presence 
# and the second column being the number of control individuals with minor allele presence for a given SNP
# the number of rows should be the number of SNPs one is analyzing 
#na: similar to na in structure, but for the case individuals 
#N = nu + na
#parameters: vector with the input in this order of the Non-Local parameters c(v, tau, k)
Initiz <- function(Collect, thin, Burnin, Gene, X.Dat, colfx, H, na, nu, N, parameters){
  
  Iter <- Collect*thin + Burnin
  samp <- length(Gene)
  NumGene = length(unique(Gene))
  OrdGene = as.numeric(as.factor(Gene))
  if (colfx == 0){p <- NumGene}else{p <- ncol(X.Dat) + NumGene}
  startcollect = 0
  
  MarkerIndexList <- list()
  for(i in 1:NumGene){
    MarkerIndexList[[i]] <-which(OrdGene==i)
  }
  
  betavec <- matrix(NA, nrow = p, ncol = Collect)
  Post.Marker <- numeric(samp)
  theta.hmat <- matrix(NA, nrow = H, ncol = Collect)
  sumcount.col <- matrix(NA, nrow = H, ncol = Collect)
  stick.col <- matrix(NA, nrow = H, ncol = Collect)
  Indic.Count <- matrix(NA, nrow = samp, ncol = Collect)
  alpha.col <- numeric(Collect)
  alpha.it <- 1
  flag <- 0 #this will count how many times we accept Alpha, the concentration parameter in the MH step#
  
  Zcolit <- rep(.5,samp)
  betavecit <- rep(1,p) 
  theta.h <- rep(1,H)
  Post.Markerit <- pbetait <- rep(.5,samp) 
  
  #Initiate DP Sticks#
  Sticks <- rep(.5, H)
  
  #If there are marker level covariates in the model we evaluate the inverse covariance matrix as well as it's sqrt#
  if (colfx != 0){
    sigmahat <- solve(diag(1, nrow = colfx) + t(X.Dat)%*%(X.Dat)) #evaluate variance for Betas#
    D = diag(eigen(sigmahat)$values) #finds diagonals 
    U = eigen(sigmahat)$vectors #finding rotation matrices of eigenvectors
    sqrsigma = U%*%sqrt(D) 
    rm(D);rm(U)
  }else{sqrsigma = NULL; sigmahat = NULL} #If there are not any covariates we do not evaluate these parameters
  
  P.Null.Mark <- mapply(function(x,y) beta(1 + x, 1 + y), N[,2], N[,1])*choose(rowSums(nu), nu[,2])*choose(rowSums(na), na[,2]) #P(X,Y|H_null) - taken from Lock and Dunson
  P.Alt.Mark <- MarginalAltLike(nu,na,N, parameters,large = 1)
  
  return(list(sigmahat, sqrsigma, Sticks, Post.Markerit,
              pbetait, betavec, Post.Marker, theta.hmat,
              sumcount.col, stick.col, alpha.col, Iter, NumGene, betavecit,
              samp, Zcolit, theta.h, alpha.it, P.Null.Mark, P.Alt.Mark, Indic.Count,
              MarkerIndexList, OrdGene))
}

####################################################
####Both of these functions evaluate ##############
######cumulative normal of the effects ############

#Input: 
#markGEff: vector of gene + XBeta effects, length = number of SNPs
# X.Dat, the covariate matrix 
# betavecit, the storage vector of updated covariate effects 
# colfx = number of columns of X.Dat (number of covariates) 
ObsBeta <- function(markGEff, X.Dat, betavecit, colfx){
  pbetait <- (c(markGEff + X.Dat%*%betavecit[1:colfx]))
  return(pbetait)
}

ObsBeta0 <- function(markGEff){
  pbetait <- (c(markGEff))
  return(pbetait)
}
####################################################


#Computes Posterior Probability#
# pbetait: output from ObsBeta or ObsBeta0 (cumulative normal of covariate+gene effects per SNP)
# Alt: Marginal Alternative Likelihoods evaluated for each SNP 
# Null: Marginal Null Likelihood evaluated for each SNP
# samp: number of SNPs 
Post.Prob <- function(pbetait, Alt, Null, samp){
  Post.Markit <- (pnorm(pbetait)* Alt)/((1-pnorm(pbetait)) *Null + pnorm(pbetait)*Alt)
  Alt.Status <- rbinom(samp, 1, prob = Post.Markit)
  
  N1 <- sum(Alt.Status)
  N0 = length(Alt.Status)-N1
  n <- length(Alt.Status)
  return(list(Post.Markit, Alt.Status, N1, N0, n))
}

###########################################################################

#Multinomial Sample Each Column#
#Now sample cluster ID#
#Multinomial Sample#
#ProbMat: output from Rcpp function "multinomprobCPP" or "multinomprobCPP0"
#ProbMat: matrix of probabilities of a gene being allocated to a certain DP atom 
# number of rows = number of atoms; number of columns = number of genes 
ClustAssignG2 <- function(ProbMat){
  indc <- numeric(ncol(ProbMat))
  for (i in 1:ncol(ProbMat)){
    prob.col <- ProbMat[,i]
    if (sum(prob.col) == 0){indc[i] <- sample(1:length(prob.col), size = 1)}else{
      indc[i] <- sample(1:length(prob.col), size = 1, prob = prob.col)
    }
  }
  sum.count <- sapply(1:nrow(ProbMat), function(x) sum(indc == x))
  return(list(indc, sum.count))
}
#output: (1) vector of length=number of genes that has atom ID for each gene
# sum.count is the total number of allocated genes per atom   

# Stick and Alpha Update #
# Log-Posterior for Alpha when we use MH-Step to sample Alpha
#log posterior of alpha#
#H: number of stick breaks 
#a, b: prior parameters for Gamma(a,b) prior
#v: vector of stick break probabilities 
post.alph <- function(H, alpha, v, a, b){
  val <- (H+a-2)*log(alpha) - (alpha/b) + (alpha-1)*sum(log(1-v)) 
  return(val)
}

# Update steps for sticks-V and Alpha
#H: number of stick breaks
#sum.count: number of genes per atom assignment 
# alpha: last alpha value sampled 
# hyp: vector of prior parameters for Gamma(a,b) prior: c(a,b)
VandAlph <- function(H, sum.count, alpha, hyp){
  sum.more <- c()
  V.new <- numeric(H)
  V.new[H] <- 1
  #Find count of cases per cluster ii greater than cluster ii
  for (ii in 2:H){
    sum.more[ii-1] <- sum(as.numeric(sum.count[(ii):H]))
    #Update V Weights#
    V.new[ii-1] <- rbeta(1 ,1 + sum.count[ii-1], alpha + sum.more[ii-1])
  }
  
  Sticks <- c(1, cumprod(1-V.new[-H]))*V.new
  #update Sticks, which is the DP weights
  #Correct the sticks that equal 1, other than the last one#
  V.new[-H] <- ifelse(V.new[-H] == 1, V.new[-H] - runif(1, 0, .00001), V.new[-H])
  
  alpha.prop <- alpha + runif(1, -1, 1)
  if(alpha.prop < 0){alpha <- alpha}else{ #accept reject step#
    new = post.alph(H=H, alpha = alpha.prop, v = V.new[-H], a = hyp[1], b = hyp[2]); 
    old = post.alph(H=H, alpha = alpha, v = V.new[-H], a = hyp[1], b = hyp[2]);
    if(new-old > log(runif(1))){
      alpha <- alpha.prop
    }
  }
  return(list(alpha, Sticks))
}
#output: accepted or present sample alpha (numeric), and sticks is a vector of weights 

## Atoms Update DP ##
#Input: 
#samp: number of SNPs; H: Number of stickbreaks; NumGene: Number of unique Genes; indc: gene-atom allocation vector of 
#length of genes; 
# Zcolit: update vector for latent Normal values Z
# X.Dat: marker-covariate matrix 
# betavecit: update vector for gene+covariate effects 
# colfx: number of covariates (excluding number of genes)
# MarkerIndicesList: a list collected from the Initialize function 
#this function is used if covariates are considered 
AtomAssign <- function(samp, H, NumGene, indc, Zcolit, X.Dat, betavecit, colfx, MarkerIndicesList){
  clustMark1 <- numeric(samp)
  clusterfreq1 <- numeric(H)
  theta.h <- numeric(H)
  for (n in 1:NumGene){
    clustMark1[MarkerIndicesList[[n]]] <- indc[n]
  }
  
  clusterfreq1 <- sapply(1:H, function(x) sum(clustMark1 == x))
  
  for (i in 1:H){
    if (clusterfreq1[i] != 0){
      mean.s.i <- (1 + sum(Zcolit[(which(clustMark1 == i))] - X.Dat[(which(clustMark1 == i)),]%*%betavecit[1:colfx]))/(1 + clusterfreq1[i])
      var.s.i <- 1/(1 + clusterfreq1[i])
      theta.h[i] <- rnorm(1, mean =  mean.s.i, sd = sqrt(var.s.i))}else{theta.h[i] <- rnorm(1, 0, 1)}
    betavecit[(which(indc ==i) + colfx)] <- theta.h[i]
  }
  
  return(list(betavecit, theta.h))
}

#Same function as AtomAssign but does not consider marker-level covariates 
AtomAssign0 <- function(samp, H, NumGene, indc, Zcolit, MarkerIndicesList){
  clustMark1 <- numeric(samp)
  clusterfreq1 <- numeric(H)
  theta.h <- numeric(H)
  betavecit <- numeric(NumGene)
  for (n in 1:NumGene){
    clustMark1[MarkerIndicesList[[n]]] <- indc[n]
  }
  
  clusterfreq1 <- sapply(1:H, function(x) sum(clustMark1 == x))
  
  for (i in 1:H){
    if (clusterfreq1[i] != 0){
      mean.s.i <- (1 + sum(Zcolit[(which(clustMark1 == i))]))/(1 + clusterfreq1[i])
      var.s.i <- 1/(1 + clusterfreq1[i])
      theta.h[i] <- rnorm(1, mean =  mean.s.i, sd = sqrt(var.s.i))}else{theta.h[i] <- rnorm(1, 0, 1)}
    betavecit[which(indc ==i)] <- theta.h[i]
  }
  
  return(list(betavecit, theta.h))
}

#Latent variable Z update 
#Input: 
#Zcolit: last iterations sample vector of Z's
#Alt.Status: output from PostProbabilities, vector of length SNPs that has indicators for SNP association (1) or not (0)
#N0: Number of SNPs with Alt.Status == 0; number of null SNPs for that iteration 
#N1: Number of SNPs with Alt.Status == 1; number of alternative SNPs for that iteration
#obs.beta: is the vector of gene effect + XBeta (constructed in the sampler)  
ZUpdate <- function(Zcolit, Alt.Status, N0, N1, obs.beta){
  avec = rep(-Inf, times = N0+N1)
  bvec=rep(0, times = N0+N1)
  avec[Alt.Status==1]=0
  bvec[Alt.Status==1]=Inf
  Zcolit <- rtruncnorm(n = N0+N1, a= avec, b = bvec, mean = obs.beta, sd = 1)
  return(Zcolit)
}

#### Chib Fixed Parameters Update ###
# Update marker-level covariates
#Input: 
#Zcolit: from above function "ZUpdate"
#markGEff: vector of gene + XBeta effects, length = number of SNPs
#sigmahat: computed in the Initialize step; the variance matrix of the covariate matrix 
#sqrsigma: square-root of the sigmahat matrix 
#betavecit: storage of last iteration's beta samples (which includes gene effects)
#X.Dat: covariate matrix 
#colfx: number of covariates (that are not genes) 
ChibUpdate <- function(Zcolit, markGEff, sigmahat, sqrsigma, betavecit, X.Dat,colfx){
  DifVec <- Zcolit - markGEff
  betahat <- sigmahat%*%(t(X.Dat)%*%(DifVec))
  
  betavecit[c(1:colfx)] <- sqrsigma%*%rnorm(length(betahat), mean = 0, sd = 1) + betahat
  
  return(betavecit)
}


#Wrap function: 
#Collect: how many samples one would like to have at the end of the analysis
#thin: thin the chain, collect every thin sample
#Burnin: warm-up period for the sampler to omit from collection 
#Gene: Vector, the length of the number of SNPs, with each entry being a numerical value 
# corresponding to the parent-gene that SNP is annotated for
# for our analyses and simulations we ran as.factor(GENE ANNOTATION VECTOR) 
# to obtain this vector. 
#X.Dat: marker-covariate matrix; number of cols = number of covariates; number of rows= number of SNPs
#colfx: number of columns of X.Dat; if colfx = 0 then the analysis runs without covariate estimation 
#H: number of stick breaks for the Dirichlet Process; we set this to 20
#nu: matrix, with the first column being the number of control individuals without minor allele presence 
# and the second column being the number of control individuals with minor allele presence for a given SNP
# the number of rows should be the number of SNPs one is analyzing 
#na: similar to na in structure, but for the case individuals 
#N = nu + na
#parameters: vector with the input in this order of the Non-Local parameters c(v, tau, k)
# hyp: vector of prior parameters for Gamma(a,b) prior: c(a,b)
DPMasterChibNonLoc <- function(Collect, Gene, Burnin, X.Dat, nu, na, N, H, colfx, thin, hyp, parameters){
  
  Initial <- Initiz(Collect, thin, Burnin, Gene, X.Dat, colfx, H, na, nu, N, parameters)  
  
  sigmahat = Initial[[1]]; sqrsigma = Initial[[2]]; Sticks = Initial[[3]];
  Post.Markit = Initial[[4]]; pbetait = Initial[[5]]; 
  
  betavec = Initial[[6]];
  
  Post.Marker = Initial[[7]]; theta.hmat = Initial[[8]]; sumcount.col = Initial[[9]];
  stick.col = Initial[[10]]; alpha.col = Initial[[11]]; Iter = Initial[[12]]; NumGene = Initial[[13]];
  betavecit = Initial[[14]];samp = Initial[[15]]; Zcolit = Initial[[16]]; theta.h = Initial[[17]];
  alpha = Initial[[18]]; H = H;  P.Null.Mark = Initial[[19]]; P.Alt.Mark = Initial[[20]]; Indic.Count = Initial[[21]]
  MarkerIndicesList = Initial[[22]]; OrdGene = Initial[[23]]
  
  rm(Initial); rm(sumcount.col); rm(theta.hmat); rm(stick.col);
  
  startcollect = 0;
  #Gibbs Sampler#
  for (nsim in 1:Iter){
    #Covariates are used#
    if(colfx != 0){
      
      MarkEff <- betavecit[colfx + OrdGene] 
      obs.beta <- ObsBeta(markGEff = MarkEff, X.Dat = X.Dat, betavecit = betavecit, colfx = colfx)
      PostProbabilities <- Post.Prob(obs.beta, P.Alt.Mark, P.Null.Mark, samp)
      
      #DP Update 
      ClustProbs = multinomprobCPP(Zcolit, Gene-1, colfx, theta.h, Sticks, betavecit, X.Dat, NumGene, H);
      ClustAssign <- ClustAssignG2(ClustProbs)
      UpdateVAlph <- VandAlph(H, ClustAssign[[2]], alpha, hyp); alpha = UpdateVAlph[[1]]; Sticks = UpdateVAlph[[2]];
      NewAtoms <- AtomAssign(samp, H, NumGene, ClustAssign[[1]], Zcolit, X.Dat, betavecit, colfx, MarkerIndicesList)
      
      #Update Betavecit with new Atoms#
      betavecit <- NewAtoms[[1]]
      theta.h <- NewAtoms[[2]]
      
      MarkEff <- betavecit[colfx + OrdGene] 
      obs.beta <- ObsBeta(MarkEff, X.Dat, betavecit, colfx)
      
      #Update Chib Parameters#
      betavecit <- ChibUpdate(Zcolit, MarkEff, sigmahat, sqrsigma, betavecit, X.Dat,colfx)
    }
    
    #Covariates were not used#
    if(colfx == 0){
      
      MarkEff <- betavecit[OrdGene] 
      obs.beta <- ObsBeta0(MarkEff)
      PostProbabilities <- Post.Prob(obs.beta, P.Alt.Mark, P.Null.Mark, samp)
      
      #DP Steps#
      ClustProbs <- multinomprob0CPP(Zcolit, Gene-1, theta.h, Sticks, NumGene, H)
      ClustAssign <- ClustAssignG2(ClustProbs);
      UpdateVAlph <- VandAlph(H, ClustAssign[[2]], alpha, hyp); alpha = UpdateVAlph[[1]]; Sticks = UpdateVAlph[[2]];
      NewAtoms <- AtomAssign0(samp, H, NumGene, ClustAssign[[1]], Zcolit, MarkerIndicesList);
      
      #Update Betavecit with new Atoms#
      betavecit <- NewAtoms[[1]]
      theta.h <- NewAtoms[[2]]
      
      
      MarkEff <- betavecit[OrdGene] 
      obs.beta <- ObsBeta0(MarkEff)
    }

    Zcolit <- ZUpdate(Zcolit, PostProbabilities[[2]], PostProbabilities[[3]], PostProbabilities[[4]], obs.beta); 
    
    
    if(nsim > Burnin & nsim%%thin == 0){
      if(startcollect%%100 == 0){
      print(sprintf("Sample Collected # %d", startcollect))}
      
      betavec[,startcollect+1] <- betavecit
      Post.Marker = Post.Marker + PostProbabilities[[1]]
      alpha.col[startcollect+1] <- UpdateVAlph[[1]]
      Indic.Count[,startcollect + 1] <- PostProbabilities[[2]]
      startcollect = startcollect + 1
    }
  }
  return(list("Beta_Mat" = betavec, "Post.Mark" = (Post.Marker)/Collect, "PostAlph" = alpha.col,
             "Indic.Count" = rowMeans(Indic.Count)))
}
#output: Beta_Mat is a matrix of Collect x Number of Covariates + number of Genes
# 1:colfx corresponds to the covariate columns, while the rest are gene effects 
#Post.Mark: posterior average of SNP probability of being associated given the data
#alpha.col: vector of samples of the concentration parameter alpha
# Indic.Count: Empirical estimate of Posterior Probability of Association
# average of Alt.Indic per SNP, used in non-local parameter determination
