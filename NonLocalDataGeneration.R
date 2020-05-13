#Data Generation for Non-Local data#
#q: random number seed
#NumGene: number of genes to simulate 
#indiv: number of maximum possible SNPs per gene to simulate
#Cats: number of marker level categorical variables 
#Cont: number of marker level continuous variables
#N.Control: number of control (phenotype negative) subjects
#N.Cases: number of cases (phenotype positive) subjects 
#setting: 
# = 0: all null gene effects 
# = 1: 20% of gene effects are +3/ 80% of gene effects are -3 (80% null, 20% positive gene eff)
# = 2: Probability of null follows 1-Beta(1,0.2) distribution (mostly null with minimal signal)
# = anything else: random normal distribution 
#k, V, tau: non-local parameters to specify for data generation 
source("SNLPMasterFile.R") # <-- Need to source this for MHSamps function

DataGenNonLoc.Spec <- function(q, NumGene, indiv, Cats, Cont, N.Control, N.Case, setting = 0, k, V, tau){
  
  set.seed(q)
  
  Ngene=NumGene
  UniqGene = c(1:Ngene)
  Gene = c() ##Gives the gene number for each marker
  for(i in 1:Ngene){
    Nmarkers = sample(c(1:indiv),1)+1 ##number of markers for each gene is taken  from 2,..,indiv with equal probability
    Gene=c(Gene,rep(UniqGene[i],Nmarkers))
  }
  
  #retain order#
  
  orderGene <- 1:length(Gene)
  ################################# have fixed effects influence SNP probability #########
  ##################### Edit this into the estimation part later ####################
  if(Cats >= 1){
    ##SNP level effects##
    #First lets do SNP categories#
    cat.probs <- rbeta(Cats, 1,1)
    cat.probs <- cat.probs/sum(cat.probs)
    cats <- sample(1:Cats, size = length(Gene), prob = cat.probs, replace = TRUE)
    #cats <- model.matrix(~factor(cats)-1) <- changed on May 13th 2020
    cats <- model.matrix(~factor(cats))
    cats <- cats[,-1]
  }else{cats = NULL}
  # the original code that is commented out did not produce a matrix of dummy variables for 
  # SNP level categories without an intercept term 
  #now simulate something continuous#
  
  if (Cont >= 1){
    m.cont <- matrix(NA, nrow = length(Gene), ncol = Cont)
    for (i in 1:Cont){
      mean.norm <- runif(1)
      m.cont[,i] <- rgamma(length(Gene), 2, 1) - rnorm(length(Gene), mean.norm, .5)
    }
  }else{m.cont = NULL}
  ###########################################################################
  #concatenate the new data with the gene level data##
  X <- model.matrix(~ factor(Gene) - 1)
  thetagene <- numeric(NumGene)
  
  #assign gene probabilities randomly
  if (setting == 0){
    thetagene <- rep(-4, NumGene)}else if(setting == 1){
      thetagene <- c(rep(3, .2*NumGene), rep(-3, .8*NumGene))}else if(setting == 2){
        thetagene <- qnorm(1 - rbeta(NumGene, 1, .2))}else{thetagene <- rnorm(NumGene, 0, 1)}
  
  if (Cont == 0 & Cats ==0){
    X.tot <- X
    X.fixed <- NULL
    betatrue <- thetagene}else{
      X.fixed <- cbind(cats, m.cont)
      #simulate from the model#
      betatruefix <- rnorm(Cont + Cats, 0, 1) 
      betatrue <- c(betatruefix, thetagene)
      X.tot <- cbind(X.fixed, X)
    }
  
  ##This can be taken from the cumulative normal dist to find our initial for prior p_b(ij)#
  trueprobs <- X.tot%*%betatrue
  
  pA_True = c(pnorm(trueprobs, 0,1))
  #SNPs probability up here#
  rm(X.tot)
  N0=N.Control;
  N1=N.Case; ##sample size for each disease status
  mats = list()
  Indicators = rep(1,length(Gene))
  
  nonLocSamp <- MHsamps(10000, v = V, tau = tau, jit= .7, k = k)
  thesamps <- nonLocSamp$samps #acceptance rate = 0.417
  
  for(i in 1:length(Gene)){ #for each gene ID (bigger list)#
    P = pA_True[i] #locate their "true linear model" probability per SNP (had to alter this)
    Indicator = rbinom(1,1,P) #for unique Gene ID i, with probability true prob, give indicators
    if(Indicator==0){ #if Gene with pA_true is 1 assign...
      p01 = rbeta(1,1,1) ###generate same parameter for groups 0,1 from uniform dist.
      #still under Gene i, SNP j, simulate binom(n0, p01)
      num0 = rbinom(1,N0,p01) #this is assuming that the proportions are the same... null hypothesis#
      num1 = rbinom(1,N1,p01)
    }
    if(Indicator==1){ #altered the 0 and 1 in the Indicator lines because the lines below this
      #indicate that we have differing probabilities for the cases and controls in regards 
      #to their minor allele counts!!! This is the alternative hypothesis. 
      alpha0 = rnorm(1)
      alpha1 = sample(thesamps, 1); 
      p0 = pnorm(alpha0-alpha1); p1 = pnorm(alpha0 + alpha1)
      ###independently generate paramater groups for 0,1
      num0 = rbinom(1,N0,p0) #with probability p0, simulate # of controls with minor allele
      num1 = rbinom(1,N1,p1)} #with probability p1, simulate # of cases with minor allele
    Indicators[i] = Indicator
    mats[[i]] = matrix(nrow=2,ncol=2,c(N0-num0,num0,N1-num1,num1)) # <- changed
    #mats[[i]] = matrix(nrow=2,ncol=2,c(num0,N0-num0,num1,N1-num1)) <- original
  }
  # swapped the columns to match the input of the model, date: May 13th 2020
  # Thanks to Y. Guo for catching this. 
  n=matrix(nrow=length(Gene),ncol = 2)
  n0=matrix(nrow=length(Gene),ncol = 2)
  n1=matrix(nrow=length(Gene),ncol = 2)
  for(i in 1:length(Gene)){
    n[i,]=rowSums(mats[[i]])
    n0[i,] = mats[[i]][,1]
    n1[i,]=mats[[i]][,2]
  }
  
  return(list("nu" = n0, "na" = n1, "N" = n,
              "Indicators" = Indicators, "TrueBeta" = betatrue,
              "MarkerGene" = Gene, "X.Dat" = X.fixed))
}
