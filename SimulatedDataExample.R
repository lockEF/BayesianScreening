source("NonLocalDataGeneration.R")

#Simulate the Data#
datafor <- DataGenNonLoc.Spec(q=10, 500, 30, 3, 2, N.Control = 100, N.Case = 100, setting = 3, k = 1, V = 4, tau = 0.75)

#with the given parameter set in this function, we are misspecifying V


#Collect: how many MCMC samples one would like to have at the end of the analysis
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
# hyp: vector of prior parameters for concentration parameter, for Gamma(a,b) prior: c(a,b)
#parameters: vector with the input in this order of the Non-Local parameters c(v, tau, k)
Trialc <- DPMasterChibNonLoc(Collect = 2000, Gene = datafor$MarkerGene, Burnin = 1000, X.Dat= datafor$X.Dat, 
                             nu = datafor$nu, na = datafor$na, N = datafor$N, H = 20, colfx = 5, 
                             thin = 5, hyp = c(1.1, 1.1), parameters = c(3, .75, 1))

summary(Trialc$Post.Mark)
hist(Trialc$Post.Mark, xlab = "Posterior Probability of Association", 
     main = "Posterior Probability of Association Histogram")

#Summary measures of Gene Effects vs. True Gene Effects#
summary(rowMeans(Trialc$Beta_Mat)[-c(1:5)]); summary(datafor$TrueBeta[-c(1:5)])

#Comparison between Truth and Cut-off Marker Posterior Probability 
table(datafor$Indicators); table(Trialc$Post.Mark >= 0.5)

#Fixed Effects Estimation as compared to Truth
rowMeans(Trialc$Beta_Mat)[1:5]; datafor$TrueBeta[1:5]

cbind(t(cbind(sapply(1:5, function(x) quantile(Trialc$Beta_Mat[x,], probs = c(.025,.5, .975))))),datafor$TrueBeta[1:5])
#last column are the true estimates
#95% credible intervals for estimates