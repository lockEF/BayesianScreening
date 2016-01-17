
###Load BayesianScreening package
library(BayesianScreening)

#################Simulate data for 500 genes:
#############################################
Ngene=500
UniqGene = c(1:Ngene)
Gene = c() ##Gives the gene number for each marke
for(i in 1:Ngene){
  Nmarkers = sample(c(1:19),1)+1 ##number of markers for each gene is taken  from 2,..,20 with equal probability
  Gene=c(Gene,rep(UniqGene[i],Nmarkers))
}

pA_True = c(rep(1,400),rep(0,100)) ## 400 genes have all null markers, 100 have all alternative markers
N0=80;
N1=80; ##sample size = 80 for each group
mats = list()
Indicators = rep(1,length(Gene))
for(i in 1:length(Gene)){
	P = pA_True[UniqGene==Gene[i]]
	Indicator = rbinom(1,1,P)
	if(Indicator==1){
    p01 = rbeta(1,1,1) ###generate same parameter for groups 0,1 from uniform dist.
		num0 = rbinom(1,N0,p01)
		num1 = rbinom(1,N1,p01)}
	if(Indicator==0){
	  p0 = rbeta(1,1,1) ###independently gennerate paramater groups for 0,1
		num0 = rbinom(1,N0,p0)
		p1 = rbeta(1,1,1)
		num1 = rbinom(1,N1,p1)}
	Indicators[i] = Indicator
mats[[i]] = matrix(nrow=2,ncol=2,c(num0,N0-num0,num1,N1-num1))
}
n=matrix(nrow=length(Gene),ncol = 2)
n0=matrix(nrow=length(Gene),ncol = 2)
n1=matrix(nrow=length(Gene),ncol = 2)
for(i in 1:length(Gene)){
	n[i,]=rowSums(mats[[i]])
	n0[i,] = mats[[i]][,1]
	n1[i,]=mats[[i]][,2]
}
####each row of n gives total number of cases in each category (e.g., minor allele/ no minor allele) for a given marker
####n0 gives number of cases in each category for group 0
####n1 gives number of cases in each categoru for group 1

########################################Get testing results
Results = MultiScale.DP.binom(n,n0,n1,UniqGene,Gene,K=20,alpha=1,NumDraws=2000)
Pg = Results[[1]] ###gene probabilities
post= Results[[2]] ####marker probabilities
error = sum(abs(post-Indicators))/length(Gene) ##posterior error
