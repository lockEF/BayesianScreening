###This file contains code to reproduce the screen for differences in methylation between GBM and LGG samples, presented in the manusctipt
###"Bayesian genome- and epigenome-wide association studies with gene level dependence"

##Illumina HumanMethylation450K data were obtained from TCGA, using the "Build archive" function
##For GBM data: https://tcga-data.nci.nih.gov/tcga/tcgaCancerDetails.jsp?diseaseType=GBM&diseaseName=Glioblastoma%20multiforme
##And LGG data: https://tcga-data.nci.nih.gov/tcga/tcgaCancerDetails.jsp?diseaseType=LGG&diseaseName=Brain%20Lower%20Grade%20Glioma
##Accessed October 6th, 2015.  
##For the code that was used for data gathering, filtering, and processing see this link: https://www.dropbox.com/s/y2swhsyi5sobj7h/DataProcessing.r?dl=0 

##the full data matrix with both GBM and atrocytoma samples can be downloaded as a .rData file at this link: https://www.dropbox.com/s/lgn45klhdoxchc2/MethMat.rData?dl=0

##Load the data:
load("MethMat.rData")
cpg.info = cpg.info[1:375862,] ##remove unlabelled cpgs at end
meth = meth[1:375862,] 

#Load BayesianScreening package
library(BayesianScreening)

#############The code below is for estimating the kernels.     
#First, estimate under different number of kernels $K=2,3,...,11
Dict = list()
Dict[[1]] = 'NA'
for(k in 2:11){Dict[[k]]= EstimateTruncatedDictionary(Xsamp,K=k, NumDraws=500)}

###Find the set of dictionary densities that maximizes log-likelihood under CV
MeanCrossVal =c()
Xsamp2= meth[sample(dim(meth)[1],1000),] ##Select another set of 1000 probes
#Run in parallel (can take 30-60 min per K)
cl <- makeCluster(8)
registerDoParallel(8)
MeanCrossVal[2:11] = foreach(k=1:10, .combine='cbind') %dopar%{
  library(gtools)
  library(msm)
  lys = c()
  for(i in 1:1000){
    j = sample(dim(meth)[2],1)
    weights = DictDensityFit(Xsamp2[i,-c(j)],mu=Dict[[k]]$Mu,Sigma=Dict[[k]]$Sigma,Concentration=Dict[[k]]$Concentration,NumDraws=500)
    y = 0
    for(ind in 1:(k+1)) y = y+weights[ind]*dtnorm(Xsamp2[i,j],Dict[[k]]$Mu[ind],Dict[[k]]$Sigma[ind],0,1)
    lys[i] = log(y)
  }
  mean(lys)
}
stopCluster(cl)

MeanCrossVal = rep(0,11)
for(k in 2:11){
lys = c()
for(i in 1:1000){
  j = sample(dim(meth)[2],1)
  weights = DictDensityFit(Xsamp2[i,-c(j)],mu=Dict[[k]]$Mu,Sigma=Dict[[k]]$Sigma,Concentration=Dict[[k]]$Concentration,NumDraws=500)
  y = 0
  for(ind in 1:(k)) y = y+weights[ind]*dtnorm(Xsamp2[i,j],Dict[[k]]$Mu[ind],Dict[[k]]$Sigma[ind],0,1)
  lys[i] = log(y)
}
MeanCrossVal[k]=mean(lys)
}
MeanCrossVal ##Maximized at K=8


###Now run testing application

Gene = as.character(cpg.info$Gene_Symbol) ####Gives gene label for each CpG
X = meth[Gene!='',] ##remove intergenic CpGs
Gene = Gene[Gene!='']
UniqGene = unique(Gene)


###
##Run testing for 1000 MCMC iterations.  This can take a while (~30 seconds per draw)
Results = MultiScale.dp.kernel(X=X, Class = class, Concentration=Dict[[8]]$Conc;, mu = Dict[[8]]$Mu, Sigma = Dict[[8]]$Sigma; alpha=1; KK=20; NumDraws=1000)
##Gene probabilities:
Pg = Results[[1]]
##histogram of gene probabilities:
hist(Pg, xlim = c(0,1),breaks=40, xlab = expression('p'[g]),main ='Histogram of gene priors')

MarkerPost = Results[[2]]
###Focus on gene 'BST2'
Indices = c(1:length(Gene))[Gene=='BST2'] ###
post.bst2 = MarkerPost[Indices]
Pg.bst2 = Pg[UniqGene=='BST2']
###Plot CpG marker locations and probabilities, with distribution of selected CpGs, for BST2
m <- cbind(c(1,4),c(1,3),c(1,2)) #for multi-panel layout
layout(m)
###first plot CpG coordinates
cpg.info[Indices,]
coords = cpg.info[Indices,4]
plot(coords,post.bst2,ylim = c(0,1), xlab = 'Chr19 position',ylab = 'Posterior probability', main = 'BST2 CpGs',cex=1.6)
points(coords[9],post.bst2[9], col = 'brown',pch=19,cex=1.6)
points(coords[7],post.bst2[7], col = 'green',pch=19,cex=1.6)
points(coords[5],post.bst2[5], col = 'purple',pch=19,cex=1.6)


####Fit densities for BST2, using its estimated prior
Fits = DictTestFixedP(X[Indices,],Class=class,mu=Dict[[8]]$Mu,Sigma=Dict[[8]]$Sigma,Concentration=Dict[[8]]$Concentration,pA=Pg.bst2,NumDraws=2000)
######Plot density fits for select CpGs
site =5
Num=Indices[[site]]
Breaks =seq(0,1,length.out=10)
ahist = hist(X[Num,Class==0],xlim=c(0,1),ylim = c(0,7.5), freq = FALSE, breaks = Breaks, col = "deepskyblue1", main = expression(paste('cg11558551, pr(',H[0],'|X)<0.001')), xlab = "Methylation",col.main='purple')
bhist = hist(X[Num,Class==1],add=T,xlim = c(0,1),freq = FALSE, breaks = Breaks,col = "tomato")
overlap = bhist
for(i in 1:length(overlap$counts)){ 
  if(ahist$density[i] > 0 & bhist$density[i] > 0){
    overlap$counts[i] = min(ahist$density[i],bhist$density[i])
  } else {
    overlap$counts[i] = 0
  }
}
plot(overlap, xlim=c(0,1), ylim=c(0,10), col='darkslateblue', add=T)
Weight0 = colMeans(Fits[[3]][site,401:2000,])
Weight1 = colMeans(Fits[[4]][site,401:2000,])
xg = seq(0,1,length.out = 200)
y0 = rep(0,length(xg))
y1 = rep(0,length(xg))
for(k in 1:K){ y0 = y0+Weight0[k]*dtnorm(xg,Dict[[K]]$Mu[k],Dict[[K]]$Sigma[k],0,1)
               y1 = y1+Weight1[k]*dtnorm(xg,Dict[[K]]$Mu[k],Dict[[K]]$Sigma[k],0,1)}
lines(xg,y0,col='blue', lwd = '3')
lines(xg,y1,col='red', lwd = '3')

site =7
Num=Indices[[site]]
Breaks =seq(0,1,length.out=10)
ahist = hist(X[Num,Class==0],xlim=c(0,1),ylim = c(0,6), freq = FALSE, breaks = Breaks, col = "deepskyblue1", main = expression(paste('cg16363586, pr(',H[0],'|X)=0.07')), xlab = "Methylation",col.main='green')
bhist = hist(X[Num,Class==1],add=T,xlim = c(0,1),freq = FALSE, breaks = Breaks,col = "tomato")
overlap = bhist
for(i in 1:length(overlap$counts)){ 
  if(ahist$density[i] > 0 & bhist$density[i] > 0){
    overlap$counts[i] = min(ahist$density[i],bhist$density[i])
  } else {
    overlap$counts[i] = 0
  }
}
plot(overlap, xlim=c(0,1), ylim=c(0,10), col='darkslateblue', add=T)
Weight0 = colMeans(Fits[[3]][site,401:2000,])
Weight1 = colMeans(Fits[[4]][site,401:2000,])
xg = seq(0,1,length.out = 200)
y0 = rep(0,length(xg))
y1 = rep(0,length(xg))
for(k in 1:K){ y0 = y0+Weight0[k]*dtnorm(xg,Dict[[K]]$Mu[k],Dict[[K]]$Sigma[k],0,1)
               y1 = y1+Weight1[k]*dtnorm(xg,Dict[[K]]$Mu[k],Dict[[K]]$Sigma[k],0,1)}
lines(xg,y0,col='blue', lwd = '3')
lines(xg,y1,col='red', lwd = '3')


site =9
Num=Indices[[site]]
Breaks =seq(0,1,length.out=10)
ahist = hist(X[Num,Class==0],xlim=c(0,1),ylim = c(0,4), freq = FALSE, breaks = Breaks, col = "deepskyblue1", main = expression(paste('cg22282590, pr(',H[0],'|X)=0.78')), xlab = "Methylation", col.main='brown')
bhist = hist(X[Num,Class==1],add=T,xlim = c(0,1),freq = FALSE, breaks = Breaks,col = "tomato")
overlap = bhist
for(i in 1:length(overlap$counts)){ 
  if(ahist$density[i] > 0 & bhist$density[i] > 0){
    overlap$counts[i] = min(ahist$density[i],bhist$density[i])
  } else {
    overlap$counts[i] = 0
  }
}
plot(overlap, xlim=c(0,1), ylim=c(0,10), col='darkslateblue', add=T)
Weight0 = colMeans(Fits[[3]][site,401:2000,])
Weight1 = colMeans(Fits[[4]][site,401:2000,])
xg = seq(0,1,length.out = 200)
y0 = rep(0,length(xg))
y1 = rep(0,length(xg))
for(k in 1:K){ y0 = y0+Weight0[k]*dtnorm(xg,Dict[[K]]$Mu[k],Dict[[K]]$Sigma[k],0,1)
               y1 = y1+Weight1[k]*dtnorm(xg,Dict[[K]]$Mu[k],Dict[[K]]$Sigma[k],0,1)}
lines(xg,y0,col='blue', lwd = '3')
lines(xg,y1,col='red', lwd = '3')


#######Plot overall distributions for significantly different CpGs

Means0.sig = rowMeans(X[MarkerPost<0.01,class==0])
Means1.sig = rowMeans(X[MarkerPost<0.01,class==1])
sds0.sig = apply(X[MarkerPost<0.01,Class==0],1,sd)
sds1.sig = apply(X[MarkerPost<0.01,Class==1],1,sd)

par(mfrow=c(1,2))
Breaks =seq(0,1,length.out=50)
ahist = hist(Means0.sig,xlim=c(0,1),ylim = c(0,4), freq = FALSE, breaks = Breaks, col = "deepskyblue1", main = "Means", xlab = "Mean")
bhist = hist(Means1.sig,add=T,xlim = c(0,1),freq = FALSE, breaks = Breaks,col = "tomato")
overlap = bhist
for(i in 1:length(overlap$counts)){ 
  if(ahist$density[i] > 0 & bhist$density[i] > 0){
    overlap$counts[i] = min(ahist$density[i],bhist$density[i])
  } else {
    overlap$counts[i] = 0
  }
}
plot(overlap, xlim=c(0,1), ylim=c(0,10), col='darkslateblue', add=T)

Breaks =seq(0,0.45,length.out=50)
ahist = hist(sds0.sig,xlim=c(0,0.4),ylim = c(0,8), freq = FALSE, breaks = Breaks, col = "deepskyblue1", main = "Standard deviations", xlab = "Standard deviation", col.main='black')
bhist = hist(sds1.sig,add=T,xlim = c(0,0.4),freq = FALSE, breaks = Breaks,col = "tomato")
overlap = bhist
for(i in 1:length(overlap$counts)){ 
  if(ahist$density[i] > 0 & bhist$density[i] > 0){
    overlap$counts[i] = min(ahist$density[i],bhist$density[i])
  } else {
    overlap$counts[i] = 0
  }
}
plot(overlap, xlim=c(0,0.4), ylim=c(0,10), col='darkslateblue', add=T)







            