\name{MultiScale.dp.kernel}
\alias{MultiScale.dp.kernel}
\title{
Fit Two-Class Model with Gene-Level Probabilities using Shared Kernels
}
\description{
Estimates posterior probability of a difference between two sample groups for each genomic variable (e.g., methylation site), using shared kernels and gene-dependent prior probabilities of equality.  
}
\usage{
MultiScale.dp.kernel(X,Class,UniqGene,Gene,mu,Sigma,Concentration=0.5,alpha=1,KK=10,NumDraws=1000)
}
\arguments{
  \item{X}{
A matrix in which rows represent variables (e.g., methylation sites) and columns represent samples.  The entries of the matrix must be continuous between 0 and 1.
}
 \item{Class}{
A vector giving a class label for each sample.
}
  \item{UniqGene}{
Vector giving the gene labels.  
}
  \item{Gene}{
Vector giving the gene label for each variable (length must be equal to the number of rows in X).
}
  \item{mu}{
Vector of kernel means.
}
  \item{Sigma}{
Vector of kernel standard deviations
}
  \item{Concentration}{
Dirichlet hyperparameter for kernel weights  
}
  \item{alpha}{
Dirichlet process concentration parameter for the gene-level probabilities 
}
  \item{KK}{
Stick-breaking threshold for Dirichlet process
}
  \item{NumDraws}{
Number of MCMC draws for posterior inference. 
}
}
\value{
%%  ~Describe the value returned
Returns an object with the following values, averaged over the MCMC iterations:
%%  If it is a LIST, use
\item{pG}{Vector giving the gene-level prior for association for each gene}
\item{posts}{Vector giving the posterior probability of association for each marker}
}
\references{
Lock, E. F. & Dunson, D. B. (2016). Bayesian genome- and epigenome-wide association studies with gene level dependence. Preprint.

Lock, E. F. & Dunson, D. B. (2015). Shared kernel Bayesian screening. Biometrika, 102 (4): 829-842.
}
\author{
Eric F. Lock
}


