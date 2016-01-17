\name{DictTestFixedP}
\alias{DictTestFixedP}
\title{
Fit Two-Class Model with Fixed Prior Probability using Shared Kernels
}
\description{
Estimates posterior probability of a difference between two sample groups for each genomic variable (e.g., methylation site), using shared kernels and a given prior probability of equality.  Returns the results for each MCMC iteration, including the kernel weights.     
}
\usage{
DictFitTestFixedP(X,Class,mu,Sigma,pA=0.5,Concentration = 0.5,NumDraws = 1000)
}
\arguments{
  \item{X}{
A matrix in which rows represent variables (e.g., methylation sites) and columns represent samples.  The entries of the matrix must be continuous between 0 and 1.
}
 \item{Class}{
A vector giving a class label for each sample.
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
  \item{NumDraws}{
Number of MCMC draws for posterior inference. 
}
}
\value{
%%  ~Describe the value returned
Returns an object with the following values:
%%  If it is a LIST, use
\item{postDraws}{M X NumDraws array, giving the posterior probability of no difference for each marker, for each MCMC draw.}
\item{tao0Draws}{M X K X NumDraws array, giving the K kernel weights for each marker, for each MCMC draw, for class 0.}
\item{tao1Draws}{M X K X NumDraws array, giving the K kernel weights for each marker, for each MCMC draw, for class 1.}
}
\references{
Lock, E. F. & Dunson, D. B. (2016). Bayesian genome- and epigenome-wide association studies with gene level dependence. Preprint.

Lock, E. F. & Dunson, D. B. (2015). Shared kernel Bayesian screening. Biometrika, 102 (4): 829-842.
}
\author{
Eric F. Lock
}


