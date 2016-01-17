\name{DictDensityFit}
\alias{DictDensityFit}
\title{
Fit Distribution from Kernels 
}
\description{
Given a vector of measurements for N samples (e.g, for a given methylation sites), and a collection of Gaussian kernels truncated between 0 and 1, estimate a weighted combination of kernels to approximate the distribution of measurements. 
}
\usage{
DictDensityFit(X,mu,Sigma,Concentration = 0.1,NumDraws = 1000)
}
\arguments{
  \item{X}{
a vector of measurements for N samples (e.g, for a given methylation sites).  The entries must be continuous between 0 and 1.
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
Returns a vector giving the weight of each kernel, averaged over the MCMC draws.
}
\references{
Lock, E. F. & Dunson, D. B. (2015). Shared kernel Bayesian screening. Biometrika, 102 (4): 829-842.
}
\author{
Eric F. Lock
}


