\name{EstimateTruncatedDictionary}
\alias{EstimateTruncatedDictionary}
\title{
Estimate Dictionary of Truncated Kernels
}
\description{
Given a matrix with methylation sites as rows and samples as columns, returns a collection of Gaussian kernels truncated between 0 and 1.   These kernels are used as a dictionary, where the distribution of values at each site  is represented a weighted combination of the kernels.  
}
\usage{
EstimateTruncatedDictionary(X, K=2, a0 = 0.5,b0 = 0.5,mu0 = 0.5,Concentration = 0.5, NumDraws = 1000)
}
\arguments{
  \item{X}{
A matrix in which rows represent variables (e.g., methylation sites) and columns represent samples.  The entries of the matrix must be continuous between 0 and 1.
}
  \item{K}{
The number of kernels
}
  \item{a0,b0}{
Gamma hyperperameters for the precision (inverse of the variance) of the kernels.  Defaults are a0=0.5, b0=0.5 (we recommend using these defaults if unsure).
}
  \item{mu0}{
Normal mean hyperparameter for the kernel means.  Default is mu0=0.5.   
}
  \item{Concentration}{
Initial value of the Dirichlet concentration parameter.  Default is 0.5.  Can be a vector of length K.  
}
  \item{NumDraws}{
Number of MCMC draws for posterior inference. 
}
}
\value{
%%  ~Describe the value returned
Returns an object with the following values, averaged over the MCMC iterations:
%%  If it is a LIST, use
\item{Mu}{Vector of length K giving the mean of each kernel}
\item{Sigma}{Vector of length K giving the standard devitation of each kernel}
\item{Concentration}{Vector of length K giving the concentration hyperparameter}
}
\references{
Lock, E. F. & Dunson, D. B. (2015). Shared kernel Bayesian screening. Biometrika, 102 (4): 829-842.
}
\author{
Eric F. Lock
}


