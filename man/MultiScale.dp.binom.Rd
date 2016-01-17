\name{MultiScale.dp.binom}
\alias{MultiScale.dp.binom}
\title{
Fit Two-Class Model with Gene-Level Probabilities for Categorical Data
}
\description{
Estimates posterior probability of a difference between two sample groups for each genomic variable (e.g., SNPs) using multinomial likelihood and the prospective Bayes factor, with gene-dependent prior probabilities of equality.  
}
\usage{
MultiScale.DP.binom(n,n0,n1,UniqGene,Gene,Concentration=c(1,1),alpha=1,K=10,NumDraws=1000)
}
\arguments{
  \item{n}{
D X Q matrix, where D is the number of variables (e.g., SNPs) and Q is the number of categories.  Each row gives total number of cases in each category for a given variable.  
}
 \item{n0}{
D X Q matrix, where each row gives number of cases in each category for group 0.
}
 \item{n1}{
D X Q matrix, where each row gives number of cases in each category for group 1.
}
  \item{UniqGene}{
Vector giving the gene labels.  
}
  \item{Gene}{
Vector giving the gene label for each variable (length must be equal to the number of rows in X).
}
  \item{Concentration}{
Dirichlet hyperparameter for the multinomial probablities (defaults to uniform)  
}
  \item{alpha}{
Dirichlet process concentration parameter for the gene-level probabilities 
}
  \item{K}{
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
\item{posts}{Vector giving the posterior probability of association for each variable}
}
\references{
Lock, E. F. & Dunson, D. B. (2016). Bayesian genome- and epigenome-wide association studies with gene level dependence. Preprint.

Balding, D. J. (2006) A tutorial on statistical methods for population association studies. Nature Reviews Genetics, 7(10), 781–79.
}
\author{
Eric F. Lock
}


