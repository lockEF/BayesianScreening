# BayesianScreening

This is an R package to compute posteriors from a gene-level hierarchical prior, as described in the manuscript "Bayesian genome- and epigenome-wide association studies with gene level dependence" [1].  Also includes functions to perform two-class testing with shared kernels for methylation array data (see [2]).  The file BinaryExample.r gives example code for simulated categorical data. The file GliomaMethylExample.r gives example code for an application to methylation screening data. 

The package can be installed in R, directly from GitHub, using the devtools library:
```
install.packages('devtools')
library(devtools)
install_github("lockEF/BayesianScreening")
```

The functions in the R script SNLPMasterFile.R can be used to perform genome-wide association studies with a gene-level prior, additional marker-level covariates, and a non-local alternative, as described in [3].  This code was developed by Adam Kaplan.  See SimulatedDataExample.R for example usage.  

[1] Lock, E. F. & Dunson, D. B. (2017). Bayesian genome- and epigenome-wide association studies with gene level dependence. Biometrics, 73 (3): 1018-1028. https://www.ncbi.nlm.nih.gov/pubmed/28083869

[2] Lock, E. F. & Dunson, D. B. (2015). Shared kernel Bayesian screening. Biometrika, 102 (4): 829-842.  https://www.ncbi.nlm.nih.gov/pubmed/27046939

[3] Adam Kaplan, Eric F. Lock, Mark Fiecas, "Bayesian GWAS with Structured and Non-Local Priors," Bioinformatics, 2019. 
