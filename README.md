
<!-- README.md is generated from README.Rmd. Please edit that file -->

# emBayes

This package incorporates our recently developed spike-and-slab quantile LASSO procedures to conduct Bayesian robust variable selection and estimation through the EM algorithm. The core module of this package is developed in C++. 

## How to install

  

## Example

    library(interep)
    data("dat")
    ## Load the environment factors, lipid factors and the response
    e=dat$e
    g=dat$z
    y=dat$y
    ## Initial value for the coefficient vector
    beta0=dat$coef
    ## True nonzero coefficients
    index=dat$index
    b = interep(e, g, y,beta0,corre="e",pmethod="mixed",lam1=dat$lam1, lam2=dat$lam2,maxits=30)
    ## Cut off the noise
    b[abs(b)<0.05]=0
    ## Compute TP and FP
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)



## Methods

This package provides implementation for methods proposed in

  
