
<!-- README.md is generated from README.Rmd. Please edit that file -->

# emBayes

This package incorporates our recently developed spike-and-slab quantile LASSO procedures to conduct Bayesian robust variable selection and estimation through the EM algorithm. The core module of this package is developed in C++. 

## How to install

  

## Example

    library(emBayes)
    data(genes)
    ## Load the clinical factors, genetic factors, response and quantile data
    clin=genes$clin
    X=genes$X
    y=genes$y
    quant=genes$quant
    ## Generate tuning vectors of desired range
    t0 <- seq(0.01,0.015,length.out=2)
    t1 <- seq(0.1,0.5,length.out=2)
    ## Perform cross-validation and obtain tuning parameters based on check loss
    CV <- cv.emBayes(y,clin,X,quant,t0,t1,k=5,func="BQLSS",error=0.01,maxiter=2000)
    s0 <- CV$CL.s0
    s1 <- CV$CL.s1
    ## Perform BQLSS under optimal tuning and calculate value of TP and FP for selecting beta
    EM <- emBayes(y,clin,X,quant,s0,s1,func="BQLSS",error=0.01,maxiter=2000)
    fit <- EM$beta
    coef <- genes$coef
    tp <- sum(fit[coef!=0]!=0)
    fp <- sum(fit[coef==0]!=0)
    list(tp=tp,fp=fp)



## Methods

This package provides implementation for methods proposed in

  
