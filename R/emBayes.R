#' @useDynLib emBayes, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' fit a model with given tuning parameters
#' 
#' This function performs penalized variable selection based on spike-and-slab quantile LASSO (BQLSS) or Bayesian linear regression with spike-and-slab LASSO prior (BLSS).
#' Typical usage is to first obtain the optimal spike scale and slab scale using cross-validation, then specify them in the emBayes function.
#' @importFrom stats coefficients lm runif coef
#' @importFrom glmnet cv.glmnet glmnet
#' @param y a vector of response variable.
#' @param clin a matrix of clinical factors. It has default value NULL.
#' @param X a matrix of genetic factors.
#' @param quant value of quantile.
#' @param s0 value of the spike scale \eqn{s_{0}}.
#' @param s1 value of the slab scale \eqn{s_{1}}.
#' @param func methods to perform variable selection. Two choices are available: "BLSS" and "BQLSS".
#' @param error cutoff value for determining convergence. The algorithm reaches convergence if the difference in the expected log-likelihood of two iterations is less than the value of error. The default value is 0.01.
#' @param maxiter the maximum number of iterations that is used in the estimation algorithm. The default value is 200.
#' @details 
#' The current version of emBayes supports two types of methods: "BLSS" and "BQLSS".
#' \itemize{
#' \item \strong{BLSS:} Bayesian LASSO with spike-and-slab prior fits a Bayesian linear regression through the EM algorithm.  
#' \item \strong{BQLSS:} spike-and-slab quantile LASSO fits a Bayesian quantile regression (based on asymmetric Laplace distribution) through the EM algorithm.
#' }
#' Users can choose the desired method by setting {func="BLSS" or "BQLSS"}.
#' @return a list with components:
#' \item{alpha}{a vector containing the estimated intercept and clinical coefficients.}
#' \item{intercept}{value of the estimated intercept.}
#' \item{clin.coe}{a vector of estimated clinical coefficients.}
#' \item{beta}{a vector of estimated beta coefficients.}
#' \item{sigma}{value of estimated asymmetric Laplace distribution scale parameter \eqn{\sigma}.}
#' \item{theta}{value of estimated probability parameter \eqn{\theta}.}
#' \item{iter}{value of number of iterations.}
#' \item{ll}{a vector of expectation of likelihood at each iteration.}
#' 
#' @examples 
#' data(genes)
#' ##load the clinical factors, genetic factors, response and quantile data
#' clin=genes$clin
#' X=genes$X
#' y=genes$y
#' quant=genes$quant
#' ##generate tuning vectors of desired range 
#' t0 <- seq(0.01,0.015,length.out=2)
#' t1 <- seq(0.1,0.5,length.out=2)
#' ##perform cross-validation and obtain tuning parameters based on check loss
#' CV <- cv.emBayes(y,clin,X,quant,t0,t1,k=5,func="BQLSS",error=0.01,maxiter=2000)
#' s0 <- CV$CL.s0
#' s1 <- CV$CL.s1
#' ##perform BQLSS under optimal tuning and calculate value of TP and FP for selecting beta 
#' EM <- emBayes(y,clin,X,quant,s0,s1,func="BQLSS",error=0.01,maxiter=2000)
#' fit <- EM$beta
#' coef <- genes$coef
#' tp <- sum(fit[coef!=0]!=0)
#' fp <- sum(fit[coef==0]!=0)
#' list(tp=tp,fp=fp)
#' 
#' @export

emBayes <- function(y,clin=NULL,X,quant,s0,s1,func,error=0.01,maxiter=200){
  
  p <- ncol(X)
  n <- nrow(X)
  inter <- rep(1,n)
  C <- cbind(inter,clin)
  q <- ncol(C)
  
  #initial
  sigma=5
  #sigma2=5
  theta=runif(1,0,1)
  #lcd=seq(0.05,0.15,length.out=10)
  #rlambda1 <- cv.regnet(X,y,response="continuous",penalty="lasso",lamb.1=lcd,folds=5,robust=TRUE)
  #lam <- as.numeric(rlambda1[1])*1
  #regbeta <- regnet(X,y,response="continuous",penalty="lasso",lamb.1=lam,robust=TRUE)
  #beta <- as.matrix(regbeta$coeff)[-1]
  lambda <- (cv.glmnet(X,y)$lambda.min)*1
  fit <- glmnet(X,y,lambda=lambda)
  beta <- coef(fit)[-1]
  
  y0 <- y-X%*%beta
  regalpha <- lm(y0 ~ clin)
  alpha <- as.numeric(coefficients(regalpha))
  
  if(func=="BQLSS"){
    ep22 <- (2/(quant*(1-quant)))
    ep1 <- (1-2*quant)/(quant*(1-quant))
    
    ll <- c()
    loglike <- logQR(y,X,C,alpha,beta,sigma,theta,s0,s1,ep1,ep22)
    lk <- loglike$logver
    vn <- loglike$vn
    vp <- loglike$vp
    Pgamma <- loglike$Pgamma
    invS <- loglike$invS
    ll <- c(ll,lk)
    
    iter <- 0
    diff <- 1
    
    while( diff > error){
      iter  <- iter+1
      EM <- EMQR(y,X,C,n,p,q,quant,alpha,beta,sigma,theta,s0,s1,Pgamma,invS,ep1,ep22,vn,vp)
      alpha <- EM$alpha
      beta <- EM$beta
      sigma <- EM$sigma
      theta <- EM$theta
      
      loglike2 <- logQR(y,X,C,alpha,beta,sigma,theta,s0,s1,ep1,ep22)
      lk2 <- loglike2$logver
      ll <- c(ll,lk2)
      diff <- abs(lk2-lk)
      vn <- loglike2$vn
      vp <- loglike2$vp
      Pgamma <- loglike2$Pgamma
      invS <- loglike2$invS
      
      lk <- lk2
    }
  }
  else{
    ll <- c()
    loglike <- logR(y,X,C,alpha,beta,sigma2,theta,s0,s1)
    lk <- loglike$logver
    ll <- c(ll,lk)
    Pgamma <- loglike$Pgamma
    invS <- loglike$invS
    
    iter <- 0
    diff <- 1
    
    while( diff > error){
      iter  <- iter+1
      EM <- EMR(y,X,C,n,p,q,alpha,beta,sigma2,theta,Pgamma,invS)
      alpha <- EM$alpha
      beta <- EM$beta
      sigma2 <- EM$sigma2
      theta <- EM$theta
      
      loglike2 <- logR(y,X,C,alpha,beta,sigma2,theta,s0,s1)
      lk2 <- loglike2$logver
      ll <- c(ll,lk2)
      diff <- abs(lk2-lk)
      Pgamma <- loglike2$Pgamma
      invS <- loglike2$invS
      
      lk <- lk2
    }
  }
  
  intercept <- alpha[1]
  clin.coe <- alpha[-1]

  return(list(alpha=alpha,intercept=intercept,clin.coe=clin.coe,beta=beta,sigma=sigma,theta=theta,iter=iter,ll=ll))
}