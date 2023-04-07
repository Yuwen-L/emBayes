#' @useDynLib SSEM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' This function performs BQLSS.
#' @importFrom stats coefficients lm runif coef
#' @importFrom glmnet cv.glmnet glmnet
#' @param y response variable
#' @param clin clinical factors
#' @param X genetic variables
#' @param quant quantile
#' @param s0 spike scale
#' @param s1 slab scale
#' @param func choose function "BLSS" or "BQLSS"
#' @param error cutoff value
#' @param maxiter max iteration number
#' @export

SSEM <- function(y,clin,X,quant,s0,s1,func,error=0.01,maxiter=2000){
  
  p <- ncol(X)
  C <- cbind(1,clin)
  q <- ncol(C)
  n <- nrow(X)
  
  #initial
  sigma=5
  sigma2=5
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

  return(list(alpha=alpha,intercept=intercept,clin.coe=clin.coe,beta=beta,sigma=sigma,theta=theta,Pgamma=Pgamma,iter=iter,ll=ll,diff=diff))
}