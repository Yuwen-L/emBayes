#' @useDynLib SSEM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' This function performs BLSS.
#' @importFrom stats coefficients lm runif
#' @importFrom regnet cv.regnet regnet
#' @param y response variable
#' @param clin clinical factors
#' @param X genetic variables
#' @param s0 spike scale
#' @param s1 slab scale
#' @param error cutoff value
#' @param maxiter max iteration number
#' @export


BLSSEM <- function(y,clin,X,s0,s1,error=0.01,maxiter=2000){
  
  p     <- ncol(X)
  C <- cbind(1,clin)
  q <- ncol(C)
  n     <- nrow(X)
  
  #initial
  sigma2=5
  theta=runif(1,0,1)
  lcd=seq(0.01,0.2,length.out=10)
  rlambda1 <- cv.regnet(X,y,response="continuous",penalty="lasso",lamb.1=lcd,folds=5,robust=TRUE)
lam <- as.numeric(rlambda1[1])*1
regbeta <- regnet(X,y,response="continuous",penalty="lasso",lamb.1=lam,robust=TRUE)
beta <- as.matrix(regbeta$coeff)[-1]
y0 <- y-X%*%beta
regalpha <- lm(y0 ~ clin)
alpha <- as.numeric(coefficients(regalpha))
  
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
  
  #y0 <- y-X%*%beta
  #reg <- lm(y0 ~ clin)
  intercept <- alpha[1]
  clin.coe <- alpha[-1]
  #intercept
    #X0 <- rep(1,n)
    #y0 <- y-X%*%beta
    #reg   <- lm(y0 ~ X0)
    #alpha <- as.numeric(coefficients(reg)[1])
  
  return(list(alpha=alpha,intercept=intercept,clin.coe=clin.coe,beta=beta,sigma2=sigma2,theta=theta,Pgamma=Pgamma,iter=iter,ll=ll,diff=diff,fl=ll[length(ll)]))
}