#' @useDynLib SSEM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' This function performs BQLSS.
#' @importFrom stats coefficients lm runif
#' @importFrom regnet cv.regnet regnet
#' @param y response variable
#' @param clin clinical factors
#' @param X genetic variables
#' @param quant quantile
#' @param s0 spike scale
#' @param s1 slab scale
#' @param error cutoff value
#' @param maxiter max iteration number
#' @export

BQLSSEM <- function(y,clin,X,quant,s0,s1,error=0.01,maxiter=2000){
  
  p <- ncol(X)
  C <- cbind(1,clin)
  q <- ncol(C)
  n <- nrow(X)
  ep22 <- (2/(quant*(1-quant)))
  ep1 <- (1-2*quant)/(quant*(1-quant))
  
  #initial
  sigma=5
  theta=runif(1,0,1)
  lcd=seq(0.05,0.15,length.out=10)
  rlambda1 <- cv.regnet(X,y,response="continuous",penalty="lasso",lamb.1=lcd,folds=5,robust=TRUE)
lam <- as.numeric(rlambda1[1])*1
regbeta <- regnet(X,y,response="continuous",penalty="lasso",lamb.1=lam,robust=TRUE)
beta <- as.matrix(regbeta$coeff)[-1]
y0 <- y-X%*%beta
regalpha <- lm(y0 ~ clin)
alpha <- as.numeric(coefficients(regalpha))

  #lambda0 <- cv.glmnet(X,y)$lambda.min
  #lambda0 <- lambda0*2
  #fit0<-glmnet(X,y,lambda=lambda0)
  #beta <- coef(fit0)[-1]
  #alpha <- coef(fit0)[1]
  
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
  
  #y0 <- y-X%*%beta
  #reg <- lm(y0 ~ clin)
  intercept <- alpha[1]
  clin.coe <- alpha[-1]
  #intercept
    #X0 <- rep(1,n)
    #y0 <- y-X%*%beta
    #reg   <- lm(y0 ~ X0)
    #alpha <- as.numeric(coefficients(reg)[1])

  return(list(alpha=alpha,intercept=intercept,clin.coe=clin.coe,beta=beta,sigma=sigma,theta=theta,Pgamma=Pgamma,iter=iter,ll=ll,diff=diff,fl=ll[length(ll)]))
}