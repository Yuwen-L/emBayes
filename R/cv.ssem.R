#' @importFrom Rcpp sourceCpp
NULL

#' k-folds cross-validation for SSEM
#' 
#' This function performs cross-validation.
#' @param y response variable
#' @param clin clinical factors
#' @param X genetic variables
#' @param quant quantile
#' @param t0 spike scale range
#' @param t1 slab scale range
#' @param k fold number
#' @param func choose function "BLSS" or "BQLSS"
#' @param error cutoff value
#' @param maxiter max iteration number
#' @details
#' When performing cross-validation for SSEM, function cv.SSEM returns two sets of optimal tuning parameters and their corresponding cross-validation error matrices. 
#' The spike scale parameter \eqn{rs0} and the slab scale parameter \eqn{rs1} are obtained based on the mean absolute error (MAE). 
#' The spike scale parameter \eqn{nrs0} and the slab scale parameter \eqn{nrs1} are obtained based on the mean squared error (MSE). 
#' Corresponding error matrices \eqn{rCV} and \eqn{nrCV} can also be obtained from the output.
#' @return a list with components:
#' \item{rs0}{the optimal spike scale under MAE.}
#' \item{rs1}{the optimal slab scale under MAE.}
#' \item{nrs0}{the optimal slab scale under MSE.}
#' \item{nrs1}{the optimal slab scale under MSE.}
#' \item{rCV}{cross-validation error matrix under MAE.}
#' \item{nrCV}{cross-validation error matrix under MSE.}
#' 
#' @export

cv.ssem <- function(y,clin,X,quant,t0,t1,k,func,error=0.01,maxiter=2000){
  l1 <- length(t0)
  l2 <- length(t1)
  n <- nrow(X)
  C <- cbind(1,clin)
  
  s <- sample(1:n,n,replace=FALSE)
  folds <- cut(s,breaks=k,labels=FALSE)
  
  CV1 <- matrix(0,l1,l2)
  CV2 <- matrix(0,l1,l2)
  for(i in 1:l1){
      for(j in 1:l2){
        rob <- rep(0,k)
        nrob <- rep(0,k)
        for(f in 1:k){
          k.sub <- which(folds==f)
          EM <- SSEM(y[-k.sub], clin[-k.sub,], X[-k.sub,],quant,t0[i],t1[j],func,error,maxiter)
          beta_k <- EM$beta
          alpha_k <- EM$alpha
          rob[f]<-sum(abs(y[k.sub]-X[k.sub,]%*%beta_k-C[k.sub,]%*%alpha_k))
          nrob[f]<-sum((y[k.sub]-X[k.sub,]%*%beta_k-C[k.sub,]%*%alpha_k)^2)
        }
        CV1[i,j] <- sum(rob)/n
        CV2[i,j] <- sum(nrob)/n
      }
  }

  indices1 <- which(CV1 == min(CV1), arr.ind=TRUE)
  indices2 <- which(CV2 == min(CV2), arr.ind=TRUE)
  return(list("rs0"=t0[indices1[1]],"rs1"=t1[indices1[2]],"rCV"=CV1,"nrs0"=t0[indices2[1]],"nrs1"=t1[indices2[2]],"nrCV"=CV2))
}