#' @importFrom Rcpp sourceCpp
NULL

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
  return(list("s0"=t0[indices1[1]],"s1"=t1[indices1[2]],"CV1"=CV1,"s02"=t0[indices2[1]],"s12"=t1[indices2[2]],"CV2"=CV2))
}