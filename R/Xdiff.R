#' The matrix of the pairwise diferences of X.
#'
#' 
#' @export
#' @description \code{Xdiff} is used to get the matrix of all the pairwise diferences of X.
#' Each row represents one difference.
#'
#' @param  X \eqn{N-by-p} matrix. The sample dataset of X. Each row represents one observation.
#'
#' @return \eqn{N(N-1)/2-by-p} matrix. The matrix of the pairwise diferences of X.
#'
#' @examples
#' library(MASS)
#' Num = 40
#' p = 3
#' q = 2
#' mu1 = rep(0,p)
#' sigma1 = diag(p)
#' mu2 = rep(0,q)
#' sigma2 = diag(q)
#' X = mvrnorm(Num,mu1,sigma1)
#' Y = mvrnorm(Num,mu2,sigma2)
#' SUM_X = rowSums(X)
#' Y[,1] = SUM_X+0.01*Y[,1]
#' Y[,2] = (SUM_X)^2+0.01*Y[,2]
#' X_diff <- Xdiff(X)
#' Y_diff <- Xdiff(Y)


Xdiff <- function(X) {
  # Author: Chuanping Yu
  # Date: Apr 17, 2018


  N <- nrow(X)
  p <- ncol(X)
  Xdiff <- matrix(0,N*(N-1)/2,p)

  ind <- 0
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      ind <- ind+1
      Xdiff[ind,] <- X[i,]-X[j,]
    }
  }

  return(Xdiff)

}
