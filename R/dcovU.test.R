#' Permutation bootstrap methods for independence test of unbiased distance covariance.
#'
#' @export
#' @importFrom energy dcovU
#' @importFrom MASS Null
#' @importFrom stats dist quantile
#' @description \code{dcovU.test} is used to do perform independence test via permutation bootstrap methods, which perform permutations for empirical distribution function.
#' The distance measure is set as Euclideam measure.
#'
#' @param X \eqn{N-by-p} matrix. The sample dataset of X. Each row represents one observation.
#' @param Y \eqn{N-by-q} matrix. The sample dataset of Y. Each row represents one observation.
#' @param R Number of replicates for independence test. Default is set to 1000.
#' @param alpha The confidence level. The default value is 0.05.
#'
#' @return A list containing the following components will be returned:
#' \describe{
#' \item{threshold}{The threshold dcovU value of given confidence interval;}
#' \item{estimate}{The estimated unbiased distance covariance square;}
#' }
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
#' \donttest{
#' dcovU.test(X,Y)
#' }


dcovU.test <- function(X,Y,R=1000,alpha=0.05) {
  # Author: Yijin Ni
  # Date: Apr 26, 2022
  
  
  # The following functions in other packages are going to be used.
  # energy::dcovU()
  # MASS:Null()
  
  method <- "Specify the number of relicates R (R > 0) ofr an independence test"
  if(!is.null(R)) {
    R <- floor(R)
    if(R < 1)
      R <- 0
    if(R > 0)
      method <- "dCov independence test (permutation test)"
  }
  else {
    R <- 0
  }
  
  # Unbiased dcov^2
  p <- ncol(X)
  q <- ncol(Y)
  Num <- nrow(X)
  
  X_0 <- X
  Y_0 <- Y
  estimate <- dcovU(X_0, Y_0)  # Unbiased dcov^2.
  
  if (R > 0)
    reps <- rep(0, R)
  
  for (iboot in 1:R) {
    idx <- sample(Num)  # Random permutation index.
    Y_perm <- Y_0[idx, ]
    reps[iboot] = dcovU(X_0, Y_perm)
  }
  threshold = quantile(reps, 1-alpha)

  return(list(threshold = threshold, estimate = estimate))
  
}
