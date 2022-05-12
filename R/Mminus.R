#' Calculate the M minus matrix.
#'
#' @export
#' @description \code{Mminus} is used to calculate the M minus matrix defined in our paper.
#'
#' @param  gM \eqn{N(N-1)/2-by-q} matrix. The g-Matrix of Y.
#' It can be obtained by function \code{gMatrix}.
#' @param Xdiff \eqn{N(N-1)/2-by-p} matrix. The matrix of the pairwise difference of X.
#' Each row represents one observation. It can be obtained by function \code{Xdiff}.
#'
#' @return The M minus matrix.
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
#' gM_Y <- gMatrix(Y)
#' g_Y <- gM_Y$gM
#' Mminus_matrix <- Mminus(g_Y,X_diff)


Mminus <- function(gM,Xdiff) {
  # Author: Chuanping Yu
  # Date: Apr 17, 2018


  p <- ncol(Xdiff)

  ind_neg <- gM<0
  gMminus <- as.matrix(gM[ind_neg,])
  Mminus0 <- -gMminus[,rep(1,p)]*Xdiff[ind_neg,]
  ind_nonzero <- rowSums(Mminus0==0)==0
  Mminus <- Mminus0[ind_nonzero,]

  return(Mminus)

}
