#' Calculate the M plus matrix.
#' 
#' @export
#' @description \code{Mplus} is used to calculate the M plus matrix defined in our paper.
#'
#' @param  gM \eqn{N(N-1)/2-by-q} matrix. The g-Matrix of Y.
#' It can be obtained by function \code{gMatrix}.
#' @param Xdiff \eqn{N(N-1)/2-by-p} matrix. The matrix of the pairwise difference of X.
#' Each row represents one observation. It can be obtained by function \code{Xdiff}.
#'
#' @return The M plus matrix.
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
#' Mplus_matrix <- Mplus(g_Y,X_diff)


Mplus <- function(gM,Xdiff) {
  # Author: Chuanping Yu
  # Date: Apr 17, 2018


  p <- ncol(Xdiff)

  ind_pos <- gM>0
  gMplus <- as.matrix(gM[ind_pos,])
  Mplus0 <- gMplus[,rep(1,p)]*Xdiff[ind_pos,]
  ind_nonzero <- rowSums(Mplus0==0)==0
  Mplus <- Mplus0[ind_nonzero,]

  return(Mplus)

}
