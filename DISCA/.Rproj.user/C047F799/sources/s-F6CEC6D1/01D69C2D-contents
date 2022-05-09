#' Calculate the \eqn{g_{ij}} values.
#'
#' @importFrom stats dist
#' @export
#' @description \code{gMatrix} is used to calulate all the \eqn{g_{ij}} values
#' for \eqn{1\le i<j\le N}, where \eqn{N} is the sample size of the data.
#'
#' @param Y \eqn{N-by-q} matrix. The sample dataset of Y. Each row represents one observation.
#' @return A list containing the following components will be returned:
#' \describe{
#' \item{gM}{\eqn{N(N-1)/2-by-1} matrix. The g-Matrix of Y.}
#' \item{sum4}{The sum of distance matrix of Y.}
#' }
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
#' gMatrix_Y <- gMatrix(Y)
#' gMatrix_X <- gMatrix(X)



gMatrix <- function(Y) {
  # Author: Chuanping Yu
  # Date: Apr 17, 2018


  N <- nrow(Y)

  dist0 <- dist(Y,p=2)  # Get distance matrix between the rows of Y.
  g0 <- as.matrix(dist0)  # Reformulate up-triangular matrix to a symmetric matrix: g0(i,j) = |Y_i-Y_j|_q
  sum2 <- matrix(-(1/(N-2))*rowSums(g0),N,1) # Negative average distance for each Y_i.
  sum2_M <- sum2[,rep(1,N)]  # The second term of the formulation of each g_ij.
  sum3 <- matrix(-(1/(N-2))*colSums(g0),1,N)
  sum3_M <- sum3[rep(1,N),]  # The third term of the formulation of each g_ij.
  sum4 <- sum(g0)/((N-1)*(N-2))
  sum4_M <- matrix(sum4,N,N)  # The fourth term of the formulation of g_ij.

  g <- g0+sum2_M+sum3_M+sum4_M
  gM <- matrix(g[lower.tri(g)],N*(N-1)/2,1)

  return(list(gM=gM,sum4=sum4*((N-1)*(N-2))))

}
