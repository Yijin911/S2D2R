#' Get the dimension reduction subspaces for both X and Y.
#'
#' @export
#' @importFrom energy dcovU
#' @importFrom MASS Null
#' @importFrom stats dist quantile
#' @description \code{DISCA} is used to do dimension reduction via distance covariance by projecting X and Y into subspaces \eqn{W_1,} and \eqn{W_2}, respectively.
#' The dimension of X and Y, \eqn{W_1,} and \eqn{W_2} can be different.
#'
#' @param X \eqn{N-by-p} matrix. The sample dataset of X. Each row represents one observation.
#' @param Y \eqn{N-by-q} matrix. The sample dataset of Y. Each row represents one observation.
#' @param e The tolerance of two values calculated by adjacent two steps by Difference of Convex Algorithm (DCA).
#' If the two values is smaller than e, we terminate the loop of DCA.
#' The default value is 10^-3.
#' @param xi One of the augmented Lagrangian parameters.
#' It should be nonnegative as well as bigger than psi.
#' The default value is 200.
#' @param psi The other augmented Lagrangian parameters.
#' The default value is 10.
#' @param rho The augmented Lagrangian parameter in subproblem.
#' The default value is 1000.
#' @param omega The increment of xi. It should be bigger than 1.
#' The default value is 10.
#' @param abs.epsilon The absolute tolerance in ADMM.
#' The default value is 10^-3.
#' @param rel.epsilon The relative tolerance in ADMM.
#' The default value is 10^-3.
#' @param alpha The confidence level. The default value is 0.05.
#' @param R Number of replicates in permutation bootstrap method for independence test.
#'
#' @return A list containing the following components will be returned:
#' \describe{
#' \item{dW1}{The dimension of subspace W1;}
#' \item{basis.W1}{The basis of subspace W1;}
#' \item{dW2}{The dimension of subspace W2;}
#' \item{basis.W2}{The basis of subspace W2}
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
#' \donttest{
#' DISCA(X,Y)
#' }


DISCA <- function(X,Y,e=10^-3,
                 xi=200,psi=10,rho=1000,omega=10,abs.epsilon=10^-3,rel.epsilon=10^-3,alpha=0.05,R=1000) {
  # Author: Chuanping Yu
  # Date: Apr 17, 2018


  # The following functions in other packages are going to be used.
  # energy::dcov()
  # MASS:Null()

  p <- ncol(X)
  q <- ncol(Y)
  Num <- nrow(X)

  X_0 <- X
  Y_0 <- Y

# ninv = qnorm(1-alpha/2)^2

  # Normalize X and Y matrix:
  X_0 <- sweep(X_0, 2, colMeans(X_0))
  Y_0 <- sweep(Y_0, 2, colMeans(Y_0))
  X_0 <- sweep(X_0, 2, sqrt(colSums(X_0^2))/sqrt(Num), FUN = "/")
  Y_0 <- sweep(Y_0, 2, sqrt(colSums(Y_0^2))/sqrt(Num), FUN = "/")
  
  X_diff <- Xdiff(X_0)  # The matrix of the pairwise difference of X.
  Y_diff <- Xdiff(Y_0)  # The matrix of the pairwise difference of Y.

  gM_Y <- gMatrix(Y_0)
  g_Y <- gM_Y$gM  # g-Matrix of Y.
# sum4_Y <- gM_Y$sum4  # sum4 of Y.

  U = matrix(0,0,p)  # each row represents a direction.
  U_n = diag(p)  # Null space of U.
  DC_X = 0
  THRES_X = numeric(0)

  X_diff_proj = X_diff
  p_1 = p

  for (DIM in 1:(p-1)){  # Check if the reduced dimension can reach DIM.
    X_diff_proj = X_diff%*%U_n
    u <- minDCOV(X_diff_proj,g_Y)
    u_ult = (1/sqrt(sum((U_n%*%u)^2)))*U_n%*%u    # New direction in U_n.
    test_sum <- dcovU.test(X_0%*%u_ult, Y, R = R, alpha = alpha)
    dc_0 = test_sum$estimate
    DC_X = c(DC_X,dc_0)
    p_1 = p_1-1;
    threshold = test_sum$threshold
    THRES_X = c(THRES_X,threshold)
    if (dc_0>threshold){
      DIM <- DIM-1  # Failure in this attempt.
      break
    }

    U <- rbind(U,t(u_ult))
    U_n <- Null(t(U))

    sprintf('Finished dimension %d of W_1\n',DIM)
  }

  if (DIM==p-1){
    test_sum <- dcovU.test(X_0%*%U_n,Y_0)
    dc_0 = test_sum$estimate
    threshold = test_sum$threshold
    DC_X = c(DC_X,dc_0)
    THRES_X = c(THRES_X,threshold)
    if (dc_0<threshold){
      DIM = p
      U_n = matrix(0,0,0)
    }
  }

  dim1 <- p-DIM

  ############################################################

  V = matrix(0,0,q)  # each row represents a direction.
  V_n = diag(q)  # Null space of V.
  DC_Y = 0
  THRES_Y = numeric(0)

  Y_diff_proj = Y_diff
  q_1 = q

  if (p-DIM==0){
    X_p = X_0
  }else{
    X_p = X_0%*%U_n
  }

  gM_Xp <- gMatrix(X_p)
  g_X <- gM_Xp$gM  # g-Matrix of Xp.
# sum4_X <- gM_Xp$sum4  # sum4 of Xp.

  for (DIM in 1:(q-1)){
    Y_diff_proj = Y_diff%*%V_n
    v = minDCOV(Y_diff_proj,g_X)
    v_ult = (1/sqrt(sum((V_n%*%v)^2)))*V_n%*%v
    test_sum <- dcovU.test(Y_0%*%v_ult, X_0, R = R, alpha = alpha)
    dc_0 = test_sum$estimate
    DC_Y = c(DC_Y,dc_0)
    q_1 = q_1-1
    threshold = test_sum$threshold
    THRES_Y = c(THRES_Y,threshold)
    if (dc_0>threshold){
      DIM <- DIM-1
      break
    }

    V <- rbind(V,t(v_ult))
    V_n <- Null(t(V))

    sprintf('Finished dimension %d of W_2\n',DIM)
  }

  if (DIM==q-1){
    test_sum <- dcovU.test(Y_0%*%v_ult, X_0, R = R, alpha = alpha)
    dc_0 = test_sum$estimate
    threshold = test_sum$threshold
    DC_Y = c(DC_Y,dc_0)
    THRES_Y = c(THRES_Y,threshold)
    if (dc_0<threshold){
      DIM = q
      V_n = matrix(0,0,0)
    }
  }

  dim2 <- q-DIM

  return(list(dW1=dim1,basis.W1=U_n,dW2=dim2,basis.W2=V_n))

}
