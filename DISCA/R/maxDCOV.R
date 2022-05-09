#' Maximize the distance covariance of \eqn{Xu} and \eqn{Y} under certain constraint.
#'
#' @export
#' @importFrom stats runif
#' @description \code{maxDCOV} is used to solve the non-convex optimization problem of
#' maximizing the distance covariance of \eqn{Xu} and \eqn{Y} under the constraint that
#' the Euclidean norm of \eqn{u} is 1.
#'
#' @param Xd \eqn{N(N-1)/2-by-p} matrix. The matrix of the pairwise difference of X.
#' Each row represents one observation. It can be obtained by function \code{Xdiff}.
#' @param gM \eqn{N(N-1)/2-by-q} matrix. The g-Matrix of Y.
#' It can be obtained by function \code{gMatrix}.
#' @param u0 The initial guess of the solution. The default value is (1,0,...,0).
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
#'
#' @return The solution \eqn{u}.
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
#' \donttest{
#' maxDCOV(X_diff,g_Y)
#' }


maxDCOV <- function(Xd,gM,u0=c(1,rep(0,ncol(Xd)-1)),e=10^-3,
                    xi=200,psi=10,rho=1000,omega=10,abs.epsilon=10^-3,rel.epsilon=10^-3) {
  # Author: Chuanping Yu
  # Date: Apr 17, 2018
  
  
  p <- ncol(Xd)
  
  M_plus <- Mminus(gM,Xd)  # The M plus matrix.
  M_minus <- Mplus(gM,Xd)  # The M minus matrix.
  
  # Normalize u0 if the l2 norm of u0 is not equal to 1.
  if (sum(u0^2)==0){
    u0 = c(1,rep(0,p-1))
  }
  if (sqrt(sum(u0^2))!=1){
    u0 <- (1/sqrt(sum(u0^2)))*u0
  }
  
  index <- 0
  
  u_pre <- u0
  u <- rep(0,p)
  
  while (min(sqrt(sum((u-u_pre)^2)),sqrt(sum((u+u_pre)^2)))/max(sqrt(sum(u^2)),1)>e){
    index <- index+1
    u <- u_pre
    sig <- sign(M_minus%*%u)
    partial <- sig
    sig_zero <- sig[sig==0]
    if (sum(sig==0)>0){
      sig_zero <- runif(length(sig_zero),min=-1,max=1)
      partial[sig==0] = sig_zero
    }
    if (sum(abs(u-rep(0,p)))==0){
      y = t(M_minus)%*%partial
    }else{
      y = ((xi-psi)/sqrt(sum(u^2)))*u+t(M_minus)%*%partial
    }
    z = M_plus%*%u;
    v = 10*rep(1,nrow(M_plus))
    for (l in 1:10^5){
      u_l = solve(xi*diag(p)+rho*(t(M_plus)%*%M_plus))%*%(y+t(M_plus)%*%(rho*z-v))
      x = (1/rho)*v+M_plus%*%u_l
      z_l = sign(x)*apply(cbind(abs(x)-1/rho,rep(0,length(x))),1,max)
      v = v+rho*(M_plus%*%u_l-z_l)
      r = M_plus%*%u_l-z_l
      s = rho*t(M_plus)%*%(z_l-z)
      if (sqrt(sum(r^2))<=sqrt(nrow(M_plus))*abs.epsilon+rel.epsilon*max(sqrt(sum((M_plus%*%u_l)^2)),sqrt(sum(z_l^2))) &&
          sqrt(sum(s^2))<=sqrt(p)*abs.epsilon+rel.epsilon*sqrt(sum(t(M_plus)%*%v)^2)){
        break
      }
      z <- z_l
    }
    u_pre <- u_l
    psi <- psi + xi*(sqrt(sum(u_pre^2))-1);
    xi <- omega*xi
  }
  
  return(u_pre)
  
}
