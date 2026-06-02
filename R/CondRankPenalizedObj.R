
CondRankPenalizedObj <- function(tau, y, X, A, wvec, internal.obj, lambda,
                                 alpha1, alpha2, Z = NULL) {
  n <- length(y)
  p <- length(tau)
  q <- length(tau)/2
  beta_vec <- tau[1:q]
  theta_vec <- tau[(q + 1):p]
  if(internal.obj!="gaussian" & internal.obj!="logistic") {
    stop("Only gaussian or logistic internal objectives are allowed")
  } else if(internal.obj=="gaussian") {
    if (is.null(Z)){
      X.tau <- as.numeric(X%*%beta_vec + X%*%theta_vec)
    } else {
      X.tau <- as.numeric(X%*%beta_vec + Z%*%theta_vec)
    }
    A.beta <- as.numeric(A%*%beta_vec)
    mu.tmp <- plogis(A.beta)

    hI <- sum(X.tau*X.tau)/(2*n) - sum(X.tau*y)/n + alpha1*sum(abs(theta_vec)) + (alpha2/2)*sum(tau*tau)
    objfn.val <- hI - lambda*log(sum(wvec*mu.tmp))
  } else if(internal.obj=="logistic") {

    objfn.val <- 0
  }
  return(objfn.val)
}
