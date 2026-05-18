RankPenalizedObj <- function(par, y, X, A, wvec, internal.obj, lambda, alpha) {
  n <- length(y)
  if(internal.obj!="gaussian" & internal.obj!="logistic") {
    stop("Only gaussian or logistic internal objectives are allowed")
  } else if(internal.obj=="gaussian") {
    X.beta <- as.numeric(X%*%par)
    A.beta <- as.numeric(A%*%par)
    mu.tmp <- plogis(A.beta)

    hI <- sum(X.beta*X.beta)/2 - sum(X.beta*y) + (alpha/2)*sum(par*par)
    objfn.val <- hI - lambda*log(sum(wvec*mu.tmp))
  } else if(internal.obj=="logistic") {
    X.beta <- X%*%par
    theta.hat <- as.numeric(X.beta)/nu
    A.beta <- c(outer(theta.hat, theta.hat, FUN="-"))
    mu.tmp <- plogis(A.beta)

    hI <- sum(y*X.beta) - sum(log(1 + exp(X.beta)))
    objfn.val <- hI/n + sum(V*mu.tmp)
  }
  return(objfn.val)
}
