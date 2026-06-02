loocv_helper <- function(y, X, external.scores, lambda,
                         alpha, Xlist, internal.obj="gaussian",
                         discrepancy="spearman", nu=NULL, maxiter=500) {

  nobs <- length(y)
  mse.val <- rep(NA, nobs)
  if(!is.null(Xlist)) {
    nlist <- length(Xlist)
  }
  for(k in 1:nobs) {
    y.train <- y[-k]
    X.train <- X[-k,,drop=FALSE]
    X.test <- matrix(X[k,], nrow=1, ncol=ncol(X))
    y.test <- y[k]
    ext.scores.train <- external.scores[-k]

    if(is.null(Xlist)) {
      rshrink <- RankingShrink(y=y.train, X=X.train, external.scores=ext.scores.train,
                               lambda = lambda, alpha = alpha,
                               internal.obj=internal.obj, discrepancy=discrepancy,
                               nu=nu, maxiter=maxiter, optimization="direct")
    } else if(!is.null(Xlist)) {
      Xlist.train <- list()
      for(h in 1:nlist) {
        Xlist.train[[h]] <- Xlist[[h]][-k,,drop=FALSE]
      }
      rshrink <- RankingShrinkMarg(y=y.train, X=X.train, Xlist=Xlist.train,
                                   external.scores=ext.scores.train,
                                   lambda = lambda, alpha = alpha,
                                   internal.obj=internal.obj, discrepancy=discrepancy,
                                   nu=nu, maxiter=maxiter, optimization="direct")
    }
    rfit <- as.numeric(X.test%*%rshrink$coef)
    mse.val[k] <- (y.test - rfit)*(y.test - rfit)
  }
  return(mean(mse.val))
}
