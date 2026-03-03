rasper <- function(y, X, external.scores, lambda,
                   alpha = 0, internal.obj="gaussian",
                   discrepancy="spearman", nu=NULL, maxiter=500,
                   tol=1e-5) {

  ## Assume, for now, that X has been standardized:
  Xnorm <- X
  na.cols <- apply(Xnorm, 2, function(x) sum(is.na(x)))
  nn <- nrow(Xnorm)
  intercept <- FALSE
  #if(sum(na.cols > 0) == 1) {
  #    Xnorm[,na.cols > 0] <- rep(1, nrow(Xnorm))
  #    intercept <- TRUE
  #} else if(sum(na.cols > 0) > 1) {
  #    stop("More than 1 column of X has no variation")
  #}

  p <- ncol(Xnorm)
  XtX <- crossprod(Xnorm, Xnorm)
  XtXp.alpha <- XtX
  diag(XtXp.alpha) <- diag(XtXp.alpha) + alpha
  Xty <- crossprod(Xnorm, y)
  if(is.null(nu)) {
    ols.beta <- as.numeric(solve(XtX, Xty))
    nu <- 0.05*sqrt(sum(ols.beta*ols.beta))
  }
  external.ranking <- rank(external.scores)
  wrank <- ConstructVmat(y=y, external.ranking=external.ranking,
                         discrepancy=discrepancy, internal.obj=internal.obj)

  Amat <- kronecker(Xnorm, rep(1/nu, nn)) - kronecker(rep(1/nu, nn), Xnorm)
  objfnvals <- rep(NA, maxiter + 1)
  beta.old <- rep(0.0, ncol(Xnorm))

  objfnvals[1] <- RankPenalizedObj(par=beta.old, y=y, X=Xnorm, A=Amat, wvec=wrank,
                                   internal.obj=internal.obj, lambda=lambda,
                                   alpha=alpha)

  for(k in 1:maxiter) {
     A.beta <- as.numeric(Amat%*%beta.old)
     VV.tmp <- wrank*plogis(A.beta)
     WW.tmp <- rep(NA, length(VV.tmp))
     small.phi <- abs(A.beta) < 1e-4
     WW.tmp[small.phi] <-  (1/4 - (A.beta[small.phi]^2)/48 + (A.beta[small.phi]^4)/480)
     WW.tmp[!small.phi] <-  1/(2*A.beta[!small.phi]) - plogis(A.beta[!small.phi])/A.beta[!small.phi]
     VVec <- VV.tmp/sum(VV.tmp)
     VWVec <- VVec*WW.tmp

     A.sandwich <- t(VWVec*Amat)%*%Amat
     AtV <- crossprod(Amat, VVec)

     X1 <- XtXp.alpha - lambda*A.sandwich
     X2 <- (lambda/2)*AtV + Xty

     beta.new <- as.numeric(solve(X1, X2))

     objfnvals[k+1] <- RankPenalizedObj(par=beta.new, y=y, X=Xnorm, A=Amat, wvec=wrank,
                                         internal.obj=internal.obj, lambda=lambda,
                                         alpha=alpha)
     par.dist <- sqrt(sum((beta.new - beta.old)^2))
     if(par.dist < tol) {
       break
     }
     beta.old <- beta.new
  }

  num.iter <- k
  wrank.normalized <- wrank/sum(wrank)
  X1 <- XtXp.alpha - (lambda/4)*A.sandwich
  df.mat <- solve(X1, XtX)
  approx.df <- sum(diag(df.mat))

  fitted.vals <- NULL
  aic <- NULL
  if(internal.obj=="gaussian") {
    fitted.vals <- as.vector(Xnorm%*%beta.new)
    rss <- mean((y - fitted.vals)^2)
    aic <- nn*log(rss) + 2*approx.df
  }
  return(list(coef=beta.new, objfnvals=objfnvals[!is.na(objfnvals)], fitted.values=fitted.vals,
              df=approx.df, aic=aic, num.iter=num.iter))
}


