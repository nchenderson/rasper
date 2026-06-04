
MaxRankCorr <- function(y, X, external.scores, nu=NULL, maxiter=100,
                        tol=1e-5) {

  nn <- nrow(X)
  external.ranking <- rank(external.scores)
  if(is.null(nu)) {
    nu <- 0.05
  }
  external.ranking <- rank(external.scores)
  wrank <- kronecker(external.ranking, rep(1, nn)) - kronecker(rep(1, nn), external.ranking) > 0
  wrank <- as.vector(wrank)

  Amat <- kronecker(X, rep(1/nu, nn)) - kronecker(rep(1/nu, nn), X)
  objfnvals <- rep(NA, maxiter + 1)
  beta.old <- rep(1.0, ncol(X))
  beta.old <- beta.old/sqrt(sum(beta.old*beta.old))

  betanormfn <- function(kappa, dvec, avec) {
    den <- (dvec/(dvec*dvec + kappa))^2
    ans <- sum(avec*avec*den) - 1
    return(ans)
  }

  objfnvals[1] <- sum(wrank*plogis(as.numeric(Amat%*%beta.old)))
  for(k in 1:maxiter) {
    A.beta <- as.numeric(Amat%*%beta.old)
    WW.tmp <- wrank*plogis(A.beta)
    VV.tmp <- rep(NA, length(WW.tmp))
    small.phi <- abs(A.beta) < 1e-4
    VV.tmp[small.phi] <-  -0.25 + (A.beta[small.phi]^2)/48 - (A.beta[small.phi]^4)/480
    VV.tmp[!small.phi] <-  1/(2*A.beta[!small.phi]) - plogis(A.beta[!small.phi])/A.beta[!small.phi]
    VV.tmp <- (-1)*VV.tmp
    WVec <- WW.tmp/sum(WW.tmp)
    VWVec <- WVec*VV.tmp

    ### Use SVD to compute update:
    svwvec <- sqrt(VWVec)
    A.tmp <- svwvec*Amat
    Y.tmp <- 0.5*sqrt(WVec/VV.tmp)

    SVD_Amat <- svd(A.tmp)
    avec <- crossprod(SVD_Amat$u, Y.tmp)

    kappa_up <- sqrt(p)*max(SVD_Amat$d^2)
    min_d <- min(SVD_Amat$d)
    kappa_down <- abs(min_d) - min_d^2

    ## Find kappa to satisfy beta^t*beta = 1 constraint
    kappa.star <- uniroot(betanormfn, lower=kappa_down, upper=kappa_up,
                          dvec=SVD_Amat$d, avec=avec)$root

    dd_kap <- SVD_Amat$d/(SVD_Amat$d^2 + kappa.star)

    beta.new <- SVD_Amat$v%*%(dd_kap*avec)
    objfnvals[k+1] <- sum(wrank*plogis(as.numeric(Amat%*%beta.new)))

    par.dist <- sqrt(sum((beta.new - beta.old)^2))
    if(par.dist < tol) {
      break
    }
    beta.old <- beta.new
  }
  return(list(coef=as.numeric(beta.new), objfnvals=objfnvals[!is.na(objfnvals)]))
}
