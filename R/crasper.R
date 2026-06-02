
crasper <- function(y, X, Z = NULL, alpha1 = 2, alpha2 = 10, lambda = 500,
                    external.scores, discrepancy = "spearman",
                    nu = NULL, tol = 1e-5, maxiter = 500){

  # Assume X is already normalized(this is for computing nu)
  nn <- nrow(X)
  p <- ncol(X)

  # Create Xa
  if(is.null(Z)){
    Xa <- cbind(X, X)
  } else {
    Xa <- cbind(X, Z)
  }
  # create left identity matrix
  IL <- cbind(diag(p), matrix(0, nrow = p, ncol = p))

  # create 2p x 2p identity matrix
  I.2p <- diag(2*p)

  # Set nu
  if(is.null(nu)) {
    XtX <- crossprod(X, X)
    Xty <- crossprod(X, y)

    ols.beta <- as.numeric(solve(XtX, Xty))
    nu <- 0.05*sqrt(sum(ols.beta*ols.beta))
  }

  # Compute A matrix
  Amat <- kronecker(X, rep(1/nu, nn)) - kronecker(rep(1/nu, nn), X)

  # Compute rankings based on external model
  external.ranking <- rank(external.scores)
  wrank <- rasper:::ConstructWmat(y=y, external.ranking=external.ranking,
                                  discrepancy=discrepancy)

  XatXa <- crossprod(Xa, Xa)
  Xaty <- crossprod(Xa, y)
  XatXap_alpha <- XatXa + nn * alpha2 * I.2p
  A_IL <- Amat %*% IL # Used for computation of Q and c

  tau.old <- rep(0, 2*p)
  theta_indicator <- c(rep(0, p), rep(1, p))

  obj_iter <- c()
  par.dist_iter <- c()

  # Did model converge before maxiter
  failed <- 0
  failcount <- 0

  for(k in 1:maxiter){
    cat(sprintf("\riteration %d/%d", k, maxiter))
    # Compute Wt and Vt
    beta.old <- tau.old[1:p]
    A.beta <- as.numeric(Amat%*%beta.old)
    WW.tmp <- wrank*plogis(A.beta)
    VV.tmp <- rep(NA, length(WW.tmp))
    small.phi <- abs(A.beta) < 1e-4
    VV.tmp[small.phi] <-  -0.25 + (A.beta[small.phi]^2)/48 - (A.beta[small.phi]^4)/480
    VV.tmp[!small.phi] <-  1/(2*A.beta[!small.phi]) - plogis(A.beta[!small.phi])/A.beta[!small.phi]
    WVec <- WW.tmp/sum(WW.tmp) # diag(Wt)
    VWVec <- WVec*VV.tmp       # diag(Vt Wt)

    QQ <- XatXap_alpha - nn*lambda*crossprod(A_IL, diag(VWVec) %*% A_IL)
    R <- chol(QQ)
    cvec <- Xaty + ((nn*lambda)/2)*t(A_IL)%*%WVec
    ytilde <- solve(R, cvec)
    #tau.new.tmp <- solve(t(R)%*%R, R%*%ytilde) # solution with lambda=0

    fit <- glmnet(x=t(R), y=ytilde,
                  family = "gaussian",
                  intercept = F,
                  standardize = F,
                  lambda = alpha1,
                  penalty.factor = theta_indicator,
                  maxit = 1000, thresh=1e-12)
    tau.new <- as.numeric(coef(fit)[-1])

    # glmnet failed to converge count
    failcount <- failcount + as.integer(fit$jerr != 0)

    ### Things that works with 4 hashmarks:
    ####A_IL_W <- Xaty + ((nn*lambda)/2)*t(A_IL)%*%VVec
    ####QQ <- crossprod(A_IL, diag(VWVec) %*% A_IL)
    ####tau.new <- solve(XatXap_alpha - nn*lambda*QQ, A_IL_W)

    if (is.null(Z)){
      obj <- RankPenalizedObjPlus(tau.new, y, X, Amat, wrank, "gaussian", lambda, alpha1, alpha2)
    } else {
      obj <- RankPenalizedObjPlus(tau.new, y, X, Amat, wrank, "gaussian", lambda, alpha1, alpha2, Z = Z)
    }

    # Store lasso loglikelihood
    obj_iter <- c(obj_iter, obj)

    # Check if tolerance is met
    par.dist <- sqrt(sum((tau.old - tau.new)^2))

    par.dist_iter <- c(par.dist_iter, par.dist)

    if(par.dist < tol) {
      break
    }

    # Update Tau
    tau.old <- tau.new

    if (k == maxiter){
      failed <- 1
    }
  }

  if (!is.null(Z)){
    X <- cbind(X, Z)
  }
  return_obj <- list(obj_iter = obj_iter, tau = tau.old, rasper_failed = failed,
                     glmnet_failcount = failcount, X = X)
  return(return_obj)
}
