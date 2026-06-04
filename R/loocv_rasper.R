loocv.rasper <- function(y, X, external.scores, lambda.seq, alpha.seq, Xlist=NULL,
                         internal.obj="gaussian", discrepancy="spearman", nu=NULL,
                         maxiter=500) {
  ## Still need to work on this for marginalized ranking parameters.
  nlam <- length(lambda.seq)
  nalpha <- length(alpha.seq)
  loo_score <- aic_vals <- rep(NA, nlam*nalpha)
  alpha_vals <- lam_vals <- rep(NA, nlam*nalpha)

  count <- 1
  if(is.null(Xlist)) {
    for(k in 1:nlam) {
      for(j in 1:nalpha) {
        loo_specific <- loocv_helper(y, X, external.scores, lambda=lambda.seq[k],
                                     alpha=alpha.seq[j], Xlist=NULL, internal.obj=internal.obj,
                                     discrepancy=discrepancy, nu=nu, maxiter=maxiter)

        fit_all <- rasper(y, X, external.scores, lambda=lambda.seq[k],
                          alpha=alpha.seq[j], internal.obj=internal.obj,
                          discrepancy=discrepancy, nu=nu, maxiter=maxiter)
        loo_score[count] <- loo_specific
        aic_vals[count] <- fit_all$aic
        alpha_vals[count] <- alpha.seq[j]
        lam_vals[count] <- lambda.seq[k]
        count <- count + 1
      }
    }
  } else if(!is.null(Xlist)) {
    for(k in 1:nlam) {
      for(j in 1:nalpha) {
        loo_specific <- loocv_helper(y, X, external.scores, lambda=lambda.seq[k],
                                     alpha=alpha.seq[j], Xlist=Xlist, internal.obj=internal.obj,
                                     discrepancy=discrepancy, nu=nu, maxiter=maxiter,
                                     optimization="direct")

        fit_all <- RankingShrinkMarg(y, X, Xlist=Xlist, external.scores, lambda=lambda.seq[k],
                                     alpha=alpha.seq[j], internal.obj=internal.obj,
                                     discrepancy=discrepancy, nu=nu, maxiter=maxiter,
                                     optimization="direct")
        loo_score[count] <- loo_specific
        aic_vals[count] <- fit_all$aic
        alpha_vals[count] <- alpha.seq[j]
        lam_vals[count] <- lambda.seq[k]
        count <- count + 1
      }
    }
  }
  loo_dframe <- data.frame(LooScore=loo_score, aic=aic_vals, lambda=lam_vals, alpha=alpha_vals)
  return(loo_dframe)
}

