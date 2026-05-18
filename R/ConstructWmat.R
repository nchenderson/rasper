ConstructWmat <- function(y, external.ranking, discrepancy) {
  ## The returned "matrix" should be n^2 x n^2
  ## This returns the values weights wij described in the rasper paper
  ## It returns the vectorized version of the matrix
  n <- length(y)
  if(discrepancy=="spearman") {
    Wmat <- rep(external.ranking, each=n)/(4*n*n)
    #Vmat <- lambda/(4*n)*rep(n - external.ranking, each=n)
  } else if(discrepancy=="kendall" & internal.obj == "auc") {
    CompareRanks <- outer(external.ranking, external.ranking, FUN="-") > 0
    CompareY <- outer(y, y) > 0
    Wmat <- (2*lambda/(n-1))*(0.5 - CompareRanks) + CompareY
  }
  return(c(Wmat))
}
