ConstructWmat <- function(y, external.ranking, discrepancy) {
  ## The returned "matrix" should be n^2 x n^2
  ## This returns the values weights wij described in the rasper paper
  ## It returns the vectorized version of the matrix
  n <- length(y)
  if(discrepancy=="spearman") {
    Wmat <- rep(external.ranking, each=n)/(4*n*n*n)
  } else if(discrepancy=="kendall") {
    CompareRanks <- kronecker(external.ranking, rep(1, n)) - kronecker(rep(1, n), external.ranking) > 0
    Wmat <- (4/(n*(n-1)))*CompareRanks
  }
  return(c(Wmat))
}
