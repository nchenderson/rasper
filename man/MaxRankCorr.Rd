\name{MaxRankCorr}
\alias{MaxRankCorr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Rank Correlation Estimator
}
\description{
Performs estimation of regression coefficients which maximize a smooth approximation
to the rank correlation coefficient between the fitted values of the regression
and a vector of external risk scores.
}
\usage{
MaxRankCorr(y, X, external.scores, nu, maxiter, tol)
}
\arguments{
  \item{y}{
The length-n vector of outcomes.
}
  \item{X}{
An n x p design matrix.
}
  \item{external.scores}{
The vector of n external risk scores
}
  \item{nu}{
Smoothing term used in smooth approximation of maximum rank correlation objective
function. Should be equal to \code{NULL} or a number greater than 0.
}
  \item{maxiter}{
The maximum number of MM algorithm iterations
}
  \item{tol}{
The stopping tolerance
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
TBA
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicholas Henderson
}
\note{
%%  ~~further notes~~
}
\seealso{
\code{\link{rasper}}
}
\examples{
n <- 100
p <- 3
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
Y <- X\%*\%c(1, 0, -1) + rnorm(n)
mu.ext <- X\%*\%c(1, 0.1, -0.3)

mcorr <- MaxRankCorr(y=Y, X=X, external.scores=mu.ext)

## Maximum rank correlation estimate of regression coefficients
mcorr$coef

## Compare fitted values with external scores
mcorr_fitted_vals <- as.numeric(X\%*\%mcorr$coef)
cor(mcorr_fitted_vals, mu.ext)

## Look at objective function across iterations
plot(mcorr$objfnvals)
summary(diff(mcorr$objfnvals))
}
\keyword{regression}
\keyword{nonparametric}
\concept{ rank association }
