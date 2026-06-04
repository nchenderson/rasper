\name{MaxRankCorr}
\alias{MaxRankCorr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum Rank Correlation
}
\description{
Performs estimation of regression coefficients which maximize a smooth approximation
of the maximum rank correlation coefficient.
}
\usage{
MaxRankCorr(y, X, external.scores, nu, maxiter, tol)
}
%- maybe also 'usage' for other objects documented here.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
n <- 100
p <- 3
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
Y <- X\%*\%c(1, 0, -1) + rnorm(n)
mu.ext <- X\%*\%c(1, 0.1, -0.3)

mcorr <- MaxRankCorr(y=Y, X=X, external.scores=mu.ext)

## Look at objective function
mcorr$coef

## Look at objective function across iterations
plot(mcorr$objfnvals)
summary(diff(mcorr$objfnvals))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
