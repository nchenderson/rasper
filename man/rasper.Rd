\name{rasper}
\alias{rasper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regression with ranking-based penalization
}
\description{
TBA
}
\usage{
rasper(y, X, external.scores, lambda, alpha = 0, internal.obj = "gaussian", discrepancy = "spearman", nu = NULL, maxiter = 500, tol = 1e-05)
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
  \item{lambda}{
The ranking penalty term.
}
  \item{alpha}{
The L2 penalty term
}
  \item{internal.obj}{
TBA
}
  \item{discrepancy}{
TBA
}
  \item{nu}{
TBA
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
nint <- 50   ## number of observations in internal data
p <- 5   ## number of covariates
beta.true <- c(2, -1, 0.5, 0, 0)
beta.ext <- c(3.5, -1.0, -1.3, -2.0, 1.0)

X <- matrix(rnorm(nint*p), nrow=nint, ncol=p)
Xnew <- matrix(rnorm(nint*p), nrow=nint, ncol=p)

## Standardize design matrices
X <- scale(X)
Xnew <- scale(Xnew)

##
mu.true <- as.numeric(X\%*\%beta.true)
mu.true.new <- as.numeric(Xnew\%*\%beta.true)
Y <- X\%*\%beta.true + rnorm(nint, sd=5)
Yprime <- X\%*\%beta.true + rnorm(nint, sd=5)

mu.ext <- as.numeric(X\%*\%beta.ext)
mu.ext.new <- as.numeric(X\%*\%beta.ext)
external.ranking <- rank(as.numeric(mu.ext))/length(as.numeric(mu.ext))
external.ranking.new <- rank(as.numeric(mu.ext.new))/length(as.numeric(mu.ext.new))

lambda.try <- 500
alpha.try <- 10
rshrink <- rasper(y=Y, X=X, external.scores=mu.ext,
                  lambda = lambda.try, alpha = alpha.try,
                  discrepancy="spearman", maxiter=500)
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
