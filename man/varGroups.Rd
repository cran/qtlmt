\name{varGroup}
\Rdversion{0.1}
\alias{varGroup}
\title{
Group variables
}
\description{
Group variables via variable selection such that the grouped variables are optimal for multivariate analysis.
}
\usage{
varGroup(x, z, zscope=NULL, k=qf(0.95,1,nrow(x)-2),
   kf=k/2, method=c("pool","best"), direction=c("both",
   "forward","backward"))
}
\arguments{
  \item{x}{
data matrix/frame. Columns are variables to select from.
}
  \item{z}{
data matrix/frame. Columns are variables with elements 1 or not 1 (any others).
}
  \item{zscope}{
which variables in \code{z} to be used for grouping.
}
  \item{k}{
entry/stay value in backward stepwise.
}
  \item{kf}{
entry/stay value in forward stepwise.
}
  \item{method}{
grouping method at each step: pool all the groups selected from each \code{zcope} variable, or choose the largest group (see examples).
}
  \item{direction}{
forward selection, backward elimination or both stepwise.
}
}

\value{
indicators of grouped variables.
}

\seealso{
\code{\link{varSelect}} and \code{\link{varStep}}
}
\examples{
data(etrait)
varGroup(traits, mdat, zscope=c(42,55), method="pool")
varGroup(traits, mdat, zscope=c(42,55), method="best")
}
