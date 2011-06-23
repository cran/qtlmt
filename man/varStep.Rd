\name{varStep}
\Rdversion{0.1}
\alias{varAdd1}
\alias{varDrop1}
\alias{varStep}
\title{
Variable selection
}
\description{
Add a variable, drop a variable, or select a subset of variables via variable selection that are optimal for multivariate analysis.
}
\usage{
varAdd1(x, grp, vin=NULL, scope=1:ncol(x), k=0)

varDrop1(x, grp, vin=1:ncol(x), k=0)

varStep(x, grp, scope, k, kf=k/2, direction=c("both",
   "forward","backward"))
}
\arguments{
  \item{x}{
data matrix/frame. Columns are variables to select from.
}
  \item{grp}{
grouping indicator of observations.
}
  \item{vin}{
which variables (i.e., columns of \code{x}) already in model.
}
  \item{scope}{
which variables (i.e., columns of \code{x}) to select from.
}
  \item{k}{
entry/stay value in backward stepwise.
}
  \item{kf}{
entry/stay value in forward stepwise.
}
  \item{direction}{
forward selection, backward elimination or both stepwise.
}
}

\value{
which variable to add (add1), which variable to drop (drop1), or a subset of variables in the final model (step).
}

\seealso{
\code{\link{varGroup}} and \code{\link{varSelect}}
}
\examples{
data(etrait)
varAdd1(traits, grp=mdat[,42], vin=10, scope=1:ncol(traits))
varStep(traits, grp=mdat[,42], k=12, scope=1:ncol(traits),
   direction="back")
}

